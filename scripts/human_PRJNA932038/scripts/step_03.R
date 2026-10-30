cfg <- load_config(); root <- project_root(); attempt <- Sys.getenv("ASC_P_ATTEMPT_DIR")
assert_packages(c("Seurat", "harmony", "clustree", "dplyr", "ggplot2", "patchwork"))
set.seed(cfg$seed)
prev <- read_previous_checkpoint(3)

# Step 03分为两个阶段：
# 1）本脚本只完成Harmony和resolution网格扫描，生成clustree后停止；
# 2）用户人工确认最终resolution后，再补充最终聚类和后续输出。
# 注意：N/T是需要保留的生物学条件，不要求两组在UMAP上完全重合。
# Harmony按sample_id减少文库/样本技术差异，但是否过度校正需结合marker判断。
resolution_grid <- c(0.2, 0.3, 0.4, 0.5, 0.6)

# 使用Step 02批准的QC + singlet对象创建Seurat对象
obj <- Seurat::CreateSeuratObject(prev$counts, project = "PRJNA932038",
                                   meta.data = prev$meta, min.cells = 0, min.features = 0)
obj <- Seurat::NormalizeData(obj, verbose = FALSE)

# 每个样本分别寻找高变基因，再选择共同的integration features，
# 避免高变基因完全被某一个样本主导。
split <- Seurat::SplitObject(obj, split.by = "sample_id")
split <- lapply(split, Seurat::FindVariableFeatures, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
features <- Seurat::SelectIntegrationFeatures(split, nfeatures = 3000)
Seurat::VariableFeatures(obj) <- features
obj <- Seurat::ScaleData(obj, features = features, verbose = FALSE)
obj <- Seurat::RunPCA(obj, features = features, npcs = cfg$integration$npcs, verbose = FALSE, seed.use = cfg$seed)

# 保存未经校正的PCA-UMAP作为参照，不能只看Harmony结果。
obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = cfg$integration$dims, reduction.name = "umap_pca", seed.use = cfg$seed, verbose = FALSE)

# Harmony校正和Harmony-UMAP
obj <- harmony::RunHarmony(obj, group.by.vars = cfg$integration$harmony_group, reduction.use = "pca",
                           dims.use = cfg$integration$dims, reduction.save = "harmony", verbose = FALSE)
obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = cfg$integration$dims, reduction.name = "umap", seed.use = cfg$seed, verbose = FALSE)
obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = cfg$integration$dims, verbose = FALSE)

# 在同一个RNA_snn图上扫描候选分辨率；此处不选择最终resolution。
obj <- Seurat::FindClusters(
  obj,
  graph.name = "RNA_snn",
  resolution = resolution_grid,
  random.seed = cfg$seed,
  verbose = FALSE
)

resolution_columns <- paste0("RNA_snn_res.", resolution_grid)
missing_resolution_columns <- setdiff(resolution_columns, colnames(obj[[]]))
if (length(missing_resolution_columns) > 0) {
  stop("FindClusters未生成预期的resolution列: ",
       paste(missing_resolution_columns, collapse = ", "))
}

# 导出每个resolution对应的cluster数量，辅助人工判断是否出现稳定平台或过度碎裂。
resolution_summary <- data.frame(
  resolution = resolution_grid,
  metadata_column = resolution_columns,
  n_clusters = vapply(
    resolution_columns,
    function(x) length(unique(obj[[x, drop = TRUE]])),
    integer(1)
  )
)
write.csv(
  resolution_summary,
  file.path(attempt, "tables", "resolution_cluster_counts.csv"),
  row.names = FALSE
)

# clustree展示相邻resolution之间各cluster如何分裂和继承。
p_resolution <- clustree::clustree(
  obj,
  prefix = "RNA_snn_res."
) +
  # ggraph默认把edge guide写成字符串"edge_colourbar"；在当前
  # ggplot2 3.5.2 + ggraph 2.2.2组合中，字符串指南无法被正确注册，
  # 会在ggsave时报告Unknown guide: edge_colourbar。
  # 显式传入guide函数可保持相同图形含义并解决兼容性问题。
  ggraph::scale_edge_colour_continuous(
    guide = ggraph::guide_edge_colourbar()
  )
save_plot_both(
  p_resolution,
  file.path(attempt, "figures", "clustree_resolution_0.2_to_0.6"),
  13,
  9
)

# Harmony前后对照图：N/T不完全重合可以是合理生物学差异，
# 审核重点是同类细胞是否仍被单一样本技术性分离。
p_pca <- Seurat::DimPlot(obj, reduction = "umap_pca", group.by = "sample_id") + ggplot2::ggtitle("Uncorrected PCA UMAP by sample")
p_hs <- Seurat::DimPlot(obj, reduction = "umap", group.by = "sample_id") + ggplot2::ggtitle("Harmony UMAP by sample")
p_hg <- Seurat::DimPlot(obj, reduction = "umap", group.by = "group") + ggplot2::ggtitle("Harmony UMAP by N/T")
save_plot_both(p_pca + p_hs + p_hg, file.path(attempt, "figures", "integration_comparison_before_resolution_selection"), 15, 6)

# ====== PHASE 2: Lock final resolution (0.4) ======
# User confirmed resolution 0.4 based on clustree plateau and UMAP review.
final_res <- 0.4
obj$seurat_clusters <- obj[["RNA_snn_res.0.4", drop = TRUE]]
SeuratObject::Idents(obj) <- "seurat_clusters"

# --- Final cluster UMAP ---
p_cl <- Seurat::DimPlot(obj, reduction = "umap", label = TRUE) + ggplot2::ggtitle("Final clusters (res=0.4, 19 clusters)")
save_plot_both(p_cl, file.path(attempt, "figures", "final_clusters_UMAP_res_0.4"), 8, 7)

# --- Cluster composition ---
composition <- as.data.frame(table(cluster = obj$seurat_clusters, sample_id = obj$sample_id), stringsAsFactors = FALSE)
write.csv(composition, file.path(attempt, "tables", "cluster_by_sample_counts.csv"), row.names = FALSE)

# --- Save final clustered object ---
checkpoint <- file.path(attempt, "objects", "PRJNA932038_global_harmony_clustered.rds")
saveRDS(obj, checkpoint)

write_waiting_review(attempt, 3, c(
  "## Resolution selection", "",
  "- Resolution 0.4 selected (19 clusters) based on clustree plateau (0.3->0.4: +1 cluster) and UMAP review.",
  "- Cluster 4 shows biological heterogeneity (3-4 UMAP sub-groups) but was kept as one unit; stromal subsetting in Step 05 will refine.",
  "- Confirm uncorrected vs Harmony UMAPs show adequate batch correction without N/T erasure.",
  "- Approve to proceed to Step 04 marker-based annotation."
))
write_step_manifest(3, attempt, checkpoint, file.path(root, "scripts", "step_03.R"),
                    extra = list(final_resolution = final_res, resolution_grid = resolution_grid))

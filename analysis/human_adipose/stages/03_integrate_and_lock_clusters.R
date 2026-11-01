cfg <- load_config(); root <- project_root(); attempt <- Sys.getenv("ASC_P_ATTEMPT_DIR")
assert_packages(c("Seurat", "harmony", "clustree", "dplyr", "ggplot2", "patchwork"))
set.seed(cfg$seed)
prev <- read_previous_checkpoint(3)

# Integrate samples and inspect the clustering-resolution grid
# Generate Harmony embeddings and clustree plots
# Final figure generation uses the selected resolution recorded in the configuration
# N/T is a biological condition; complete overlap in UMAP is not required
# Harmony corrects by sample_id; inspect markers for possible overcorrection
resolution_grid <- c(0.2, 0.3, 0.4, 0.5, 0.6)

# Create the Seurat object from the QC-filtered singlet checkpoint
obj <- Seurat::CreateSeuratObject(prev$counts, project = "PRJNA932038",
                                   meta.data = prev$meta, min.cells = 0, min.features = 0)
obj <- Seurat::NormalizeData(obj, verbose = FALSE)

# Select variable genes per sample and derive shared integration features
# This limits domination by a single sample
split <- Seurat::SplitObject(obj, split.by = "sample_id")
split <- lapply(split, Seurat::FindVariableFeatures, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
features <- Seurat::SelectIntegrationFeatures(split, nfeatures = 3000)
Seurat::VariableFeatures(obj) <- features
obj <- Seurat::ScaleData(obj, features = features, verbose = FALSE)
obj <- Seurat::RunPCA(obj, features = features, npcs = cfg$integration$npcs, verbose = FALSE, seed.use = cfg$seed)

# Retain an uncorrected PCA-based UMAP for comparison
obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = cfg$integration$dims, reduction.name = "umap_pca", seed.use = cfg$seed, verbose = FALSE)

# Run Harmony and compute the Harmony-based UMAP
obj <- harmony::RunHarmony(obj, group.by.vars = cfg$integration$harmony_group, reduction.use = "pca",
                           dims.use = cfg$integration$dims, reduction.save = "harmony", verbose = FALSE)
obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = cfg$integration$dims, reduction.name = "umap", seed.use = cfg$seed, verbose = FALSE)
obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = cfg$integration$dims, verbose = FALSE)

# Scan candidate resolutions on the same RNA_snn graph
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

# Export cluster counts across resolutions
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

# Visualize cluster transitions across adjacent resolutions
p_resolution <- clustree::clustree(
  obj,
  prefix = "RNA_snn_res."
) +
  # Use an explicit edge-guide function for compatibility
  # The ggplot2 3.5.2 / ggraph 2.2.2 combination may not register
  # the string guide name edge_colourbar during ggsave
  # Passing the guide function preserves the intended legend
  ggraph::scale_edge_colour_continuous(
    guide = ggraph::guide_edge_colourbar()
  )
save_plot_both(
  p_resolution,
  file.path(attempt, "figures", "clustree_resolution_0.2_to_0.6"),
  13,
  9
)

# Compare embeddings before and after Harmony correction
# Inspect whether the same cell type remains separated by sample
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
write_step_manifest(3, attempt, checkpoint, file.path(root, "stages", "03_integrate_and_lock_clusters.R"),
                    extra = list(final_resolution = final_res, resolution_grid = resolution_grid))

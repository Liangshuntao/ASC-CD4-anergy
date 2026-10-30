############################################################
## PRJNA932038 人乳腺脂肪组织单核RNA测序（snRNA-seq）
## Step 01：数据读取、逐样本质量控制和逐样本scDblFinder双细胞检测
##
## 重要说明
## - 原始Cell Ranger矩阵和上一版分析结果均只读，不会被覆盖。
## - 双细胞检测按照每个10x文库/样本分别运行，不能把6个样本合并后再检测。
## - 默认运行到singlet（单细胞核）检查点后停止，等待人工核对。
## - 只有人工核对tables/01_scDblFinder_doublet_summary_by_sample.csv后，
##   才能把RUN_DOWNSTREAM改为TRUE继续归一化、PCA和聚类。
## - 本脚本下游部分仍是“未经批次校正的探索性分析”，不能直接用于正式注释结论。
############################################################

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
})

## 1. 检查scDblFinder所需的Bioconductor包；这里只检查，不在运行时自动安装
required_bioc <- c("SingleCellExperiment", "SummarizedExperiment", "scDblFinder")
missing_bioc <- required_bioc[
  !vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  stop(
    "Missing required packages: ", paste(missing_bioc, collapse = ", "),
    ". Install them with BiocManager::install(c(",
    paste(sprintf("'%s'", missing_bioc), collapse = ", "), "))."
  )
}

## 2. 用户可调整参数 -----------------------------------------------------------
SEED <- 1234

## FALSE：只运行到singlet检查点（当前推荐设置）
## TRUE：人工批准doublet结果后才允许继续下游分析
RUN_DOWNSTREAM <- FALSE

## 在核对singlet-only clustree、marker一致性和每位患者的cluster支持度前，
## 不预先指定最终resolution。NA表示“尚未人工选择”。
FINAL_RESOLUTION <- NA_real_

project_dir <- "D:/Users/lstsw/Documents/ASC_CD4/final_result"
data_dir <- file.path(project_dir, "final_result")

## 3. 输入和输出路径 -----------------------------------------------------------
## 使用全新的输出目录，防止覆盖上一版SIMPLE分析结果
out_dir <- file.path(project_dir, "ASC_P_human_validation_simple_scDblFinder")
dirs <- file.path(out_dir, c("figures", "figures/QC_by_sample",
                            "figures/doublet", "figures/UMAP",
                            "objects", "tables", "logs"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

## 如果已经生成过singlet检查点，则立即停止，避免误覆盖已有结果
singlet_checkpoint <- file.path(
  out_dir, "objects/PRJNA932038_merged_after_samplewise_QC_singlets_only.rds"
)
if (file.exists(singlet_checkpoint)) {
  stop(
    "Singlet checkpoint already exists and will not be overwritten: ",
    singlet_checkpoint,
    "\nRename the existing output directory before starting a new attempt."
  )
}

## 4. 六个样本的患者、N/T分组、SRA编号和Cell Ranger矩阵路径
## N：normal breast adipose tissue；T：tumor-adjacent adipose tissue
sample_info <- tibble::tribble(
  ~sample_id, ~patient_id, ~group, ~sra_run,       ~matrix_path,
  "S53_N",    "S53",       "N",    "SRR23526669", file.path(data_dir, "s1_SRR23526669", "filtered_feature_bc_matrix"),
  "S53_T",    "S53",       "T",    "SRR23526668", file.path(data_dir, "s2_SRR23526668", "filtered_feature_bc_matrix"),
  "S76_N",    "S76",       "N",    "SRR23526667", file.path(data_dir, "s3_SRR23526667", "filtered_feature_bc_matrix"),
  "S76_T",    "S76",       "T",    "SRR23526666", file.path(data_dir, "s4_SRR23526666", "filtered_feature_bc_matrix"),
  "S81_N",    "S81",       "N",    "SRR23526665", file.path(data_dir, "s5_SRR23526665", "filtered_feature_bc_matrix"),
  "S81_T",    "S81",       "T",    "SRR23526664", file.path(data_dir, "s6_SRR23526664", "filtered_feature_bc_matrix")
)

## 5. 保留上一版人工确定的逐样本QC阈值
## 注意：这些阈值仍需通过QC图和敏感性分析进行人工审计，本次不暗中修改
sample_qc_cutoffs <- tibble::tribble(
  ~sample_id, ~min_features, ~max_features, ~min_counts, ~max_counts, ~max_mt, ~max_HB,
  "S53_N",    400,           10000,         500,         65000,       5,       1,
  "S53_T",    400,           9000,          500,         52000,       5,       1,
  "S76_N",    450,           7500,          550,         40000,       5,       1,
  "S76_T",    420,           7000,          540,         40000,       5,       1,
  "S81_N",    600,           6000,          700,         21000,       5,       1,
  "S81_T",    400,           4500,          500,         18000,       5,       1
)

stopifnot(
  nrow(sample_info) == 6,
  !anyDuplicated(sample_info$sample_id),
  setequal(sample_info$sample_id, sample_qc_cutoffs$sample_id)
)

missing_paths <- sample_info$matrix_path[!dir.exists(sample_info$matrix_path)]
if (length(missing_paths) > 0) {
  stop("Missing matrix directories:\n", paste(missing_paths, collapse = "\n"))
}

write.csv(sample_info, file.path(out_dir, "tables/sample_info.csv"), row.names = FALSE)
write.csv(sample_qc_cutoffs,
          file.path(out_dir, "tables/01_sample_QC_cutoffs_used.csv"),
          row.names = FALSE)

## 6. 读取10x矩阵；若存在多种feature类型，只保留Gene Expression
read_gene_expression <- function(path) {
  x <- Read10X(data.dir = path, unique.features = TRUE)
  if (is.list(x)) {
    if (!"Gene Expression" %in% names(x)) {
      stop("Read10X returned multiple assays but no 'Gene Expression': ", path)
    }
    x <- x[["Gene Expression"]]
  }
  x
}

## 兼容Seurat v4和v5的原始counts提取函数；scDblFinder必须使用原始counts
get_counts_matrix <- function(obj) {
  if (inherits(obj[["RNA"]], "Assay5")) {
    LayerData(obj, assay = "RNA", layer = "counts")
  } else {
    GetAssayData(obj, assay = "RNA", slot = "counts")
  }
}

## 每张核查图同时保存PDF和300 dpi PNG
save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300)
}

## 7. 逐样本创建Seurat对象并计算QC指标 -----------------------------------------
scRNAlist_raw <- list()
pre_qc_summary <- list()
qc_quantiles <- list()

for (i in seq_len(nrow(sample_info))) {
  info <- sample_info[i, ]
  sid <- info$sample_id
  message("Importing ", sid, " (", info$sra_run, ")")

  counts <- read_gene_expression(info$matrix_path)
  pre_qc_summary[[sid]] <- data.frame(
    sample_id = sid,
    patient_id = info$patient_id,
    group = info$group,
    sra_run = info$sra_run,
    n_genes_matrix = nrow(counts),
    n_cells_matrix = ncol(counts)
  )

  ## 为忠实保留上一版SIMPLE脚本，这里继续使用min.cells=3、min.features=300。
  ## 因此该对象应称为“scDblFinder前对象”，并非完全未经筛选的原始对象。
  obj <- CreateSeuratObject(
    counts = counts, project = sid, min.cells = 3, min.features = 300
  )
  obj$sample_id <- sid
  obj$patient_id <- info$patient_id
  obj$group <- info$group
  obj$sra_run <- info$sra_run
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  ## 计算血红蛋白相关基因比例，用于辅助识别红细胞污染
  hb_genes <- intersect(
    c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"),
    rownames(obj)
  )
  obj[["percent.HB"]] <- if (length(hb_genes) > 0) {
    PercentageFeatureSet(obj, features = hb_genes)
  } else {
    rep(0, ncol(obj))
  }

  md <- obj[[]]
  qc_quantiles[[sid]] <- data.frame(
    sample_id = sid,
    n_cells_pre_doublet = ncol(obj),
    nFeature_q01 = unname(quantile(md$nFeature_RNA, 0.01)),
    nFeature_median = median(md$nFeature_RNA),
    nFeature_q99 = unname(quantile(md$nFeature_RNA, 0.99)),
    nCount_q01 = unname(quantile(md$nCount_RNA, 0.01)),
    nCount_median = median(md$nCount_RNA),
    nCount_q99 = unname(quantile(md$nCount_RNA, 0.99)),
    percentMT_median = median(md$percent.mt),
    percentMT_q99 = unname(quantile(md$percent.mt, 0.99)),
    percentHB_median = median(md$percent.HB),
    percentHB_q99 = unname(quantile(md$percent.HB, 0.99))
  )

  ## 绘制每个样本的nFeature、nCount、线粒体比例和血红蛋白比例
  qc_long <- tidyr::pivot_longer(
    cbind(md, sample_id = sid),
    cols = c(nFeature_RNA, nCount_RNA, percent.mt, percent.HB),
    names_to = "metric", values_to = "value"
  )
  p_qc <- ggplot(qc_long, aes(sample_id, value)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.2) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    theme_classic() + labs(title = paste0(sid, " pre-doublet QC"), x = NULL)
  save_plot_both(
    p_qc,
    file.path(out_dir, "figures/QC_by_sample", paste0(sid, "_QC")),
    8, 7
  )

  scRNAlist_raw[[sid]] <- obj
}

write.csv(bind_rows(pre_qc_summary),
          file.path(out_dir, "tables/01_pre_qc_summary.csv"), row.names = FALSE)
write.csv(bind_rows(qc_quantiles),
          file.path(out_dir, "tables/01_QC_quantiles_by_sample.csv"), row.names = FALSE)

## 8. 逐样本QC过滤并运行scDblFinder -------------------------------------------
## 基本QC在前，scDblFinder在后；每个样本独立运行
scRNAlist_qc <- list()
scRNAlist_singlet <- list()
qc_summary <- list()
doublet_summary <- list()

for (sid in sample_info$sample_id) {
  message("QC + scDblFinder: ", sid)
  obj <- scRNAlist_raw[[sid]]
  cutoff <- sample_qc_cutoffs[sample_qc_cutoffs$sample_id == sid, ]
  stopifnot(nrow(cutoff) == 1)

  n_before <- ncol(obj)
  obj_qc <- subset(
    obj,
    subset = nFeature_RNA > cutoff$min_features &
      nFeature_RNA < cutoff$max_features &
      nCount_RNA > cutoff$min_counts &
      nCount_RNA < cutoff$max_counts &
      percent.mt < cutoff$max_mt &
      percent.HB < cutoff$max_HB
  )
  n_after_qc <- ncol(obj_qc)

  qc_summary[[sid]] <- data.frame(
    sample_id = sid, n_before_QC = n_before, n_after_QC = n_after_qc,
    n_removed_by_QC = n_before - n_after_qc,
    percent_removed_by_QC = round(100 * (n_before - n_after_qc) / n_before, 2)
  )

  ## 使用QC后的原始counts构建SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = get_counts_matrix(obj_qc)),
    colData = S4Vectors::DataFrame(obj_qc[[]])
  )
  set.seed(SEED)
  ## 当前每次循环仅包含一个10x文库，因此无需再指定samples列
  sce <- scDblFinder::scDblFinder(sce)

  dbl_md <- as.data.frame(SummarizedExperiment::colData(sce))
  required_cols <- c("scDblFinder.class", "scDblFinder.score")
  if (!all(required_cols %in% colnames(dbl_md))) {
    stop("scDblFinder output columns missing for ", sid)
  }
  ## 必须保证scDblFinder输出顺序与Seurat细胞顺序完全一致
  if (!identical(rownames(dbl_md), colnames(obj_qc))) {
    stop("Cell order mismatch after scDblFinder for ", sid)
  }

  obj_qc$scDblFinder.class <- as.character(dbl_md$scDblFinder.class)
  obj_qc$scDblFinder.score <- as.numeric(dbl_md$scDblFinder.score)
  n_doublet <- sum(obj_qc$scDblFinder.class == "doublet")
  n_singlet <- sum(obj_qc$scDblFinder.class == "singlet")
  if (n_doublet + n_singlet != ncol(obj_qc)) {
    stop("Unexpected scDblFinder classes for ", sid)
  }

  doublet_summary[[sid]] <- data.frame(
    sample_id = sid, n_after_QC = n_after_qc,
    n_singlet = n_singlet, n_doublet = n_doublet,
    doublet_rate_percent = round(100 * n_doublet / n_after_qc, 2)
  )

  cell_table <- data.frame(
    cell_barcode = colnames(obj_qc), sample_id = sid,
    scDblFinder_class = obj_qc$scDblFinder.class,
    scDblFinder_score = obj_qc$scDblFinder.score,
    nFeature_RNA = obj_qc$nFeature_RNA,
    nCount_RNA = obj_qc$nCount_RNA,
    percent.mt = obj_qc$percent.mt,
    percent.HB = obj_qc$percent.HB
  )
  write.csv(cell_table,
            file.path(out_dir, "tables", paste0("01_scDblFinder_classification_", sid, ".csv")),
            row.names = FALSE)

  p_dbl <- ggplot(obj_qc[[]],
                  aes(scDblFinder.class, scDblFinder.score,
                      fill = scDblFinder.class)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.2) +
    theme_classic() + theme(legend.position = "none") +
    labs(title = paste0(sid, " scDblFinder"), x = NULL)
  save_plot_both(
    p_dbl,
    file.path(out_dir, "figures/doublet", paste0("01_scDblFinder_score_", sid)),
    5, 5
  )

  ## 同时保存含doublet标签的QC对象和仅保留singlet的对象
  scRNAlist_qc[[sid]] <- obj_qc
  scRNAlist_singlet[[sid]] <- subset(
    obj_qc, subset = scDblFinder.class == "singlet"
  )
}

qc_summary <- bind_rows(qc_summary)
doublet_summary <- bind_rows(doublet_summary)
write.csv(qc_summary,
          file.path(out_dir, "tables/01_QC_filter_summary_by_sample.csv"),
          row.names = FALSE)
write.csv(doublet_summary,
          file.path(out_dir, "tables/01_scDblFinder_doublet_summary_by_sample.csv"),
          row.names = FALSE)

## 9. 保存逐样本检查点 ---------------------------------------------------------
## 第一个RDS保留doublet标签，便于审计；第二个RDS只包含singlet
saveRDS(scRNAlist_qc,
        file.path(out_dir, "objects/PRJNA932038_samplewise_after_QC_before_doublet_filter.rds"))
saveRDS(scRNAlist_singlet,
        file.path(out_dir, "objects/PRJNA932038_samplewise_after_QC_singlets_only.rds"))

## 10. 只合并singlet对象；后续禁止重新使用scRNAlist_qc覆盖它
scRNA1 <- merge(
  x = scRNAlist_singlet[[1]],
  y = scRNAlist_singlet[2:length(scRNAlist_singlet)],
  add.cell.ids = names(scRNAlist_singlet),
  project = "PRJNA932038_singlets"
)
## 合并后的强制防错检查：不得残留doublet，barcode不得重复
if (!all(scRNA1$scDblFinder.class == "singlet")) {
  stop("Merged checkpoint unexpectedly contains non-singlet cells.")
}
if (anyDuplicated(colnames(scRNA1))) {
  stop("Duplicated cell barcodes after merge.")
}

write.csv(as.data.frame(table(scRNA1$sample_id)),
          file.path(out_dir, "tables/01_cell_counts_after_QC_singlet_by_sample.csv"),
          row.names = FALSE)
write.csv(as.data.frame(table(scRNA1$group)),
          file.path(out_dir, "tables/01_cell_counts_after_QC_singlet_by_group.csv"),
          row.names = FALSE)
saveRDS(scRNA1, singlet_checkpoint)

## 记录R和包版本，方便复现
writeLines(capture.output(sessionInfo()), file.path(out_dir, "logs/sessionInfo.txt"))
writeLines(
  c(
    "WAITING_FOR_HUMAN_APPROVAL",
    "Review tables/01_scDblFinder_doublet_summary_by_sample.csv",
    "Review all figures/doublet/ plots",
    "Do not enable RUN_DOWNSTREAM until the doublet rates are approved."
  ),
  file.path(out_dir, "REVIEW_SINGLET_CHECKPOINT.md")
)

## 11. 人工审批闸门 -------------------------------------------------------------
## 默认在这里正常退出，不进行归一化、PCA、UMAP或聚类
if (!RUN_DOWNSTREAM) {
  message("WAITING_FOR_HUMAN_APPROVAL")
  message("Singlet checkpoint: ", singlet_checkpoint)
  quit(save = "no", status = 0)
}

## 12. 可选的探索性下游分析 -----------------------------------------------------
## 只有批准doublet结果并手动设置RUN_DOWNSTREAM <- TRUE后才会运行。
## 本区块故意不自动执行批次校正；未经校正的PCA/UMAP必须与patient/sample
## 校正结果并列比较后，才能批准正式注释对象。
if (inherits(scRNA1[["RNA"]], "Assay5")) {
  scRNA1 <- JoinLayers(scRNA1)
}
scRNA1 <- NormalizeData(scRNA1, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst", nfeatures = 3000)
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1), npcs = 50)

## 沿用上一版的前40个PC；仍需结合ElbowPlot进行人工复核
pc_use <- 1:40
scRNA1 <- FindNeighbors(scRNA1, reduction = "pca", dims = pc_use)
scRNA1 <- FindClusters(scRNA1, resolution = seq(0.2, 1.0, by = 0.1))

if (!requireNamespace("clustree", quietly = TRUE)) {
  stop("Package 'clustree' is required for resolution review.")
}
p_tree <- clustree::clustree(scRNA1, prefix = "RNA_snn_res.")
save_plot_both(
  p_tree,
  file.path(out_dir, "figures/01_singlet_clustree_resolution_0.2_to_1.0"),
  12, 9
)

## 生成resolution网格和clustree后再次停止，防止未经人工核对直接选择0.3
if (is.na(FINAL_RESOLUTION)) {
  saveRDS(scRNA1,
          file.path(out_dir, "objects/PRJNA932038_singlets_resolution_grid_no_final_choice.rds"))
  stop(
    "Resolution grid completed. FINAL_RESOLUTION is intentionally unset. ",
    "Review clustree, markers and per-patient support before choosing it."
  )
}

resolution_column <- paste0("RNA_snn_res.", FINAL_RESOLUTION)
if (!resolution_column %in% colnames(scRNA1[[]])) {
  stop("Requested resolution column does not exist: ", resolution_column)
}
## 只有人工填写了网格中真实存在的FINAL_RESOLUTION，才设置最终cluster身份
Idents(scRNA1) <- resolution_column
scRNA1$seurat_clusters <- scRNA1[[resolution_column, drop = TRUE]]
scRNA1 <- RunUMAP(scRNA1, reduction = "pca", dims = pc_use, seed.use = SEED)

p_cluster <- DimPlot(scRNA1, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
p_sample <- DimPlot(scRNA1, group.by = "sample_id")
p_group <- DimPlot(scRNA1, group.by = "group")
save_plot_both(p_cluster, file.path(out_dir, "figures/UMAP/01_UMAP_clusters"), 8, 7)
save_plot_both(p_sample, file.path(out_dir, "figures/UMAP/01_UMAP_sample"), 8, 7)
save_plot_both(p_group, file.path(out_dir, "figures/UMAP/01_UMAP_group"), 8, 7)

saveRDS(
  scRNA1,
  file.path(out_dir, "objects/PRJNA932038_singlets_uncorrected_exploratory_dimreduced.rds")
)

# Step 04：主要细胞类型注释候选材料生成
#
# 设计原则：
# 1) 这一步只生成“人工注释材料”，不自动写死最终 celltype_major。
# 2) 注释以 canonical marker 为主，SingleR/HPCA 只作为辅助参考。
# 3) 乳腺癌脂肪组织中，Human Adipose Atlas / adipose marker 只能帮助判断
#    “像哪类脂肪组织相关细胞”，不能单独证明 CAF 亚型或功能状态。
# 4) 所有比例与样本组成按 sample/patient 汇总，避免把细胞当作生物学重复。

cfg <- load_config()
root <- project_root()
attempt <- Sys.getenv("ASC_P_ATTEMPT_DIR")

assert_packages(c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble",
  "ggplot2", "patchwork", "SingleR", "celldex", "SummarizedExperiment", "future", "future.apply"
))
source(file.path(root, "R", "analysis_helpers.R"))
set.seed(cfg$seed)

# Enable parallelization for FindAllMarkers and other Seurat ops
future::plan("multicore", workers = 4)
options(future.globals.maxSize = 4000 * 1024^2)

# celldex/ExperimentHub 在 Windows 上可能默认写入 AppData cache；
# CodeBuddy 沙盒/权限环境下容易失败。这里统一改到项目内 cache。
bioc_cache <- file.path(root, "R", "bioc_cache")
dir.create(bioc_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(EXPERIMENT_HUB_CACHE = bioc_cache)
Sys.setenv(ANNOTATION_HUB_CACHE = bioc_cache)
if (requireNamespace("ExperimentHub", quietly = TRUE)) {
  try(ExperimentHub::setExperimentHubOption("CACHE", bioc_cache), silent = TRUE)
}
if (requireNamespace("AnnotationHub", quietly = TRUE)) {
  try(AnnotationHub::setAnnotationHubOption("CACHE", bioc_cache), silent = TRUE)
}

# -----------------------------
# 0. 兼容性与输入预检
# -----------------------------

package_check <- data.frame(
  package = c("Seurat", "SeuratObject", "SingleR", "celldex",
              "SummarizedExperiment", "Matrix", "dplyr", "tidyr",
              "ggplot2", "patchwork"),
  version = vapply(c("Seurat", "SeuratObject", "SingleR", "celldex",
                     "SummarizedExperiment", "Matrix", "dplyr", "tidyr",
                     "ggplot2", "patchwork"),
                   function(p) as.character(utils::packageVersion(p)), character(1)),
  library_path = vapply(c("Seurat", "SeuratObject", "SingleR", "celldex",
                          "SummarizedExperiment", "Matrix", "dplyr", "tidyr",
                          "ggplot2", "patchwork"),
                        function(p) dirname(system.file(package = p)), character(1)),
  stringsAsFactors = FALSE
)
write.csv(package_check, file.path(attempt, "tables", "step04_package_versions.csv"), row.names = FALSE)

compatibility_checks <- data.frame(
  check = c(
    "SeuratObject::LayerData available",
    "Seurat::DimPlot available",
    "Seurat::DotPlot available",
    "SingleR::SingleR available",
    "celldex::HumanPrimaryCellAtlasData available"
  ),
  pass = c(
    exists("LayerData", where = asNamespace("SeuratObject"), inherits = FALSE),
    exists("DimPlot", where = asNamespace("Seurat"), inherits = FALSE),
    exists("DotPlot", where = asNamespace("Seurat"), inherits = FALSE),
    exists("SingleR", where = asNamespace("SingleR"), inherits = FALSE),
    exists("HumanPrimaryCellAtlasData", where = asNamespace("celldex"), inherits = FALSE)
  ),
  stringsAsFactors = FALSE
)
write.csv(compatibility_checks, file.path(attempt, "tables", "step04_compatibility_checks.csv"), row.names = FALSE)
if (!all(compatibility_checks$pass)) {
  stop("Step04 compatibility precheck failed; see step04_compatibility_checks.csv")
}

obj <- read_previous_checkpoint(4)
Seurat::DefaultAssay(obj) <- "RNA"

required_meta <- c("sample_id", "patient_id", "group", "seurat_clusters")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta)) {
  stop("Step04 input metadata missing: ", paste(missing_meta, collapse = ", "))
}
if (!"umap" %in% Seurat::Reductions(obj)) stop("Step04 input lacks Harmony UMAP reduction named `umap`")

Seurat::Idents(obj) <- "seurat_clusters"
cluster_levels <- levels(Seurat::Idents(obj))

# Seurat v5 兼容：显式使用 layer='data'，避免旧 slot 写法。
rna_data <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
if (!inherits(rna_data, "Matrix")) {
  rna_data <- methods::as(rna_data, "dgCMatrix")
}

input_audit <- data.frame(
  metric = c("n_cells", "n_genes", "n_samples", "n_patients", "n_clusters",
             "assay", "data_layer_class", "reductions"),
  value = c(
    ncol(obj), nrow(obj),
    length(unique(obj$sample_id)), length(unique(obj$patient_id)),
    length(cluster_levels),
    Seurat::DefaultAssay(obj), class(rna_data)[1],
    paste(Seurat::Reductions(obj), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(input_audit, file.path(attempt, "tables", "step04_input_audit.csv"), row.names = FALSE)

# -----------------------------
# 1. marker 面板
# -----------------------------

# 依据 breast_cancer_adipose_subcluster_annotation_guide.md 完善的大类 marker。
# 注意：这里用于候选判读，不自动替代人工注释。
annotation_markers <- list(
  adipocyte_like = c("ADIPOQ", "PLIN1", "CIDEA", "FABP4", "LPL", "PPARG"),
  stromal_ASPC_FAP = c("PDGFRA", "DCN", "LUM", "COL1A1", "COL3A1", "THY1", "PI16", "FBN1", "PRG4", "PAMR1"),
  pericyte = c("RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3"),
  smooth_muscle_myofibroblast = c("ACTA2", "TAGLN", "MYL9", "MYH11"),
  endothelial_blood = c("PECAM1", "VWF", "KDR", "CLDN5", "CDH5"),
  endothelial_lymphatic = c("PROX1", "LYVE1", "PDPN", "FLT4"),
  macrophage_myeloid = c("LYZ", "C1QA", "C1QB", "CD68", "MRC1", "CSF1R", "S100A8", "S100A9", "FCGR3A", "CD14"),
  T_NK = c("CD3D", "CD3E", "TRAC", "IL7R", "CD4", "CD8A", "NKG7", "GNLY", "KLRD1"),
  B_plasma = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1"),
  epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  mast = c("TPSAB1", "TPSB2", "KIT", "CPA3"),
  proliferating = c("MKI67", "TOP2A", "STMN1", "UBE2C")
)

marker_catalog <- tibble::tibble(
  marker_set = rep(names(annotation_markers), lengths(annotation_markers)),
  gene = unlist(annotation_markers, use.names = FALSE),
  present = unlist(annotation_markers, use.names = FALSE) %in% rownames(obj)
)
write.csv(marker_catalog, file.path(attempt, "tables", "annotation_marker_catalog.csv"), row.names = FALSE)

present_genes_local <- function(genes) intersect(genes, rownames(obj))

# -----------------------------
# 2. cluster marker score（先算 score，后生成 markers 摘要）
# -----------------------------

score_cols <- character()
for (nm in names(annotation_markers)) {
  genes <- present_genes_local(annotation_markers[[nm]])
  col <- paste0("marker_score_", make.names(nm))
  if (length(genes)) {
    mat <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data", features = genes)
    obj[[col]] <- Matrix::colMeans(mat)
  } else {
    obj[[col]] <- NA_real_
  }
  score_cols <- c(score_cols, col)
}

cluster_scores <- obj@meta.data |>
  tibble::rownames_to_column("cell") |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(score_cols), ~median(.x, na.rm = TRUE)),
    cells = dplyr::n(),
    n_samples = dplyr::n_distinct(sample_id),
    n_patients = dplyr::n_distinct(patient_id),
    .groups = "drop"
  )

score_matrix <- as.matrix(cluster_scores[, score_cols, drop = FALSE])
score_matrix[!is.finite(score_matrix)] <- NA_real_
marker_primary_idx <- max.col(score_matrix, ties.method = "first")
score_sorted <- t(apply(score_matrix, 1, function(x) sort(x, decreasing = TRUE, na.last = TRUE)))
cluster_scores$marker_primary <- names(annotation_markers)[marker_primary_idx]
cluster_scores$marker_second <- names(annotation_markers)[match(colnames(score_matrix)[apply(score_matrix, 1, function(x) order(x, decreasing = TRUE, na.last = TRUE)[2])], score_cols)]
cluster_scores$marker_score_top <- score_sorted[, 1]
cluster_scores$marker_score_second <- score_sorted[, 2]
cluster_scores$marker_score_margin <- cluster_scores$marker_score_top - cluster_scores$marker_score_second
write.csv(cluster_scores, file.path(attempt, "tables", "cluster_marker_scores.csv"), row.names = FALSE)

# --- Generate per-cluster marker summary table (pct.1/pct.2 for best marker set) ---
markers <- data.frame(
  p_val = numeric(), avg_log2FC = numeric(), pct.1 = numeric(),
  pct.2 = numeric(), cluster = character(), gene = character(),
  p_val_adj = numeric(), stringsAsFactors = FALSE
)
for (cl in cluster_levels) {
  cl_idx <- which(cluster_levels == cl)
  top_set <- names(annotation_markers)[which.max(score_matrix[cl_idx, ])]
  top_genes <- present_genes_local(annotation_markers[[top_set]])
  if (length(top_genes) > 0) {
    cells_in <- which(obj$seurat_clusters == cl)
    cells_out <- which(obj$seurat_clusters != cl)
    pct_in <- Matrix::rowMeans(rna_data[top_genes, cells_in, drop = FALSE] > 0)
    pct_out <- Matrix::rowMeans(rna_data[top_genes, cells_out, drop = FALSE] > 0)
    markers <- rbind(markers, data.frame(
      p_val = NA, avg_log2FC = NA, pct.1 = pct_in, pct.2 = pct_out,
      cluster = cl, gene = top_genes, p_val_adj = NA,
      stringsAsFactors = FALSE
    ))
  }
}
write.table(markers, file.path(attempt, "tables", "cluster_markers.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# 3. 样本/分组组成：用于识别样本主导或 N/T 偏倚
# -----------------------------

cluster_sample_counts <- as.data.frame(table(cluster = obj$seurat_clusters, sample_id = obj$sample_id),
                                       stringsAsFactors = FALSE)
cluster_group_counts <- as.data.frame(table(cluster = obj$seurat_clusters, group = obj$group),
                                      stringsAsFactors = FALSE)
write.csv(cluster_sample_counts, file.path(attempt, "tables", "cluster_sample_composition.csv"), row.names = FALSE)
write.csv(cluster_group_counts, file.path(attempt, "tables", "cluster_group_composition.csv"), row.names = FALSE)

cluster_support <- obj@meta.data |>
  tibble::rownames_to_column("cell") |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::summarise(
    cells = dplyr::n(),
    n_samples = dplyr::n_distinct(sample_id),
    n_patients = dplyr::n_distinct(patient_id),
    dominant_sample = names(sort(table(sample_id), decreasing = TRUE))[1],
    dominant_sample_fraction = max(table(sample_id)) / dplyr::n(),
    N_fraction = mean(group == "N"),
    T_fraction = mean(group == "T"),
    .groups = "drop"
  )
write.csv(cluster_support, file.path(attempt, "tables", "cluster_support_summary.csv"), row.names = FALSE)

# -----------------------------
# 4. SingleR 辅助注释
# -----------------------------

singleR_cluster <- data.frame(
  cluster = as.character(cluster_levels),
  SingleR_label = NA_character_,
  SingleR_pruned_label = NA_character_,
  SingleR_delta_next = NA_real_,
  stringsAsFactors = FALSE
)
singleR_note <- paste0(
  "SingleR not completed. ExperimentHub/AnnotationHub cache set to project-local path: ",
  bioc_cache
)

tryCatch({
  ref <- celldex::HumanPrimaryCellAtlasData()
  bp <- if (requireNamespace("BiocParallel", quietly = TRUE)) BiocParallel::SerialParam() else NULL
  if (is.null(bp)) {
    pred <- SingleR::SingleR(
      test = rna_data,
      ref = ref,
      labels = ref$label.main,
      clusters = obj$seurat_clusters
    )
  } else {
    pred <- SingleR::SingleR(
      test = rna_data,
      ref = ref,
      labels = ref$label.main,
      clusters = obj$seurat_clusters,
      BPPARAM = bp
    )
  }
  singleR_cluster <- data.frame(
    cluster = rownames(pred),
    SingleR_label = pred$labels,
    SingleR_pruned_label = pred$pruned.labels,
    SingleR_delta_next = pred$delta.next,
    stringsAsFactors = FALSE
  )
  singleR_note <- "SingleR HPCA completed in cluster mode; use as auxiliary evidence only."
}, error = function(e) {
  singleR_note <<- paste("SingleR unavailable or failed:", conditionMessage(e))
})
writeLines(singleR_note, file.path(attempt, "reports", "SingleR_status.txt"), useBytes = TRUE)
write.csv(singleR_cluster, file.path(attempt, "tables", "SingleR_cluster_predictions.csv"), row.names = FALSE)

# -----------------------------
# 5. 人工注释审核表与低置信提示
# -----------------------------

review <- cluster_scores |>
  dplyr::select(seurat_clusters, cells, n_samples, n_patients,
                marker_primary, marker_second, marker_score_top,
                marker_score_second, marker_score_margin,
                dplyr::all_of(score_cols)) |>
  dplyr::rename(cluster = seurat_clusters) |>
  dplyr::mutate(cluster = as.character(cluster)) |>
  dplyr::left_join(cluster_support |> dplyr::mutate(cluster = as.character(seurat_clusters)) |>
                     dplyr::select(cluster, dominant_sample, dominant_sample_fraction, N_fraction, T_fraction),
                   by = "cluster") |>
  dplyr::left_join(singleR_cluster, by = "cluster") |>
  dplyr::mutate(
    confidence_flags = paste(
      ifelse(marker_score_margin < 0.02, "low_marker_margin", ""),
      ifelse(dominant_sample_fraction > 0.60, "sample_dominated", ""),
      ifelse(n_samples < length(unique(obj$sample_id)), "not_all_samples", ""),
      ifelse(is.na(SingleR_label), "no_SingleR", ""),
      sep = ";"
    ),
    confidence_flags = gsub("(^;+|;+$)", "", gsub(";+", ";", confidence_flags)),
    confidence_flags = ifelse(nchar(confidence_flags) == 0, "none", confidence_flags),
    # 用户人工填写；approve_step.R 会要求 celltype_major 非空后才能批准 Step04。
    celltype_major = "",
    adipose_context_label = "",
    reviewer_notes = ""
  ) |>
  dplyr::arrange(as.integer(cluster))

write.csv(review, file.path(attempt, "tables", "cluster_annotation_review.csv"), row.names = FALSE)

low_confidence <- review |>
  dplyr::filter(confidence_flags != "none" | is.na(SingleR_label))
write.csv(low_confidence, file.path(attempt, "tables", "low_confidence_or_mixed_clusters.csv"), row.names = FALSE)

# -----------------------------
# 5b. 三重交叉验证共识注释（marker + SingleR + PTPRC/CD45）
#     替代 Azimuth/Human Adipose Atlas 参考映射，基于鉴定指南的 marker 体系
# -----------------------------

# --- 计算 PTPRC (CD45) 免疫状态 ---
if ("PTPRC" %in% rownames(obj)) {
  ptprc_dat <- Seurat::FetchData(obj, vars = c("seurat_clusters", "PTPRC"))
  ptprc_cluster <- ptprc_dat |>
    dplyr::group_by(seurat_clusters) |>
    dplyr::summarise(
      PTPRC_pct = mean(PTPRC > 0) * 100,
      PTPRC_mean = mean(PTPRC),
      .groups = "drop"
    )
} else {
  ptprc_cluster <- data.frame(seurat_clusters = cluster_levels, PTPRC_pct = NA_real_, PTPRC_mean = NA_real_)
}

# --- 共识决策引擎 ---
# 规则优先级：PTPRC > marker_score_margin > SingleR > marker_primary
consensus <- review |>
  dplyr::left_join(ptprc_cluster |> dplyr::mutate(seurat_clusters = as.character(seurat_clusters)),
                   by = c("cluster" = "seurat_clusters")) |>
  dplyr::mutate(
    # Step 1: PTPRC immune gate — CD45+/PTPRC+ > 20% pct → immune lineage
    immune_gate = !is.na(PTPRC_pct) & PTPRC_pct > 20,

    # Step 2: Consensus celltype based on marker + SingleR + PTPRC
    consensus_celltype = dplyr::case_when(
      # --- Immune lineages (PTPRC gate) ---
      immune_gate & marker_primary == "mast"                    ~ "mast",
      immune_gate & grepl("T_cell|T cell|NK", SingleR_label)    ~ "T_NK",
      immune_gate & grepl("B_cell|B cell", SingleR_label)       ~ "B_plasma",
      immune_gate & grepl("Macrophage|Monocyte", SingleR_label) ~ "macrophage_myeloid",
      immune_gate & grepl("proliferat", tolower(marker_primary))~ "proliferating_immune",
      immune_gate                                                ~ "macrophage_myeloid",

      # --- Non-immune: use marker primary if margin is clear ---
      !immune_gate & marker_score_margin >= 0.15 &
        marker_primary %in% c("stromal_ASPC_FAP", "endothelial_blood",
                               "endothelial_lymphatic", "smooth_muscle_myofibroblast",
                               "pericyte", "adipocyte_like", "epithelial_tumor",
                               "B_plasma", "T_NK")
        ~ marker_primary,

      # --- Adipocyte with low margin → check SingleR ---
      !immune_gate & marker_primary == "adipocyte_like" &
        grepl("Fibroblast|Tissue_stem", SingleR_label) ~ "stromal_ASPC_FAP",

      # --- Fallback ---
      TRUE ~ marker_primary
    ),

    # Map consensus to human-readable broad labels (compatible with Step 05 filters)
    celltype_major_auto = dplyr::recode(consensus_celltype,
      "stromal_ASPC_FAP"             = "stromal/FAP",
      "adipocyte_like"               = "adipocyte_like",
      "pericyte"                     = "mural/smooth muscle",
      "smooth_muscle_myofibroblast"  = "mural/smooth muscle",
      "endothelial_blood"            = "endothelial",
      "endothelial_lymphatic"        = "endothelial",
      "macrophage_myeloid"           = "macrophage/myeloid",
      "T_NK"                         = "T/NK cells",
      "B_plasma"                     = "B/plasma cells",
      "epithelial_tumor"             = "epithelial",
      "mast"                         = "macrophage/myeloid",
      "proliferating"                = "proliferating",
      "proliferating_immune"         = "proliferating",
      .default = "unknown"
    ),

    # Fine annotations for stromal/adipose clusters (Layer 2 of annotation guide)
    stromal_subtype = dplyr::case_when(
      celltype_major_auto == "stromal/FAP" &
        marker_score_stromal_ASPC_FAP >= 0.3 &
        marker_score_pericyte < 0.2 &
        marker_score_smooth_muscle_myofibroblast < 0.2
        ~ "ASPC/FAP-like (PDGFRA+ COL1A1+ DCN+)",

      celltype_major_auto == "stromal/FAP" &
        marker_score_stromal_ASPC_FAP >= 0.3 &
        (marker_score_pericyte >= 0.2 | marker_score_smooth_muscle_myofibroblast >= 0.2)
        ~ "mixed stromal/FAP-mural",

      celltype_major_auto == "mural/smooth muscle" &
        marker_score_smooth_muscle_myofibroblast > marker_score_pericyte
        ~ "smooth muscle / myofibroblast (ACTA2+ TAGLN+)",

      celltype_major_auto == "mural/smooth muscle" &
        marker_score_pericyte >= marker_score_smooth_muscle_myofibroblast
        ~ "pericyte-like (RGS5+ PDGFRB+ CSPG4+)",

      celltype_major_auto == "adipocyte_like" &
        marker_score_adipocyte_like >= 0.5
        ~ "adipocyte-like (ADIPOQ+ PLIN1+ FABP4+)",

      celltype_major_auto == "endothelial" &
        marker_score_endothelial_blood > marker_score_endothelial_lymphatic
        ~ "blood endothelial (PECAM1+ VWF+ CDH5+)",

      celltype_major_auto == "endothelial" &
        marker_score_endothelial_lymphatic >= marker_score_endothelial_blood
        ~ "lymphatic endothelial (PROX1+ LYVE1+ PDPN+)",

      TRUE ~ ""
    )
  )

# Write automated consensus table (for comparison with human review)
consensus_out <- consensus |>
  dplyr::select(cluster, cells, marker_primary, marker_score_margin,
                PTPRC_pct, SingleR_label, immune_gate,
                consensus_celltype, celltype_major_auto, stromal_subtype,
                confidence_flags)
write.csv(consensus_out, file.path(attempt, "tables", "consensus_annotation_auto.csv"), row.names = FALSE)

# Write pre-filled review table (human can override)
review_filled <- review |>
  dplyr::left_join(consensus |> dplyr::select(cluster, celltype_major_auto, stromal_subtype),
                   by = "cluster") |>
  dplyr::mutate(
    celltype_major = celltype_major_auto,
    adipose_context_label = stromal_subtype,
    reviewer_notes = "AUTO: consensus (marker+SingleR+PTPRC). Human review required before approval."
  ) |>
  dplyr::select(-celltype_major_auto, -stromal_subtype)
write.csv(review_filled, file.path(attempt, "tables", "cluster_annotation_review_autofilled.csv"), row.names = FALSE)

# ----- Consensus summary report -----
consensus_summary <- consensus |>
  dplyr::count(celltype_major_auto, name = "n_clusters") |>
  dplyr::arrange(dplyr::desc(n_clusters))
write.csv(consensus_summary, file.path(attempt, "tables", "consensus_celltype_summary.csv"), row.names = FALSE)

# -----------------------------
# 6. 可视化输出
# -----------------------------

dot_genes <- unique(unlist(annotation_markers, use.names = FALSE))
dot_genes <- present_genes_local(dot_genes)
if (length(dot_genes)) {
  p_dot <- Seurat::DotPlot(obj, features = dot_genes, group.by = "seurat_clusters") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
  save_plot_both(p_dot, file.path(attempt, "figures", "global_marker_dotplot"), 15, 18, dpi = 300)
}

p_umap_cluster <- Seurat::DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggplot2::ggtitle("Step04 input clusters")
save_plot_both(p_umap_cluster, file.path(attempt, "figures", "global_UMAP_clusters"), 8, 7, dpi = 300)

p_umap_sample <- Seurat::DimPlot(obj, reduction = "umap", group.by = "sample_id") +
  ggplot2::ggtitle("Step04 UMAP by sample")
p_umap_group <- Seurat::DimPlot(obj, reduction = "umap", group.by = "group") +
  ggplot2::ggtitle("Step04 UMAP by N/T")
save_plot_both(p_umap_sample + p_umap_group, file.path(attempt, "figures", "global_UMAP_sample_group"), 13, 6, dpi = 300)

feature_genes <- present_genes_local(c(
  "ADIPOQ", "PLIN1", "PDGFRA", "DCN", "COL1A1", "RGS5", "ACTA2",
  "PECAM1", "PROX1", "LYZ", "CD3D", "MS4A1", "EPCAM"
))
if (length(feature_genes)) {
  p_feat <- Seurat::FeaturePlot(
    obj,
    features = feature_genes,
    reduction = "umap",
    ncol = 4,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95"
  )
  save_plot_both(p_feat, file.path(attempt, "figures", "canonical_marker_featureplots"), 14, 12, dpi = 300)
}

# -----------------------------
# 7. 保存检查点与人工审核报告
# -----------------------------

checkpoint <- file.path(attempt, "objects", "PRJNA932038_global_preannotation.rds")
saveRDS(obj, checkpoint)

review_lines <- c(
  "## Purpose",
  "",
  "- This step generated global major-cell-type annotation materials only.",
  "- `celltype_major` was intentionally left blank for human review.",
  "- SingleR/HPCA is auxiliary; canonical markers and sample support remain primary.",
  "- Human/adipose atlas-style labels should be treated as candidate context labels, not definitive CAF or functional-state calls.",
  "",
  "## Required human actions before approval",
  "",
  "1. Review `figures/global_marker_dotplot.pdf`, `figures/canonical_marker_featureplots.pdf`, and UMAPs.",
  "2. Review `tables/cluster_markers.tsv`, `tables/cluster_marker_scores.csv`, and `tables/SingleR_cluster_predictions.csv`.",
  "3. Fill every `celltype_major` in `tables/cluster_annotation_review.csv`.",
  "4. Optional: fill `adipose_context_label` for stromal/adipose/vascular clusters.",
  "5. Pay special attention to clusters flagged in `tables/low_confidence_or_mixed_clusters.csv`, especially sample-dominated clusters such as S81_N-enriched clusters.",
  "",
  "## Suggested broad labels",
  "",
  "- adipocyte_like",
  "- stromal_ASPC_FAP",
  "- pericyte",
  "- smooth_muscle_myofibroblast",
  "- endothelial_blood",
  "- endothelial_lymphatic",
  "- macrophage_myeloid",
  "- T_NK",
  "- B_plasma",
  "- epithelial_tumor",
  "- mast",
  "- proliferating",
  "- doublet_ambiguous",
  "- unknown",
  "",
  "## Status"
)
write_waiting_review(attempt, 4, review_lines)

write_step_manifest(
  4, attempt, checkpoint, file.path(root, "scripts", "step_04.R"),
  extra = list(
    annotation_strategy = "marker_primary_SingleR_auxiliary_human_review_required",
    singleR_status = singleR_note,
    marker_sets = names(annotation_markers)
  )
)

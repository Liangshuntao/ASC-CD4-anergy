#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else getwd()
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE) else NA_character_

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(scales)
  library(digest)
  library(yaml)
})

project_path <- function(...) file.path(root, ...)
out_dir <- project_path("review", "current_results_review_20260715", "reviewer_figure_ASC_P")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step04_path <- project_path("steps", "step04", "attempt_011", "objects", "PRJNA932038_global_annotation_v2_candidate.rds")
step04_map_path <- project_path("steps", "step04", "attempt_011", "tables", "annotation_v2_cluster_mapping.csv")
prefix <- "Reviewer_Figure_HumanBreastAdipose_ASC_P_like_with_panelE_REVISED_FINAL"
output_stems <- file.path(out_dir, c(
  "Figure_A_human_breast_adipose_fixed_annotation_with_panelE_REVISED_FINAL",
  "Figure_B_human_breast_adipose_ASC_P_3of6_candidate_map_REVISED_FINAL",
  "Figure_C_human_breast_adipose_lineage_DotPlot_with_panelE_REVISED_FINAL",
  "Figure_D_human_breast_adipose_ASC_P_3of6_four_group_DotPlot_REVISED_FINAL",
  "Figure_E_human_breast_adipose_ASC_P_ENG_THY1_within_3of6_REVISED_FINAL",
  prefix
))
summary_path <- file.path(out_dir, "HumanBreastAdipose_ASC_P_3of6_candidate_summary_REVISED.tsv")
core_summary_path <- file.path(out_dir, "HumanBreastAdipose_ASC_P_ENG_THY1_within_3of6_summary_REVISED.tsv")
candidate_table_path <- file.path(out_dir, "HumanBreastAdipose_ASC_P_candidate_nuclei_REVISED.tsv")
sample_tissue_path <- file.path(out_dir, "HumanBreastAdipose_ASC_P_candidate_by_sample_and_tissue_REVISED.tsv")
denominator_audit_path <- file.path(out_dir, "HumanBreastAdipose_ASC_P_denominator_audit_REVISED.tsv")
marker_audit_path <- file.path(out_dir, "HumanBreastAdipose_ASC_P_marker_gene_availability_REVISED.tsv")
plot_data_paths <- file.path(out_dir, c(
  "Figure_A_annotation_embedding_data_REVISED.tsv",
  "Figure_B_3of6_candidate_embedding_data_REVISED.tsv",
  "Figure_C_lineage_marker_summary_REVISED.tsv",
  "Figure_D_3of6_six_gene_boundary_summary_REVISED.tsv",
  "Figure_E_ENG_THY1_within_3of6_embedding_data_REVISED.tsv"
))
legend_path <- file.path(out_dir, paste0(prefix, "_legend.md"))
report_path <- file.path(out_dir, "PRJNA932038_ASC_P_like_with_panelE_REVISED_RUN_REPORT.md")
session_path <- file.path(out_dir, "PRJNA932038_ASC_P_like_with_panelE_REVISED_sessionInfo.txt")
manifest_path <- file.path(out_dir, "reviewer_figure_with_panelE_REVISED_manifest.yml")

all_output_paths <- c(
  paste0(output_stems, ".pdf"), paste0(output_stems, ".png"), paste0(output_stems, ".tiff"),
  summary_path, core_summary_path, candidate_table_path, sample_tissue_path,
  denominator_audit_path, marker_audit_path, plot_data_paths, legend_path,
  report_path, session_path, manifest_path
)

write_tsv_utf8 <- function(x, path) {
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "", fileEncoding = "UTF-8")
}

save_plot <- function(plot, stem, width, height) {
  paths <- paste0(stem, c(".pdf", ".png", ".tiff"))
  ggplot2::ggsave(paths[[1]], plot, width = width, height = height,
                  device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf",
                  bg = "white", limitsize = FALSE)
  ggplot2::ggsave(paths[[2]], plot, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(paths[[3]], plot, width = width, height = height, dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE)
  if (!all(file.exists(paths)) || any(file.info(paths)$size <= 0)) stop("One or more plot outputs are missing or empty")
  paths
}

resolve_feature <- function(requested, features) {
  exact <- which(features == requested)
  if (length(exact) == 1L) return(list(exact = TRUE, case_insensitive = TRUE, resolved = features[exact]))
  ci <- which(toupper(features) == toupper(requested))
  if (length(ci) == 1L) return(list(exact = FALSE, case_insensitive = TRUE, resolved = features[ci]))
  list(exact = FALSE, case_insensitive = FALSE, resolved = NA_character_)
}

summary_rows <- function(groups, genes, group_levels, data_layer, counts_layer) {
  out <- list()
  for (group_name in group_levels) {
    group_cells <- names(groups)[!is.na(groups) & groups == group_name]
    if (!length(group_cells)) {
      if (identical(group_name, "Fibroblast/CAF ASC_P-like candidates (>=3/6)")) {
        for (gene in genes) out[[length(out) + 1L]] <- data.frame(
          group = group_name, gene = gene, nucleus_number = 0L,
          percent_detected_raw_counts = 0, average_normalized_RNA = NA_real_, stringsAsFactors = FALSE
        )
      }
      next
    }
    for (gene in genes) {
      if (!gene %in% rownames(data_layer) || !gene %in% rownames(counts_layer)) next
      raw <- as.numeric(counts_layer[gene, group_cells, drop = TRUE])
      norm <- as.numeric(data_layer[gene, group_cells, drop = TRUE])
      out[[length(out) + 1L]] <- data.frame(
        group = group_name, gene = gene, nucleus_number = length(group_cells),
        percent_detected_raw_counts = 100 * mean(raw > 0),
        average_normalized_RNA = mean(norm), stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) stop("No group-by-gene summaries generated")
  do.call(rbind, out)
}

main <- function() {
  if (any(file.exists(all_output_paths)) && Sys.getenv("REVIEW_FIGURE_OVERWRITE") != "1") {
    stop("Output exists; set REVIEW_FIGURE_OVERWRITE=1 to redraw: ", paste(all_output_paths[file.exists(all_output_paths)], collapse = ", "))
  }
  if (!file.exists(step04_path)) stop("Missing required Step04 object: ", step04_path)
  if (!file.exists(step04_map_path)) stop("Missing required formal annotation map: ", step04_map_path)

  obj <- readRDS(step04_path)
  meta <- obj[[]]
  required_meta <- c("seurat_clusters", "celltype_v2", "major_celltype_v2", "sample_id", "patient_id", "group")
  if (!all(required_meta %in% colnames(meta))) stop("Missing required audited metadata fields: ", paste(setdiff(required_meta, colnames(meta)), collapse = ", "))

  sample_col <- "sample_id"
  patient_col <- "patient_id"
  tissue_col <- "group"
  tissue_label_map <- c(N = "distant-normal breast adipose tissue", T = "tumor-adjacent breast adipose tissue")
  for (field in c(sample_col, patient_col, tissue_col)) {
    values <- as.character(meta[[field]])
    if (anyNA(values) || any(!nzchar(values))) stop("Audited metadata field contains missing/empty values: ", field)
  }
  tissue_values <- as.character(meta[[tissue_col]])
  names(tissue_values) <- rownames(meta)
  if (!setequal(unique(tissue_values), names(tissue_label_map))) stop("Unexpected group values; expected exactly N and T, observed: ", paste(sort(unique(tissue_values)), collapse = ", "))
  sample_values <- as.character(meta[[sample_col]]); names(sample_values) <- rownames(meta)
  patient_values <- as.character(meta[[patient_col]]); names(patient_values) <- rownames(meta)
  sample_map <- unique(data.frame(sample = sample_values, patient = patient_values, tissue = tissue_values, stringsAsFactors = FALSE))
  if (anyDuplicated(sample_map$sample)) stop("A sample maps to multiple patient or tissue labels")

  if (!"RNA" %in% names(obj@assays)) stop("Step04 object lacks RNA assay")
  rna_layers <- SeuratObject::Layers(obj[["RNA"]])
  if (!all(c("counts", "data") %in% rna_layers)) stop("RNA counts and data layers are both required")

  mapping <- utils::read.csv(step04_map_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_mapping <- c("cluster", "celltype_v2", "major_celltype_v2", "n_cells")
  if (!all(required_mapping %in% colnames(mapping))) stop("Formal map lacks required columns: ", paste(setdiff(required_mapping, colnames(mapping)), collapse = ", "))
  mapping$cluster <- as.character(mapping$cluster)
  if (anyDuplicated(mapping$cluster)) stop("Formal map contains duplicated clusters")
  meta$cluster <- as.character(meta$seurat_clusters)
  if (!all(unique(meta$cluster) %in% mapping$cluster)) stop("Formal map does not cover all object clusters")
  meta$celltype_v2 <- unname(setNames(mapping$celltype_v2, mapping$cluster)[meta$cluster])
  meta$major_celltype_v2 <- unname(setNames(mapping$major_celltype_v2, mapping$cluster)[meta$cluster])
  actual_cluster_n <- as.integer(table(factor(meta$cluster, levels = mapping$cluster)))
  names(actual_cluster_n) <- mapping$cluster
  if (any(as.integer(mapping$n_cells) != actual_cluster_n[mapping$cluster])) stop("Formal map n_cells does not match actual object cluster counts")

  fibroblast_caf <- rownames(meta)[meta$major_celltype_v2 == "Fibroblast/CAF"]
  if (!length(fibroblast_caf)) stop("No formal Fibroblast/CAF denominator")

  asc_genes <- c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM")
  lineage_sets <- list(
    CAF = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "DPT", "SFRP2", "PDGFRA", "COL10A1", "COMP"),
    Endothelial = c("PECAM1", "VWF", "CLDN5", "EMCN", "KDR", "PLVAP", "ESAM", "RAMP2"),
    Epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT14", "KRT5", "KRT17", "MUC1"),
    Immune = c("PTPRC", "CD3D", "CD3E", "TRAC", "CD79A", "MS4A1", "NKG7", "GNLY", "LST1", "TYROBP"),
    Mural = c("PDGFRB", "RGS5", "CSPG4", "ACTA2", "TAGLN", "MYH11", "NOTCH3", "MCAM", "CAVIN3"),
    Adipocyte = c("ADIPOQ", "PLIN1", "CIDEC", "FABP4", "LEP")
  )
  requested_markers <- unique(c(unlist(lineage_sets, use.names = FALSE), asc_genes))
  feature_names <- rownames(obj[["RNA"]])
  marker_resolution <- lapply(requested_markers, resolve_feature, features = feature_names)
  marker_audit <- data.frame(
    requested_gene = requested_markers,
    present_exact = vapply(marker_resolution, `[[`, logical(1), "exact"),
    present_case_insensitive = vapply(marker_resolution, `[[`, logical(1), "case_insensitive"),
    resolved_feature = vapply(marker_resolution, `[[`, character(1), "resolved"),
    used_in_panel = requested_markers %in% unlist(lineage_sets, use.names = FALSE) | requested_markers %in% asc_genes,
    stringsAsFactors = FALSE
  )
  if (anyNA(marker_audit$resolved_feature[marker_audit$requested_gene %in% asc_genes])) stop("Six-gene panel is incomplete")
  resolved_asc <- setNames(marker_audit$resolved_feature[match(asc_genes, marker_audit$requested_gene)], asc_genes)
  plot_genes <- unique(na.omit(marker_audit$resolved_feature))
  data_layer <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data", features = plot_genes)
  counts_layer <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts", features = plot_genes)
  for (gene in asc_genes) {
    actual <- unname(resolved_asc[[gene]])
    if (actual != gene) {
      rownames(data_layer)[rownames(data_layer) == actual] <- gene
      rownames(counts_layer)[rownames(counts_layer) == actual] <- gene
    }
  }
  asc_counts <- counts_layer[asc_genes, fibroblast_caf, drop = FALSE]
  asc_detected <- asc_counts > 0
  n_detected_6 <- Matrix::colSums(asc_detected)
  candidate_3of6 <- names(n_detected_6)[n_detected_6 >= 3L]
  core_cells <- candidate_3of6[asc_detected["ENG", candidate_3of6] & asc_detected["THY1", candidate_3of6]
  ]
  n_fibroblast_caf <- length(fibroblast_caf)
  n_candidate_3of6 <- length(candidate_3of6)
  n_core <- length(core_cells)
  pct_candidate <- 100 * n_candidate_3of6 / n_fibroblast_caf
  pct_core <- if (n_candidate_3of6) 100 * n_core / n_candidate_3of6 else NA_real_

  emb_reduction <- if ("umap" %in% names(obj@reductions)) "umap" else if ("umap_pca" %in% names(obj@reductions)) "umap_pca" else stop("No approved UMAP reduction")
  emb <- SeuratObject::Embeddings(obj, reduction = emb_reduction)[, 1:2, drop = FALSE]
  if (!setequal(rownames(emb), rownames(meta))) stop("UMAP and metadata cell sets differ")
  if (!all(c(fibroblast_caf, candidate_3of6, core_cells) %in% rownames(emb))) stop("Fibroblast/CAF or candidate cells are missing from UMAP")
  reduction_label <- if (emb_reduction == "umap") "UMAP" else "PCA UMAP"
  plot_meta <- meta[rownames(emb), , drop = FALSE]

  candidate_meta <- meta[candidate_3of6, , drop = FALSE]
  candidate_table <- data.frame(
    nucleus_barcode = candidate_3of6,
    sample = sample_values[candidate_3of6], patient = patient_values[candidate_3of6],
    tissue_condition = unname(tissue_label_map[tissue_values[candidate_3of6]]),
    tissue_raw_label = tissue_values[candidate_3of6], cluster = candidate_meta$cluster,
    formal_celltype = candidate_meta$celltype_v2, formal_major_celltype = candidate_meta$major_celltype_v2,
    stringsAsFactors = FALSE
  )
  for (gene in asc_genes) candidate_table[[paste0(gene, "_raw_detected")]] <- as.integer(asc_detected[gene, candidate_3of6])
  candidate_table$n_detected_6 <- as.integer(n_detected_6[candidate_3of6])
  candidate_table$candidate_3of6 <- rep(TRUE, n_candidate_3of6)
  candidate_table$ENG_THY1_codetected_within_3of6 <- candidate_table$ENG_raw_detected == 1L & candidate_table$THY1_raw_detected == 1L
  write_tsv_utf8(candidate_table, candidate_table_path)

  fib_meta <- meta[fibroblast_caf, , drop = FALSE]
  sample_tissue_keys <- unique(data.frame(sample = sample_values[fibroblast_caf], patient = patient_values[fibroblast_caf], tissue_raw_label = tissue_values[fibroblast_caf], stringsAsFactors = FALSE))
  sample_tissue_rows <- lapply(seq_len(nrow(sample_tissue_keys)), function(i) {
    s <- sample_tissue_keys$sample[i]; t <- sample_tissue_keys$tissue_raw_label[i]
    denom <- fibroblast_caf[sample_values[fibroblast_caf] == s & tissue_values[fibroblast_caf] == t]
    cand <- intersect(candidate_3of6, denom); core <- intersect(core_cells, denom)
    data.frame(sample = s, patient = sample_tissue_keys$patient[i], tissue_raw_label = t, tissue_condition = unname(tissue_label_map[t]),
               n_fibroblast_caf = length(denom), n_candidate_3of6 = length(cand), pct_candidate_3of6 = 100 * length(cand) / length(denom),
               n_eng_thy1_within_3of6 = length(core), pct_eng_thy1_within_3of6 = if (length(cand)) 100 * length(core) / length(cand) else NA_real_, stringsAsFactors = FALSE)
  })
  sample_tissue_summary <- do.call(rbind, sample_tissue_rows)
  write_tsv_utf8(sample_tissue_summary, sample_tissue_path)
  candidate_by_sample <- table(sample_values[candidate_3of6])
  core_by_sample <- table(sample_values[core_cells])
  max_candidate_share <- if (n_candidate_3of6) max(candidate_by_sample) / n_candidate_3of6 else NA_real_
  max_core_share <- if (n_core) max(core_by_sample) / n_core else NA_real_

  summary_table <- data.frame(
    metric = c("Fibroblast/CAF denominator", "Partial ASC_P-like candidate nuclei (>=3/6)", "ENG+THY1 co-detected within >=3/6 candidates"),
    n = c(n_fibroblast_caf, n_candidate_3of6, n_core), denominator = c(NA_integer_, n_fibroblast_caf, n_candidate_3of6),
    percent = c(NA_real_, pct_candidate, pct_core), stringsAsFactors = FALSE
  )
  for (t in names(tissue_label_map)) {
    denom <- fibroblast_caf[tissue_values[fibroblast_caf] == t]
    cand <- intersect(candidate_3of6, denom); core <- intersect(core_cells, denom)
    summary_table <- rbind(summary_table, data.frame(metric = paste0("Tissue source: ", tissue_label_map[t]), n = length(cand), denominator = length(denom), percent = 100 * length(cand) / length(denom), stringsAsFactors = FALSE))
  }
  write_tsv_utf8(summary_table, summary_path)
  write_tsv_utf8(data.frame(metric = c("candidate_3of6", "ENG_THY1_codetected_within_3of6"), n = c(n_candidate_3of6, n_core), denominator = c(n_fibroblast_caf, n_candidate_3of6), percent = c(pct_candidate, pct_core), stringsAsFactors = FALSE), core_summary_path)

  candidate_cluster_n <- table(factor(meta$cluster[match(candidate_3of6, rownames(meta))], levels = mapping$cluster))
  core_cluster_n <- table(factor(meta$cluster[match(core_cells, rownames(meta))], levels = mapping$cluster))
  denominator_audit <- mapping |>
    filter(major_celltype_v2 == "Fibroblast/CAF") |>
    transmute(cluster = as.character(cluster), formal_celltype = as.character(celltype_v2), formal_major_celltype = as.character(major_celltype_v2), mapping_nuclei = as.integer(n_cells), actual_nuclei = as.integer(actual_cluster_n[cluster]), n_candidate_3of6 = as.integer(candidate_cluster_n[cluster]), n_eng_thy1_within_3of6 = as.integer(core_cluster_n[cluster]))
  if (any(denominator_audit$mapping_nuclei != denominator_audit$actual_nuclei)) stop("Cluster denominator audit mismatch")
  write_tsv_utf8(denominator_audit, denominator_audit_path)
  write_tsv_utf8(marker_audit, marker_audit_path)

  cluster_levels <- mapping$cluster[order(suppressWarnings(as.numeric(mapping$cluster)))]
  annotation_labels <- setNames(paste0("Cluster ", mapping$cluster, ": ", mapping$celltype_v2), mapping$cluster)
  annotation_df <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], cluster = plot_meta$cluster, annotation = factor(annotation_labels[plot_meta$cluster], levels = annotation_labels[cluster_levels]), stringsAsFactors = FALSE)
  annotation_centres <- annotation_df |> group_by(cluster) |> summarise(dim1 = median(dim1), dim2 = median(dim2), .groups = "drop")
  cluster_palette <- setNames(grDevices::hcl.colors(length(cluster_levels), "Dynamic"), annotation_labels[cluster_levels])
  pA <- ggplot(annotation_df, aes(dim1, dim2, colour = annotation)) + geom_point(size = 0.16, alpha = 0.68) + geom_label(data = annotation_centres, aes(dim1, dim2, label = cluster), inherit.aes = FALSE, colour = "black", fill = "white", linewidth = 0.2, size = 3.8, fontface = "bold", alpha = 0.82) + coord_equal() + scale_colour_manual(values = cluster_palette, drop = FALSE) + labs(title = "Broad annotation of normal and\ncancer-adjacent breast adipose tissue", x = paste0(reduction_label, "-1"), y = paste0(reduction_label, "-2"), colour = NULL) + theme_classic(base_size = 10) + theme(aspect.ratio = 1, legend.position = "right", legend.text = element_text(size = 5.6), legend.margin = margin(0, 0, 0, 0), plot.title = element_text(size = 9.5, lineheight = 0.95)) + guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1.8)))

  candidate_df <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], candidate = rownames(emb) %in% candidate_3of6, stringsAsFactors = FALSE)
  pB <- ggplot(candidate_df, aes(dim1, dim2)) + geom_point(colour = "grey85", size = 0.14, alpha = 0.55) + geom_point(data = candidate_df[candidate_df$candidate, , drop = FALSE], colour = "#D55E00", size = 0.42, alpha = 0.95) + coord_equal() + labs(title = "Partial ASC_P-like candidates in the\nFibroblast/CAF compartment", subtitle = sprintf(">=3/6 genes detected: n = %d (%.3f%% of Fibroblast/CAF)", n_candidate_3of6, pct_candidate), x = paste0(reduction_label, "-1"), y = paste0(reduction_label, "-2")) + theme_classic(base_size = 10) + theme(plot.title = element_text(size = 9.5, lineheight = 0.95), plot.subtitle = element_text(size = 8.5, lineheight = 0.95))

  lineage_order <- unique(unlist(lineage_sets, use.names = FALSE))
  lineage_groups <- plot_meta$cluster; names(lineage_groups) <- rownames(plot_meta)
  lineage_summary <- summary_rows(lineage_groups, lineage_order, cluster_levels, data_layer, counts_layer)
  lineage_summary$lineage <- unname(setNames(rep(names(lineage_sets), lengths(lineage_sets)), lineage_order)[lineage_summary$gene])
  lineage_summary$gene <- factor(lineage_summary$gene, levels = lineage_order); lineage_summary$group <- factor(lineage_summary$group, levels = rev(cluster_levels))
  pC <- ggplot(lineage_summary, aes(gene, group, size = percent_detected_raw_counts, colour = average_normalized_RNA)) + geom_point() + scale_size(range = c(0.35, 6.5), name = "Nuclei detected (%)") + scale_colour_viridis_c(option = "C", name = "Mean normalized expression") + facet_grid(. ~ lineage, scales = "free_x", space = "free_x") + labs(title = "Lineage-marker expression across annotated nuclei", x = NULL, y = "Cluster") + theme_classic(base_size = 8.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), strip.background = element_rect(fill = "grey94", colour = NA), legend.position = "bottom")

  group <- rep("Other nuclei", nrow(plot_meta)); names(group) <- rownames(plot_meta)
  group[plot_meta$major_celltype_v2 == "Fibroblast/CAF"] <- "Other Fibroblast/CAF"
  group[plot_meta$major_celltype_v2 %in% c("Endothelial", "Endothelial-like / unresolved", "Lymphatic endothelial")] <- "Vascular endothelial"
  group[plot_meta$major_celltype_v2 == "Mural"] <- "Mural/pericyte"
  group[candidate_3of6] <- "Fibroblast/CAF ASC_P-like candidates (>=3/6)"
  group_levels <- c("Mural/pericyte", "Vascular endothelial", "Other Fibroblast/CAF", "Fibroblast/CAF ASC_P-like candidates (>=3/6)")
  required_groups <- setdiff(group_levels, "Fibroblast/CAF ASC_P-like candidates (>=3/6)")
  if (!all(required_groups %in% unique(group))) stop("Panel D missing required non-candidate groups")
  group <- factor(group, levels = group_levels); group_values <- as.character(group); names(group_values) <- rownames(plot_meta)
  four_summary <- summary_rows(group_values, asc_genes, group_levels, data_layer, counts_layer)
  four_summary$group <- factor(four_summary$group, levels = group_levels); four_summary$gene <- factor(four_summary$gene, levels = asc_genes)
  four_summary$mean_scaled_expression <- ave(four_summary$average_normalized_RNA, four_summary$gene, FUN = function(x) { z <- as.numeric(scale(x)); z[!is.finite(z)] <- 0; pmax(-2, pmin(2, z)) })
  pD <- ggplot(four_summary, aes(gene, group, size = percent_detected_raw_counts, colour = mean_scaled_expression)) + geom_point() + scale_size(range = c(0.5, 8.5), name = "Nuclei detected (%)") + scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "Mean scaled expression") + labs(title = "ASC_P-associated gene expression across boundary compartments", x = NULL, y = NULL) + theme_classic(base_size = 9.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

  panel_e_df <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], panel_e_group = "Other nuclei", stringsAsFactors = FALSE)
  panel_e_df$panel_e_group[rownames(emb) %in% candidate_3of6] <- ">=3/6 candidate without ENG+THY1 co-detection"
  panel_e_df$panel_e_group[rownames(emb) %in% core_cells] <- "ENG+THY1 co-detected within >=3/6"
  panel_e_df$panel_e_group <- factor(panel_e_df$panel_e_group, levels = c("Other nuclei", ">=3/6 candidate without ENG+THY1 co-detection", "ENG+THY1 co-detected within >=3/6"))
  pE <- ggplot(panel_e_df, aes(dim1, dim2, colour = panel_e_group)) + geom_point(size = 0.14, alpha = 0.5) + coord_equal() + scale_colour_manual(values = c("Other nuclei" = "grey88", ">=3/6 candidate without ENG+THY1 co-detection" = "#F4A261", "ENG+THY1 co-detected within >=3/6" = "#9E2A2B"), labels = c("Other nuclei", ">=3/6 candidate\nwithout ENG+THY1", "ENG+THY1 core"), drop = FALSE, name = NULL) + labs(title = "ENG+THY1 within >=3/6\ncandidates", subtitle = sprintf("n = %d of %d candidates (%.2f%%)", n_core, n_candidate_3of6, pct_core), x = paste0(reduction_label, "-1"), y = paste0(reduction_label, "-2")) + theme_classic(base_size = 10) + theme(legend.position = "bottom", legend.text = element_text(size = 6), plot.title = element_text(size = 9.5, lineheight = 0.95), plot.subtitle = element_text(size = 8.5, lineheight = 0.95))

  figure_outputs <- c(save_plot(pA, output_stems[[1]], 10, 11), save_plot(pB, output_stems[[2]], 10, 8), save_plot(pC, output_stems[[3]], 16, 7), save_plot(pD, output_stems[[4]], 11, 5), save_plot(pE, output_stems[[5]], 10, 8))
  bottom_row <- pD | pE
  composite <- ((pA | pB) / pC / bottom_row) + plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 14))) + plot_layout(heights = c(1.15, 0.85, 0.80))
  figure_outputs <- c(figure_outputs, save_plot(composite, output_stems[[6]], 18, 16))
  write_tsv_utf8(annotation_df, plot_data_paths[[1]]); write_tsv_utf8(candidate_df, plot_data_paths[[2]]); write_tsv_utf8(lineage_summary, plot_data_paths[[3]]); write_tsv_utf8(four_summary, plot_data_paths[[4]]); write_tsv_utf8(panel_e_df, plot_data_paths[[5]])

  legend <- paste0("# Figure legend\n\n", "**Partial ASC_P-like candidate fibroblast nuclei in normal and cancer-adjacent human breast adipose tissue.**\n\n", sprintf("(A) Approved Step04 annotation. (B) Fibroblast/CAF nuclei with raw RNA counts > 0 for at least three of six genes (ENG, THY1, NT5E, PDGFRA, VCAM1 and MCAM); n = %d (%.3f%% of the approved Fibroblast/CAF denominator, n = %d). (C) Lineage-marker expression across annotated nuclei. (D) Six-gene context across >=3/6 candidates, other Fibroblast/CAF, vascular endothelial and mural/pericyte compartments. (E) ENG and THY1 co-detection within the >=3/6 candidates; n = %d of %d (%.3f%%). Panel E is nested and not independent validation. T denotes tumor-adjacent breast adipose tissue and N denotes distant-normal breast adipose tissue. These are transcriptional candidates, not a definitively established human ASC_P population.", n_candidate_3of6, pct_candidate, n_fibroblast_caf, n_core, n_candidate_3of6, pct_core))
  writeLines(legend, legend_path, useBytes = TRUE)

  report_lines <- c("# PRJNA932038 ASC_P-like reviewer figure with Panel E: revised run report", "", paste0("Status: SUCCESS; generated at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z")), "", "## Audited input", paste0("- Step04 object: ", step04_path), paste0("- Formal map: ", step04_map_path), "- Tissue field: group; explicit mapping N = distant-normal breast adipose tissue, T = tumor-adjacent breast adipose tissue.", "- No reclustering, no new UMAP, no original RDS modification, and no Step07 result used.", "- Layout-only revision: Panel A legend moved to the right, Panel D and Panel E placed on one row, and composite export at 18 x 16 inches.", "", "## Counts", sprintf("- Fibroblast/CAF denominator: %d nuclei.", n_fibroblast_caf), sprintf("- >=3/6 candidates: %d (%.3f%%).", n_candidate_3of6, pct_candidate), sprintf("- ENG+THY1 within >=3/6: %d (%.3f%% of candidates).", n_core, pct_core), sprintf("- Candidate sample levels: %d; candidate patient levels: %d.", length(unique(sample_values[candidate_3of6])), length(unique(patient_values[candidate_3of6]))), sprintf("- Maximum candidate share from one sample: %.3f%%; maximum ENG+THY1 core share from one sample: %.3f%%.", 100 * max_candidate_share, 100 * max_core_share), "- Full sample-by-tissue counts are in the accompanying TSV.", "", "## Evidence boundary", "- Panel D is a context summary and Panel E is a nested core; neither is independent validation.", "- The result supports a human breast adipose transcriptomic correlate, not an identical human ASC_P population, protein-defined CD105+/CD90+ identity, intratumoral validation, or conserved function.")
  writeLines(report_lines, report_path, useBytes = TRUE)
  capture.output(sessionInfo(), file = session_path)
  manifest <- list(created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), status = "SUCCESS", input_step04 = "steps/step04/attempt_011/objects/PRJNA932038_global_annotation_v2_candidate.rds", input_step04_sha256 = digest::digest(file = step04_path, algo = "sha256", serialize = FALSE), input_step04_mapping = "steps/step04/attempt_011/tables/annotation_v2_cluster_mapping.csv", input_step04_mapping_sha256 = digest::digest(file = step04_map_path, algo = "sha256", serialize = FALSE), script_sha256 = if (!is.na(script_path) && file.exists(script_path)) digest::digest(file = script_path, algo = "sha256", serialize = FALSE) else NA_character_, candidate_definition = "Formal Fibroblast/CAF denominator; raw RNA counts > 0 for >=3 of six fixed genes", nested_core_definition = "Within >=3/6 candidates; raw RNA counts > 0 for ENG and THY1", n_fibroblast_caf = n_fibroblast_caf, n_candidate_3of6 = n_candidate_3of6, pct_candidate_3of6 = pct_candidate, n_eng_thy1_within_3of6 = n_core, pct_eng_thy1_within_3of6 = pct_core, sample_field = sample_col, patient_field = patient_col, tissue_field = tissue_col, tissue_label_map = list(N = unname(tissue_label_map[["N"]]), T = unname(tissue_label_map[["T"]])), umap_reduction = emb_reduction, no_step07_result_used = TRUE, outputs = sub(paste0(root, "/"), "", c(figure_outputs, summary_path, core_summary_path, candidate_table_path, sample_tissue_path, denominator_audit_path, marker_audit_path, plot_data_paths, legend_path, report_path, session_path)))
  yaml::write_yaml(manifest, manifest_path)
  message("Revised Panel E reviewer figure written to ", out_dir)
}

tryCatch(
  main(),
  error = function(e) {
    blocked_lines <- c("# PRJNA932038 ASC_P-like reviewer figure: BLOCKED", "", paste0("Status: BLOCKED; ", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z")), "", "## Error", conditionMessage(e), "", "Formal outputs were not accepted as SUCCESS.")
    if (Sys.getenv("REVIEW_FIGURE_OVERWRITE") == "1" || !file.exists(report_path)) writeLines(blocked_lines, report_path, useBytes = TRUE)
    stop(e)
  }
)

# Final reviewer-facing NormalTotal ASC_P-like candidate figure (stage 03b).
# Plotting-only revision: no reclustering, normalization, integration,
# differential-expression testing, trajectory analysis, or new annotation.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so --file is available", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)

locate_project_root <- function(path) {
  candidates <- character()
  if (nzchar(path)) {
    candidates <- c(candidates, path)
    current <- dirname(path)
    while (!identical(current, dirname(current))) {
      candidates <- c(candidates, current)
      current <- dirname(current)
    }
    candidates <- c(candidates, getwd())
  }
  candidates <- unique(candidates)
  for (candidate in candidates) {
    candidate_path <- tryCatch(normalizePath(candidate, winslash = "/", mustWork = FALSE), error = function(e) candidate)
    repeat {
      if (file.exists(file.path(candidate_path, "config", "config_GSE161529_ASC_P.R")) &&
          file.exists(file.path(candidate_path, "scripts", "utils_GSE161529_ASC_P.R"))) {
        return(normalizePath(candidate_path, winslash = "/", mustWork = TRUE))
      }
      parent <- dirname(candidate_path)
      if (identical(parent, candidate_path)) break
      candidate_path <- parent
    }
  }

  workspace_root <- normalizePath(file.path("D:/Users/lstsw/Documents/Adipose_derived stromal_metabolic niches"), winslash = "/", mustWork = FALSE)
  candidate_project <- file.path(workspace_root, "adipose_stromal_metabolic_niche_lit", "public_data_analysis", "GSE161529_ASC_P_validation")
  if (file.exists(file.path(candidate_project, "config", "config_GSE161529_ASC_P.R")) &&
      file.exists(file.path(candidate_project, "scripts", "utils_GSE161529_ASC_P.R"))) {
    return(normalizePath(candidate_project, winslash = "/", mustWork = TRUE))
  }

  stop("Could not locate the GSE161529 ASC_P project root from script path: ", path, call. = FALSE)
}

project_root <- locate_project_root(script_path)
Sys.setenv(GSE161529_ASC_P_PROJECT_ROOT = project_root)
source(file.path(project_root, "scripts", "utils_GSE161529_ASC_P.R"))
source(file.path(project_root, "config", "config_GSE161529_ASC_P.R"), local = .GlobalEnv)
cfg <- get("cfg", envir = .GlobalEnv)
cfg$overwrite <- TRUE
ensure_packages(c("SeuratObject", "Matrix", "ggplot2", "patchwork", "scales"))
set.seed(cfg$seed)

input_annotated <- file.path(cfg$updated_object_dir, "SeuratObject_NormTotal_annotated.rds")
input_original <- file.path(cfg$input_dir, "SeuratObject_NormTotal.rds")
input_path <- if (file.exists(input_annotated)) input_annotated else input_original
if (!file.exists(input_path)) {
  stop("Neither annotated nor original NormalTotal RDS is available", call. = FALSE)
}

annotation_path <- file.path(cfg$table_dir, "NormTotal_final_cluster_annotation.tsv")
if (!file.exists(annotation_path)) {
  stop("Missing fixed NormalTotal annotation table: ", annotation_path, call. = FALSE)
}

reviewer_figure_dir <- file.path(cfg$figure_dir, "reviewer_ASC_P")
dir.create(reviewer_figure_dir, recursive = TRUE, showWarnings = FALSE)

state <- begin_step(
  cfg,
  "03_NormTotal_ASC_P_candidate_reviewer_figures_FINAL",
  c(input_path, annotation_path)
)

make_summary <- function(data_mat, counts_mat, cells, group, gene_set) {
  genes <- intersect(gene_set, intersect(rownames(data_mat), rownames(counts_mat)))
  if (!length(genes)) stop("No requested genes are available in both RNA layers", call. = FALSE)

  rows <- list()
  for (group_i in levels(group)) {
    selected <- cells[!is.na(group) & group == group_i]
    if (!length(selected)) next

    for (gene_i in genes) {
      raw <- as.numeric(counts_mat[gene_i, selected, drop = TRUE])
      norm <- as.numeric(data_mat[gene_i, selected, drop = TRUE])
      rows[[length(rows) + 1L]] <- data.frame(
        group = group_i,
        gene = gene_i,
        cell_number = length(selected),
        percent_expressed_counts = 100 * mean(raw > 0),
        average_normalized_RNA = mean(norm),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) stop("No group-by-gene summaries could be calculated", call. = FALSE)
  do.call(rbind, rows)
}

save_final_figure <- function(plot, stem, width = 18, height = 16.5) {
  pdf_path <- paste0(stem, ".pdf")
  png_path <- paste0(stem, ".png")
  tiff_path <- paste0(stem, ".tiff")
  message("Saving figure outputs to: ", pdf_path, ", ", png_path, ", ", tiff_path)

  assert_new_path(pdf_path, cfg$overwrite)
  assert_new_path(png_path, cfg$overwrite)
  assert_new_path(tiff_path, cfg$overwrite)

  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"

  ggplot2::ggsave(
    filename = pdf_path, plot = plot, width = width, height = height,
    units = "in", device = pdf_device, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = png_path, plot = plot, width = width, height = height,
    units = "in", device = "png", dpi = 600, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = tiff_path, plot = plot, width = width, height = height,
    units = "in", device = "tiff", dpi = 600, compression = "lzw",
    bg = "white", limitsize = FALSE
  )

  c(pdf_path, png_path, tiff_path)
}

# -----------------------------------------------------------------------------
# Load the existing object and fixed broad annotations
# -----------------------------------------------------------------------------

obj <- load_seurat_object(input_path)
meta <- safe_meta_data(obj)
cluster_col <- resolve_meta_column(meta, cfg$cluster_candidates, TRUE, "NormalTotal author cluster")
patient_col <- resolve_meta_column(meta, cfg$patient_candidates, TRUE, "NormalTotal patient")

if (!"RNA" %in% safe_assay_names(obj)) stop("NormalTotal lacks the RNA assay", call. = FALSE)
counts <- get_rna_matrix(obj, "counts")
data_mat <- get_rna_matrix(obj, "data")

cells <- Reduce(intersect, list(rownames(meta), colnames(counts), colnames(data_mat)))
if (!length(cells)) stop("No cells shared by metadata, RNA counts, and RNA data", call. = FALSE)

meta <- meta[cells, , drop = FALSE]
counts <- counts[, cells, drop = FALSE]
data_mat <- data_mat[, cells, drop = FALSE]

fixed_annotation <- utils::read.delim(
  annotation_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (!all(c("cluster", "celltype_final") %in% colnames(fixed_annotation))) {
  stop("Fixed annotation table lacks cluster/celltype_final", call. = FALSE)
}

fixed_annotation$cluster <- as.character(fixed_annotation$cluster)
fixed_annotation <- fixed_annotation[match(as.character(0:7), fixed_annotation$cluster), , drop = FALSE]
if (anyNA(fixed_annotation$cluster) || !identical(fixed_annotation$cluster, as.character(0:7))) {
  stop("Fixed annotation table must contain author clusters 0-7 exactly once", call. = FALSE)
}

cluster_values <- as.character(meta[[cluster_col]])
if (!all(cluster_values %in% as.character(0:7))) {
  stop("Unexpected cluster values: ",
       paste(sort(setdiff(unique(cluster_values), as.character(0:7))), collapse = ", "),
       call. = FALSE)
}

annotation_index <- match(cluster_values, fixed_annotation$cluster)
meta$celltype_broad_final <- fixed_annotation$celltype_final[annotation_index]

# -----------------------------------------------------------------------------
# Recompute the fixed >=3/6 candidate definition from raw RNA counts
# -----------------------------------------------------------------------------

asc_genes <- c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM")
lineage_sets <- list(
  "Luminal progenitor" = c("EPCAM", "KRT19", "KIT", "ALDH1A3", "LTF"),
  "Fibroblast" = c("COL1A1", "COL1A2", "DCN", "LUM", "DPT", "PDGFRA"),
  "Basal/myoepithelial" = c("KRT5", "KRT14", "KRT17", "TP63", "DST"),
  "Mature luminal" = c("KRT8", "KRT18", "ESR1", "FOXA1", "TFF1", "AGR2"),
  "Vascular endothelial" = c("PECAM1", "VWF", "CLDN5", "KDR", "PLVAP"),
  "Mural/pericyte" = c("RGS5", "PDGFRB", "CSPG4", "NOTCH3", "MYH11", "MCAM"),
  "Immune" = c("PTPRC", "CD3D", "MS4A1", "NKG7", "LST1")
)

available_both <- intersect(rownames(counts), rownames(data_mat))
missing_asc <- setdiff(asc_genes, available_both)
if (length(missing_asc)) {
  stop("One or more ASC_P genes are absent from RNA counts/data: ",
       paste(missing_asc, collapse = ", "), call. = FALSE)
}

missing_lineage <- setdiff(unique(unlist(lineage_sets, use.names = FALSE)), available_both)
if (length(missing_lineage)) {
  message("Lineage markers absent from both RNA layers and omitted: ",
          paste(missing_lineage, collapse = ", "))
}

asc_detect <- as.matrix(counts[asc_genes, cells, drop = FALSE]) > 0
meta$ASC_P_n_detected <- as.integer(colSums(asc_detect))
meta$ASC_P_like_candidate_ge3 <- cluster_values == "1" & meta$ASC_P_n_detected >= 3L

cluster1_number <- sum(cluster_values == "1")
if (!cluster1_number) stop("Author cluster 1 contains no cells", call. = FALSE)
candidate_number <- sum(meta$ASC_P_like_candidate_ge3)
candidate_percent <- 100 * candidate_number / cluster1_number

if (candidate_number != 545L || round(candidate_percent, 2) != 4.90) {
  stop(
    "Candidate audit failed. Expected 545 cells and 4.90%; recomputed ",
    candidate_number, " cells and ", sprintf("%.2f", candidate_percent),
    "%. Do not change the candidate definition; inspect the input object/layers.",
    call. = FALSE
  )
}

candidate_subtitle <- sprintf(
  "Cluster 1 fibroblasts with \u22653/6 genes detected: n = %s (%.2f%%)",
  format(candidate_number, big.mark = ","), candidate_percent
)

# -----------------------------------------------------------------------------
# Gene-level summaries for panels C and D
# -----------------------------------------------------------------------------

lineage_gene_order <- unique(unlist(lineage_sets, use.names = FALSE))
cluster_factor <- factor(cluster_values, levels = as.character(0:7))
lineage_summary <- make_summary(
  data_mat, counts, cells, cluster_factor, lineage_gene_order
)

gene_to_lineage <- setNames(
  rep(names(lineage_sets), lengths(lineage_sets)),
  unlist(lineage_sets, use.names = FALSE)
)
lineage_summary$lineage <- factor(
  unname(gene_to_lineage[as.character(lineage_summary$gene)]),
  levels = names(lineage_sets)
)
lineage_summary$gene <- factor(lineage_summary$gene, levels = lineage_gene_order)
lineage_summary$group <- factor(lineage_summary$group, levels = as.character(0:7))

group_top_to_bottom <- c(
  "Cluster 1 ASC_P-like candidates (\u22653/6)",
  "Cluster 1 other fibroblasts",
  "Cluster 4 vascular endothelial",
  "Cluster 5 mural/pericyte"
)

plot_group <- rep(NA_character_, nrow(meta))
plot_group[cluster_values == "1" & meta$ASC_P_like_candidate_ge3] <- group_top_to_bottom[[1]]
plot_group[cluster_values == "1" & !meta$ASC_P_like_candidate_ge3] <- group_top_to_bottom[[2]]
plot_group[cluster_values == "4"] <- group_top_to_bottom[[3]]
plot_group[cluster_values == "5"] <- group_top_to_bottom[[4]]

four_group <- factor(plot_group, levels = group_top_to_bottom)
four_summary <- make_summary(data_mat, counts, cells, four_group, asc_genes)
four_summary$group <- factor(four_summary$group, levels = group_top_to_bottom)
four_summary$gene <- factor(four_summary$gene, levels = asc_genes)
four_summary$mean_normalized_z_by_gene <- ave(
  four_summary$average_normalized_RNA,
  four_summary$gene,
  FUN = function(x) {
    z <- as.numeric(scale(x))
    z[!is.finite(z)] <- 0
    pmax(-2, pmin(2, z))
  }
)

# -----------------------------------------------------------------------------
# Reuse the author embedding
# -----------------------------------------------------------------------------

reduction <- if ("tsne" %in% safe_reduction_names(obj)) {
  "tsne"
} else if ("umap" %in% safe_reduction_names(obj)) {
  "umap"
} else {
  stop("Neither tsne nor umap is available", call. = FALSE)
}

reduction_label <- if (reduction == "tsne") "t-SNE" else "UMAP"
emb <- safe_embedding(obj, reduction)
emb <- emb[rownames(emb) %in% rownames(meta), , drop = FALSE]
plot_meta <- meta[rownames(emb), , drop = FALSE]

# -----------------------------------------------------------------------------
# Panel A: author clusters 0-7 with fixed broad cell-type labels
# -----------------------------------------------------------------------------

celltype_palette <- c(
  "Luminal progenitor" = "#009E73",
  "Fibroblast" = "#D55E00",
  "Basal/myoepithelial" = "#CC79A7",
  "Mature luminal" = "#56B4E9",
  "Vascular endothelial" = "#0072B2",
  "Mural/pericyte-VSMC" = "#E69F00",
  "Immune" = "#999999",
  "Unresolved PTPRC+ENG cluster" = "#7B61A8"
)

legend_celltype <- fixed_annotation$celltype_final
legend_celltype[fixed_annotation$cluster == "7"] <- "Unresolved PTPRC\u207aENG\u207a"
legend_labels <- setNames(
  paste0("Cluster ", fixed_annotation$cluster, ": ", legend_celltype),
  fixed_annotation$cluster
)

annotation_df <- data.frame(
  dim1 = emb[, 1],
  dim2 = emb[, 2],
  cluster = factor(as.character(plot_meta[[cluster_col]]), levels = as.character(0:7)),
  annotation = factor(
    legend_labels[as.character(plot_meta[[cluster_col]])],
    levels = legend_labels[as.character(0:7)]
  ),
  stringsAsFactors = FALSE
)

annotation_palette <- setNames(
  unname(celltype_palette[fixed_annotation$celltype_final]),
  legend_labels[fixed_annotation$cluster]
)

annotation_centres <- stats::aggregate(
  cbind(dim1, dim2) ~ cluster, data = annotation_df, FUN = stats::median
)

p_annotation <- ggplot2::ggplot(
  annotation_df, ggplot2::aes(dim1, dim2, colour = annotation)
) +
  ggplot2::geom_point(size = 0.16, alpha = 0.75) +
  ggplot2::geom_label(
    data = annotation_centres,
    ggplot2::aes(dim1, dim2, label = cluster),
    inherit.aes = FALSE,
    colour = "black", fill = "white", linewidth = 0.22,
    size = 5, fontface = "bold", alpha = 0.82
  ) +
  ggplot2::coord_equal(expand = TRUE) +
  ggplot2::scale_colour_manual(values = annotation_palette, drop = FALSE) +
  ggplot2::labs(
    title = "Broad cell-type annotation of normal human breast tissue",
    x = paste0(reduction_label, "-1"),
    y = paste0(reduction_label, "-2"),
    colour = NULL
  ) +
  ggplot2::theme_classic(base_size = 10, base_family = "sans") +
  ggplot2::theme(
    aspect.ratio = 1,
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 7),
    plot.title = ggplot2::element_text(face = "bold", size = 11)
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 2.2))
  )

# -----------------------------------------------------------------------------
# Panel B: fixed >=3/6 candidates within cluster 1
# -----------------------------------------------------------------------------

candidate_df <- data.frame(
  dim1 = emb[, 1],
  dim2 = emb[, 2],
  candidate = as.logical(plot_meta$ASC_P_like_candidate_ge3)
)

p_candidate_map <- ggplot2::ggplot(candidate_df, ggplot2::aes(dim1, dim2)) +
  ggplot2::geom_point(colour = "#D9D9D9", size = 0.14, alpha = 0.55) +
  ggplot2::geom_point(
    data = candidate_df[candidate_df$candidate, , drop = FALSE],
    colour = "#D55E00", size = 0.34, alpha = 0.92
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Partial ASC_P-like candidate fibroblasts in normal breast tissue",
    subtitle = candidate_subtitle,
    x = paste0(reduction_label, "-1"),
    y = paste0(reduction_label, "-2")
  ) +
  ggplot2::theme_classic(base_size = 10, base_family = "sans") +
  ggplot2::theme(
    aspect.ratio = 1,
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    plot.subtitle = ggplot2::element_text(size = 9)
  )

# -----------------------------------------------------------------------------
# Panel C: gene-level lineage evidence; cluster 0 is displayed at the top
# -----------------------------------------------------------------------------

p_lineage <- ggplot2::ggplot(
  lineage_summary,
  ggplot2::aes(gene, group, size = percent_expressed_counts,
               colour = average_normalized_RNA)
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = rev(as.character(0:7))) +
  ggplot2::scale_size(range = c(0.4, 7), name = "Cells detected (%)") +
  ggplot2::scale_colour_viridis_c(
    option = "C", name = "Mean normalized expression"
  ) +
  ggplot2::facet_grid(. ~ lineage, scales = "free_x", space = "free_x") +
  ggplot2::labs(
    title = "Lineage-marker expression across normal breast cell clusters",
    x = NULL, y = "Cluster"
  ) +
  ggplot2::theme_classic(base_size = 9, base_family = "sans") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7.5),
    axis.text.y = ggplot2::element_text(size = 8.5),
    strip.background = ggplot2::element_rect(fill = "grey94", colour = NA),
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold", size = 11)
  )

# -----------------------------------------------------------------------------
# Panel D: four lineage-context groups; candidate group is displayed at the top
# -----------------------------------------------------------------------------

p_four <- ggplot2::ggplot(
  four_summary,
  ggplot2::aes(gene, group, size = percent_expressed_counts,
               colour = mean_normalized_z_by_gene)
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = rev(group_top_to_bottom)) +
  ggplot2::scale_size(range = c(0.5, 9), name = "Cells detected (%)") +
  ggplot2::scale_colour_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-2, 2), oob = scales::squish,
    name = "Mean scaled expression"
  ) +
  ggplot2::labs(
    title = "Partial ASC_P-associated gene expression in candidate fibroblasts",
    x = NULL, y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10, base_family = "sans") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8.5),
    axis.text.y = ggplot2::element_text(size = 8.5),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    plot.margin = ggplot2::margin(4, 6, 2, 6)
  )

# -----------------------------------------------------------------------------
# Separate A-D panels and patient summary panel (new outputs with _FINAL suffix)
# -----------------------------------------------------------------------------

patient_values <- as.character(plot_meta[[patient_col]])
patient_summary <- do.call(rbind, lapply(sort(unique(patient_values)), function(patient_id) {
  idx <- patient_values == patient_id & cluster_values == "1"
  data.frame(
    patient = patient_id,
    cluster1_cell_number = sum(idx),
    n_ge_3 = sum(meta$ASC_P_like_candidate_ge3[idx]),
    pct_ge_3_of_cluster1 = 100 * safe_divide(sum(meta$ASC_P_like_candidate_ge3[idx]), sum(idx)),
    stringsAsFactors = FALSE
  )
}))
patient_summary <- patient_summary[order(patient_summary$pct_ge_3_of_cluster1, decreasing = TRUE), , drop = FALSE]
patient_summary$patient <- factor(patient_summary$patient, levels = patient_summary$patient)

p_patient <- ggplot2::ggplot(
  patient_summary,
  ggplot2::aes(patient, pct_ge_3_of_cluster1)
) +
  ggplot2::geom_col(fill = "#D55E00", width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(n_ge_3, "/", cluster1_cell_number)),
    vjust = -0.25,
    size = 2.5
  ) +
  ggplot2::labs(
    title = "Cluster 1 ASC_P-like candidate distribution by patient",
    x = NULL,
    y = "Candidate cells in cluster 1 (%)"
  ) +
  ggplot2::theme_classic(base_size = 9, base_family = "sans") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  ) +
  ggplot2::expand_limits(y = max(patient_summary$pct_ge_3_of_cluster1, na.rm = TRUE) * 1.12)

panel_outputs <- c(
  save_final_figure(p_annotation, file.path(reviewer_figure_dir, "Figure_A_NormalTotal_fixed_annotation_FINAL"), width = 10, height = 11),
  save_final_figure(p_candidate_map, file.path(reviewer_figure_dir, "Figure_C_NormalTotal_ASC_P_like_candidate_map_FINAL"), width = 10, height = 8),
  save_final_figure(p_lineage, file.path(reviewer_figure_dir, "Figure_B_NormalTotal_gene_level_lineage_DotPlot_FINAL"), width = 16, height = 6),
  save_final_figure(p_four, file.path(reviewer_figure_dir, "Figure_D_NormalTotal_ASC_P_four_group_DotPlot_FINAL"), width = 9, height = 6),
  save_final_figure(p_patient, file.path(reviewer_figure_dir, "Supplementary_NormalTotal_ASC_P_candidate_by_patient_FINAL"), width = 10, height = 6)
)

copy_output_stem <- function(src_stem, dst_stem) {
  for (ext in c("pdf", "png", "tiff")) {
    src_path <- paste0(src_stem, ".", ext)
    dst_path <- paste0(dst_stem, ".", ext)
    if (file.exists(src_path)) file.copy(src_path, dst_path, overwrite = TRUE)
  }
}

copy_output_stem(
  file.path(reviewer_figure_dir, "Figure_A_NormalTotal_fixed_annotation_FINAL"),
  file.path(reviewer_figure_dir, "Figure_A_NormalTotal_fixed_annotation")
)
copy_output_stem(
  file.path(reviewer_figure_dir, "Figure_B_NormalTotal_gene_level_lineage_DotPlot_FINAL"),
  file.path(reviewer_figure_dir, "Figure_B_NormalTotal_gene_level_lineage_DotPlot")
)
copy_output_stem(
  file.path(reviewer_figure_dir, "Figure_C_NormalTotal_ASC_P_like_candidate_map_FINAL"),
  file.path(reviewer_figure_dir, "Figure_C_NormalTotal_ASC_P_like_candidate_map")
)
copy_output_stem(
  file.path(reviewer_figure_dir, "Figure_D_NormalTotal_ASC_P_four_group_DotPlot_FINAL"),
  file.path(reviewer_figure_dir, "Figure_D_NormalTotal_ASC_P_four_group_DotPlot")
)
copy_output_stem(
  file.path(reviewer_figure_dir, "Supplementary_NormalTotal_ASC_P_candidate_by_patient_FINAL"),
  file.path(reviewer_figure_dir, "Supplementary_NormalTotal_ASC_P_candidate_by_patient")
)
copy_output_stem(
  file.path(reviewer_figure_dir, "Reviewer_Figure_NormalTotal_ASC_P_like_FINAL"),
  file.path(reviewer_figure_dir, "Reviewer_Figure_NormalTotal_ASC_P_like")
)

# -----------------------------------------------------------------------------
# Assemble and save only the final A-D reviewer figure
# -----------------------------------------------------------------------------

top_row <- p_annotation | p_candidate_map
composite <- (top_row / p_lineage / p_four) +
  patchwork::plot_layout(heights = c(1.0, 1.0, 0.65)) +
  patchwork::plot_annotation(tag_levels = "A")

composite <- composite &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 13),
    plot.tag.position = c(0, 1)
  )

final_stem <- file.path(
  reviewer_figure_dir,
  "Reviewer_Figure_NormalTotal_ASC_P_like_FINAL"
)
figure_outputs <- save_final_figure(composite, final_stem)

if (!all(vapply(c(panel_outputs, figure_outputs), is_valid_file, logical(1), min_bytes = 1000))) {
  stop("At least one final reviewer figure is missing or empty", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Dynamic final figure legend
# -----------------------------------------------------------------------------

legend_path <- file.path(
  reviewer_figure_dir,
  "Reviewer_Figure_NormalTotal_ASC_P_like_legend_FINAL.md"
)

legend_lines <- c(
  "# Figure legend", "",
  "**Identification of partial ASC_P-like candidate fibroblasts in normal human breast tissue.**",
  "",
  paste0(
    "(A) ", reduction_label,
    " visualization of broad cell-type annotations based on canonical lineage-marker expression. ",
    "(B) Localization of partial ASC_P-like candidate cells, defined as cluster 1 fibroblasts ",
    "with detectable expression of at least three of six ASC_P-associated genes ",
    "(ENG, THY1, NT5E, PDGFRA, VCAM1, and MCAM). A total of ",
    format(candidate_number, big.mark = ","),
    " candidate cells were identified, representing ",
    sprintf("%.2f", candidate_percent),
    "% of cluster 1 fibroblasts. ",
    "(C) Dot plot showing canonical lineage-marker expression across the eight author-defined clusters. ",
    "Dot size represents the percentage of cells with detectable RNA counts, and color represents ",
    "mean normalized expression. ",
    "(D) Expression of the six ASC_P-associated genes in candidate cells, other cluster 1 fibroblasts, ",
    "vascular endothelial cells, and mural/pericyte cells. Candidate fibroblasts preferentially expressed ",
    "THY1, NT5E, PDGFRA, and VCAM1, with moderate ENG and limited MCAM expression, supporting a partial ",
    "rather than complete ASC_P-like transcriptional phenotype. Because candidate cells were selected ",
    "using the same six-gene panel, panel D describes their expression composition and lineage context ",
    "rather than providing independent validation."
  )
)

write_lines_safe(legend_lines, legend_path, cfg$overwrite)

if (!is_valid_file(legend_path, min_bytes = 100)) {
  stop("Final figure legend is missing or empty", call. = FALSE)
}

outputs <- c(figure_outputs, legend_path)
end_step(
  cfg, state, outputs, status = "SUCCESS",
  notes = c(
    "Plotting-only final revision; no cells, labels, or expression values were changed.",
    paste0("Candidate audit passed: ", candidate_number, " cells; ",
           sprintf("%.2f", candidate_percent), "% of cluster 1."),
    "Only the final PDF, PNG, TIFF, and dynamic figure legend were generated."
  )
)

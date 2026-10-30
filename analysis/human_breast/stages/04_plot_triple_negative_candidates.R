# Explicit execution context supplied by run_breast_analysis.R.
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.code_root <- normalizePath(file.path(dirname(.script), ".."), winslash = "/", mustWork = TRUE)
.work_root <- Sys.getenv("ASC_BREAST_WORK_ROOT")
.data_root <- Sys.getenv("ASC_BREAST_DATA_ROOT")
if (!nzchar(.work_root) || !nzchar(.data_root)) stop("Use run_breast_analysis.R with explicit input/output arguments")
# TNBC ASC_P-like reviewer Figure A-D (revised executable version).
# Requires an approved cluster map; never uses the draft map implicitly.
# No QC, normalization, integration, reclustering, or inferential testing.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so --file is available", call. = FALSE)
project_root <- .work_root
source(file.path(.code_root, "lib", "common.R"))
cfg <- source_project_config()
if (identical(Sys.getenv("GSE161529_OVERWRITE"), "1")) cfg$overwrite <- TRUE
ensure_packages(c("SeuratObject", "Matrix", "digest", "jsonlite", "ggplot2", "patchwork", "scales"))
set.seed(cfg$seed)

object_id <- "TNBC"
input_path <- file.path(cfg$input_dir, "SeuratObject_TNBC.rds")
map_path <- file.path(cfg$table_dir, "TNBC_cluster_annotation_map_approved.tsv")
if (!file.exists(input_path)) stop("Missing TNBC RDS: ", input_path, call. = FALSE)
if (!file.exists(map_path)) stop("Missing approved TNBC map. First review TNBC_formal_cluster_annotation.tsv and save TNBC_cluster_annotation_map_approved.tsv: ", map_path, call. = FALSE)
state <- begin_step(cfg, "03b_TNBC_ASC_P_reviewer_figure_FINAL", c(input_path, map_path))

message("[03b TNBC] Stage 1/5: loading the original TNBC object and approved map")
obj <- load_seurat_object(input_path)
meta <- safe_meta_data(obj)
cluster_col <- resolve_meta_column(meta, cfg$cluster_candidates, TRUE, "TNBC author cluster")
if (!"RNA" %in% safe_assay_names(obj)) stop("TNBC lacks the RNA assay", call. = FALSE)
counts <- get_rna_matrix(obj, "counts"); data_mat <- get_rna_matrix(obj, "data")
cells <- Reduce(intersect, list(rownames(meta), colnames(counts), colnames(data_mat)))
if (!length(cells)) stop("No cells are shared by metadata, RNA counts, and RNA data", call. = FALSE)
meta <- meta[cells, , drop = FALSE]; counts <- counts[, cells, drop = FALSE]; data_mat <- data_mat[, cells, drop = FALSE]
cluster_values <- trimws(as.character(meta[[cluster_col]]))
if (anyNA(cluster_values) || any(!nzchar(cluster_values))) stop("TNBC cluster metadata contains missing/empty values", call. = FALSE)
sort_cluster_ids <- function(x) {
  numeric_x <- suppressWarnings(as.numeric(x))
  if (length(x) && all(is.finite(numeric_x))) x[order(numeric_x)] else sort(x)
}
clusters <- sort_cluster_ids(unique(cluster_values))

map <- utils::read.delim(map_path, stringsAsFactors = FALSE, check.names = FALSE)
required_map <- c("cluster", "celltype_final", "compartment", "candidate_denominator")
if (!all(required_map %in% colnames(map))) stop("Approved map must contain: ", paste(required_map, collapse = ", "), call. = FALSE)
map$cluster <- trimws(as.character(map$cluster))
map$celltype_final <- trimws(as.character(map$celltype_final))
map$compartment <- tolower(trimws(as.character(map$compartment)))
if (anyNA(map$cluster) || any(!nzchar(map$cluster))) stop("Approved map contains missing/empty cluster IDs", call. = FALSE)
if (anyDuplicated(map$cluster)) stop("Approved map contains duplicated clusters", call. = FALSE)
if (nrow(map) != length(clusters) || !setequal(map$cluster, clusters)) {
  stop("Approved map must contain every TNBC author cluster exactly once and no extra clusters", call. = FALSE)
}
if (anyNA(map$celltype_final) || any(!nzchar(map$celltype_final)) ||
    anyNA(map$compartment) || any(!nzchar(map$compartment))) {
  stop("Approved map contains missing/empty celltype_final or compartment labels", call. = FALSE)
}
map <- map[match(clusters, map$cluster), , drop = FALSE]
denom_text <- tolower(trimws(as.character(map$candidate_denominator)))
true_values <- c("true", "1", "yes", "y")
false_values <- c("false", "0", "no", "n")
if (anyNA(denom_text) || any(!denom_text %in% c(true_values, false_values))) {
  stop("candidate_denominator must contain only TRUE/FALSE (or 1/0, yes/no, y/n)", call. = FALSE)
}
denom <- denom_text %in% true_values
if (!any(denom)) stop("Approved map has no candidate_denominator=TRUE cluster", call. = FALSE)
denominator_clusters <- map$cluster[denom]
if (!identical(denominator_clusters, "5")) {
  stop(
    "TNBC ASC_P-like candidate denominator must be Cluster 5 only; approved map currently specifies: ",
    paste(denominator_clusters, collapse = ", "),
    call. = FALSE
  )
}
allowed_denominator_compartments <- c("caf", "fibroblast", "fibroblast/caf", "stromal_fibroblast")
if (any(denom & !map$compartment %in% allowed_denominator_compartments)) {
  stop("candidate_denominator=TRUE must be restricted to an approved CAF/fibroblast compartment", call. = FALSE)
}
non_caf_compartments <- c("endothelial", "vascular endothelial", "mural", "pericyte", "mural/pericyte", "immune", "epithelial", "tumor", "malignant")
if (any(denom & map$compartment %in% non_caf_compartments)) {
  stop("candidate_denominator=TRUE includes a non-CAF compartment; review the approved map", call. = FALSE)
}

asc_genes <- c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM")
missing_asc <- setdiff(asc_genes, intersect(rownames(counts), rownames(data_mat)))
if (length(missing_asc)) stop("ASC_P genes missing from one or both RNA layers: ", paste(missing_asc, collapse = ", "), call. = FALSE)
asc_detect <- as.matrix(counts[asc_genes, cells, drop = FALSE]) > 0
meta$ASC_P_n_detected <- as.integer(colSums(asc_detect))
cell_is_denominator <- denom[match(cluster_values, map$cluster)]
meta$ASC_P_like_candidate_ge3 <- cell_is_denominator & meta$ASC_P_n_detected >= 3L
n_candidates <- sum(meta$ASC_P_like_candidate_ge3); n_denom <- sum(cell_is_denominator)
if (!n_candidates || !n_denom) stop("No TNBC candidates or candidate denominator cells detected", call. = FALSE)
pct_candidates <- 100 * n_candidates / n_denom
message(sprintf("[03b TNBC] Stage 2/5: candidate definition complete (%d candidates; %.2f%% of approved Cluster 5 CAF compartment)", n_candidates, pct_candidates))

group_levels <- c("CAF ASC_P-like candidates (>=3/6)", "Other CAF", "Vascular endothelial")
comp <- as.character(map$compartment[match(cluster_values, map$cluster)])
plot_group <- rep("Other cells", length(cells))
plot_group[meta$ASC_P_like_candidate_ge3] <- group_levels[1]
plot_group[!meta$ASC_P_like_candidate_ge3 & cell_is_denominator] <- group_levels[2]
plot_group[comp %in% c("endothelial", "vascular endothelial")] <- group_levels[3]
if (!all(group_levels %in% plot_group)) {
  stop("One or more required Panel D groups are absent: ", paste(setdiff(group_levels, unique(plot_group)), collapse = ", "), call. = FALSE)
}
meta$ASC_P_plot_group <- plot_group

lineage_sets <- list(
  CAF = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "DPT", "SFRP2", "PDGFRA", "COL10A1", "COMP"),
  Mural = c("PDGFRB", "RGS5", "CSPG4", "ACTA2", "TAGLN", "MYH11", "NOTCH3", "MCAM", "CAVIN3"),
  Endothelial = c("PECAM1", "VWF", "CLDN5", "EMCN", "KDR", "PLVAP", "ESAM", "RAMP2"),
  Immune = c("PTPRC", "CD3D", "CD3E", "TRAC", "CD79A", "MS4A1", "NKG7", "GNLY", "LST1", "TYROBP"),
  Epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT14", "KRT5", "KRT17", "MUC1")
)
summary_rows <- function(groups, genes, group_levels = unique(groups)) {
  out <- list()
  for (g in group_levels) {
    g_cells <- cells[!is.na(groups) & groups == g]
    if (!length(g_cells)) next
    for (gene in genes) {
      if (!gene %in% rownames(counts) || !gene %in% rownames(data_mat)) next
      raw <- as.numeric(counts[gene, g_cells, drop = TRUE]); norm <- as.numeric(data_mat[gene, g_cells, drop = TRUE])
      out[[length(out) + 1L]] <- data.frame(group = g, gene = gene, cell_number = length(g_cells), percent_expressed_counts = 100 * mean(raw > 0), average_normalized_RNA = mean(norm), stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) stop("No group-by-gene summaries could be calculated", call. = FALSE)
  do.call(rbind, out)
}
lineage_order <- unique(unlist(lineage_sets, use.names = FALSE))
message("[03b TNBC] Stage 3/5: calculating lineage and six-gene summary tables")
lineage_summary <- summary_rows(cluster_values, lineage_order, group_levels = clusters)
lineage_summary$lineage <- unname(setNames(rep(names(lineage_sets), lengths(lineage_sets)), unlist(lineage_sets, use.names = FALSE))[as.character(lineage_summary$gene)])
lineage_summary$gene <- factor(lineage_summary$gene, levels = lineage_order)
lineage_summary$group <- factor(lineage_summary$group, levels = rev(clusters))
four_summary <- summary_rows(plot_group, asc_genes, group_levels = group_levels)
four_summary$group <- factor(four_summary$group, levels = rev(group_levels))
four_summary$gene <- factor(four_summary$gene, levels = asc_genes)
four_summary$mean_scaled_expression <- ave(four_summary$average_normalized_RNA, four_summary$gene, FUN = function(x) { z <- as.numeric(scale(x)); z[!is.finite(z)] <- 0; pmax(-2, pmin(2, z)) })

reduction <- if ("tsne" %in% safe_reduction_names(obj)) "tsne" else if ("umap" %in% safe_reduction_names(obj)) "umap" else stop("TNBC has neither tSNE nor UMAP")
reduction_label <- if (reduction == "tsne") "t-SNE" else "UMAP"
message(sprintf("[03b TNBC] Stage 4/5: constructing Figure A-D plot objects using %s", reduction_label))
emb <- safe_embedding(obj, reduction); emb <- emb[rownames(emb) %in% cells, , drop = FALSE]
if (!nrow(emb)) stop("No embedding cells overlap the aligned metadata/RNA cells", call. = FALSE)
plot_meta <- meta[rownames(emb), , drop = FALSE]; plot_clusters <- cluster_values[match(rownames(emb), cells)]
labels <- setNames(paste0("Cluster ", map$cluster, ": ", map$celltype_final), map$cluster)
annotation_df <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], annotation = factor(labels[plot_clusters], levels = labels[clusters]), cluster = plot_clusters, stringsAsFactors = FALSE)
annotation_centres <- stats::aggregate(cbind(dim1, dim2) ~ cluster, annotation_df, stats::median)
palette <- setNames(grDevices::hcl.colors(length(clusters), "Dark 3"), labels[clusters])
pA <- ggplot2::ggplot(annotation_df, ggplot2::aes(dim1, dim2, colour = annotation)) + ggplot2::geom_point(size = 0.15, alpha = 0.7) + ggplot2::geom_label(data = annotation_centres, ggplot2::aes(dim1, dim2, label = cluster), inherit.aes = FALSE, colour = "black", fill = "white", linewidth = 0.22, size = 5, fontface = "bold", alpha = 0.8) + ggplot2::coord_equal() + ggplot2::scale_colour_manual(values = palette, drop = FALSE) + ggplot2::labs(title = "Broad cell-type annotation of TNBC breast tumor tissue", x = paste0(reduction_label, "-1"), y = paste0(reduction_label, "-2"), colour = NULL) + ggplot2::theme_classic(base_size = 10) + ggplot2::theme(aspect.ratio = 1, legend.position = "bottom", legend.text = ggplot2::element_text(size = 7)) + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 2)))

candidate_df <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], candidate = plot_meta$ASC_P_like_candidate_ge3)
pB <- ggplot2::ggplot(candidate_df, ggplot2::aes(dim1, dim2)) + ggplot2::geom_point(colour = "grey85", size = 0.14, alpha = 0.55) + ggplot2::geom_point(data = candidate_df[candidate_df$candidate, , drop = FALSE], colour = "#D55E00", size = 0.34, alpha = 0.92) + ggplot2::coord_equal() + ggplot2::labs(title = "Partial ASC_P-like candidate fibroblasts in TNBC breast tumors", subtitle = sprintf("Cluster 5 CAFs with >=3/6 genes detected: n = %d (%.2f%%)", n_candidates, pct_candidates), x = paste0(reduction_label, "-1"), y = paste0(reduction_label, "-2")) + ggplot2::theme_classic(base_size = 10)
pC <- ggplot2::ggplot(lineage_summary, ggplot2::aes(gene, group, size = percent_expressed_counts, colour = average_normalized_RNA)) + ggplot2::geom_point() + ggplot2::scale_size(range = c(0.4, 7), name = "Cells detected (%)") + ggplot2::scale_colour_viridis_c(option = "C", name = "Mean normalized expression") + ggplot2::facet_grid(. ~ lineage, scales = "free_x", space = "free_x") + ggplot2::labs(title = "Lineage-marker expression across cell clusters from TNBC breast tumors", x = NULL, y = "Cluster") + ggplot2::theme_classic(base_size = 9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8), strip.background = ggplot2::element_rect(fill = "grey94", colour = NA), legend.position = "bottom")
pD <- ggplot2::ggplot(four_summary, ggplot2::aes(gene, group, size = percent_expressed_counts, colour = mean_scaled_expression)) + ggplot2::geom_point() + ggplot2::scale_size(range = c(0.5, 9), name = "Cells detected (%)") + ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "Mean scaled expression") + ggplot2::labs(title = "Partial ASC_P-associated gene expression in TNBC CAF and vascular compartments", x = NULL, y = NULL) + ggplot2::theme_classic(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = ggplot2::margin(5, 5, 5, 5))

save_plot <- function(plot, stem, width, height) {
  paths <- paste0(stem, c(".pdf", ".png", ".tiff"))
  for (path in paths) assert_new_path(path, cfg$overwrite)
  ggplot2::ggsave(paths[[1]], plot, width = width, height = height, device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf", bg = "white", limitsize = FALSE)
  ggplot2::ggsave(paths[[2]], plot, width = width, height = height, dpi = 600, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(paths[[3]], plot, width = width, height = height, dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE)
  paths
}
out_dir <- file.path(cfg$figure_dir, "reviewer_ASC_P")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("[03b TNBC] Stage 5/5: exporting Figure A-D and reviewer composite")
# Keep panel positions and tags identical to the reviewed NormalTotal mother figure:
# A, fixed annotation; B, candidate localization; C, lineage boundary DotPlot; D, six-gene context DotPlot.
figure_outputs <- c(save_plot(pA, file.path(out_dir, "Figure_A_TNBC_fixed_annotation_FINAL"), 10, 11), save_plot(pB, file.path(out_dir, "Figure_B_TNBC_ASC_P_like_candidate_map_FINAL"), 10, 8), save_plot(pC, file.path(out_dir, "Figure_C_TNBC_gene_level_lineage_DotPlot_FINAL"), 16, 7), save_plot(pD, file.path(out_dir, "Figure_D_TNBC_ASC_P_three_group_DotPlot_FINAL"), 9, 6))
composite <- (pA | pB) / pC / pD + patchwork::plot_annotation(tag_levels = "A", theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))) + patchwork::plot_layout(heights = c(1, 1, 0.65))
figure_outputs <- c(figure_outputs, save_plot(composite, file.path(out_dir, "Reviewer_Figure_TNBC_ASC_P_like_FINAL"), 18, 19))
legend_path <- file.path(out_dir, "Reviewer_Figure_TNBC_ASC_P_like_legend_FINAL.md")
legend_text <- sprintf(paste0("# Figure legend\n\n**Partial ASC_P-like candidate fibroblasts in TNBC human breast tumors.**\n\n(A) %s visualization using the approved TNBC cluster map. (B) Candidate CAFs were defined using RNA counts > 0 for at least three of six genes (ENG, THY1, NT5E, PDGFRA, VCAM1, MCAM) within the approved Cluster 5 CAF compartment; n = %d (%.2f%%). (C) Gene-level lineage-marker expression across author clusters. Dot size represents cells detected and colour represents mean normalized RNA. (D) Six-gene expression across candidate CAFs, other CAFs, and vascular endothelial cells. TNBC has no clean independent mural/pericyte author cluster, so a mural group is not forced. These are partial ASC_P-like transcriptional candidates, not a definitive human ASC_P population. Because candidate cells were selected using the same six-gene panel, panel D describes their expression composition and lineage context rather than providing independent validation."), reduction_label, n_candidates, pct_candidates)
write_lines_safe(legend_text, legend_path, cfg$overwrite)
figure_outputs <- c(figure_outputs, legend_path)
image_outputs <- setdiff(figure_outputs, legend_path)
if (!all(vapply(image_outputs, is_valid_file, logical(1), min_bytes = 1000)) || !is_valid_file(legend_path, min_bytes = 100)) stop("One or more TNBC Figure A-D outputs are missing or empty", call. = FALSE)
end_step(cfg, state, figure_outputs, status = "SUCCESS", notes = c(sprintf("Candidate count: %d; %.2f%% of approved Cluster 5 CAF compartment.", n_candidates, pct_candidates), "Approved TNBC mapping and six-gene raw-count detection were used; no original RDS was modified. TNBC has no mural/pericyte author cluster; Panel D uses three groups."))


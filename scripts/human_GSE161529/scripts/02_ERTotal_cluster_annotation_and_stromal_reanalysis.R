# ERTotal cluster annotation draft and ASC_P candidate precheck.
# Read-only descriptive analysis: reuses the author object, RNA counts/data,
# author cluster, and existing t-SNE/UMAP. No QC, normalization, integration,
# dimensional reduction, clustering, or inferential patient-level testing.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so --file is available", call. = FALSE)
project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "scripts", "utils_GSE161529_ASC_P.R"))
cfg <- source_project_config()
ensure_packages(c("SeuratObject", "Matrix", "digest", "jsonlite", "ggplot2", "pheatmap"))
set.seed(cfg$seed)

object_id <- "ERTotal"
input_path <- file.path(cfg$input_dir, "SeuratObject_ERTotal.rds")
marker_path <- file.path(cfg$table_dir, "ERTotal_cluster_FindMarkers.tsv")
if (!file.exists(input_path)) stop("Missing ERTotal RDS: ", input_path, call. = FALSE)
if (!file.exists(marker_path)) stop("Run scripts/01c_ERTotal_FindMarkers.R first; missing: ", marker_path, call. = FALSE)
state <- begin_step(cfg, "02_ERTotal_cluster_annotation_and_stromal_reanalysis", c(input_path, marker_path))

obj <- load_seurat_object(input_path)
meta <- safe_meta_data(obj)
cluster_col <- resolve_meta_column(meta, cfg$cluster_candidates, TRUE, "ERTotal author cluster")
patient_col <- resolve_meta_column(meta, cfg$patient_candidates, TRUE, "ERTotal patient")
sample_col <- resolve_meta_column(meta, cfg$sample_candidates, TRUE, "ERTotal sample")
if (!"RNA" %in% safe_assay_names(obj)) stop("ERTotal lacks RNA assay", call. = FALSE)
counts <- get_rna_matrix(obj, "counts")
data_mat <- get_rna_matrix(obj, "data")
cells <- Reduce(intersect, list(rownames(meta), colnames(counts), colnames(data_mat)))
if (!length(cells)) stop("No cells shared by metadata, RNA counts, and RNA data", call. = FALSE)
meta <- meta[cells, , drop = FALSE]
counts <- counts[, cells, drop = FALSE]
data_mat <- data_mat[, cells, drop = FALSE]
cluster_values <- as.character(meta[[cluster_col]])
clusters <- sort(unique(cluster_values))

markers <- utils::read.delim(marker_path, stringsAsFactors = FALSE, check.names = FALSE)
required_marker_cols <- c("author_cluster", "gene", "p_val_adj", "avg_log2FC")
if (!all(required_marker_cols %in% colnames(markers))) stop("ERTotal FindMarkers table lacks required columns", call. = FALSE)

marker_sets <- list(
  CAF_fibroblast = cfg$fibroblast_markers,
  mural_pericyte = cfg$mural_markers,
  endothelial = cfg$endothelial_markers,
  immune = c("PTPRC", "CD3D", "CD3E", "TRAC", "CD8A", "CD4", "MS4A1", "CD79A", "CD68", "CD14", "TYROBP", "LST1", "FCER1G", "C1QC", "NKG7", "GNLY", "MZB1", "JCHAIN"),
  epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT14", "KRT5", "KRT17", "KRT7", "MUC1", "KRT1")
)
score_rows <- list()
for (cl in clusters) {
  cl_cells <- cells[cluster_values == cl]
  for (lineage in names(marker_sets)) {
    genes <- intersect(marker_sets[[lineage]], rownames(data_mat))
    if (!length(genes)) next
    block <- data_mat[genes, cl_cells, drop = FALSE]
    score_rows[[length(score_rows) + 1L]] <- data.frame(
      object_id = object_id, cluster = cl, lineage = lineage,
      marker_gene_number = length(genes), percent_expressed = 100 * mean(as.matrix(block) > 0),
      mean_normalized_RNA = mean(as.matrix(block)), cell_number = length(cl_cells),
      patient_number = length(unique(as.character(meta[cl_cells, patient_col]))), stringsAsFactors = FALSE
    )
  }
}
scores <- do.call(rbind, score_rows)

# Draft mapping is deliberately conservative and requires manual review.
# It is a compartment hypothesis, not an approved biological annotation.
wide <- reshape(scores[, c("cluster", "lineage", "percent_expressed", "mean_normalized_RNA")],
                timevar = "lineage", idvar = "cluster", direction = "wide")
getv <- function(cl, lineage, metric = "mean_normalized_RNA") {
  col <- paste0(metric, ".", lineage)
  if (!col %in% colnames(wide)) return(NA_real_)
  x <- wide[wide$cluster == cl, col]
  if (!length(x) || !is.finite(x[[1]])) NA_real_ else x[[1]]
}
draft_compartment <- vapply(clusters, function(cl) {
  vals <- vapply(names(marker_sets), function(nm) getv(cl, nm), numeric(1))
  if (all(!is.finite(vals))) return("unknown")
  top <- names(which.max(replace(vals, !is.finite(vals), -Inf)))
  if (top == "epithelial_tumor") "epithelial/tumor" else if (top == "immune") "immune" else if (top == "endothelial") "endothelial" else if (top == "mural_pericyte") "mural" else if (top == "CAF_fibroblast") "CAF" else "unknown"
}, character(1))

top_evidence <- vapply(clusters, function(cl) {
  z <- markers[as.character(markers$author_cluster) == cl, , drop = FALSE]
  z <- z[order(z$p_val_adj, -z$avg_log2FC), , drop = FALSE]
  paste(utils::head(z$gene, 10L), collapse = ";")
}, character(1))
mixed_flag <- vapply(clusters, function(cl) {
  z <- scores[scores$cluster == cl, , drop = FALSE]
  vals <- z$mean_normalized_RNA[match(c("CAF_fibroblast", "mural_pericyte", "endothelial", "immune", "epithelial_tumor"), z$lineage)]
  sum(is.finite(vals) & vals > 0.10) >= 2L
}, logical(1))
cluster_cells <- data.frame(cluster = clusters, cell_number = vapply(clusters, function(cl) sum(cluster_values == cl), integer(1)), stringsAsFactors = FALSE)
cluster_cells$patient_number <- vapply(clusters, function(cl) length(unique(as.character(meta[cluster_values == cl, patient_col]))), integer(1))
map <- data.frame(
  object_id = object_id, cluster = clusters, draft_compartment = draft_compartment,
  cell_number = cluster_cells$cell_number, patient_number = cluster_cells$patient_number,
  mixed_flag = mixed_flag, candidate_denominator_draft = draft_compartment == "CAF",
  evidence_top10_FindMarkers = top_evidence, mapping_status = "draft_requires_manual_review",
  stringsAsFactors = FALSE
)

map_path <- file.path(cfg$table_dir, "ERTotal_cluster_annotation_map_draft.tsv")
score_path <- file.path(cfg$table_dir, "ERTotal_cluster_marker_scores.tsv")
top_path <- file.path(cfg$table_dir, "ERTotal_cluster_top15_markers.tsv")
write_tsv_safe(map, map_path, cfg$overwrite)
write_tsv_safe(scores, score_path, cfg$overwrite)
top <- do.call(rbind, lapply(split(markers, markers$author_cluster), function(x) utils::head(x[order(x$p_val_adj, -x$avg_log2FC), , drop = FALSE], 15L)))
write_tsv_safe(top, top_path, cfg$overwrite)

# Draft-only object copy; never use this as the approved mapping input.
meta$celltype_broad_draft <- map$draft_compartment[match(cluster_values, map$cluster)]
meta$ERTotal_mapping_mixed_flag <- mixed_flag[match(cluster_values, map$cluster)]
obj@meta.data <- meta
draft_obj_path <- file.path(cfg$updated_object_dir, "SeuratObject_ERTotal_annotation_draft.rds")
save_rds_safe(obj, draft_obj_path, cfg$overwrite)

# Existing embedding with numbered author clusters and draft compartment legend.
reduction <- if ("tsne" %in% safe_reduction_names(obj)) "tsne" else if ("umap" %in% safe_reduction_names(obj)) "umap" else stop("ERTotal has neither tSNE nor UMAP")
emb <- safe_embedding(obj, reduction); emb <- emb[intersect(rownames(emb), cells), , drop = FALSE]
plot_df <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], cluster = factor(cluster_values[match(rownames(emb), cells)], levels = clusters), stringsAsFactors = FALSE)
p <- ggplot2::ggplot(plot_df, ggplot2::aes(dim1, dim2, colour = cluster)) + ggplot2::geom_point(size = 0.12, alpha = 0.7) + ggplot2::coord_equal() + ggplot2::theme_classic(base_size = 10) + ggplot2::labs(title = "ERTotal author clusters: draft compartment review", subtitle = "Mapping is provisional; approve cluster boundaries before candidate analysis", x = paste0(toupper(reduction), "-1"), y = paste0(toupper(reduction), "-2"), colour = "Author cluster")
fig_stem <- file.path(cfg$figure_dir, "identity_boundary", "ERTotal_cluster_annotation_draft")
fig_pdf <- paste0(fig_stem, ".pdf"); fig_png <- paste0(fig_stem, ".png")
assert_new_path(fig_pdf, cfg$overwrite); assert_new_path(fig_png, cfg$overwrite)
ggplot2::ggsave(fig_pdf, p, width = 10, height = 8, device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf", bg = "white")
ggplot2::ggsave(fig_png, p, width = 10, height = 8, dpi = 300, bg = "white")

report_path <- file.path(cfg$report_dir, "02_ERTotal_CLUSTER_ANNOTATION_DRAFT.md")
report <- c("# ERTotal cluster annotation draft", "", "This is a provisional compartment mapping for review; it is not an approved biological annotation.", "No QC, normalization, integration, dimensional reduction, clustering, or inferential patient-level testing was performed.", "", paste0("Input SHA256: `", sha256_file(input_path), "`"), paste0("Cluster column: `", cluster_col, "`; clusters: ", paste(clusters, collapse = ", ")), paste0("Patient column: `", patient_col, "`; patients: ", length(unique(meta[[patient_col]]))), paste0("FindMarkers source: `", normalizePath(marker_path, winslash = "/"), "`"), "", "Review each cluster using top FindMarkers, marker scores, patient support, and mixed-lineage flags.", "Do not run the Figure A-D script until an approved map has been created from the draft table.")
write_lines_safe(report, report_path, cfg$overwrite)
end_step(cfg, state, c(map_path, score_path, top_path, draft_obj_path, fig_pdf, fig_png, report_path), status = "SUCCESS", notes = c("Draft mapping only; manual approval required.", "Original ERTotal RDS was not modified."))

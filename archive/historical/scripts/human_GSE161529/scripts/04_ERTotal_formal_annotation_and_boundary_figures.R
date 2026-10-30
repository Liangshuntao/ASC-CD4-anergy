# Apply reviewed ERTotal mapping and create formal annotation boundary figures.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "scripts", "utils_GSE161529_ASC_P.R"))
cfg <- source_project_config()
if (identical(Sys.getenv("GSE161529_OVERWRITE"), "1")) cfg$overwrite <- TRUE
ensure_packages(c("SeuratObject", "Matrix", "digest", "jsonlite", "ggplot2"))
set.seed(cfg$seed)
input_path <- file.path(cfg$input_dir, "SeuratObject_ERTotal.rds")
map_path <- file.path(cfg$table_dir, "ERTotal_formal_cluster_annotation.tsv")
if (!file.exists(input_path)) stop("Missing ERTotal RDS", call. = FALSE)
if (!file.exists(map_path)) stop("Missing fixed annotation TSV", call. = FALSE)
state <- begin_step(cfg, "04_ERTotal_formal_annotation_and_boundary_figures", c(input_path, map_path))
obj <- load_seurat_object(input_path)
meta <- safe_meta_data(obj)
cluster_col <- resolve_meta_column(meta, cfg$cluster_candidates, TRUE, "ERTotal cluster")
patient_col <- resolve_meta_column(meta, cfg$patient_candidates, TRUE, "ERTotal patient")
counts <- get_rna_matrix(obj, "counts"); data_mat <- get_rna_matrix(obj, "data")
cells <- Reduce(intersect, list(rownames(meta), colnames(counts), colnames(data_mat)))
meta <- meta[cells, , drop = FALSE]; counts <- counts[, cells, drop = FALSE]; data_mat <- data_mat[, cells, drop = FALSE]
cluster_values <- as.character(meta[[cluster_col]]); clusters <- sort(unique(cluster_values))
map <- utils::read.delim(map_path, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("analysis_id", "author_cluster", "cell_number", "cell_fraction_percent", "patient_number", "broad_lineage", "formal_cell_type", "malignant_status", "confidence", "positive_marker_evidence", "mixed_cluster_flag", "recommended_next_action")
if (!all(need %in% colnames(map))) stop("Fixed mapping lacks required columns", call. = FALSE)
map$author_cluster <- as.character(map$author_cluster); map <- map[match(clusters, map$author_cluster), , drop = FALSE]
if (anyNA(map$author_cluster) || anyDuplicated(map$author_cluster) || !identical(map$author_cluster, clusters)) stop("Mapping does not match ERTotal clusters exactly", call. = FALSE)
meta$ERTotal_celltype_final <- map$formal_cell_type[match(cluster_values, map$author_cluster)]
meta$ERTotal_broad_lineage <- map$broad_lineage[match(cluster_values, map$author_cluster)]
meta$ERTotal_malignant_status <- map$malignant_status[match(cluster_values, map$author_cluster)]
meta$ERTotal_annotation_confidence <- map$confidence[match(cluster_values, map$author_cluster)]
meta$ERTotal_mixed_cluster_flag <- tolower(as.character(map$mixed_cluster_flag[match(cluster_values, map$author_cluster)])) %in% c("true", "1", "yes")
obj@meta.data <- meta
annotated_path <- file.path(cfg$updated_object_dir, "SeuratObject_ERTotal_formally_annotated.rds")
save_rds_safe(obj, annotated_path, cfg$overwrite)
evidence_path <- file.path(cfg$table_dir, "ERTotal_cluster_annotation_evidence.tsv")
write_tsv_safe(map, evidence_path, cfg$overwrite)

labels <- setNames(paste0("Cluster ", map$author_cluster, ": ", map$formal_cell_type), map$author_cluster)
palette <- setNames(grDevices::hcl.colors(length(clusters), "Dark 3"), labels[clusters])
plot_embedding <- function(reduction, stem, title) {
  emb <- safe_embedding(obj, reduction); emb <- emb[intersect(rownames(emb), cells), , drop = FALSE]
  d <- data.frame(x = emb[, 1], y = emb[, 2], label = factor(labels[cluster_values[match(rownames(emb), cells)]], levels = labels[clusters]), cluster = cluster_values[match(rownames(emb), cells)])
  cen <- stats::aggregate(cbind(x, y) ~ cluster, d, stats::median)
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, colour = label)) + ggplot2::geom_point(size = 0.14, alpha = 0.7) + ggplot2::geom_label(data = cen, ggplot2::aes(x, y, label = cluster), inherit.aes = FALSE, colour = "black", fill = "white", linewidth = 0.2, size = 4, fontface = "bold", alpha = 0.8) + ggplot2::coord_equal() + ggplot2::scale_colour_manual(values = palette, drop = FALSE) + ggplot2::labs(title = title, x = paste0(toupper(reduction), "-1"), y = paste0(toupper(reduction), "-2"), colour = NULL) + ggplot2::theme_classic(base_size = 10) + ggplot2::theme(aspect.ratio = 1, legend.position = "bottom", legend.text = ggplot2::element_text(size = 7)) + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 2)))
  pdf <- file.path(cfg$figure_dir, "identity_boundary", paste0(stem, ".pdf")); png <- file.path(cfg$figure_dir, "identity_boundary", paste0(stem, ".png")); assert_new_path(pdf, cfg$overwrite); assert_new_path(png, cfg$overwrite)
  ggplot2::ggsave(pdf, p, width = 11, height = 10, device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf", bg = "white"); ggplot2::ggsave(png, p, width = 11, height = 10, dpi = 300, bg = "white"); c(pdf, png)
}
figure_outputs <- character(); reds <- safe_reduction_names(obj)
if ("tsne" %in% reds) figure_outputs <- c(figure_outputs, plot_embedding("tsne", "ERTotal_final_annotation_tSNE", "ERTotal formal cluster annotation (t-SNE)"))
if ("umap" %in% reds) figure_outputs <- c(figure_outputs, plot_embedding("umap", "ERTotal_final_annotation_UMAP", "ERTotal formal cluster annotation (UMAP)"))

lineage_sets <- list(Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "ESR1", "PGR", "FOXA1", "GATA3", "MUC1", "KRT5", "KRT14", "KRT17", "TP63"), CAF = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL11A1", "LRRC15", "SFRP2", "CTHRC1", "POSTN", "MMP11", "FAP"), Endothelial = c("PECAM1", "VWF", "CLDN5", "CDH5", "ACKR1", "EMCN", "PLVAP", "FLT1", "PTPRB", "CLEC14A", "ENG"), Mural = c("RGS5", "CSPG4", "MCAM", "PDGFRB", "ACTA2", "MYH11", "NOTCH3", "RERGL", "COX4I2", "CARMN"), Immune = c("PTPRC", "CD3D", "CD3E", "TRAC", "CD8A", "IL7R", "NKG7", "GNLY", "C1QA", "C1QB", "C1QC", "TREM2", "LYZ", "JCHAIN", "MZB1", "CD79A"))
genes <- intersect(unique(unlist(lineage_sets, use.names = FALSE)), intersect(rownames(counts), rownames(data_mat)))
rows <- list()
for (cl in clusters) for (gene in genes) { z <- cells[cluster_values == cl]; rows[[length(rows) + 1L]] <- data.frame(cluster = cl, gene = gene, percent_expressed_counts = 100 * mean(as.numeric(counts[gene, z]) > 0), mean_normalized_RNA = mean(as.numeric(data_mat[gene, z]))) }
dot <- do.call(rbind, rows); dot$cluster <- factor(dot$cluster, levels = rev(clusters)); dot$gene <- factor(dot$gene, levels = genes); dot$lineage <- unname(setNames(rep(names(lineage_sets), lengths(lineage_sets)), unlist(lineage_sets, use.names = FALSE))[as.character(dot$gene)])
p_dot <- ggplot2::ggplot(dot, ggplot2::aes(gene, cluster, size = percent_expressed_counts, colour = mean_normalized_RNA)) + ggplot2::geom_point() + ggplot2::scale_size(range = c(0.4, 7), name = "Cells detected (%)") + ggplot2::scale_colour_viridis_c(option = "C", name = "Mean normalized expression") + ggplot2::facet_grid(. ~ lineage, scales = "free_x", space = "free_x") + ggplot2::labs(title = "ERTotal lineage-marker expression across formal clusters", x = NULL, y = "Author cluster") + ggplot2::theme_classic(base_size = 9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7), strip.background = ggplot2::element_rect(fill = "grey94", colour = NA), legend.position = "bottom")
dot_pdf <- file.path(cfg$figure_dir, "identity_boundary", "ERTotal_lineage_marker_DotPlot.pdf"); dot_png <- file.path(cfg$figure_dir, "identity_boundary", "ERTotal_lineage_marker_DotPlot.png"); assert_new_path(dot_pdf, cfg$overwrite); assert_new_path(dot_png, cfg$overwrite); ggplot2::ggsave(dot_pdf, p_dot, width = 16, height = 8, bg = "white"); ggplot2::ggsave(dot_png, p_dot, width = 16, height = 8, dpi = 300, bg = "white"); figure_outputs <- c(figure_outputs, dot_pdf, dot_png)

report_path <- file.path(cfg$report_dir, "ERTotal_FORMAL_ANNOTATION_APPLY_REPORT.md")
report <- c("# ERTotal formal annotation application report", "", "Reviewed fixed cluster mapping was written to a new Seurat object copy.", "The original ERTotal RDS was not modified. No ERTotalSub or ASC_P-like analysis was performed.", "", paste0("Input SHA256: `", sha256_file(input_path), "`"), paste0("Mapping SHA256: `", sha256_file(map_path), "`"), paste0("Clusters: ", paste(clusters, collapse = ", "), "; cells: ", nrow(meta), "; patients: ", length(unique(as.character(meta[[patient_col]])))), paste0("Annotated copy: `", normalizePath(annotated_path, winslash = "/"), "`"), "", "Boundary reminders: Cluster 1 remains mixed T/NK; Clusters 5/6 retain provisional malignant-status language; Cluster 9 osteoclast-like is a transcriptional-state label; Cluster 10 is small and requires cautious patient-level interpretation.")
write_lines_safe(report, report_path, cfg$overwrite)
outputs <- c(annotated_path, evidence_path, figure_outputs, report_path)
end_step(cfg, state, outputs, status = "SUCCESS", notes = c("Fixed ERTotal mapping applied to a new object copy.", "No ERTotalSub or ASC_P-like analysis performed."))

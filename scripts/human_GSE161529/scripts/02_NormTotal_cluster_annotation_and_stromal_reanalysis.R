source(file.path(normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE), "scripts", "utils_GSE161529_ASC_P.R"))
cfg <- source_project_config()
if (identical(Sys.getenv("GSE161529_OVERWRITE"), "1")) cfg$overwrite <- TRUE
ensure_packages(c("SeuratObject", "Matrix", "digest", "jsonlite", "ggplot2", "pheatmap"))
set.seed(cfg$seed)

input_path <- file.path(cfg$input_dir, "SeuratObject_NormTotal.rds")
if (!file.exists(input_path)) stop("Missing NormalTotal RDS: ", input_path, call. = FALSE)
state <- begin_step(cfg, "02_NormTotal_cluster_annotation_and_stromal_reanalysis", input_path)

obj <- load_seurat_object(input_path)
meta <- safe_meta_data(obj)
cluster_col <- resolve_meta_column(meta, cfg$cluster_candidates, TRUE, "NormalTotal cluster")
patient_col <- resolve_meta_column(meta, cfg$patient_candidates, TRUE, "NormalTotal patient")
sample_col <- resolve_meta_column(meta, c("orig.ident", "sample", "sample_id", "group"), TRUE, "NormalTotal sample")
if (!"RNA" %in% safe_assay_names(obj)) stop("RNA assay missing", call. = FALSE)
clusters <- sort(unique(as.character(meta[[cluster_col]])))

marker_sets <- list(
  basal_myoepithelial = c("KRT5", "KRT14", "KRT17", "TP63", "ACTA2", "TAGLN", "MYLK", "DST", "SNAI2"),
  luminal_progenitor = c("EPCAM", "KRT19", "KIT", "ELF5", "TNFRSF11A", "ALDH1A3", "LTF"),
  mature_luminal = c("KRT8", "KRT18", "ESR1", "PGR", "FOXA1", "GATA3", "TFF1", "AGR2"),
  fibroblast_stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "DPT", "PDGFRA", "SFRP2", "COL14A1", "CD34"),
  vascular_endothelial = c("PECAM1", "VWF", "CLDN5", "KDR", "ESAM", "EMCN", "PLVAP", "RAMP2", "ENG"),
  lymphatic_endothelial = c("PDPN", "CCL21", "PROX1", "FLT4", "LYVE1"),
  mural_perivascular = c("RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "MYH11", "ACTA2", "TAGLN", "DES"),
  immune = c("PTPRC", "CD3D", "CD3E", "TRBC1", "IL7R", "CD4", "CD8A", "CD79A", "MS4A1", "CD37", "NKG7", "GNLY", "LST1", "TYROBP", "FCER1G", "CTSS", "CD68", "CD14")
)
asc_genes <- c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM")
all_genes <- unique(c(unlist(marker_sets, use.names = FALSE), asc_genes))
present_data <- intersect(all_genes, rownames(get_rna_matrix(obj, "data")))
present_counts <- intersect(all_genes, rownames(get_rna_matrix(obj, "counts")))
missing_genes <- setdiff(all_genes, union(present_data, present_counts))
data_mat <- extract_gene_matrix(obj, present_data, "data")
counts_mat <- extract_gene_matrix(obj, present_counts, "counts")
cell_ids <- intersect(rownames(meta), intersect(colnames(data_mat), colnames(counts_mat)))
meta <- meta[cell_ids, , drop = FALSE]
data_mat <- data_mat[, cell_ids, drop = FALSE]
counts_mat <- counts_mat[, cell_ids, drop = FALSE]

cluster_summary <- function(mat, genes, layer_name) {
  genes <- intersect(genes, rownames(mat))
  if (!length(genes)) return(data.frame())
  out <- list()
  for (cl in clusters) {
    cells <- rownames(meta)[as.character(meta[[cluster_col]]) == cl]
    block <- mat[genes, cells, drop = FALSE]
    out[[length(out) + 1L]] <- data.frame(
      analysis_id = "NormalTotal", cluster = cl, module = layer_name,
      gene_number_present = length(genes), percent_expressed = mean(as.matrix(block) > 0) * 100,
      mean_normalized_RNA = mean(as.matrix(block)),
      cell_number = length(cells), patient_number = length(unique(as.character(meta[cells, patient_col]))),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

module_scores <- do.call(rbind, lapply(names(marker_sets), function(nm) cluster_summary(data_mat, marker_sets[[nm]], nm)))
marker_availability <- data.frame(
  gene = all_genes, data_present = all_genes %in% present_data, counts_present = all_genes %in% present_counts,
  stringsAsFactors = FALSE
)
write_tsv_safe(marker_availability, file.path(cfg$table_dir, "NormTotal_marker_availability.tsv"), cfg$overwrite)

# Cluster-level cell numbers and patient/sample composition.
cluster_cells <- as.data.frame(table(cluster = as.character(meta[[cluster_col]])), stringsAsFactors = FALSE)
colnames(cluster_cells)[2] <- "cell_number"
cluster_cells$percentage <- 100 * cluster_cells$cell_number / nrow(meta)
cluster_cells$patient_number <- vapply(cluster_cells$cluster, function(cl) length(unique(as.character(meta[as.character(meta[[cluster_col]]) == cl, patient_col]))), integer(1))
write_tsv_safe(cluster_cells, file.path(cfg$table_dir, "NormTotal_cluster_cell_numbers.tsv"), cfg$overwrite)
sample_cluster <- as.data.frame(table(sample = as.character(meta[[sample_col]]), cluster = as.character(meta[[cluster_col]])), stringsAsFactors = FALSE)
colnames(sample_cluster)[3] <- "cell_number"
write_tsv_safe(sample_cluster, file.path(cfg$table_dir, "NormTotal_sample_by_cluster.tsv"), cfg$overwrite)

# Evidence-based broad annotations. Cluster 0/3 and 2 are decided using module scores;
# cluster 7 remains unresolved unless cell-level evidence is clearly non-mixed.
score_wide <- reshape(module_scores[, c("cluster", "module", "percent_expressed")],
                      timevar = "module", idvar = "cluster", direction = "wide")
get_score <- function(cl, module) {
  col <- paste0("percent_expressed.", module)
  if (!col %in% colnames(score_wide)) return(0)
  val <- score_wide[score_wide$cluster == cl, col]
  if (!length(val) || !is.finite(val[[1]])) 0 else val[[1]]
}
ann <- data.frame(cluster = clusters, stringsAsFactors = FALSE)
ann$celltype_final <- vapply(clusters, function(cl) {
  if (cl == "0") return(if (get_score(cl, "luminal_progenitor") >= get_score(cl, "mature_luminal")) "Luminal progenitor" else "Luminal epithelial 1")
  if (cl == "1") return("Fibroblast")
  if (cl == "2") return(if (get_score(cl, "basal_myoepithelial") >= get_score(cl, "mural_perivascular")) "Basal/myoepithelial" else "Mural/pericyte-VSMC")
  if (cl == "3") return(if (get_score(cl, "mature_luminal") >= get_score(cl, "luminal_progenitor")) "Mature luminal" else "Luminal epithelial 2")
  if (cl == "4") return("Vascular endothelial")
  if (cl == "5") return("Mural/pericyte-VSMC")
  if (cl == "6") return("Immune")
  "Unresolved PTPRC+ENG cluster"
}, character(1))
ann$confidence <- c("moderate", "high", "moderate", "moderate", "high", "high", "high", "low")
ann$positive_markers <- c(
  "EPCAM;KRT19;KRT18;KRT8", "DPT;LUM;DCN;COL1A1;COL1A2;PDGFRA;SFRP2;COL3A1",
  "ACTA2;TAGLN;KRT5/KRT14/TP63 review", "KRT8;KRT18;KRT19;EPCAM",
  "PLVAP;CLDN5;PECAM1;VWF;ESAM;KDR;EMCN", "RGS5;PDGFRB;TAGLN;MYH11;MCAM;NOTCH3",
  "PTPRC;immune lineage panel", "PTPRC;ENG"
)
ann$exclusion_markers <- c(
  "stromal;immune;endothelial", "EPCAM/KRT;PTPRC;PECAM1", "RGS5/PDGFRB/MYH11 versus basal panel",
  "stromal;immune;endothelial", "KRT/EPCAM;PTPRC", "EPCAM/KRT;PTPRC;PECAM1/VWF",
  "EPCAM/KRT;PECAM1/VWF", "lymphatic/endothelial and doublet review"
)
ann$patient_reproducibility <- vapply(ann$cluster, function(cl) cluster_cells$patient_number[match(cl, cluster_cells$cluster)], integer(1))
ann$main_evidence <- c(
  "Author epithelial exclusion plus luminal module comparison", "FindMarkers and fibroblast module dominance",
  "Basal/myoepithelial versus mural module comparison", "Author epithelial exclusion plus luminal module comparison",
  "Endothelial module dominance", "Mural module dominance", "PTPRC marker and immune module dominance",
  "PTPRC+ENG mixed boundary; cell-level co-detection required"
)
ann$remaining_uncertainty <- c(
  "Progenitor versus mature luminal distinction is provisional", "None major at broad level",
  "ACTA2/TAGLN can reflect mural cells; basal markers must be checked", "Progenitor versus mature luminal distinction is provisional",
  "ENG/MCAM do not override canonical endothelial program", "Potential endothelial overlap is secondary", "Fine immune subtype not resolved", "Do not force lymphatic, immune, or doublet label without cell-level evidence"
)
write_tsv_safe(ann, file.path(cfg$table_dir, "NormTotal_final_cluster_annotation.tsv"), cfg$overwrite)
evidence <- merge(ann, module_scores, by = "cluster", all.x = TRUE)
write_tsv_safe(evidence, file.path(cfg$table_dir, "NormTotal_annotation_evidence.tsv"), cfg$overwrite)

# Cell-level ASC_P detection from raw RNA counts (>0), including all signature genes.
asc_present <- intersect(asc_genes, rownames(counts_mat))
asc_detect <- as.matrix(counts_mat[asc_present, , drop = FALSE]) > 0
asc_n <- colSums(asc_detect)
asc_all6 <- if (length(asc_present) == length(asc_genes)) colSums(asc_detect) == length(asc_genes) else rep(FALSE, length(cell_ids))
cell_ann <- ann$celltype_final[match(as.character(meta[[cluster_col]]), ann$cluster)]
cell_df <- data.frame(cell = cell_ids, patient = as.character(meta[[patient_col]]), sample = as.character(meta[[sample_col]]), cluster = as.character(meta[[cluster_col]]), celltype_final = cell_ann,
                      n_ASC_P_detected = as.integer(asc_n), ASC_P_all6_detected = asc_all6,
                      nCount_RNA = if ("nCount_RNA" %in% names(meta)) meta$nCount_RNA else NA_real_,
                      nFeature_RNA = if ("nFeature_RNA" %in% names(meta)) meta$nFeature_RNA else NA_real_,
                      stringsAsFactors = FALSE)
rownames(cell_df) <- cell_ids
cell_df$ASC_P_candidate_stromal <- !cell_df$celltype_final %in% c("Immune", "Vascular endothelial", "Unresolved PTPRC+ENG cluster", "Luminal progenitor", "Mature luminal", "Luminal epithelial 1", "Luminal epithelial 2", "Basal/myoepithelial") & cell_df$n_ASC_P_detected >= 2L
write_tsv_safe(cell_df[cell_df$ASC_P_candidate_stromal, , drop = FALSE], file.path(cfg$table_dir, "NormTotal_ASC_candidate_cells.tsv"), cfg$overwrite)

# Patient-level support and donor composition; patients, not cells, are the reproducibility unit.
patient_cluster <- do.call(rbind, lapply(sort(unique(cell_df$patient)), function(pt) {
  pt_idx <- cell_df$patient == pt
  pt_total <- sum(pt_idx)
  do.call(rbind, lapply(clusters, function(cl) {
    idx <- pt_idx & cell_df$cluster == cl
    data.frame(patient = pt, cluster = cl, cell_number = sum(idx), patient_total_cells = pt_total,
               percent_within_patient = if (pt_total) 100 * sum(idx) / pt_total else NA_real_,
               n_ge_2 = sum(cell_df$n_ASC_P_detected[idx] >= 2), n_ge_3 = sum(cell_df$n_ASC_P_detected[idx] >= 3),
               n_eq_6 = sum(cell_df$ASC_P_all6_detected[idx]), stringsAsFactors = FALSE)
  }))
}))
write_tsv_safe(patient_cluster, file.path(cfg$table_dir, "NormTotal_patient_by_cluster.tsv"), cfg$overwrite)
candidate_background <- cell_df$celltype_final %in% c("Fibroblast", "Mural/pericyte-VSMC")
candidate_patient <- do.call(rbind, lapply(sort(unique(cell_df$patient)), function(pt) {
  idx <- cell_df$patient == pt
  bg <- idx & candidate_background
  cand <- idx & cell_df$ASC_P_candidate_stromal
  data.frame(patient = pt, stromal_mural_background_cells = sum(bg), candidate_cells_ge_2 = sum(cand),
             candidate_percent_of_stromal_mural = if (sum(bg)) 100 * sum(cand) / sum(bg) else NA_real_,
             candidate_n_ge_3 = sum(cand & cell_df$n_ASC_P_detected >= 3), candidate_n_eq_6 = sum(cand & cell_df$ASC_P_all6_detected),
             stringsAsFactors = FALSE)
}))
write_tsv_safe(candidate_patient, file.path(cfg$table_dir, "NormTotal_ASC_candidate_by_patient.tsv"), cfg$overwrite)

asc_cluster <- do.call(rbind, lapply(clusters, function(cl) {
  idx <- cell_df$cluster == cl
  data.frame(cluster = cl, cell_number = sum(idx), n_ge_2 = sum(cell_df$n_ASC_P_detected[idx] >= 2),
             n_ge_3 = sum(cell_df$n_ASC_P_detected[idx] >= 3), n_eq_6 = sum(cell_df$ASC_P_all6_detected[idx]),
             pct_ge_2 = mean(cell_df$n_ASC_P_detected[idx] >= 2) * 100, pct_ge_3 = mean(cell_df$n_ASC_P_detected[idx] >= 3) * 100,
             pct_eq_6 = mean(cell_df$ASC_P_all6_detected[idx]) * 100, stringsAsFactors = FALSE)
}))
asc_celltype <- do.call(rbind, lapply(unique(cell_df$celltype_final), function(ct) {
  idx <- cell_df$celltype_final == ct
  data.frame(celltype_final = ct, cell_number = sum(idx), n_ge_2 = sum(cell_df$n_ASC_P_detected[idx] >= 2),
             n_ge_3 = sum(cell_df$n_ASC_P_detected[idx] >= 3), n_eq_6 = sum(cell_df$ASC_P_all6_detected[idx]),
             pct_ge_2 = mean(cell_df$n_ASC_P_detected[idx] >= 2) * 100, pct_ge_3 = mean(cell_df$n_ASC_P_detected[idx] >= 3) * 100,
             pct_eq_6 = mean(cell_df$ASC_P_all6_detected[idx]) * 100, stringsAsFactors = FALSE)
}))
write_tsv_safe(asc_cluster, file.path(cfg$table_dir, "NormTotal_ASC_signature_detection_by_cluster.tsv"), cfg$overwrite)
write_tsv_safe(asc_celltype, file.path(cfg$table_dir, "NormTotal_ASC_signature_detection_by_celltype.tsv"), cfg$overwrite)

# Cluster 7 diagnostic: counts-based co-detection and available QC/doublet metadata.
cluster7 <- cell_df$cluster == "7"
cluster7_cells <- rownames(cell_df)[cluster7]
immune_genes <- intersect(marker_sets$immune, rownames(counts_mat)); endo_genes <- intersect(marker_sets$vascular_endothelial, rownames(counts_mat)); lymph_genes <- intersect(marker_sets$lymphatic_endothelial, rownames(counts_mat))
immune_detect <- colSums(as.matrix(counts_mat[immune_genes, , drop = FALSE]) > 0)
endo_detect <- colSums(as.matrix(counts_mat[endo_genes, , drop = FALSE]) > 0)
lymph_detect <- colSums(as.matrix(counts_mat[lymph_genes, , drop = FALSE]) > 0)
diag7 <- data.frame(
  metric = c("cluster7_cell_number", "cluster7_immune_marker_positive", "cluster7_endothelial_marker_positive", "cluster7_lymphatic_marker_positive", "cluster7_PTPRC_and_endothelial_same_cell", "cluster7_PTPRC_and_lymphatic_same_cell", "cluster7_mean_nCount_RNA", "cluster7_mean_nFeature_RNA", "cluster7_percent_mt_available", "cluster7_doublet_metadata_available"),
  value = c(sum(cluster7), sum(immune_detect[cluster7] > 0), sum(endo_detect[cluster7] > 0), sum(lymph_detect[cluster7] > 0), sum((counts_mat["PTPRC", cluster7_cells] > 0) & (endo_detect[cluster7] > 0), na.rm = TRUE), sum((counts_mat["PTPRC", cluster7_cells] > 0) & (lymph_detect[cluster7] > 0), na.rm = TRUE), mean(cell_df$nCount_RNA[cluster7], na.rm = TRUE), mean(cell_df$nFeature_RNA[cluster7], na.rm = TRUE), as.integer(any(grepl("percent.mt", names(meta), ignore.case = TRUE))), as.integer(any(grepl("doublet", names(meta), ignore.case = TRUE)))),
  stringsAsFactors = FALSE
)
write_tsv_safe(diag7, file.path(cfg$table_dir, "NormTotal_cluster7_diagnostic_summary.tsv"), cfg$overwrite)
cluster7_provider <- patient_cluster[patient_cluster$cluster == "7", c("patient", "cell_number", "patient_total_cells", "percent_within_patient"), drop = FALSE]
cluster7_provider$percent_of_cluster7 <- if (sum(cluster7) > 0) 100 * cluster7_provider$cell_number / sum(cluster7) else NA_real_
write_tsv_safe(cluster7_provider, file.path(cfg$table_dir, "NormTotal_cluster7_by_patient.tsv"), cfg$overwrite)

# ASC gene-level summaries by cluster and broad cell type.
asc_gene_rows <- list()
for (g in asc_present) for (cl in clusters) {
  cells <- cell_ids[as.character(meta[cell_ids, cluster_col]) == cl]
  asc_gene_rows[[length(asc_gene_rows) + 1L]] <- data.frame(group_type = "cluster", group = cl, gene = g, percent_expressed = mean(as.matrix(counts_mat[g, cells, drop = FALSE]) > 0) * 100, mean_normalized_RNA = mean(data_mat[g, cells]), stringsAsFactors = FALSE)
}
for (g in asc_present) for (ct in unique(cell_df$celltype_final)) {
  cells <- cell_df$cell[cell_df$celltype_final == ct]
  asc_gene_rows[[length(asc_gene_rows) + 1L]] <- data.frame(group_type = "celltype", group = ct, gene = g, percent_expressed = mean(as.matrix(counts_mat[g, cells, drop = FALSE]) > 0) * 100, mean_normalized_RNA = mean(data_mat[g, cells]), stringsAsFactors = FALSE)
}
asc_gene_summary <- do.call(rbind, asc_gene_rows)

# Per-gene localization table for the requested exclusion/support modules.
localization_rows <- list()
localization_sets <- c(marker_sets, list(ASC_P = asc_genes))
for (module_name in names(localization_sets)) {
  genes <- intersect(localization_sets[[module_name]], intersect(rownames(data_mat), rownames(counts_mat)))
  if (!length(genes)) next
  for (g in genes) for (cl in clusters) {
    cells <- cell_ids[as.character(meta[cell_ids, cluster_col]) == cl]
    localization_rows[[length(localization_rows) + 1L]] <- data.frame(module = module_name, gene = g, cluster = cl,
      percent_expressed_counts = 100 * mean(as.matrix(counts_mat[g, cells, drop = FALSE]) > 0),
      mean_normalized_RNA = mean(data_mat[g, cells]), cell_number = length(cells), patient_number = length(unique(meta[cells, patient_col])),
      stringsAsFactors = FALSE)
  }
}
write_tsv_safe(do.call(rbind, localization_rows), file.path(cfg$table_dir, "NormTotal_ASC_localization_by_cluster.tsv"), cfg$overwrite)

# Annotated copy with new metadata only; original RDS remains untouched.
meta$celltype_broad <- cell_ann
meta$celltype_broad_confidence <- ann$confidence[match(cell_ann, ann$celltype_final)]
meta$celltype_broad_evidence <- ann$main_evidence[match(cell_ann, ann$celltype_final)]
meta$ASC_P_n_detected <- as.integer(asc_n[match(rownames(meta), cell_ids)])
meta$ASC_P_all6_detected <- asc_all6[match(rownames(meta), cell_ids)]
meta$ASC_P_candidate_stromal <- cell_df$ASC_P_candidate_stromal[match(rownames(meta), cell_df$cell)]
obj@meta.data <- meta
annotated_path <- file.path(cfg$updated_object_dir, "SeuratObject_NormTotal_annotated.rds")
save_rds_safe(obj, annotated_path, cfg$overwrite, compress = TRUE)

# Plot helpers.
emb <- safe_embedding(obj, "tsne")
emb <- emb[intersect(rownames(emb), cell_ids), , drop = FALSE]
plot_groups <- function(groups, title, path_base, palette = NULL) {
  d <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], group = as.character(groups[match(rownames(emb), names(groups))]), stringsAsFactors = FALSE)
  if (anyNA(d$group)) d$group[is.na(d$group)] <- "NA"
  p <- ggplot2::ggplot(d, ggplot2::aes(dim1, dim2, colour = group)) + ggplot2::geom_point(size = 0.12, alpha = 0.65) + ggplot2::coord_equal() + ggplot2::theme_classic() + ggplot2::labs(title = title, x = "tSNE-1", y = "tSNE-2", colour = NULL)
  pdf <- paste0(path_base, ".pdf"); png <- paste0(path_base, ".png"); assert_new_path(pdf, cfg$overwrite); assert_new_path(png, cfg$overwrite)
  ggplot2::ggsave(pdf, p, width = 9, height = 7); ggplot2::ggsave(png, p, width = 9, height = 7, dpi = 300); c(pdf, png)
}
plot_groups(setNames(as.character(meta[[cluster_col]]), rownames(meta)), "NormalTotal author clusters", file.path(cfg$figure_dir, "audit", "NormTotal_original_clusters"))
plot_groups(setNames(meta$celltype_broad, rownames(meta)), "NormalTotal final broad celltype", file.path(cfg$figure_dir, "identity_boundary", "NormTotal_final_celltype_DimPlot"))

dot <- module_scores
dot$module <- factor(dot$module, levels = names(marker_sets))
p_dot <- ggplot2::ggplot(dot, ggplot2::aes(module, cluster, size = percent_expressed, colour = mean_normalized_RNA)) + ggplot2::geom_point() + ggplot2::scale_size(range = c(0.5, 8)) + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + ggplot2::labs(title = "NormalTotal lineage marker DotPlot", x = NULL, y = "Author cluster", colour = "Mean RNA", size = "% detected")
dot_pdf <- file.path(cfg$figure_dir, "identity_boundary", "NormTotal_lineage_marker_DotPlot.pdf"); dot_png <- file.path(cfg$figure_dir, "identity_boundary", "NormTotal_lineage_marker_DotPlot.png"); assert_new_path(dot_pdf, cfg$overwrite); assert_new_path(dot_png, cfg$overwrite); ggplot2::ggsave(dot_pdf, p_dot, width = 10, height = 7); ggplot2::ggsave(dot_png, p_dot, width = 10, height = 7, dpi = 300)

make_gene_dot <- function(genes, path_base, title) {
  genes <- intersect(genes, intersect(rownames(data_mat), rownames(counts_mat)))
  if (!length(genes)) {
    write_lines_safe(c(paste0("# ", title), "", "No requested marker genes were present in both RNA data and counts layers."), paste0(path_base, ".md"), cfg$overwrite)
    return(invisible(NULL))
  }
  rows <- list()
  for (g in genes) for (cl in clusters) {
    cells <- cell_ids[as.character(meta[cell_ids, cluster_col]) == cl]
    rows[[length(rows) + 1L]] <- data.frame(group = cl, gene = g,
      percent_expressed = mean(as.matrix(counts_mat[g, cells, drop = FALSE]) > 0) * 100,
      mean_normalized_RNA = mean(data_mat[g, cells]), stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows)
  p <- ggplot2::ggplot(d, ggplot2::aes(gene, group, size = percent_expressed, colour = mean_normalized_RNA)) + ggplot2::geom_point() + ggplot2::scale_size(range = c(0.5, 8)) + ggplot2::theme_classic() + ggplot2::labs(title = title, x = NULL, y = "Author cluster", colour = "Mean RNA", size = "% detected")
  pdf <- paste0(path_base, ".pdf"); png <- paste0(path_base, ".png"); assert_new_path(pdf, cfg$overwrite); assert_new_path(png, cfg$overwrite); ggplot2::ggsave(pdf, p, width = 10, height = 6); ggplot2::ggsave(png, p, width = 10, height = 6, dpi = 300)
}
make_gene_dot(c("KIT", "ELF5", "TNFRSF11A", "ALDH1A3", "LTF", "KRT19", "ESR1", "PGR", "FOXA1", "GATA3", "TFF1", "AGR2", "KRT8", "KRT18"), file.path(cfg$figure_dir, "identity_boundary", "NormTotal_luminal_0_vs_3_validation"), "NormalTotal luminal validation")
make_gene_dot(c("KRT5", "KRT14", "KRT17", "TP63", "MYLK", "DST", "EPCAM", "RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "MYH11"), file.path(cfg$figure_dir, "identity_boundary", "NormTotal_cluster2_validation"), "NormalTotal cluster 2 basal/myoepithelial versus mural")
make_gene_dot(c("PDPN", "CCL21", "PROX1", "FLT4", "LYVE1", "PECAM1", "VWF", "CLDN5", "KDR", "ESAM", "EMCN", "PLVAP", "PTPRC", "CD3D", "NKG7", "CD79A", "MS4A1", "LST1", "TYROBP", "CD68", "CD14"), file.path(cfg$figure_dir, "identity_boundary", "NormTotal_cluster6_7_validation"), "NormalTotal cluster 6/7 diagnostic")

# ASC signature plots and candidate map.
asc_dot <- asc_gene_summary[asc_gene_summary$group_type == "cluster", , drop = FALSE]
p_asc_dot <- ggplot2::ggplot(asc_dot, ggplot2::aes(gene, group, size = percent_expressed, colour = mean_normalized_RNA)) + ggplot2::geom_point() + ggplot2::scale_size(range = c(0.5, 8)) + ggplot2::theme_classic() + ggplot2::labs(title = "NormalTotal ASC_P six-gene detection", x = NULL, y = "Author cluster", colour = "Mean RNA", size = "% detected")
adp <- file.path(cfg$figure_dir, "global_localization", "NormTotal_ASC_signature_DotPlot.pdf"); adn <- file.path(cfg$figure_dir, "global_localization", "NormTotal_ASC_signature_DotPlot.png"); assert_new_path(adp, cfg$overwrite); assert_new_path(adn, cfg$overwrite); ggplot2::ggsave(adp, p_asc_dot, width = 10, height = 6); ggplot2::ggsave(adn, p_asc_dot, width = 10, height = 6, dpi = 300)
feature_rows <- list(); for (g in asc_present) feature_rows[[length(feature_rows) + 1L]] <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], value = as.numeric(counts_mat[g, rownames(emb)] > 0), gene = g)
feature_df <- do.call(rbind, feature_rows)
p_feat <- ggplot2::ggplot(feature_df, ggplot2::aes(dim1, dim2, colour = value)) + ggplot2::geom_point(size = 0.12) + ggplot2::scale_colour_gradient(low = "grey90", high = "red") + ggplot2::coord_equal() + ggplot2::facet_wrap(~gene, ncol = 3) + ggplot2::theme_classic() + ggplot2::labs(title = "NormalTotal ASC_P genes by tSNE", x = "tSNE-1", y = "tSNE-2", colour = "Detected")
afp <- file.path(cfg$figure_dir, "global_localization", "NormTotal_ASC_signature_FeaturePlot.pdf"); afn <- file.path(cfg$figure_dir, "global_localization", "NormTotal_ASC_signature_FeaturePlot.png"); assert_new_path(afp, cfg$overwrite); assert_new_path(afn, cfg$overwrite); ggplot2::ggsave(afp, p_feat, width = 10, height = 8); ggplot2::ggsave(afn, p_feat, width = 10, height = 8, dpi = 300)
p_hist <- ggplot2::ggplot(cell_df, ggplot2::aes(n_ASC_P_detected, fill = celltype_final)) + ggplot2::geom_histogram(binwidth = 1, boundary = -0.5, position = "dodge") + ggplot2::theme_classic() + ggplot2::labs(title = "NormalTotal ASC_P genes detected per cell", x = "Number of ASC_P genes detected", y = "Cells", fill = NULL)
hp <- file.path(cfg$figure_dir, "global_localization", "NormTotal_ASC_genes_detected_per_cell.pdf"); hn <- file.path(cfg$figure_dir, "global_localization", "NormTotal_ASC_genes_detected_per_cell.png"); assert_new_path(hp, cfg$overwrite); assert_new_path(hn, cfg$overwrite); ggplot2::ggsave(hp, p_hist, width = 10, height = 7); ggplot2::ggsave(hn, p_hist, width = 10, height = 7, dpi = 300)
candidate_group <- ifelse(cell_df$ASC_P_candidate_stromal, "ASC_P-like stromal candidate", "Other")
plot_groups(setNames(candidate_group, cell_df$cell), "NormalTotal ASC_P-like stromal candidate map", file.path(cfg$figure_dir, "candidates", "NormTotal_ASC_candidate_stromal_map"))

# Object audit and report.
audit <- data.frame(object_id = "NormTotal", input_path = normalizePath(input_path, winslash = "/"), input_sha256 = sha256_file(input_path), cells = nrow(meta), genes_RNA = nrow(get_rna_matrix(obj, "data")), clusters = length(clusters), cluster_column = cluster_col, patient_column = patient_col, sample_column = sample_col, reductions = paste(safe_reduction_names(obj), collapse = ";"), metadata_columns = paste(colnames(meta), collapse = ";"), percent_mt_available = any(grepl("percent.mt", colnames(meta), ignore.case = TRUE)), doublet_metadata_available = any(grepl("doublet", colnames(meta), ignore.case = TRUE)), stringsAsFactors = FALSE)
write_tsv_safe(audit, file.path(cfg$table_dir, "NormTotal_object_audit.tsv"), cfg$overwrite)
fib_asc <- asc_celltype[asc_celltype$celltype_final == "Fibroblast", , drop = FALSE]
mural_asc <- asc_celltype[asc_celltype$celltype_final == "Mural/pericyte-VSMC", , drop = FALSE]
cluster4_asc6 <- asc_cluster$n_eq_6[asc_cluster$cluster == "4"]
cluster7_ptprc_endo <- diag7$value[diag7$metric == "cluster7_PTPRC_and_endothelial_same_cell"]
cluster7_ptprc_lymph <- diag7$value[diag7$metric == "cluster7_PTPRC_and_lymphatic_same_cell"]
candidate_total <- sum(cell_df$ASC_P_candidate_stromal)
candidate_patients <- length(unique(cell_df$patient[cell_df$ASC_P_candidate_stromal]))
report_lines <- c(
  "# NORMTOTAL cluster annotation and ASC_P preliminary stromal screen", "",
  "## Data and method", "",
  paste0("- Input: `", normalizePath(input_path, winslash = "/"), "`"), paste0("- Input SHA256: `", sha256_file(input_path), "`"),
  paste0("- Cells: ", nrow(meta), "; author clusters: ", paste(clusters, collapse = ", "), "; patients: ", length(unique(meta[[patient_col]]))),
  "- Existing author reductions and clusters were reused. No QC, normalization, integration, reduction, or clustering was rerun.",
  "- New annotations are stored only in an annotated copy; the author RDS is unchanged.", "",
  "## Cluster conclusions", "",
  paste0("- Cluster 0: ", ann$celltype_final[ann$cluster == "0"], "; luminal progenitor versus mature luminal evidence remains provisional."),
  "- Cluster 1: Fibroblast; the strongest normal stromal/ASC_P-related compartment.",
  paste0("- Cluster 2: ", ann$celltype_final[ann$cluster == "2"], "; ACTA2/TAGLN requires basal versus mural comparison."),
  paste0("- Cluster 3: ", ann$celltype_final[ann$cluster == "3"], "; luminal distinction remains provisional."),
  "- Cluster 4: Vascular endothelial; ENG/MCAM do not override the canonical endothelial program.",
  "- Cluster 5: Mural/pericyte-VSMC; a key ASC_P-related comparison compartment, not automatically CAF.",
  "- Cluster 6: Immune; PTPRC-positive and supported across all patients.",
  "- Cluster 7: Unresolved PTPRC+ENG cluster pending cell-level co-detection/doublet evidence.", "",
  paste0("- Cluster 2 basal/myoepithelial module detection was ", round(get_score("2", "basal_myoepithelial"), 2), "% versus ", round(get_score("2", "mural_perivascular"), 2), "% for mural/perivascular; the basal assignment is therefore supported but remains marked moderate because ACTA2/TAGLN are shared."),
  paste0("- Cluster 7 counts-based diagnostic: ", cluster7_ptprc_endo, " cells co-detected PTPRC and a vascular marker, and ", cluster7_ptprc_lymph, " co-detected PTPRC and a lymphatic marker; percent.mt and doublet metadata were unavailable."),
  "## ASC_P interpretation", "",
  "ASC_P candidates are restricted to non-immune, non-epithelial, non-endothelial stromal/mural backgrounds and are reported as transcriptional candidates only.",
  "The six-gene detection tables use raw RNA counts > 0 and report n>=2, n>=3, and n=6; they do not claim functional equivalence to mouse or human CD105+/CD90+ ASC.", "",
  paste0("- Fibroblast ASC_P screen: ", fib_asc$n_ge_2, " cells with >=2 genes, ", fib_asc$n_ge_3, " with >=3, ", fib_asc$n_eq_6, " with all 6; mural/pericyte-VSMC: ", mural_asc$n_ge_2, " / ", mural_asc$n_ge_3, " / ", mural_asc$n_eq_6, "."),
  paste0("- Candidate stromal cells (>=2 raw-count genes after immune/epithelial/endothelial exclusion): ", candidate_total, " cells across ", candidate_patients, " patients; endothelial cluster 4 had ", cluster4_asc6, " six-of-six cells but was excluded as an endothelial compartment."),
  "## Secondary clustering decision", "",
  "No secondary stromal/mural clustering was performed in this run. The existing evidence supports first completing broad annotation and ASC_P localization; a second clustering decision should follow inspection of cluster 1/5 heterogeneity and patient-level support.", "",
  "## Outputs", "",
  paste0("- Annotated copy: `", normalizePath(annotated_path, winslash = "/"), "`"),
  "- Tables: NormTotal_object_audit.tsv; NormTotal_cluster_cell_numbers.tsv; NormTotal_sample_by_cluster.tsv; NormTotal_patient_by_cluster.tsv; NormTotal_marker_availability.tsv; NormTotal_annotation_evidence.tsv; NormTotal_final_cluster_annotation.tsv; NormTotal_cluster7_diagnostic_summary.tsv; NormTotal_cluster7_by_patient.tsv; NormTotal_ASC_signature_detection_by_cluster.tsv; NormTotal_ASC_signature_detection_by_celltype.tsv; NormTotal_ASC_candidate_cells.tsv; NormTotal_ASC_candidate_by_patient.tsv; NormTotal_ASC_localization_by_cluster.tsv.",
  "- Figures: original clusters, lineage DotPlot, luminal validation, cluster 2 validation, cluster 6/7 validation, final celltype map, ASC_P DotPlot/FeaturePlot, per-cell detection, and candidate stromal map."
)
report_out <- file.path(cfg$report_dir, "NORMTOTAL_CLUSTER_ANNOTATION_REPORT.md")
write_lines_safe(report_lines, report_out, cfg$overwrite)

end_step(cfg, state, c(file.path(cfg$table_dir, c("NormTotal_object_audit.tsv", "NormTotal_cluster_cell_numbers.tsv", "NormTotal_sample_by_cluster.tsv", "NormTotal_patient_by_cluster.tsv", "NormTotal_marker_availability.tsv", "NormTotal_annotation_evidence.tsv", "NormTotal_final_cluster_annotation.tsv", "NormTotal_cluster7_diagnostic_summary.tsv", "NormTotal_cluster7_by_patient.tsv", "NormTotal_ASC_signature_detection_by_cluster.tsv", "NormTotal_ASC_signature_detection_by_celltype.tsv", "NormTotal_ASC_candidate_cells.tsv", "NormTotal_ASC_candidate_by_patient.tsv", "NormTotal_ASC_localization_by_cluster.tsv")), annotated_path, report_out), status = "SUCCESS", notes = c("NormalTotal annotation and ASC_P screen completed without modifying author RDS."))

cfg <- load_config(); root <- project_root(); attempt <- Sys.getenv("ASC_P_ATTEMPT_DIR")
assert_packages(c("Matrix", "dplyr", "ggplot2", "patchwork", "data.table"))
set.seed(cfg$seed)

# --- Load all samples using data.table::fread for matrix market files ---
counts_list <- list()
meta_list <- list()
mt_genes <- NULL

for (i in seq_len(nrow(cfg$samples))) {
  s <- cfg$samples[i, ]
  d <- file.path(cfg$input_root, s$folder, "filtered_feature_bc_matrix")

  # Read features
  feat <- data.table::fread(file.path(d, "features.tsv.gz"), header = FALSE,
                             col.names = c("gene_id", "gene_name", "feature_type"))
  gene_names <- feat$gene_name
  if (is.null(mt_genes)) mt_genes <- grep("^MT-", gene_names, value = TRUE)

  # Read barcodes
  bcs <- data.table::fread(file.path(d, "barcodes.tsv.gz"), header = FALSE,
                            col.names = "barcode")
  new_names <- paste0(s$sample_id, "_", bcs$barcode)

  # Read matrix market
  mat <- Matrix::readMM(file.path(d, "matrix.mtx.gz"))
  rownames(mat) <- gene_names
  colnames(mat) <- new_names
  mat <- as(mat, "dgCMatrix")

  counts_list[[i]] <- mat

  # Build metadata
  nc <- length(new_names)
  meta_list[[i]] <- data.frame(
    orig.ident   = rep(s$sample_id, nc),
    sample_id    = rep(s$sample_id, nc),
    patient_id   = rep(s$patient_id, nc),
    group        = rep(s$group, nc),
    run_accession = rep(s$run_accession, nc),
    stringsAsFactors = FALSE
  )
  rownames(meta_list[[i]]) <- new_names
}

# --- Merge matrices on common genes ---
common_genes <- Reduce(intersect, lapply(counts_list, rownames))
all_counts <- do.call(cbind, lapply(counts_list, function(m) m[common_genes, , drop = FALSE]))
all_meta <- do.call(rbind, meta_list)
stopifnot(identical(colnames(all_counts), rownames(all_meta)))
stopifnot(ncol(all_counts) == 63450L)
stopifnot(!anyDuplicated(colnames(all_counts)))

# --- Compute QC metrics manually ---
all_meta$nCount_RNA <- Matrix::colSums(all_counts)
all_meta$nFeature_RNA <- Matrix::colSums(all_counts > 0)

# percent.mt
mt_idx <- which(common_genes %in% mt_genes)
if (length(mt_idx) > 0) {
  mt_counts <- Matrix::colSums(all_counts[mt_idx, , drop = FALSE])
  all_meta$percent.mt <- mt_counts / all_meta$nCount_RNA * 100
} else {
  all_meta$percent.mt <- 0
}

# --- Output cell count table ---
counts_table <- as.data.frame(table(all_meta$sample_id), stringsAsFactors = FALSE)
names(counts_table) <- c("sample_id", "cells")
write.csv(counts_table, file.path(attempt, "tables", "cells_before_qc.csv"), row.names = FALSE)

# --- QC violin plots ---
qc_long <- data.frame(
  sample_id = all_meta$sample_id,
  nFeature_RNA = all_meta$nFeature_RNA,
  nCount_RNA   = all_meta$nCount_RNA,
  percent_mt   = all_meta$percent.mt,
  stringsAsFactors = FALSE
)

p1a <- ggplot2::ggplot(qc_long, ggplot2::aes(x = sample_id, y = nFeature_RNA, fill = sample_id)) +
  ggplot2::geom_violin(scale = "width") +
  ggplot2::theme_bw() + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(x = "", y = "nFeature_RNA")

p1b <- ggplot2::ggplot(qc_long, ggplot2::aes(x = sample_id, y = nCount_RNA, fill = sample_id)) +
  ggplot2::geom_violin(scale = "width") +
  ggplot2::theme_bw() + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(x = "", y = "nCount_RNA")

p1c <- ggplot2::ggplot(qc_long, ggplot2::aes(x = sample_id, y = percent_mt, fill = sample_id)) +
  ggplot2::geom_violin(scale = "width") +
  ggplot2::theme_bw() + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(x = "", y = "percent.mt")

p1 <- patchwork::wrap_plots(p1a, p1b, p1c, ncol = 1)
save_plot_both(p1, file.path(attempt, "figures", "QC_violin_before_filtering"), 11, 10)

# --- QC scatter plot ---
p2 <- ggplot2::ggplot(all_meta, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, colour = sample_id)) +
  ggplot2::geom_point(size = 0.15, alpha = 0.25) +
  ggplot2::scale_x_log10() + ggplot2::theme_bw()
save_plot_both(p2, file.path(attempt, "figures", "QC_scatter_nCount_nFeature_before"), 8, 6)

# --- Save checkpoint: list with counts matrix + metadata ---
checkpoint_obj <- list(counts = all_counts, meta = all_meta, mt_genes = mt_genes)
checkpoint <- file.path(attempt, "objects", "PRJNA932038_preQC.rds")
saveRDS(checkpoint_obj, checkpoint)

write_waiting_review(attempt, 1, c("## Required checks", "",
  "- Confirm all six samples and unique barcodes (63,450 nuclei).",
  "- Inspect per-sample nFeature/nCount/percent.mt distributions.",
  "- Confirm no sample path or metadata swap.",
  "- NOTE: Step 01 uses base Matrix + ggplot2 instead of Seurat merge to avoid S4SXP bug."))
write_step_manifest(1, attempt, checkpoint, file.path(root, "scripts", "step_01.R"))

cfg <- load_config(); root <- project_root(); attempt <- Sys.getenv("ASC_P_ATTEMPT_DIR")
assert_packages(c("Matrix", "dplyr", "ggplot2", "patchwork", "SummarizedExperiment"))
set.seed(cfg$seed)

# Step 01 checkpoint is list(counts=dgCMatrix, meta=data.frame, mt_genes=character)
prev <- read_previous_checkpoint(2)

# ------ Per-sample QC thresholds (on metadata from Step 01) ------
meta <- prev$meta
meta$qc_primary <- FALSE
meta$qc_article <- FALSE
threshold_rows <- list()
for (sid in unique(meta$sample_id)) {
  ix <- which(meta$sample_id == sid)
  nf <- meta$nFeature_RNA[ix]; nc <- meta$nCount_RNA[ix]
  nf_upper <- median(nf) + cfg$qc$mad_multiplier * mad(nf)
  nc_upper <- median(nc) + cfg$qc$mad_multiplier * mad(nc)
  article_nf <- mean(nf) + c(-1, 1) * cfg$qc$article_sd_multiplier * sd(nf)
  article_nc <- mean(nc) + c(-1, 1) * cfg$qc$article_sd_multiplier * sd(nc)
  meta$qc_primary[ix] <- nf >= cfg$qc$min_features & nf <= nf_upper & nc > 0 & nc <= nc_upper & meta$percent.mt[ix] < cfg$qc$max_percent_mt
  meta$qc_article[ix] <- nf >= article_nf[1] & nf <= article_nf[2] & nc >= article_nc[1] & nc <= article_nc[2] & meta$percent.mt[ix] < cfg$qc$max_percent_mt
  threshold_rows[[sid]] <- data.frame(sample_id = sid, nFeature_upper = nf_upper, nCount_upper = nc_upper,
    article_nFeature_low = article_nf[1], article_nFeature_high = article_nf[2],
    article_nCount_low = article_nc[1], article_nCount_high = article_nc[2])
}

# ------ scDblFinder per-sample via dgCMatrix SCE (bypasses Seurat entirely) ------
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scDblFinder))
set.seed(cfg$seed)
doublet_class <- rep("singlet", ncol(prev$counts))
doublet_score <- rep(0, ncol(prev$counts))
names(doublet_class) <- names(doublet_score) <- colnames(prev$counts)
for (sid in unique(meta$sample_id)) {
  cells <- colnames(prev$counts)[meta$sample_id == sid]
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = prev$counts[, cells, drop = FALSE]))
  sce <- scDblFinder::scDblFinder(sce)
  doublet_class[cells] <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
  doublet_score[cells] <- as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score)
}
meta$scDblFinder.class <- doublet_class
meta$scDblFinder.score <- doublet_score
meta$qc_keep <- meta$qc_primary & meta$scDblFinder.class == "singlet"

# ------ QC summary table ------
summary <- meta |>
  tibble::rownames_to_column("cell") |>
  dplyr::group_by(sample_id, patient_id, group) |>
  dplyr::summarise(before = dplyr::n(), primary_pass = sum(qc_primary), article_pass = sum(qc_article),
                   singlets = sum(scDblFinder.class == "singlet"), final_keep = sum(qc_keep), .groups = "drop")
write.csv(summary, file.path(attempt, "tables", "qc_retention_by_sample.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(threshold_rows), file.path(attempt, "tables", "qc_thresholds_by_sample.csv"), row.names = FALSE)

# ------ QC violin after filtering (ggplot2, no Seurat) ------
qc_meta <- meta[meta$qc_keep, , drop = FALSE]
qc_long <- data.frame(
  sample_id = qc_meta$sample_id,
  nFeature_RNA = qc_meta$nFeature_RNA,
  nCount_RNA   = qc_meta$nCount_RNA,
  percent_mt   = qc_meta$percent.mt,
  stringsAsFactors = FALSE
)
p1a <- ggplot2::ggplot(qc_long, ggplot2::aes(x = sample_id, y = nFeature_RNA, fill = sample_id)) +
  ggplot2::geom_violin(scale = "width") + ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(x = "", y = "nFeature_RNA", title = "After primary QC + scDblFinder")
p1b <- ggplot2::ggplot(qc_long, ggplot2::aes(x = sample_id, y = nCount_RNA, fill = sample_id)) +
  ggplot2::geom_violin(scale = "width") + ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(x = "", y = "nCount_RNA")
p1c <- ggplot2::ggplot(qc_long, ggplot2::aes(x = sample_id, y = percent_mt, fill = sample_id)) +
  ggplot2::geom_violin(scale = "width") + ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(x = "", y = "percent.mt")
p1 <- patchwork::wrap_plots(p1a, p1b, p1c, ncol = 1)
save_plot_both(p1, file.path(attempt, "figures", "QC_violin_after_filtering"), 11, 10)

# ------ S81_N QC review ------
writeLines(c("# S81_N QC review", "", paste(capture.output(print(summary[summary$sample_id == "S81_N", ])), collapse = "\n"),
             "", "Cell Ranger reported 42% reads in cells. Do not remove S81_N alone; Step 08 removes the entire S81 pair for sensitivity."),
           file.path(attempt, "reports", "S81_N_QC_REVIEW.md"))

# ------ Checkpoint: filtered counts matrix + metadata ------
qc_counts <- prev$counts[, meta$qc_keep, drop = FALSE]
qc_meta <- meta[meta$qc_keep, , drop = FALSE]
checkpoint_obj <- list(counts = qc_counts, meta = qc_meta, thresholds = threshold_rows, mt_genes = prev$mt_genes)
checkpoint <- file.path(attempt, "objects", "PRJNA932038_QC_singlets.rds")
saveRDS(checkpoint_obj, checkpoint)

write_waiting_review(attempt, 2, c("## Required checks", "", "- Review sample-specific thresholds and retention.",
  "- Review scDblFinder rates.", "- Review S81_N separately.",
  "- Approve primary QC or request a new attempt; do not edit this attempt in place.",
  "- NOTE: All Seurat calls bypassed; scDblFinder run per-sample via matrix->SCE; plots via ggplot2."))
write_step_manifest(2, attempt, checkpoint, file.path(root, "stages", "02_filter_qc_and_doublets.R"))


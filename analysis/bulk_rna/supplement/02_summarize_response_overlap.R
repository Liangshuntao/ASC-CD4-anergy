# Explicit file interface. No automatic installation or data download.
.cli <- commandArgs(trailingOnly = TRUE)
.value <- function(k) {
  i <- which(.cli == k)
  if (length(i) != 1L || i == length(.cli)) stop("Required argument: ", k)
  .cli[i + 1L]
}
if ("--help" %in% .cli) {
  cat("Usage: Rscript SCRIPT --input EXPRESSION.csv --out OUTPUT_DIRECTORY --did-table PRIMARY_DID.csv [--dry-run]\n")
  quit(status = 0L)
}
.input <- normalizePath(.value("--input"), winslash = "/", mustWork = TRUE)
.output <- normalizePath(.value("--out"), winslash = "/", mustWork = FALSE)
if (identical(tolower(dirname(.input)), tolower(.output))) stop("Use separate input and output directories")
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(openxlsx)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
project_root <- .output
base_output_dir <- file.path(.output, "within_subset")
output_dir <- file.path(.output, "submission")
input_file <- .input
p_result_file <- file.path(base_output_dir, "ASC_P_Post_vs_Pre_limma_full.tsv")
n_result_file <- file.path(base_output_dir, "ASC_N_Post_vs_Pre_limma_full.tsv")
did_file <- normalizePath(.value("--did-table"), winslash = "/", mustWork = TRUE)

required_inputs <- c(input_file, p_result_file, n_result_file, did_file)
if (any(!file.exists(required_inputs))) {
  stop("Missing required inputs: ", paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
}

fdr_cutoff <- 0.05
lfc_cutoff <- 0.5

raw <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
sample_cols <- grep("^Normalized read count ", names(raw), value = TRUE)
sample_names <- sub("^Normalized read count ", "", sample_cols)
expected_samples <- c(
  paste0("ASC_N_", 1:3), paste0("ASC_P_", 1:3),
  paste0("after_ASC_N_", 1:3), paste0("after_ASC_P_", 1:3)
)
if (!setequal(sample_names, expected_samples) || length(sample_names) != 12L) {
  stop("Expected 12 ASC normalized-count sample columns were not found exactly")
}

normalized_counts <- as.matrix(raw[, sample_cols, drop = FALSE])
storage.mode(normalized_counts) <- "numeric"
if (anyNA(normalized_counts) || any(!is.finite(normalized_counts)) || any(normalized_counts < 0)) {
  stop("Normalized expression matrix contains invalid values")
}
if (!any(abs(normalized_counts - round(normalized_counts)) > 1e-8)) {
  stop("Input unexpectedly appears integer-like; re-audit method selection before proceeding")
}

gene_symbol <- trimws(gsub("'", "", raw[["Gene Symbol"]], fixed = TRUE))
gene_key <- ifelse(is.na(gene_symbol) | gene_symbol == "", as.character(raw[["Gene ID"]]), gene_symbol)
gene_id <- as.character(raw[["Gene ID"]])

# The current project has no duplicate cleaned symbols. The aggregation below
# keeps reruns deterministic if duplicates are introduced later.
count_df <- data.frame(Gene = gene_key, Gene_ID = gene_id, normalized_counts, check.names = FALSE)
gene_ids <- aggregate(Gene_ID ~ Gene, data = count_df, FUN = function(x) paste(unique(x), collapse = ";"))
count_gene <- aggregate(. ~ Gene, data = count_df[, c("Gene", sample_cols), drop = FALSE], FUN = median, na.rm = TRUE)
count_mat <- as.matrix(count_gene[, -1, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- count_gene$Gene
colnames(count_mat) <- sample_names

cols <- list(
  P_Pre = paste0("ASC_P_", 1:3),
  P_Post = paste0("after_ASC_P_", 1:3),
  N_Pre = paste0("ASC_N_", 1:3),
  N_Post = paste0("after_ASC_N_", 1:3)
)

expression_summary <- data.frame(
  Gene = rownames(count_mat),
  Gene_ID = gene_ids$Gene_ID[match(rownames(count_mat), gene_ids$Gene)],
  ASC_P_mean_normalized_count_Pre = rowMeans(count_mat[, cols$P_Pre, drop = FALSE]),
  ASC_P_mean_normalized_count_Post = rowMeans(count_mat[, cols$P_Post, drop = FALSE]),
  ASC_N_mean_normalized_count_Pre = rowMeans(count_mat[, cols$N_Pre, drop = FALSE]),
  ASC_N_mean_normalized_count_Post = rowMeans(count_mat[, cols$N_Post, drop = FALSE]),
  stringsAsFactors = FALSE
)
expression_summary$ASC_P_mean_normalized_count_contrast_all <- rowMeans(count_mat[, c(cols$P_Pre, cols$P_Post), drop = FALSE])
expression_summary$ASC_N_mean_normalized_count_contrast_all <- rowMeans(count_mat[, c(cols$N_Pre, cols$N_Post), drop = FALSE])

read_result <- function(path, label) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("Gene", "logFC", "AveExpr", "P.Value", "FDR", "method", "primary_design")
  if (!all(required %in% names(x))) stop(label, " result table is missing required columns")
  if (nrow(x) != nrow(count_mat)) stop(label, " result row count does not match expression universe")
  if (!identical(unique(x$method), "limma")) stop(label, " result method is not uniformly limma")
  x$direction_recomputed <- ifelse(
    x$FDR < fdr_cutoff & x$logFC >= lfc_cutoff, "Up",
    ifelse(x$FDR < fdr_cutoff & x$logFC <= -lfc_cutoff, "Down", "NS")
  )
  if ("direction" %in% names(x) && any(x$direction != x$direction_recomputed)) {
    stop(label, " stored direction differs from threshold-recomputed direction")
  }
  x
}

res_p <- read_result(p_result_file, "ASC_P")
res_n <- read_result(n_result_file, "ASC_N")

make_submission_table <- function(res, subset_name) {
  idx <- match(res$Gene, expression_summary$Gene)
  if (anyNA(idx)) stop("Expression summaries failed to match ", subset_name, " results")
  if (subset_name == "ASC_P") {
    mean_pre <- expression_summary$ASC_P_mean_normalized_count_Pre[idx]
    mean_post <- expression_summary$ASC_P_mean_normalized_count_Post[idx]
    mean_all <- expression_summary$ASC_P_mean_normalized_count_contrast_all[idx]
  } else {
    mean_pre <- expression_summary$ASC_N_mean_normalized_count_Pre[idx]
    mean_post <- expression_summary$ASC_N_mean_normalized_count_Post[idx]
    mean_all <- expression_summary$ASC_N_mean_normalized_count_contrast_all[idx]
  }
  data.frame(
    Gene = res$Gene,
    Gene_ID = expression_summary$Gene_ID[idx],
    log2FC = res$logFC,
    p_value = res$P.Value,
    FDR = res$FDR,
    mean_normalized_count_contrast_all = mean_all,
    mean_normalized_count_Pre = mean_pre,
    mean_normalized_count_Post = mean_post,
    AveExpr_log2_all12 = res$AveExpr,
    direction = res$direction_recomputed,
    significant = res$direction_recomputed != "NS",
    contrast = paste(subset_name, "Post vs Pre"),
    method = "limma",
    input_scale = "log2(normalized read count + 1)",
    design = "unpaired subtype-by-timepoint model",
    stringsAsFactors = FALSE
  )
}

table_p <- make_submission_table(res_p, "ASC_P")
table_n <- make_submission_table(res_n, "ASC_N")
sig_p <- table_p[table_p$significant, , drop = FALSE]
sig_n <- table_n[table_n$significant, , drop = FALSE]

paired <- merge(
  table_p[, c("Gene", "Gene_ID", "log2FC", "FDR", "direction")],
  table_n[, c("Gene", "log2FC", "FDR", "direction")],
  by = "Gene", suffixes = c("_ASC_P", "_ASC_N"), sort = FALSE
)

paired$response_class <- with(paired, ifelse(
  direction_ASC_P == "Up" & direction_ASC_N == "Up", "shared_up",
  ifelse(direction_ASC_P == "Down" & direction_ASC_N == "Down", "shared_down",
  ifelse(direction_ASC_P == "Up" & direction_ASC_N == "NS", "ASC_P_only_up",
  ifelse(direction_ASC_P == "Down" & direction_ASC_N == "NS", "ASC_P_only_down",
  ifelse(direction_ASC_N == "Up" & direction_ASC_P == "NS", "ASC_N_only_up",
  ifelse(direction_ASC_N == "Down" & direction_ASC_P == "NS", "ASC_N_only_down",
  ifelse(direction_ASC_P == "Up" & direction_ASC_N == "Down", "discordant_P_up_N_down",
  ifelse(direction_ASC_P == "Down" & direction_ASC_N == "Up", "discordant_P_down_N_up",
         "nonresponsive")))))))))

set_names <- c(
  "shared_up", "ASC_P_only_up", "ASC_N_only_up",
  "shared_down", "ASC_P_only_down", "ASC_N_only_down",
  "discordant_P_up_N_down", "discordant_P_down_N_up"
)
gene_sets <- setNames(lapply(set_names, function(z) paired$Gene[paired$response_class == z]), set_names)

p_up <- paired$Gene[paired$direction_ASC_P == "Up"]
n_up <- paired$Gene[paired$direction_ASC_N == "Up"]
p_down <- paired$Gene[paired$direction_ASC_P == "Down"]
n_down <- paired$Gene[paired$direction_ASC_N == "Down"]

summary_counts <- data.frame(
  category = c(
    "ASC_P_significant_total", "ASC_N_significant_total",
    "ASC_P_up", "ASC_N_up", "shared_up", "ASC_P_only_up_strict_N_NS", "ASC_N_only_up_strict_P_NS",
    "ASC_P_down", "ASC_N_down", "shared_down", "ASC_P_only_down_strict_N_NS", "ASC_N_only_down_strict_P_NS",
    "discordant_P_up_N_down", "discordant_P_down_N_up",
    "Venn_ASC_P_unique_up_includes_discordant", "Venn_ASC_N_unique_up_includes_discordant",
    "Venn_ASC_P_unique_down_includes_discordant", "Venn_ASC_N_unique_down_includes_discordant"
  ),
  count = c(
    nrow(sig_p), nrow(sig_n),
    length(p_up), length(n_up), length(gene_sets$shared_up), length(gene_sets$ASC_P_only_up), length(gene_sets$ASC_N_only_up),
    length(p_down), length(n_down), length(gene_sets$shared_down), length(gene_sets$ASC_P_only_down), length(gene_sets$ASC_N_only_down),
    length(gene_sets$discordant_P_up_N_down), length(gene_sets$discordant_P_down_N_up),
    length(setdiff(p_up, n_up)), length(setdiff(n_up, p_up)),
    length(setdiff(p_down, n_down)), length(setdiff(n_down, p_down))
  ),
  stringsAsFactors = FALSE
)

did <- read.csv(did_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("Gene", "logFC", "P.Value", "FDR") %in% names(did))) stop("DID table lacks required columns")
did_sig <- did[did$FDR < fdr_cutoff & abs(did$logFC) >= lfc_cutoff, , drop = FALSE]
p_response_genes <- sig_p$Gene
did_and_p <- intersect(did_sig$Gene, p_response_genes)
did_not_p <- setdiff(did_sig$Gene, p_response_genes)
p_not_did <- setdiff(p_response_genes, did_sig$Gene)
did_context <- data.frame(
  metric = c(
    "ASC_P_Post_vs_Pre_significant", "ASC_N_Post_vs_Pre_significant", "DID_significant",
    "DID_and_ASC_P_response", "DID_not_ASC_P_response", "ASC_P_response_not_DID",
    "fraction_ASC_P_response_also_DID", "fraction_DID_also_ASC_P_response"
  ),
  value = c(
    length(p_response_genes), nrow(sig_n), nrow(did_sig), length(did_and_p), length(did_not_p), length(p_not_did),
    length(did_and_p) / length(p_response_genes), length(did_and_p) / nrow(did_sig)
  ),
  stringsAsFactors = FALSE
)

# Hallmark over-representation analysis for the six strict shared/unique sets.
# This follows the project's msigdbr Hallmark source while using Fisher/ORA,
# which is the appropriate test for discrete post hoc gene lists.
hallmark_raw <- suppressWarnings(msigdbr(species = "Mus musculus", category = "H"))
hallmark_map <- unique(hallmark_raw[, c("gs_name", "gene_symbol")])
hallmark_list <- split(hallmark_map$gene_symbol, hallmark_map$gs_name)
universe <- unique(paired$Gene)

ora_one <- function(genes, category) {
  genes <- intersect(unique(genes), universe)
  N <- length(universe)
  n <- length(genes)
  rows <- lapply(names(hallmark_list), function(pathway) {
    pathway_genes <- intersect(unique(hallmark_list[[pathway]]), universe)
    M <- length(pathway_genes)
    overlap <- intersect(genes, pathway_genes)
    k <- length(overlap)
    if (M < 10L || M > 500L || n == 0L) return(NULL)
    d <- N - M - n + k
    odds <- if (d < 0L) NA_real_ else unname(fisher.test(matrix(c(k, n - k, M - k, d), nrow = 2), alternative = "greater")$estimate)
    data.frame(
      category = category,
      pathway = pathway,
      gene_set_size = n,
      pathway_size_in_universe = M,
      overlap_count = k,
      GeneRatio = paste0(k, "/", n),
      BgRatio = paste0(M, "/", N),
      fold_enrichment = if (k == 0L) 0 else (k / n) / (M / N),
      odds_ratio = odds,
      p_value = phyper(k - 1L, M, N - M, n, lower.tail = FALSE),
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || !nrow(out)) return(data.frame())
  out$FDR <- p.adjust(out$p_value, method = "BH")
  out$significant <- out$FDR < 0.05
  out[order(out$FDR, -out$fold_enrichment, out$pathway), , drop = FALSE]
}

enrichment_set_names <- c("shared_up", "shared_down", "ASC_P_only_up", "ASC_P_only_down", "ASC_N_only_up", "ASC_N_only_down")
enrichment_results <- setNames(lapply(enrichment_set_names, function(z) ora_one(gene_sets[[z]], z)), enrichment_set_names)
enrichment_all <- do.call(rbind, enrichment_results)
rownames(enrichment_all) <- NULL

cat("Preflight passed\n")
cat("Project root:", project_root, "\n")
cat("Expression universe:", length(universe), "genes\n")
cat("ASC_P/ASC_N significant:", nrow(sig_p), "/", nrow(sig_n), "\n")
cat("Strict unique up P/N:", length(gene_sets$ASC_P_only_up), "/", length(gene_sets$ASC_N_only_up), "\n")
cat("Strict unique down P/N:", length(gene_sets$ASC_P_only_down), "/", length(gene_sets$ASC_N_only_down), "\n")
cat("Discordant P-up/N-down and P-down/N-up:", length(gene_sets$discordant_P_up_N_down), "/", length(gene_sets$discordant_P_down_N_up), "\n")
cat("DID/P-response overlap:", length(did_and_p), "of", length(p_response_genes), "P-response genes and", nrow(did_sig), "DID genes\n")
cat("Hallmark sets:", length(hallmark_list), "\n")
if (dry_run) quit(save = "no", status = 0)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, name) {
  write.table(x, file.path(output_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

style_workbook <- function(wb, sheets) {
  for (sheet in sheets) {
    freezePane(wb, sheet, firstRow = TRUE)
    setColWidths(wb, sheet, cols = 1:30, widths = "auto")
  }
}

write_deg_workbook <- function(full, sig, path, subset_name) {
  wb <- createWorkbook()
  addWorksheet(wb, "README")
  note <- data.frame(
    item = c("contrast", "method", "input", "multiple_testing", "significance", "baseMean_note", "replacement_expression_columns", "design_note"),
    value = c(
      paste(subset_name, "Post vs Pre"), "limma", "fractional normalized read counts transformed as log2(count+1)",
      "Benjamini-Hochberg FDR within this contrast", "FDR < 0.05 and |log2FC| >= 0.5",
      "DESeq2 baseMean is unavailable because raw integer counts are not present and DESeq2 was not used",
      "mean_normalized_count_contrast_all, mean_normalized_count_Pre, mean_normalized_count_Post, and AveExpr_log2_all12",
      "unpaired model; donor/batch metadata were unavailable"
    ), stringsAsFactors = FALSE
  )
  writeData(wb, "README", note)
  addWorksheet(wb, "full_results"); writeData(wb, "full_results", full, withFilter = TRUE)
  addWorksheet(wb, "significant_DEGs"); writeData(wb, "significant_DEGs", sig, withFilter = TRUE)
  style_workbook(wb, c("README", "full_results", "significant_DEGs"))
  saveWorkbook(wb, path, overwrite = TRUE)
}

write_deg_workbook(table_p, sig_p, file.path(output_dir, "Table_Sx_ASC_P_Post_vs_Pre_limma.xlsx"), "ASC_P")
write_deg_workbook(table_n, sig_n, file.path(output_dir, "Table_Sy_ASC_N_Post_vs_Pre_limma.xlsx"), "ASC_N")
write_tsv(table_p, "Table_Sx_ASC_P_Post_vs_Pre_full.tsv")
write_tsv(sig_p, "Table_Sx_ASC_P_Post_vs_Pre_significant.tsv")
write_tsv(table_n, "Table_Sy_ASC_N_Post_vs_Pre_full.tsv")
write_tsv(sig_n, "Table_Sy_ASC_N_Post_vs_Pre_significant.tsv")

gene_set_detail <- function(genes, category) {
  out <- paired[match(sort(unique(genes)), paired$Gene), , drop = FALSE]
  out$response_class <- category
  out
}

wb_sets <- createWorkbook()
addWorksheet(wb_sets, "README")
set_notes <- data.frame(
  item = c("threshold", "shared_definition", "strict_unique_definition", "discordant_definition", "Venn_note"),
  value = c(
    "FDR < 0.05 and |log2FC| >= 0.5",
    "significant in both subsets with the same direction",
    "significant in one subset while the other subset is NS",
    "significant in both subsets with opposite directions",
    "Two-set Venn unique regions are mathematical set differences and therefore include discordant genes; strict unique and discordant genes are separated in this workbook"
  ), stringsAsFactors = FALSE
)
writeData(wb_sets, "README", set_notes)
addWorksheet(wb_sets, "all_gene_classification"); writeData(wb_sets, "all_gene_classification", paired, withFilter = TRUE)
for (sheet in set_names) {
  addWorksheet(wb_sets, sheet)
  writeData(wb_sets, sheet, gene_set_detail(gene_sets[[sheet]], sheet), withFilter = TRUE)
}
addWorksheet(wb_sets, "summary_counts"); writeData(wb_sets, "summary_counts", summary_counts)
style_workbook(wb_sets, names(wb_sets))
saveWorkbook(wb_sets, file.path(output_dir, "Table_Sz_shared_unique_gene_sets.xlsx"), overwrite = TRUE)
write_tsv(summary_counts, "Table_Sz_shared_unique_summary_counts.tsv")

wb_did <- createWorkbook()
addWorksheet(wb_did, "summary"); writeData(wb_did, "summary", did_context)
addWorksheet(wb_did, "DID_significant"); writeData(wb_did, "DID_significant", did_sig, withFilter = TRUE)
addWorksheet(wb_did, "DID_and_ASC_P_response"); writeData(wb_did, "DID_and_ASC_P_response", did_sig[match(did_and_p, did_sig$Gene), , drop = FALSE], withFilter = TRUE)
addWorksheet(wb_did, "DID_not_ASC_P_response"); writeData(wb_did, "DID_not_ASC_P_response", did_sig[match(did_not_p, did_sig$Gene), , drop = FALSE], withFilter = TRUE)
addWorksheet(wb_did, "ASC_P_response_not_DID"); writeData(wb_did, "ASC_P_response_not_DID", table_p[match(p_not_did, table_p$Gene), , drop = FALSE], withFilter = TRUE)
style_workbook(wb_did, names(wb_did))
saveWorkbook(wb_did, file.path(output_dir, "Table_Sw_DID_context.xlsx"), overwrite = TRUE)
write_tsv(did_context, "Table_Sw_DID_context_summary.tsv")

wb_enr <- createWorkbook()
addWorksheet(wb_enr, "README")
enr_notes <- data.frame(
  item = c("database", "species_mapping", "method", "universe", "multiple_testing", "tested_categories"),
  value = c(
    "MSigDB Hallmark via msigdbr", "Mus musculus ortholog mapping, matching the existing project convention",
    "one-sided Fisher/hypergeometric over-representation analysis",
    paste(length(universe), "genes tested in both limma contrasts"),
    "Benjamini-Hochberg correction across 50 Hallmark pathways separately within each response category",
    paste(enrichment_set_names, collapse = "; ")
  ), stringsAsFactors = FALSE
)
writeData(wb_enr, "README", enr_notes)
addWorksheet(wb_enr, "all_results"); writeData(wb_enr, "all_results", enrichment_all, withFilter = TRUE)
for (sheet in enrichment_set_names) {
  addWorksheet(wb_enr, sheet)
  writeData(wb_enr, sheet, enrichment_results[[sheet]], withFilter = TRUE)
}
style_workbook(wb_enr, names(wb_enr))
saveWorkbook(wb_enr, file.path(output_dir, "Table_Sv_shared_unique_Hallmark_enrichment.xlsx"), overwrite = TRUE)
write_tsv(enrichment_all, "Table_Sv_shared_unique_Hallmark_enrichment_all.tsv")

plot_overlap <- function(p_set, n_set, title_text, pdf_name, png_name) {
  theta <- seq(0, 2 * pi, length.out = 400)
  circles <- rbind(
    data.frame(x = -0.72 + cos(theta), y = sin(theta), set = "ASC_P"),
    data.frame(x = 0.72 + cos(theta), y = sin(theta), set = "ASC_N")
  )
  counts <- data.frame(
    x = c(-1.08, 0, 1.08), y = c(0, 0, 0),
    label = c(length(setdiff(p_set, n_set)), length(intersect(p_set, n_set)), length(setdiff(n_set, p_set)))
  )
  p <- ggplot(circles, aes(x, y, group = set, fill = set)) +
    geom_polygon(alpha = 0.28, color = NA) +
    geom_path(aes(color = set), linewidth = 0.6) +
    geom_text(data = counts, aes(x, y, label = label), inherit.aes = FALSE, size = 4, fontface = "bold") +
    annotate("text", x = -1.15, y = 1.22, label = "ASC_P", size = 4, fontface = "bold") +
    annotate("text", x = 1.15, y = 1.22, label = "ASC_N", size = 4, fontface = "bold") +
    annotate("text", x = 0, y = -1.32, label = "BH FDR < 0.05; |log2FC| >= 0.5", size = 3) +
    scale_fill_manual(values = c(ASC_P = "#D95F02", ASC_N = "#1B9E77")) +
    scale_color_manual(values = c(ASC_P = "#D95F02", ASC_N = "#1B9E77")) +
    coord_fixed(xlim = c(-2.05, 2.05), ylim = c(-1.55, 1.55), clip = "off") +
    labs(title = title_text) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none",
          plot.margin = margin(8, 12, 8, 12))
  ggsave(file.path(output_dir, pdf_name), p, width = 6.0, height = 4.6, device = cairo_pdf)
  ggsave(file.path(output_dir, png_name), p, width = 6.0, height = 4.6, dpi = 300, bg = "white")
}

plot_overlap(p_up, n_up, "Shared and unique upregulated responses",
             "FigSx_ASC_P_ASC_N_overlap_up.pdf", "FigSx_ASC_P_ASC_N_overlap_up.png")
plot_overlap(p_down, n_down, "Shared and unique downregulated responses",
             "FigSx_ASC_P_ASC_N_overlap_down.pdf", "FigSx_ASC_P_ASC_N_overlap_down.png")

sig_enrichment <- enrichment_all[enrichment_all$FDR < 0.05 & enrichment_all$overlap_count > 0, , drop = FALSE]
if (nrow(sig_enrichment)) {
  top_enrichment <- do.call(rbind, lapply(enrichment_set_names, function(z) {
    x <- sig_enrichment[sig_enrichment$category == z, , drop = FALSE]
    head(x[order(x$FDR, -x$fold_enrichment), , drop = FALSE], 5)
  }))
  top_enrichment$pathway_label <- gsub("^HALLMARK_", "", top_enrichment$pathway)
  top_enrichment$pathway_label <- gsub("_", " ", top_enrichment$pathway_label)
  top_enrichment$category <- factor(top_enrichment$category, levels = enrichment_set_names)
  top_enrichment$plot_label <- paste(top_enrichment$category, top_enrichment$pathway_label, sep = "___")
  top_enrichment$plot_label <- factor(top_enrichment$plot_label, levels = rev(top_enrichment$plot_label))
  palette <- c(
    shared_up = "#D95F02", shared_down = "#1B9E77",
    ASC_P_only_up = "#E69F00", ASC_P_only_down = "#56B4E9",
    ASC_N_only_up = "#CC79A7", ASC_N_only_down = "#0072B2"
  )
  category_labels <- c(
    shared_up = "Shared up", shared_down = "Shared down",
    ASC_P_only_up = "ASC_P-only up", ASC_P_only_down = "ASC_P-only down",
    ASC_N_only_up = "ASC_N-only up", ASC_N_only_down = "ASC_N-only down"
  )
  p_enr <- ggplot(top_enrichment, aes(x = -log10(FDR), y = plot_label, fill = category)) +
    geom_col(width = 0.72) +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_fill_manual(values = palette, drop = FALSE) +
    facet_wrap(~category, scales = "free", ncol = 2,
               labeller = as_labeller(category_labels)) +
    labs(
      x = expression(-log[10](FDR)), y = NULL,
      title = "Hallmark enrichment of shared and subset-specific responses",
      caption = "No Hallmark pathway met FDR < 0.05 for ASC_P-only up or ASC_N-only up."
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.caption = element_text(hjust = 0, size = 8),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  ggsave(file.path(output_dir, "FigSy_Hallmark_enrichment_shared_unique.pdf"), p_enr, width = 10, height = 9, device = cairo_pdf)
  ggsave(file.path(output_dir, "FigSy_Hallmark_enrichment_shared_unique.png"), p_enr, width = 10, height = 9, dpi = 300, bg = "white")
} else {
  warning("No significant Hallmark enrichment results; enrichment figure was not generated")
}

top_terms <- function(category, n = 3L) {
  x <- enrichment_results[[category]]
  x <- x[x$FDR < 0.05 & x$overlap_count > 0, , drop = FALSE]
  if (!nrow(x)) return("no Hallmark pathways at FDR < 0.05")
  terms <- head(x$pathway, n)
  terms <- gsub("^HALLMARK_", "", terms)
  terms <- tolower(gsub("_", " ", terms))
  paste(terms, collapse = ", ")
}

pct_p_response_did <- 100 * length(did_and_p) / length(p_response_genes)
pct_did_p_response <- 100 * length(did_and_p) / nrow(did_sig)
response_md <- c(
  "# Complete response package for Reviewer 2 Comment 6",
  "",
  "## What was added",
  "",
  "- Submission-ready ASC_P Post-vs-Pre and ASC_N Post-vs-Pre full/significant DEG tables.",
  "- Direction-specific overlap plots and an exhaustive shared/strict-unique/discordant gene classification.",
  "- DID-context analysis quantifying how much of the complete ASC_P response is represented by the interaction term.",
  "- Hallmark over-representation analysis for six shared/strict-unique response categories.",
  "",
  "## Key counts",
  "",
  sprintf("- ASC_P Post vs Pre: %d significant genes (%d up, %d down).", nrow(sig_p), length(p_up), length(p_down)),
  sprintf("- ASC_N Post vs Pre: %d significant genes (%d up, %d down).", nrow(sig_n), length(n_up), length(n_down)),
  sprintf("- Shared concordant responses: %d up and %d down.", length(gene_sets$shared_up), length(gene_sets$shared_down)),
  sprintf("- Strict subset-only responses (other subset NS): ASC_P %d up/%d down; ASC_N %d up/%d down.",
          length(gene_sets$ASC_P_only_up), length(gene_sets$ASC_P_only_down), length(gene_sets$ASC_N_only_up), length(gene_sets$ASC_N_only_down)),
  sprintf("- Discordant responses: %d ASC_P-up/ASC_N-down and %d ASC_P-down/ASC_N-up genes.",
          length(gene_sets$discordant_P_up_N_down), length(gene_sets$discordant_P_down_N_up)),
  sprintf("- DID context: %d significant DID genes; %d overlap the ASC_P Post-vs-Pre response. Thus %.1f%% of the complete ASC_P response is also DID-significant, while %.1f%% of DID genes are significant in the ASC_P within-subset response.",
          nrow(did_sig), length(did_and_p), pct_p_response_did, pct_did_p_response),
  "",
  "## Hallmark enrichment highlights",
  "",
  paste0("- Shared up: ", top_terms("shared_up"), "."),
  paste0("- Shared down: ", top_terms("shared_down"), "."),
  paste0("- ASC_P-only up: ", top_terms("ASC_P_only_up"), "."),
  paste0("- ASC_P-only down: ", top_terms("ASC_P_only_down"), "."),
  paste0("- ASC_N-only up: ", top_terms("ASC_N_only_up"), "."),
  paste0("- ASC_N-only down: ", top_terms("ASC_N_only_down"), "."),
  "",
  "## Results insertion",
  "",
  sprintf("To place the difference-in-differences interaction signature in the context of the complete coculture response, we additionally analyzed Post versus Pre expression changes within each ASC subset. ASC_P and ASC_N exhibited %d and %d significant differentially expressed genes, respectively (BH FDR < 0.05 and |log2FC| >= 0.5). The two subsets shared %d concordantly upregulated and %d concordantly downregulated genes. Requiring the other subset to be nonsignificant identified %d ASC_P-only upregulated, %d ASC_P-only downregulated, %d ASC_N-only upregulated, and %d ASC_N-only downregulated genes; %d genes showed significant opposite-direction responses. Of %d significant DID genes, %d overlapped the ASC_P within-subset response, corresponding to %.1f%% of all ASC_P-responsive genes. Thus, the DID signature represents a selective component rather than the entirety of the transcriptional response to coculture. Hallmark enrichment further resolved pathway programs shared by both subsets and those restricted to ASC_P or ASC_N (Figure Sy and Table Sv).",
          nrow(sig_p), nrow(sig_n), length(gene_sets$shared_up), length(gene_sets$shared_down),
          length(gene_sets$ASC_P_only_up), length(gene_sets$ASC_P_only_down),
          length(gene_sets$ASC_N_only_up), length(gene_sets$ASC_N_only_down),
          length(gene_sets$discordant_P_up_N_down) + length(gene_sets$discordant_P_down_N_up),
          nrow(did_sig), length(did_and_p), pct_p_response_did),
  "",
  "## Methods insertion",
  "",
  "The supplied fractional normalized read counts were transformed as log2(normalized count + 1) and analyzed with limma using the established unpaired subtype-by-timepoint model. Separate contrasts tested ASC_P Post versus Pre and ASC_N Post versus Pre. P values were adjusted by the Benjamini-Hochberg method, and responsive genes were defined as FDR < 0.05 and |log2FC| >= 0.5. Shared genes were significant in both subsets with concordant direction; strict subset-only genes were significant in one subset while nonsignificant in the other; significant opposite-direction genes were classified separately. Hallmark over-representation was tested with a one-sided Fisher/hypergeometric test using all 9,255 genes tested in both contrasts as the background, followed by BH correction separately within each response category.",
  "",
  "## Response-letter insertion",
  "",
  sprintf("We thank the reviewer for this important suggestion. We agree that the DID interaction term alone does not show the complete transcriptional response of each ASC subset. We therefore added full ASC_P Post-vs-Pre and ASC_N Post-vs-Pre differential-expression results (Tables Sx and Sy), direction-specific overlap analyses separating shared, strict subset-only, and discordant responses (Figure Sx and Table Sz), and Hallmark enrichment of the shared and subset-specific gene sets (Figure Sy and Table Sv). ASC_P and ASC_N showed %d and %d significant Post-vs-Pre genes, with %d shared upregulated and %d shared downregulated genes. Importantly, only %d of %d ASC_P-responsive genes were also DID-significant (%.1f%%), showing that the interaction signature is a selective component of the broader coculture response. The DID analysis is retained because it specifically identifies ASC_P-selective changes beyond those occurring in ASC_N, while the new analyses provide the complete shared and subset-specific response context requested by the reviewer.",
          nrow(sig_p), nrow(sig_n), length(gene_sets$shared_up), length(gene_sets$shared_down),
          length(did_and_p), length(p_response_genes), pct_p_response_did),
  "",
  "## Expression-summary note",
  "",
  "The project contains fractional normalized counts rather than raw integer counts. Therefore DESeq2 was not used and a DESeq2 `baseMean` cannot be reported. The submission tables provide transparent method-appropriate alternatives: mean normalized count across all six samples in each within-subset contrast, separate Pre and Post mean normalized counts, and limma AveExpr on the log2 scale."
)
writeLines(response_md, file.path(output_dir, "REVIEWER2_COMPLETE_RESPONSE_PACKAGE.md"), useBytes = TRUE)

generated <- sort(list.files(output_dir, full.names = FALSE))
completion <- c(
  "# Complete Reviewer 2 response analysis: completion report",
  "",
  "## Design and inputs",
  "",
  "- Input: `20260126 bulk ASC.csv` (fractional normalized counts; 9,255 genes; 12 samples).",
  "- Statistical unit: sample/replicate column, assumed biological as in the original project; independent donor identity is not documented.",
  "- Groups: ASC_P Pre/Post and ASC_N Pre/Post, n = 3 per group.",
  "- Model: existing unpaired limma subtype-by-timepoint model; no batch/donor covariate was available.",
  "- Threshold: BH FDR < 0.05 and |log2FC| >= 0.5.",
  "- Enrichment: Mouse-symbol Hallmark ortholog sets via msigdbr; one-sided Fisher/hypergeometric ORA; 9,255 tested genes as universe; BH correction within each of six categories.",
  "- This is a discovery/revision analysis with no independent validation cohort.",
  "",
  "## Important interpretation",
  "",
  sprintf("- Strict unique counts exclude %d opposite-direction genes, which are reported separately.", length(gene_sets$discordant_P_up_N_down) + length(gene_sets$discordant_P_down_N_up)),
  sprintf("- %d/%d (%.1f%%) ASC_P-responsive genes are DID-significant; the DID is a selective subset of the complete response.", length(did_and_p), length(p_response_genes), pct_p_response_did),
  "- DESeq2 baseMean is unavailable and was not fabricated; method-appropriate normalized-count means and limma AveExpr are provided.",
  "",
  "## Generated files",
  "",
  paste0("- `", generated, "`"),
  "",
  "## Rerun commands",
  "",
  "```powershell",
  "& 'E:\\R\\R-4.3.3\\bin\\x64\\Rscript.exe' --vanilla 'D:\\R\\20260126 ASC bulk  SCI\\scripts\\reviewer2_ASC_bulk_complete_response.R' --dry-run",
  "& 'E:\\R\\R-4.3.3\\bin\\x64\\Rscript.exe' --vanilla 'D:\\R\\20260126 ASC bulk  SCI\\scripts\\reviewer2_ASC_bulk_complete_response.R'",
  "```"
)
writeLines(completion, file.path(output_dir, "CODEX_COMPLETE_RESPONSE_COMPLETION_REPORT.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

required_outputs <- c(
  "Table_Sx_ASC_P_Post_vs_Pre_limma.xlsx", "Table_Sy_ASC_N_Post_vs_Pre_limma.xlsx",
  "Table_Sz_shared_unique_gene_sets.xlsx", "Table_Sw_DID_context.xlsx",
  "Table_Sv_shared_unique_Hallmark_enrichment.xlsx",
  "FigSx_ASC_P_ASC_N_overlap_up.pdf", "FigSx_ASC_P_ASC_N_overlap_up.png",
  "FigSx_ASC_P_ASC_N_overlap_down.pdf", "FigSx_ASC_P_ASC_N_overlap_down.png",
  "REVIEWER2_COMPLETE_RESPONSE_PACKAGE.md", "CODEX_COMPLETE_RESPONSE_COMPLETION_REPORT.md"
)
if (nrow(sig_enrichment)) {
  required_outputs <- c(required_outputs, "FigSy_Hallmark_enrichment_shared_unique.pdf", "FigSy_Hallmark_enrichment_shared_unique.png")
}
missing_outputs <- required_outputs[!file.exists(file.path(output_dir, required_outputs))]
if (length(missing_outputs)) stop("Missing required outputs: ", paste(missing_outputs, collapse = ", "))
if (any(file.info(file.path(output_dir, required_outputs))$size == 0)) stop("One or more required outputs are empty")

cat("Complete Reviewer 2 response package generated at:\n", output_dir, "\n")
print(summary_counts, row.names = FALSE)
print(did_context, row.names = FALSE)

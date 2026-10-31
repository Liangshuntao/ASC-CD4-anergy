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
  library(limma)
  library(ggplot2)
  library(openxlsx)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
project_root <- .output
input_file <- .input
output_dir <- file.path(.output, "within_subset")

required_packages <- c("limma", "ggplot2", "openxlsx", "msigdbr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

fdr_cutoff <- 0.05
lfc_cutoff <- 0.5
primary_design_label <- "unpaired limma 2x2 interaction model"

raw <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required_id_cols <- c("Gene ID", "Gene Symbol")
if (!all(required_id_cols %in% names(raw))) stop("Missing gene identifier columns")

sample_cols <- grep("^Normalized read count ", names(raw), value = TRUE)
sample_names <- sub("^Normalized read count ", "", sample_cols)
expected_samples <- c(
  paste0("ASC_N_", 1:3), paste0("ASC_P_", 1:3),
  paste0("after_ASC_N_", 1:3), paste0("after_ASC_P_", 1:3)
)
if (!setequal(sample_names, expected_samples) || length(sample_names) != 12L) {
  stop("The 12 expected ASC sample columns were not found exactly")
}

normalized_counts <- as.matrix(raw[, sample_cols, drop = FALSE])
storage.mode(normalized_counts) <- "numeric"
if (anyNA(normalized_counts) || any(!is.finite(normalized_counts))) stop("Expression matrix contains missing/nonfinite values")
if (any(normalized_counts < 0)) stop("Expression matrix contains negative values")
if (!any(abs(normalized_counts - round(normalized_counts)) > 1e-8)) {
  warning("Input appears integer-like; this script intentionally preserves the established normalized-count limma workflow")
}

gene_symbol <- trimws(gsub("'", "", raw[["Gene Symbol"]], fixed = TRUE))
gene_key <- ifelse(is.na(gene_symbol) | gene_symbol == "", as.character(raw[["Gene ID"]]), gene_symbol)
log_expr <- log2(normalized_counts + 1)

# Match the existing script's symbol-level median aggregation. No duplicates are
# present in the current file, but retaining this step makes reruns deterministic.
aggregate_df <- data.frame(Gene = gene_key, log_expr, check.names = FALSE)
expr_gene <- aggregate(. ~ Gene, data = aggregate_df, FUN = median, na.rm = TRUE)
expr_mat <- as.matrix(expr_gene[, -1, drop = FALSE])
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr_gene$Gene
colnames(expr_mat) <- sample_names

sample_info <- data.frame(
  sample = sample_names,
  timepoint = ifelse(grepl("^after_", sample_names), "Post", "Pre"),
  subtype = ifelse(grepl("ASC_P_", sample_names), "ASC_P", "ASC_N"),
  replicate_suffix = sub(".*_([0-9]+)$", "\\1", sample_names),
  stringsAsFactors = FALSE
)
sample_info$timepoint <- factor(sample_info$timepoint, levels = c("Pre", "Post"))
sample_info$subtype <- factor(sample_info$subtype, levels = c("ASC_N", "ASC_P"))
rownames(sample_info) <- sample_info$sample
sample_info <- sample_info[colnames(expr_mat), , drop = FALSE]

group_counts <- with(sample_info, table(subtype, timepoint))
if (any(group_counts != 3L)) stop("Expected exactly three samples per subtype/timepoint group")

design <- model.matrix(~ 0 + subtype * timepoint, data = sample_info)
colnames(design) <- make.names(colnames(design))
if (qr(design)$rank != ncol(design)) stop("Primary design matrix is not full rank")

contrast_matrix <- makeContrasts(
  ASC_P_Post_vs_Pre = timepointPost + subtypeASC_P.timepointPost,
  ASC_N_Post_vs_Pre = timepointPost,
  levels = design
)

pca <- prcomp(t(expr_mat), scale. = FALSE)
pca_variance <- 100 * summary(pca)$importance[2, 1:2]

cat("Preflight passed\n")
cat("Project root:", project_root, "\n")
cat("Input:", input_file, "\n")
cat("Genes:", nrow(expr_mat), " Samples:", ncol(expr_mat), "\n")
cat("Design rank:", qr(design)$rank, "/", ncol(design), "\n")
cat("Primary method:", primary_design_label, "\n")
cat("Threshold: BH FDR <", fdr_cutoff, "and |log2FC| >=", lfc_cutoff, "\n")
if (dry_run) quit(save = "no", status = 0)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(sample_info, file.path(output_dir, "sample_metadata_inferred.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

fit <- lmFit(expr_mat, design)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

extract_result <- function(coef_name, contrast_label) {
  out <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  out <- data.frame(Gene = rownames(out), out, row.names = NULL, check.names = FALSE)
  out$FDR <- p.adjust(out$P.Value, method = "BH")
  out$log2FoldChange <- out$logFC
  out$pvalue <- out$P.Value
  out$padj <- out$FDR
  out$direction <- ifelse(
    out$FDR < fdr_cutoff & out$logFC >= lfc_cutoff, "Up",
    ifelse(out$FDR < fdr_cutoff & out$logFC <= -lfc_cutoff, "Down", "NS")
  )
  out$significant <- out$direction != "NS"
  out$contrast <- contrast_label
  out$method <- "limma"
  out$input_scale <- "log2(normalized read count + 1)"
  out$primary_design <- primary_design_label
  out
}

res_p <- extract_result("ASC_P_Post_vs_Pre", "ASC_P Post vs Pre")
res_n <- extract_result("ASC_N_Post_vs_Pre", "ASC_N Post vs Pre")
sig_p <- res_p[res_p$significant, , drop = FALSE]
sig_n <- res_n[res_n$significant, , drop = FALSE]

write_tsv <- function(x, filename) {
  write.table(x, file.path(output_dir, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# Canonical method-correct files.
write_tsv(res_p, "ASC_P_Post_vs_Pre_limma_full.tsv")
write_tsv(sig_p, "ASC_P_Post_vs_Pre_limma_sig.tsv")
write_tsv(res_n, "ASC_N_Post_vs_Pre_limma_full.tsv")
write_tsv(sig_n, "ASC_N_Post_vs_Pre_limma_sig.tsv")

# Only method-correct limma filenames are exported; legacy duplicate aliases omitted.

sets <- list(
  ASC_P_up = res_p$Gene[res_p$direction == "Up"],
  ASC_N_up = res_n$Gene[res_n$direction == "Up"],
  ASC_P_down = res_p$Gene[res_p$direction == "Down"],
  ASC_N_down = res_n$Gene[res_n$direction == "Down"]
)
sets$shared_up <- intersect(sets$ASC_P_up, sets$ASC_N_up)
sets$ASC_P_only_up <- setdiff(sets$ASC_P_up, sets$ASC_N_up)
sets$ASC_N_only_up <- setdiff(sets$ASC_N_up, sets$ASC_P_up)
sets$shared_down <- intersect(sets$ASC_P_down, sets$ASC_N_down)
sets$ASC_P_only_down <- setdiff(sets$ASC_P_down, sets$ASC_N_down)
sets$ASC_N_only_down <- setdiff(sets$ASC_N_down, sets$ASC_P_down)

summary_counts <- data.frame(
  category = c(
    "ASC_P_significant_total", "ASC_N_significant_total",
    "ASC_P_up", "ASC_N_up", "shared_up", "ASC_P_only_up", "ASC_N_only_up",
    "ASC_P_down", "ASC_N_down", "shared_down", "ASC_P_only_down", "ASC_N_only_down"
  ),
  count = c(
    nrow(sig_p), nrow(sig_n),
    length(sets$ASC_P_up), length(sets$ASC_N_up), length(sets$shared_up),
    length(sets$ASC_P_only_up), length(sets$ASC_N_only_up),
    length(sets$ASC_P_down), length(sets$ASC_N_down), length(sets$shared_down),
    length(sets$ASC_P_only_down), length(sets$ASC_N_only_down)
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary_counts, "ASC_P_ASC_N_shared_unique_summary_counts.tsv")

gene_set_table <- function(genes, label) {
  base <- data.frame(Gene = sort(unique(genes)), stringsAsFactors = FALSE)
  if (!nrow(base)) {
    return(data.frame(Gene = character(), set = character(), ASC_P_logFC = numeric(),
                      ASC_P_FDR = numeric(), ASC_N_logFC = numeric(), ASC_N_FDR = numeric()))
  }
  p_stats <- res_p[, c("Gene", "logFC", "FDR")]
  names(p_stats)[-1] <- c("ASC_P_logFC", "ASC_P_FDR")
  n_stats <- res_n[, c("Gene", "logFC", "FDR")]
  names(n_stats)[-1] <- c("ASC_N_logFC", "ASC_N_FDR")
  out <- merge(base, p_stats, by = "Gene", all.x = TRUE, sort = FALSE)
  out <- merge(out, n_stats, by = "Gene", all.x = TRUE, sort = FALSE)
  out$set <- label
  out[, c("Gene", "set", "ASC_P_logFC", "ASC_P_FDR", "ASC_N_logFC", "ASC_N_FDR")]
}

wb <- createWorkbook()
addWorksheet(wb, "ASC_P_Post_vs_Pre_full"); writeData(wb, "ASC_P_Post_vs_Pre_full", res_p)
addWorksheet(wb, "ASC_N_Post_vs_Pre_full"); writeData(wb, "ASC_N_Post_vs_Pre_full", res_n)
addWorksheet(wb, "ASC_P_sig"); writeData(wb, "ASC_P_sig", sig_p)
addWorksheet(wb, "ASC_N_sig"); writeData(wb, "ASC_N_sig", sig_n)
for (sheet in c("shared_up", "ASC_P_only_up", "ASC_N_only_up", "shared_down", "ASC_P_only_down", "ASC_N_only_down")) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, gene_set_table(sets[[sheet]], sheet))
}
addWorksheet(wb, "summary_counts"); writeData(wb, "summary_counts", summary_counts)
saveWorkbook(wb, file.path(output_dir, "ASC_P_ASC_N_shared_unique_gene_sets.xlsx"), overwrite = TRUE)

plot_overlap <- function(p_set, n_set, shared_set, direction_label, pdf_name, png_name) {
  theta <- seq(0, 2 * pi, length.out = 400)
  circles <- rbind(
    data.frame(x = -0.72 + cos(theta), y = sin(theta), set = "ASC_P"),
    data.frame(x = 0.72 + cos(theta), y = sin(theta), set = "ASC_N")
  )
  counts <- data.frame(
    x = c(-1.08, 0, 1.08), y = c(0, 0, 0),
    label = c(length(setdiff(p_set, n_set)), length(shared_set), length(setdiff(n_set, p_set)))
  )
  p <- ggplot(circles, aes(x, y, group = set, fill = set)) +
    geom_polygon(alpha = 0.28, color = NA) +
    geom_path(aes(color = set), linewidth = 0.6) +
    geom_text(data = counts, aes(x, y, label = label), inherit.aes = FALSE, size = 4, fontface = "bold") +
    annotate("text", x = -1.15, y = 1.22, label = "ASC_P", size = 4, fontface = "bold") +
    annotate("text", x = 1.15, y = 1.22, label = "ASC_N", size = 4, fontface = "bold") +
    annotate("text", x = 0, y = -1.32,
             label = "BH FDR < 0.05; |log2FC| >= 0.5", size = 3) +
    scale_fill_manual(values = c(ASC_P = "#D95F02", ASC_N = "#1B9E77")) +
    scale_color_manual(values = c(ASC_P = "#D95F02", ASC_N = "#1B9E77")) +
    coord_fixed(xlim = c(-2.05, 2.05), ylim = c(-1.55, 1.55), clip = "off") +
    labs(title = paste("Shared and unique", tolower(direction_label), "responses")) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none",
          plot.margin = margin(8, 12, 8, 12))
  ggsave(file.path(output_dir, pdf_name), p, width = 6.0, height = 4.6, device = cairo_pdf)
  ggsave(file.path(output_dir, png_name), p, width = 6.0, height = 4.6, dpi = 300, bg = "white")
}

plot_overlap(sets$ASC_P_up, sets$ASC_N_up, sets$shared_up, "Upregulated",
             "FigS_ASC_P_ASC_N_overlap_up.pdf", "FigS_ASC_P_ASC_N_overlap_up.png")
plot_overlap(sets$ASC_P_down, sets$ASC_N_down, sets$shared_down, "Downregulated",
             "FigS_ASC_P_ASC_N_overlap_down.pdf", "FigS_ASC_P_ASC_N_overlap_down.png")

qc_lines <- c(
  paste("genes", nrow(expr_mat), sep = "\t"),
  paste("samples", ncol(expr_mat), sep = "\t"),
  paste("fractional_values", sum(abs(normalized_counts - round(normalized_counts)) > 1e-8), sep = "\t"),
  paste("missing_values", sum(is.na(normalized_counts)), sep = "\t"),
  paste("negative_values", sum(normalized_counts < 0), sep = "\t"),
  paste("design_rank", paste0(qr(design)$rank, "/", ncol(design)), sep = "\t"),
  paste("PC1_variance_percent", round(pca_variance[[1]], 2), sep = "\t"),
  paste("PC2_variance_percent", round(pca_variance[[2]], 2), sep = "\t")
)
writeLines(qc_lines, file.path(output_dir, "analysis_QC_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

required_outputs <- c(
  "ASC_P_Post_vs_Pre_limma_full.tsv", "ASC_P_Post_vs_Pre_limma_sig.tsv",
  "ASC_N_Post_vs_Pre_limma_full.tsv", "ASC_N_Post_vs_Pre_limma_sig.tsv",
  "ASC_P_ASC_N_shared_unique_gene_sets.xlsx",
  "FigS_ASC_P_ASC_N_overlap_up.pdf", "FigS_ASC_P_ASC_N_overlap_up.png",
  "FigS_ASC_P_ASC_N_overlap_down.pdf", "FigS_ASC_P_ASC_N_overlap_down.png"
)
missing_outputs <- required_outputs[!file.exists(file.path(output_dir, required_outputs))]
if (length(missing_outputs)) stop("Missing required outputs: ", paste(missing_outputs, collapse = ", "))
if (any(file.info(file.path(output_dir, required_outputs))$size == 0)) stop("One or more required outputs are empty")

cat("Analysis complete. Outputs:", output_dir, "\n")
print(summary_counts, row.names = FALSE)

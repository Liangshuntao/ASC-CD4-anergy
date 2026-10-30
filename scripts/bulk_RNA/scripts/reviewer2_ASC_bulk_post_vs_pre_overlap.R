#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(openxlsx)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
root_arg <- grep("^--project-root=", args, value = TRUE)

if (length(root_arg)) {
  project_root <- sub("^--project-root=", "", root_arg[[1]])
} else {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
  project_root <- dirname(dirname(script_path))
}

project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
input_file <- file.path(project_root, "20260126 bulk ASC.csv")
output_dir <- file.path(project_root, "revision_outputs", "reviewer2_transcriptomic_signature")

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

# Exact legacy filenames requested by the revision brief. Their method column
# explicitly records limma; DESeq2 is not run on fractional normalized values.
write_tsv(res_p, "ASC_P_Post_vs_Pre_DESeq2_full.tsv")
write_tsv(sig_p, "ASC_P_Post_vs_Pre_DESeq2_sig.tsv")
write_tsv(res_n, "ASC_N_Post_vs_Pre_DESeq2_full.tsv")
write_tsv(sig_n, "ASC_N_Post_vs_Pre_DESeq2_sig.tsv")

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

count_of <- function(category) summary_counts$count[match(category, summary_counts$category)]
summary_md <- c(
  "# Reviewer 2 transcriptomic signature summary",
  "",
  "## Key results",
  "",
  sprintf("- ASC_P Post vs Pre: %d significant DEGs.", count_of("ASC_P_significant_total")),
  sprintf("- ASC_N Post vs Pre: %d significant DEGs.", count_of("ASC_N_significant_total")),
  sprintf("- Upregulated: %d shared, %d ASC_P-only, and %d ASC_N-only genes.",
          count_of("shared_up"), count_of("ASC_P_only_up"), count_of("ASC_N_only_up")),
  sprintf("- Downregulated: %d shared, %d ASC_P-only, and %d ASC_N-only genes.",
          count_of("shared_down"), count_of("ASC_P_only_down"), count_of("ASC_N_only_down")),
  "- Responsive genes were defined as BH FDR < 0.05 and |log2FC| >= 0.5, matching the established Figure 8 threshold.",
  "- The primary model was unpaired because no donor/subject metadata were available to confirm that replicate suffixes denote matched biological samples.",
  "",
  "These complementary contrasts address Reviewer 2's concern by showing the complete transcriptional response to coculture within each ASC subset, while the direction-specific overlap separates common responses from subset-restricted responses. The original difference-in-differences interaction remains the appropriate analysis for ASC_P-specific change beyond the response observed in ASC_N.",
  "",
  "## Suggested response-letter citation",
  "",
  "We suggest citing the overlap plots as **Figure Sx: Shared and unique transcriptional responses of ASC subsets to coculture** and the workbook/full result tables as **Table Sx: Differential expression results for ASC_P and ASC_N Post-vs-Pre coculture responses**.",
  "",
  "## Results insertion",
  "",
  sprintf("In addition to the difference-in-differences interaction analysis, we examined the complete coculture response within each ASC subset. ASC_P and ASC_N showed %d and %d significant Post-vs-Pre differentially expressed genes, respectively (BH FDR < 0.05 and |log2FC| >= 0.5). Direction-specific overlap analysis identified %d shared upregulated and %d shared downregulated genes, together with subset-restricted responses (%d ASC_P-only upregulated, %d ASC_N-only upregulated, %d ASC_P-only downregulated, and %d ASC_N-only downregulated genes). These results place the ASC_P-specific interaction signature within the broader shared and subset-specific transcriptional response to coculture.",
          count_of("ASC_P_significant_total"), count_of("ASC_N_significant_total"),
          count_of("shared_up"), count_of("shared_down"), count_of("ASC_P_only_up"),
          count_of("ASC_N_only_up"), count_of("ASC_P_only_down"), count_of("ASC_N_only_down")),
  "",
  "## Methods insertion",
  "",
  "Normalized read counts were transformed as log2(normalized count + 1), consistent with the original ASC bulk RNA-seq workflow, and analyzed using limma with a 2 x 2 subtype-by-timepoint model. Contrasts tested ASC_P Post versus ASC_P Pre and ASC_N Post versus ASC_N Pre. P values were adjusted by the Benjamini-Hochberg method; responsive genes were defined by FDR < 0.05 and |log2 fold change| >= 0.5. Shared upregulated or downregulated genes met the threshold with concordant direction in both contrasts; subset-only genes met the threshold in one subset but were not significant in the same direction in the other subset. Samples were treated as unpaired because donor/subject identifiers were not available in the project metadata.",
  "",
  "## Response letter insertion",
  "",
  sprintf("We thank the reviewer for this constructive suggestion. We agree that the difference-in-differences interaction term alone does not display the complete transcriptional response of each ASC subset to coculture. We therefore added separate ASC_P Post-vs-Pre and ASC_N Post-vs-Pre analyses and a direction-specific shared/unique overlap analysis. These analyses identified %d significant DEGs in ASC_P and %d in ASC_N, including %d shared upregulated and %d shared downregulated genes. The overlap plots are provided as Figure Sx, and the full differential-expression results and shared/unique gene sets are provided as Table Sx and the accompanying supplemental files. We retained the original difference-in-differences analysis because it specifically tests the ASC_P-selective interaction response, whereas the new analyses provide the broader context of shared and subset-restricted coculture responses.",
          count_of("ASC_P_significant_total"), count_of("ASC_N_significant_total"),
          count_of("shared_up"), count_of("shared_down")),
  "",
  "## Method note",
  "",
  "The input contains fractional normalized read counts, so DESeq2 was not used. The exact `*_DESeq2_*` filenames requested in the task brief are retained as compatibility aliases, but every table explicitly identifies the method as limma. Canonical `*_limma_*` files are also supplied.",
  "",
  "Pathway enrichment was not run for the six overlap gene sets. The existing project uses competitive CAMERA testing on the full expression matrix and a model contrast; applying it directly to post hoc shared/unique gene lists would not answer the requested set-level enrichment question. Introducing an over-representation framework solely for this supplement would violate the instruction to preserve the established enrichment method."
)
writeLines(summary_md, file.path(output_dir, "REVIEWER2_TRANSCRIPTOMIC_SIGNATURE_SUMMARY.md"), useBytes = TRUE)

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

generated <- sort(list.files(output_dir, full.names = FALSE))
completion_md <- c(
  "# Codex completion report",
  "",
  "## Scripts inspected",
  "",
  "- `20260126 bulk ASC.R`",
  "- KS template: `KS_bulkRNA_MultiGroup_DEGs.R` and its multi-group demonstration script",
  "",
  "## New script",
  "",
  "- `scripts/reviewer2_ASC_bulk_post_vs_pre_overlap.R`",
  "",
  "## Input",
  "",
  "- `20260126 bulk ASC.csv` (9,255 genes; 12 samples; fractional normalized read counts)",
  "- Sample metadata were inferred from the 12 expression-column names because no separate metadata file was present.",
  "",
  "## Analysis decisions",
  "",
  "- Method: limma on log2(normalized read count + 1), matching the established Figure 8 workflow.",
  "- Contrasts: ASC_P Post vs Pre and ASC_N Post vs Pre.",
  "- Threshold: BH FDR < 0.05 and |log2FC| >= 0.5.",
  "- Statistical unit: biological sample/replicate column (n = 3 per subtype-timepoint group).",
  "- Design: unpaired. Suffixes 1-3 suggest possible matching, but no donor/subject metadata confirm pairing; the original Figure 8 limma model is also unpaired.",
  "- Batch/confounders: none included because no batch metadata were available. PCA QC is reported, but cannot establish absence of batch effects.",
  "- Multiple testing: Benjamini-Hochberg correction within each contrast.",
  "- Validation: structural matrix checks, exact sample-set checks, replicate-count checks, full-rank design check, output existence checks, and workbook/figure inspection.",
  "- Pathway enrichment: skipped for the six overlap sets because the project's CAMERA method is contrast-based and is not directly applicable to post hoc gene-list ORA.",
  "",
  "## Generated outputs",
  "",
  paste0("- `", generated, "`"),
  "",
  "## Warnings and limitations",
  "",
  "- DESeq2 is inappropriate for the fractional normalized matrix. Required `*_DESeq2_*` paths are compatibility aliases containing limma results and an explicit `method` column; canonical limma files are supplied.",
  "- The analysis is based on three replicate columns per group and has no independent validation cohort.",
  "- Pairing and batch structure should be revisited if donor/subject or experimental-batch metadata become available.",
  "",
  "## Exact rerun commands",
  "",
  "```powershell",
  "& 'E:\\R\\R-4.3.3\\bin\\x64\\Rscript.exe' --vanilla 'D:\\R\\20260126 ASC bulk  SCI\\scripts\\reviewer2_ASC_bulk_post_vs_pre_overlap.R' --dry-run",
  "& 'E:\\R\\R-4.3.3\\bin\\x64\\Rscript.exe' --vanilla 'D:\\R\\20260126 ASC bulk  SCI\\scripts\\reviewer2_ASC_bulk_post_vs_pre_overlap.R'",
  "```"
)
writeLines(completion_md, file.path(output_dir, "CODEX_COMPLETION_REPORT.md"), useBytes = TRUE)

required_outputs <- c(
  "ASC_P_Post_vs_Pre_DESeq2_full.tsv", "ASC_P_Post_vs_Pre_DESeq2_sig.tsv",
  "ASC_N_Post_vs_Pre_DESeq2_full.tsv", "ASC_N_Post_vs_Pre_DESeq2_sig.tsv",
  "ASC_P_ASC_N_shared_unique_gene_sets.xlsx",
  "FigS_ASC_P_ASC_N_overlap_up.pdf", "FigS_ASC_P_ASC_N_overlap_up.png",
  "FigS_ASC_P_ASC_N_overlap_down.pdf", "FigS_ASC_P_ASC_N_overlap_down.png",
  "REVIEWER2_TRANSCRIPTOMIC_SIGNATURE_SUMMARY.md", "CODEX_COMPLETION_REPORT.md"
)
missing_outputs <- required_outputs[!file.exists(file.path(output_dir, required_outputs))]
if (length(missing_outputs)) stop("Missing required outputs: ", paste(missing_outputs, collapse = ", "))
if (any(file.info(file.path(output_dir, required_outputs))$size == 0)) stop("One or more required outputs are empty")

cat("Analysis complete. Outputs:", output_dir, "\n")
print(summary_counts, row.names = FALSE)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
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

base_dir <- file.path(project_root, "revision_outputs", "reviewer2_transcriptomic_signature")
output_dir <- file.path(base_dir, "submission_ready")
input_file <- file.path(project_root, "20260126 bulk ASC.csv")
p_result_file <- file.path(output_dir, "Table_Sx_ASC_P_Post_vs_Pre_full.tsv")
n_result_file <- file.path(output_dir, "Table_Sy_ASC_N_Post_vs_Pre_full.tsv")

required_inputs <- c(input_file, p_result_file, n_result_file)
if (any(!file.exists(required_inputs))) {
  stop("Missing required inputs: ", paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
}

fdr_cutoff <- 0.05
lfc_cutoff <- 0.5
expected_counts <- c(ASC_P = 4126L, ASC_N = 5201L)

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

mat <- as.matrix(raw[, sample_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
if (anyNA(mat) || any(!is.finite(mat)) || any(mat < 0)) stop("Expression matrix contains invalid values")
if (!any(abs(mat - round(mat)) > 1e-8)) stop("Input unexpectedly appears integer-like; re-audit method selection")

gene_symbol <- trimws(gsub("'", "", raw[["Gene Symbol"]], fixed = TRUE))
gene_key <- ifelse(is.na(gene_symbol) | gene_symbol == "", as.character(raw[["Gene ID"]]), gene_symbol)
log_expr <- log2(mat + 1)
expr_df <- data.frame(Gene = gene_key, log_expr, check.names = FALSE)
expr_gene <- aggregate(. ~ Gene, data = expr_df, FUN = median, na.rm = TRUE)
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
if (qr(design)$rank != ncol(design)) stop("Design matrix is not full rank")
contrast_matrix <- makeContrasts(
  ASC_P_Post_vs_Pre = timepointPost + subtypeASC_P.timepointPost,
  ASC_N_Post_vs_Pre = timepointPost,
  levels = design
)

read_deg <- function(path, label) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("Gene", "log2FC", "p_value", "FDR", "direction", "method")
  if (!all(required %in% names(x))) stop(label, " DEG table lacks required columns")
  if (nrow(x) != nrow(expr_mat)) stop(label, " DEG table row count differs from expression universe")
  if (!identical(unique(x$method), "limma")) stop(label, " DEG table is not uniformly limma")
  recomputed <- ifelse(
    x$FDR < fdr_cutoff & x$log2FC >= lfc_cutoff, "Up",
    ifelse(x$FDR < fdr_cutoff & x$log2FC <= -lfc_cutoff, "Down", "NS")
  )
  if (any(recomputed != x$direction)) stop(label, " stored direction differs from threshold-recomputed direction")
  x$direction <- recomputed
  x
}

deg_p <- read_deg(p_result_file, "ASC_P")
deg_n <- read_deg(n_result_file, "ASC_N")
actual_counts <- c(ASC_P = sum(deg_p$direction != "NS"), ASC_N = sum(deg_n$direction != "NS"))
if (!identical(unname(actual_counts), unname(expected_counts))) {
  stop("Unexpected significant DEG counts: ", paste(names(actual_counts), actual_counts, collapse = "; "))
}

# Full ranked Hallmark analysis using the existing project CAMERA framework.
hallmark_raw <- suppressWarnings(msigdbr(species = "Mus musculus", category = "H"))
hallmark_map <- unique(hallmark_raw[, c("gs_name", "gene_symbol")])
hallmark_list <- split(hallmark_map$gene_symbol, hallmark_map$gs_name)
hallmark_index <- ids2indices(hallmark_list, rownames(expr_mat))

camera_result <- function(contrast_name, label) {
  out <- camera(
    y = expr_mat,
    index = hallmark_index,
    design = design,
    contrast = contrast_matrix[, contrast_name]
  )
  out <- data.frame(Pathway = rownames(out), out, row.names = NULL, check.names = FALSE)
  out$FDR <- p.adjust(out$PValue, method = "BH")
  out$contrast <- label
  out$signed_log10_FDR <- ifelse(out$Direction == "Up", 1, -1) * (-log10(pmax(out$FDR, .Machine$double.xmin)))
  out[order(out$FDR, out$Pathway), , drop = FALSE]
}

camera_p <- camera_result("ASC_P_Post_vs_Pre", "ASC_P Post vs Pre")
camera_n <- camera_result("ASC_N_Post_vs_Pre", "ASC_N Post vs Pre")

camera_compare <- merge(
  camera_p[, c("Pathway", "NGenes", "Direction", "PValue", "FDR", "signed_log10_FDR")],
  camera_n[, c("Pathway", "NGenes", "Direction", "PValue", "FDR", "signed_log10_FDR")],
  by = "Pathway", suffixes = c("_ASC_P", "_ASC_N"), sort = FALSE
)
camera_compare$classification <- with(camera_compare, ifelse(
  FDR_ASC_P < fdr_cutoff & FDR_ASC_N < fdr_cutoff & Direction_ASC_P == "Up" & Direction_ASC_N == "Up", "shared_up",
  ifelse(FDR_ASC_P < fdr_cutoff & FDR_ASC_N < fdr_cutoff & Direction_ASC_P == "Down" & Direction_ASC_N == "Down", "shared_down",
  ifelse(FDR_ASC_P < fdr_cutoff & FDR_ASC_N < fdr_cutoff & Direction_ASC_P != Direction_ASC_N, "discordant",
  ifelse(FDR_ASC_P < fdr_cutoff & FDR_ASC_N >= fdr_cutoff, "ASC_P_only",
  ifelse(FDR_ASC_N < fdr_cutoff & FDR_ASC_P >= fdr_cutoff, "ASC_N_only", "not_significant"))))))
camera_compare$min_FDR <- pmin(camera_compare$FDR_ASC_P, camera_compare$FDR_ASC_N)
camera_compare <- camera_compare[order(camera_compare$min_FDR, camera_compare$Pathway), , drop = FALSE]

camera_summary <- as.data.frame(table(camera_compare$classification), stringsAsFactors = FALSE)
names(camera_summary) <- c("classification", "pathway_count")

cat("Preflight passed\n")
cat("Project root:", project_root, "\n")
cat("Expression universe:", nrow(expr_mat), "genes; samples:", ncol(expr_mat), "\n")
cat("Design rank:", qr(design)$rank, "/", ncol(design), "\n")
cat("ASC_P/ASC_N significant DEG counts:", actual_counts[["ASC_P"]], "/", actual_counts[["ASC_N"]], "\n")
cat("Hallmark pathways:", nrow(camera_compare), "\n")
cat("CAMERA significant ASC_P/ASC_N:", sum(camera_p$FDR < 0.05), "/", sum(camera_n$FDR < 0.05), "\n")
print(camera_summary, row.names = FALSE)
if (dry_run) quit(save = "no", status = 0)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_volcano <- function(deg, subset_name, pdf_name, png_name) {
  up_n <- sum(deg$direction == "Up")
  down_n <- sum(deg$direction == "Down")
  top_up <- head(deg$Gene[deg$direction == "Up"][order(deg$FDR[deg$direction == "Up"])], 8)
  top_down <- head(deg$Gene[deg$direction == "Down"][order(deg$FDR[deg$direction == "Down"])], 8)
  labels <- unique(c(top_up, top_down))
  plot_df <- deg
  plot_df$y <- -log10(pmax(plot_df$FDR, .Machine$double.xmin))
  plot_df$direction <- factor(plot_df$direction, levels = c("Up", "Down", "NS"))
  p <- ggplot(plot_df, aes(x = log2FC, y = y)) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1.6) +
    scale_color_manual(
      values = c(Up = "#D95F02", Down = "#1B9E77", NS = "grey85"),
      labels = c(Up = paste0("Up (", up_n, ")"), Down = paste0("Down (", down_n, ")"), NS = "NS"),
      drop = FALSE
    ) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(
      data = plot_df[plot_df$Gene %in% labels, , drop = FALSE],
      aes(label = Gene), size = 3.5, fontface = "bold",
      box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf, min.segment.length = 0
    ) +
    labs(
      title = paste0(subset_name, ": Post vs Pre coculture response"),
      subtitle = "limma; n = 3 per timepoint; BH FDR < 0.05 and |log2FC| >= 0.5",
      x = "log2 fold change (Post / Pre)",
      y = expression(-log[10](FDR)),
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  ggsave(file.path(output_dir, pdf_name), p, width = 10, height = 8, device = cairo_pdf)
  ggsave(file.path(output_dir, png_name), p, width = 10, height = 8, dpi = 300, bg = "white")
}

plot_volcano(
  deg_p, "ASC_P",
  "FigSxA_ASC_P_Post_vs_Pre_volcano.pdf", "FigSxA_ASC_P_Post_vs_Pre_volcano.png"
)
plot_volcano(
  deg_n, "ASC_N",
  "FigSxB_ASC_N_Post_vs_Pre_volcano.pdf", "FigSxB_ASC_N_Post_vs_Pre_volcano.png"
)

# Two-column CAMERA heatmap across all Hallmark pathways.
camera_long <- rbind(
  data.frame(
    Pathway = camera_compare$Pathway, Contrast = "ASC_P Post vs Pre",
    Direction = camera_compare$Direction_ASC_P, FDR = camera_compare$FDR_ASC_P,
    score = pmax(-20, pmin(20, camera_compare$signed_log10_FDR_ASC_P)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Pathway = camera_compare$Pathway, Contrast = "ASC_N Post vs Pre",
    Direction = camera_compare$Direction_ASC_N, FDR = camera_compare$FDR_ASC_N,
    score = pmax(-20, pmin(20, camera_compare$signed_log10_FDR_ASC_N)),
    stringsAsFactors = FALSE
  )
)
pathway_order <- camera_compare$Pathway[order(camera_compare$min_FDR, decreasing = TRUE)]
camera_long$Pathway <- factor(camera_long$Pathway, levels = pathway_order)
camera_long$Contrast <- factor(camera_long$Contrast, levels = c("ASC_P Post vs Pre", "ASC_N Post vs Pre"))
camera_long$label <- ifelse(camera_long$FDR < 0.05, "*", "")

p_camera <- ggplot(camera_long, aes(x = Contrast, y = Pathway, fill = score)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 3.2, fontface = "bold") +
  scale_fill_gradient2(
    low = "#1B9E77", mid = "white", high = "#D95F02", midpoint = 0,
    limits = c(-20, 20), breaks = c(-20, -10, 0, 10, 20),
    name = "Signed\n-log10(FDR)"
  ) +
  labs(
    title = "Hallmark CAMERA: complete coculture response of ASC subsets",
    subtitle = "Orange: upregulated; green: downregulated; * BH FDR < 0.05; color capped at |20|",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )
ggsave(file.path(output_dir, "FigSxE_ASC_P_ASC_N_Hallmark_CAMERA_comparison.pdf"), p_camera, width = 7.5, height = 11, device = cairo_pdf)
ggsave(file.path(output_dir, "FigSxE_ASC_P_ASC_N_Hallmark_CAMERA_comparison.png"), p_camera, width = 7.5, height = 11, dpi = 300, bg = "white")

# Shared pathway panel: pathways significant in both contrasts with the same direction.
shared_camera <- camera_compare[camera_compare$classification %in% c("shared_up", "shared_down"), , drop = FALSE]
if (nrow(shared_camera)) {
  shared_camera$shared_score <- ifelse(
    shared_camera$classification == "shared_up", 1, -1
  ) * pmin(-log10(pmax(shared_camera$FDR_ASC_P, .Machine$double.xmin)),
           -log10(pmax(shared_camera$FDR_ASC_N, .Machine$double.xmin)))
  shared_camera <- shared_camera[order(abs(shared_camera$shared_score), decreasing = TRUE), , drop = FALSE]
  shared_plot <- head(shared_camera, 15)
  shared_plot$Pathway_label <- gsub("^HALLMARK_", "", shared_plot$Pathway)
  shared_plot$Pathway_label <- gsub("_", " ", shared_plot$Pathway_label)
  shared_plot$Pathway_label <- factor(shared_plot$Pathway_label, levels = rev(shared_plot$Pathway_label))
  p_shared <- ggplot(shared_plot, aes(x = shared_score, y = Pathway_label, fill = classification)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c(shared_up = "#D95F02", shared_down = "#1B9E77"),
                      labels = c(shared_up = "Shared up", shared_down = "Shared down"), name = NULL) +
    labs(
      title = "Shared Hallmark response of ASC_P and ASC_N",
      subtitle = "Conservative score uses the weaker -log10(FDR) of the two contrasts",
      x = "Signed shared -log10(FDR)", y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  ggsave(file.path(output_dir, "FigSxF_Shared_Hallmark_CAMERA.pdf"), p_shared, width = 7.5, height = 6.5, device = cairo_pdf)
  ggsave(file.path(output_dir, "FigSxF_Shared_Hallmark_CAMERA.png"), p_shared, width = 7.5, height = 6.5, dpi = 300, bg = "white")
}

wb <- createWorkbook()
addWorksheet(wb, "README")
readme <- data.frame(
  item = c("method", "input", "design", "contrasts", "gene_sets", "multiple_testing", "direction", "shared_definition"),
  value = c(
    "limma CAMERA competitive gene-set testing on the full expression matrix",
    "fractional normalized counts transformed as log2(count+1)",
    "unpaired subtype-by-timepoint model; n=3 per subtype/timepoint group; no donor/batch covariate available",
    "ASC_P Post vs Pre; ASC_N Post vs Pre",
    "MSigDB Hallmark via msigdbr, Mus musculus ortholog mapping matching the existing project",
    "Benjamini-Hochberg FDR separately within each contrast across 50 Hallmark pathways",
    "CAMERA Direction Up/Down refers to Post relative to Pre",
    "FDR < 0.05 in both contrasts with the same Direction"
  ), stringsAsFactors = FALSE
)
writeData(wb, "README", readme)
addWorksheet(wb, "ASC_P_CAMERA"); writeData(wb, "ASC_P_CAMERA", camera_p, withFilter = TRUE)
addWorksheet(wb, "ASC_N_CAMERA"); writeData(wb, "ASC_N_CAMERA", camera_n, withFilter = TRUE)
addWorksheet(wb, "pathway_comparison"); writeData(wb, "pathway_comparison", camera_compare, withFilter = TRUE)
addWorksheet(wb, "classification_summary"); writeData(wb, "classification_summary", camera_summary)
for (sheet in names(wb)) {
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1:30, widths = "auto")
}
saveWorkbook(wb, file.path(output_dir, "Table_St_ASC_P_ASC_N_Hallmark_CAMERA.xlsx"), overwrite = TRUE)
write.table(camera_p, file.path(output_dir, "Table_St_ASC_P_Hallmark_CAMERA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(camera_n, file.path(output_dir, "Table_St_ASC_N_Hallmark_CAMERA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(camera_compare, file.path(output_dir, "Table_St_ASC_P_ASC_N_Hallmark_CAMERA_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

top_shared <- shared_camera
top_shared <- top_shared[order(top_shared$min_FDR), , drop = FALSE]
top_shared_names <- head(gsub("_", " ", sub("^HALLMARK_", "", top_shared$Pathway)), 8)
report <- c(
  "# ASC Post-vs-Pre volcano and CAMERA run report",
  "",
  "## Analysis",
  "",
  sprintf("- ASC_P Post vs Pre: %d significant genes (%d up; %d down).", actual_counts[["ASC_P"]], sum(deg_p$direction == "Up"), sum(deg_p$direction == "Down")),
  sprintf("- ASC_N Post vs Pre: %d significant genes (%d up; %d down).", actual_counts[["ASC_N"]], sum(deg_n$direction == "Up"), sum(deg_n$direction == "Down")),
  sprintf("- CAMERA significant pathways: ASC_P %d; ASC_N %d.", sum(camera_p$FDR < 0.05), sum(camera_n$FDR < 0.05)),
  paste0("- Pathway classification: ", paste(paste0(camera_summary$classification, "=", camera_summary$pathway_count), collapse = "; "), "."),
  paste0("- Top shared same-direction Hallmarks: ", if (length(top_shared_names)) paste(top_shared_names, collapse = ", ") else "none at FDR < 0.05", "."),
  "",
  "## Statistical notes",
  "",
  "- Biological/statistical unit: sample/replicate column, assumed biological as in the existing project; donor identity is not documented.",
  "- Model: established unpaired limma subtype-by-timepoint model; no batch/donor covariate was available.",
  "- DEG threshold: BH FDR < 0.05 and |log2FC| >= 0.5.",
  "- CAMERA uses the full ranked expression evidence and is complementary to the shared/unique DEG-list ORA already provided.",
  "",
  "## Outputs",
  "",
  "- Figure SxA: ASC_P Post-vs-Pre volcano.",
  "- Figure SxB: ASC_N Post-vs-Pre volcano.",
  "- Figure SxE: two-contrast Hallmark CAMERA heatmap.",
  "- Figure SxF: Hallmark pathways significant in both subsets with the same direction.",
  "- Table St: complete CAMERA results and pathway classification.",
  "",
  "## Rerun commands",
  "",
  "```powershell",
  "& 'E:\\R\\R-4.3.3\\bin\\x64\\Rscript.exe' --vanilla 'D:\\R\\20260126 ASC bulk  SCI\\scripts\\reviewer2_ASC_bulk_volcano_CAMERA.R' --dry-run",
  "& 'E:\\R\\R-4.3.3\\bin\\x64\\Rscript.exe' --vanilla 'D:\\R\\20260126 ASC bulk  SCI\\scripts\\reviewer2_ASC_bulk_volcano_CAMERA.R'",
  "```"
)
writeLines(report, file.path(output_dir, "ASC_POST_PRE_VOLCANO_CAMERA_RUN_REPORT.md"), useBytes = TRUE)

required_outputs <- c(
  "FigSxA_ASC_P_Post_vs_Pre_volcano.pdf", "FigSxA_ASC_P_Post_vs_Pre_volcano.png",
  "FigSxB_ASC_N_Post_vs_Pre_volcano.pdf", "FigSxB_ASC_N_Post_vs_Pre_volcano.png",
  "FigSxE_ASC_P_ASC_N_Hallmark_CAMERA_comparison.pdf", "FigSxE_ASC_P_ASC_N_Hallmark_CAMERA_comparison.png",
  "Table_St_ASC_P_ASC_N_Hallmark_CAMERA.xlsx", "Table_St_ASC_P_Hallmark_CAMERA.tsv",
  "Table_St_ASC_N_Hallmark_CAMERA.tsv", "Table_St_ASC_P_ASC_N_Hallmark_CAMERA_comparison.tsv",
  "ASC_POST_PRE_VOLCANO_CAMERA_RUN_REPORT.md"
)
if (nrow(shared_camera)) {
  required_outputs <- c(required_outputs, "FigSxF_Shared_Hallmark_CAMERA.pdf", "FigSxF_Shared_Hallmark_CAMERA.png")
}
missing_outputs <- required_outputs[!file.exists(file.path(output_dir, required_outputs))]
if (length(missing_outputs)) stop("Missing required outputs: ", paste(missing_outputs, collapse = ", "))
if (any(file.info(file.path(output_dir, required_outputs))$size == 0)) stop("One or more required outputs are empty")

cat("Volcano and CAMERA package generated at:\n", output_dir, "\n")
print(camera_summary, row.names = FALSE)

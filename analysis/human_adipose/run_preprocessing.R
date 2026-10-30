#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
value <- function(k, default = NULL) {
  i <- which(args == k)
  if (!length(i)) return(default)
  if (length(i) != 1L || i == length(args) || startsWith(args[i + 1L], "--")) stop("Invalid argument: ", k)
  args[i + 1L]
}
if (!length(args) || "--help" %in% args) {
  cat("Usage: Rscript run_preprocessing.R --stage read|qc|integrate --out NEW_OUTPUT [--counts-dir SIX_SAMPLES | --previous CHECKPOINT.rds --reviewed-input] [--preflight]\n")
  quit(status = if ("--help" %in% args) 0L else 2L)
}
stage <- value("--stage")
if (is.null(stage) || !stage %in% c("read", "qc", "integrate")) stop("Select read, qc or integrate")
if (is.null(value("--out"))) stop("--out is required")
out <- normalizePath(value("--out"), winslash = "/", mustWork = FALSE)
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty output directory")
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE))
if (stage == "read") {
  if (is.null(value("--counts-dir"))) stop("--counts-dir is required")
  input <- normalizePath(value("--counts-dir"), winslash = "/", mustWork = TRUE)
  Sys.setenv(ASC_ADIPOSE_COUNTS_ROOT = input)
  if (tolower(input) == tolower(out)) stop("Keep input and output separate")
} else {
  if (is.null(value("--previous")) || !"--reviewed-input" %in% args)
    stop("Supply a reviewed checkpoint with --previous and --reviewed-input")
  input <- normalizePath(value("--previous"), winslash = "/", mustWork = TRUE)
  if (tolower(dirname(input)) == tolower(out)) stop("Keep input and output separate")
}
load_config <- function() {
  e <- new.env(parent = baseenv())
  sys.source(file.path(root, "config", "preprocessing_parameters.R"), envir = e)
  e$CONFIG
}
project_root <- function() root
assert_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install dependencies separately: ", paste(missing, collapse = ", "))
}
dependencies <- switch(stage,
  read = c("Matrix", "dplyr", "ggplot2", "patchwork", "data.table"),
  qc = c("Matrix", "dplyr", "ggplot2", "patchwork", "SingleCellExperiment", "SummarizedExperiment", "scDblFinder", "tibble"),
  integrate = c("Seurat", "harmony", "clustree", "ggraph", "dplyr", "ggplot2", "patchwork"))
assert_packages(dependencies)
cfg_check <- load_config()
if (stage == "read") {
  paths <- file.path(input, cfg_check$samples$folder, "filtered_feature_bc_matrix")
  if (!all(dir.exists(paths))) stop("Missing one or more of the six count-matrix directories")
}
if ("--preflight" %in% args) { message("Preflight PASS; no counts read or analysis executed"); quit(status = 0L) }
for (pkg in dependencies) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
for (subdir in c("figures", "tables", "objects", "reports", "logs")) dir.create(file.path(out, subdir), recursive = TRUE, showWarnings = FALSE)
Sys.setenv(ASC_P_ATTEMPT_DIR = out)
read_previous_checkpoint <- function(step) readRDS(input)
save_plot_both <- function(plot, stem, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}
write_waiting_review <- function(attempt_dir, step, lines) {
  writeLines(c(paste("# Stage", step, "review required"), lines, "No approval has been generated."),
             file.path(attempt_dir, "reports", "review_required.md"))
}
write_step_manifest <- function(step, attempt_dir, checkpoint, script_path, extra = list()) {
  capture.output(list(step = step, checkpoint = checkpoint, source = script_path, review_required = TRUE, extra = extra),
                 file = file.path(attempt_dir, "reports", "replay_provenance.txt"))
}
scripts <- c(read = "01_read_count_matrices.R", qc = "02_filter_qc_and_doublets.R", integrate = "03_integrate_and_lock_clusters.R")
sys.source(file.path(root, "stages", scripts[[stage]]), envir = new.env(parent = globalenv()))
capture.output(sessionInfo(), file = file.path(out, "logs", "sessionInfo.txt"))

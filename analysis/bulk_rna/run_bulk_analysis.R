#!/usr/bin/env Rscript
# Interface only: original computations are in provenance-linked stage files.
args <- commandArgs(trailingOnly = TRUE)
usage <- function() cat(paste(
  "Usage: Rscript analysis/bulk_rna/run_bulk_analysis.R",
  "--input EXPRESSION.csv --out NEW_OUTPUT_DIRECTORY",
  "[--stage model|all] [--preflight]\n"))
if ("--help" %in% args || !length(args)) {
  usage()
  quit(status = if ("--help" %in% args) 0L else 2L)
}
value <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (!length(idx)) return(default)
  if (length(idx) != 1L || idx == length(args) || startsWith(args[idx + 1L], "--"))
    stop("Supply exactly one value for ", flag)
  args[idx + 1L]
}
stage <- value("--stage", "model")
if (!stage %in% c("model", "all")) stop("--stage must be model or all")
input_arg <- value("--input")
out_arg <- value("--out")
if (is.null(input_arg) || is.null(out_arg)) stop("--input and --out are required")
input_file <- normalizePath(input_arg, winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(out_arg, winslash = "/", mustWork = FALSE)
if (identical(tolower(dirname(input_file)), tolower(out_dir)))
  stop("Output directory must differ from the input directory")
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
  stop("Output directory must be new or empty: ", out_dir)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this entry point using Rscript")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
dependencies <- c("tidyverse", "limma", "tibble", "dplyr")
if (stage == "all") dependencies <- unique(c(dependencies, "pheatmap", "msigdbr", "ggplot2", "ggrepel"))
missing <- dependencies[!vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install dependencies separately before running: ", paste(missing, collapse = ", "))
header <- names(read.csv(input_file, nrows = 0L, check.names = FALSE))
samples <- c(paste0("ASC_N_", 1:3), paste0("ASC_P_", 1:3),
             paste0("after_ASC_N_", 1:3), paste0("after_ASC_P_", 1:3))
required <- c("Gene ID", "Gene Symbol", paste("Normalized read count", samples))
if (!setequal(header, required)) stop("Input header does not match the documented 12-sample normalized-count file")
if ("--preflight" %in% args) {
  cat("Preflight PASS: input header, dependencies and output isolation checked. No analysis executed.\n")
  quit(status = 0L)
}
stages <- c("01_prepare_and_fit_interaction_model.R", "02_plot_interaction_heatmap.R",
            "03_test_hallmark_camera.R", "04_plot_interaction_volcano.R",
            "05_fit_cell_cycle_sensitivity.R", "06_plot_sensitivity_results.R",
            "07_plot_pca_variance.R")
if (stage == "model") stages <- stages[1L]
stopifnot(all(file.exists(file.path(script_dir, "stages", stages))))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
original_wd <- getwd()
setwd(out_dir) # Explicitly configured output, never an inherited workstation path.
analysis_env <- new.env(parent = globalenv())
analysis_env$input_file <- input_file
tryCatch({
  for (name in stages) {
    message("Running ", name)
    sys.source(file.path(script_dir, "stages", name), envir = analysis_env)
  }
  capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
}, finally = setwd(original_wd))
message("Completed configured stages in: ", out_dir)

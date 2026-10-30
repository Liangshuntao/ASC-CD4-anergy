#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
value <- function(k, default = NULL) {
  i <- which(args == k)
  if (!length(i)) return(default)
  if (length(i) != 1L || i == length(args) || startsWith(args[i + 1L], "--")) stop("Invalid value for ", k)
  args[i + 1L]
}
if ("--help" %in% args || !length(args)) {
  cat("Usage: Rscript run_clinical_analysis.R --stage tcga|metabric --out NEW_OUTPUT [--tcga-object PREPARED.rds --subtypes ORIGINAL.tsv | --metabric-dir ORIGINAL_EXPORT] [--preflight]\n")
  quit(status = if ("--help" %in% args) 0L else 2L)
}
stage <- value("--stage")
if (!stage %in% c("tcga", "metabric")) stop("Select --stage tcga or metabric")
if (is.null(value("--out"))) stop("--out is required")
out <- normalizePath(value("--out"), winslash = "/", mustWork = FALSE)
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new/empty output directory")
packages <- c("GSVA", "survival", "survminer", "dplyr", "ggpubr", "ggplot2", "patchwork")
if (stage == "tcga") packages <- c(packages, "SummarizedExperiment", "cowplot", "pheatmap", "clusterProfiler", "msigdbr", "limma", "tidyr", "enrichplot")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install dependencies separately: ", paste(missing, collapse = ", "))
env <- new.env(parent = globalenv())
if (stage == "tcga") {
  if (is.null(value("--tcga-object")) || is.null(value("--subtypes"))) stop("Supply --tcga-object and --subtypes")
  object <- normalizePath(value("--tcga-object"), mustWork = TRUE)
  env$subtype_file <- normalizePath(value("--subtypes"), mustWork = TRUE)
  if (tolower(out) %in% tolower(c(dirname(object), dirname(env$subtype_file)))) stop("Keep output separate from input directories")
  scripts <- c("01_score_tcga_cohort.R", "02_analyze_tcga_survival_and_correlations.R",
               "03_analyze_tcga_pathway_landscape.R", "06_analyze_custom_composite_score.R", "07_analyze_target_expression.R")
} else {
  if (is.null(value("--metabric-dir"))) stop("Supply --metabric-dir")
  env$metabric_dir <- normalizePath(value("--metabric-dir"), winslash = "/", mustWork = TRUE)
  needed <- file.path(env$metabric_dir, c("data_clinical_patient.txt", "data_clinical_sample.txt", "data_mrna_illumina_microarray.txt"))
  if (!all(file.exists(needed))) stop("METABRIC export lacks required files")
  if (tolower(out) == tolower(env$metabric_dir)) stop("Keep output separate from input export")
  scripts <- c("04_prepare_metabric_scores.R", "05_analyze_metabric_disease_specific_survival.R")
}
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE))
stopifnot(all(file.exists(file.path(root, "stages", scripts))))
if ("--preflight" %in% args) {
  message("Preflight PASS: paths, required files and dependencies; no scoring or models executed")
  quit(status = 0L)
}
for (pkg in packages) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
if (stage == "tcga") env$brca_data <- readRDS(object)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
oldwd <- getwd(); setwd(out)
tryCatch({
  for (script in scripts) sys.source(file.path(root, "stages", script), envir = env)
  capture.output(sessionInfo(), file = "sessionInfo.txt")
}, finally = {graphics.off(); setwd(oldwd)})

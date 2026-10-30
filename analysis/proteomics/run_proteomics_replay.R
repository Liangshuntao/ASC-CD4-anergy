#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || "--help" %in% args) {
  cat("Usage: Rscript run_proteomics_replay.R --checkpoint STATE.rds --out NEW_OUTPUT [--stage model|all] [--asc-pathways ORIGINAL.csv] [--preflight]\n")
  quit(status = if ("--help" %in% args) 0L else 2L)
}
value <- function(k, default = NULL) {
  i <- which(args == k)
  if (!length(i)) return(default)
  if (length(i) != 1L || i == length(args) || startsWith(args[i + 1L], "--")) stop("Invalid argument: ", k)
  args[i + 1L]
}
if (is.null(value("--checkpoint")) || is.null(value("--out"))) stop("--checkpoint and --out are required")
input <- normalizePath(value("--checkpoint"), winslash = "/", mustWork = TRUE)
out <- normalizePath(value("--out"), winslash = "/", mustWork = FALSE)
if (tolower(dirname(input)) == tolower(out)) stop("Keep output separate from checkpoint input")
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty output directory")
mode <- value("--stage", "model")
if (!mode %in% c("model", "all")) stop("--stage must be model or all")
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE))
dependencies <- c("limma", "dplyr", "tibble")
if (mode == "all") dependencies <- c(dependencies, "ggplot2", "pheatmap")
missing <- dependencies[!vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install dependencies separately: ", paste(missing, collapse = ", "))
state <- readRDS(input)
required <- c("expr_imp_asc", "expr_imp_cd4", "sample_info_asc", "sample_info_cd4",
              "ann_asc", "ann_cd4", "asc_gsva", "asc_gsva_Hallmark", "asc_gsva_kegg",
              "cd4_gsva", "expr_asc_gene_median", "expr_cd4_for_gsva", "annot_col")
if (!is.list(state) || !all(required %in% names(state))) stop("Checkpoint lacks required published-state objects")
stopifnot(ncol(state$expr_imp_asc) == 6L, ncol(state$expr_imp_cd4) == 6L,
          identical(colnames(state$expr_imp_asc), as.character(state$sample_info_asc$sample)),
          identical(colnames(state$expr_imp_cd4), as.character(state$sample_info_cd4$sample)))
pathways <- normalizePath(value("--asc-pathways", file.path(root, "metadata", "stromal_selected_pathways.csv")), mustWork = TRUE)
if ("--preflight" %in% args) {
  message("Checkpoint/preflight PASS; no imputation or scientific analysis executed")
  quit(status = 0L)
}
for (pkg in dependencies) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
stages <- sort(list.files(file.path(root, "stages"), pattern = "\\.R$", full.names = TRUE))
if (mode == "model") stages <- stages[1L]
env <- list2env(state, parent = globalenv())
env$asc_pathways_file <- pathways
dir.create(out, recursive = TRUE, showWarnings = FALSE)
oldwd <- getwd(); setwd(out)
tryCatch({
  for (stage in stages) sys.source(stage, envir = env)
  capture.output(sessionInfo(), file = "sessionInfo.txt")
}, finally = { graphics.off(); setwd(oldwd) })

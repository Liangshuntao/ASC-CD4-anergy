#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
value <- function(k) {
  i <- which(args == k)
  if (length(i) != 1L || i == length(args) || startsWith(args[i + 1L], "--")) stop("Required: ", k)
  args[i + 1L]
}
if (!length(args) || "--help" %in% args) {
  cat("Usage: Rscript run_breast_analysis.R --data-dir AUTHOR_OBJECTS --normal-object REVIEWED_NORMAL.rds --formal-object REVIEWED_TNK.rds --out NEW_OUTPUT [--preflight]\n")
  quit(status = if ("--help" %in% args) 0L else 2L)
}
data_dir <- normalizePath(value("--data-dir"), winslash = "/", mustWork = TRUE)
normal <- normalizePath(value("--normal-object"), winslash = "/", mustWork = TRUE)
formal <- normalizePath(value("--formal-object"), winslash = "/", mustWork = TRUE)
out <- normalizePath(value("--out"), winslash = "/", mustWork = FALSE)
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty output directory")
if (tolower(out) %in% tolower(c(data_dir, dirname(normal), dirname(formal)))) stop("Keep output separate from inputs")
required_objects <- file.path(data_dir, paste0("SeuratObject_", c("NormTotal", "ERTotal", "HER2", "TNBC"), ".rds"))
if (!all(file.exists(required_objects))) stop("Missing one or more of the four author Seurat objects")
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE))
metadata <- list.files(file.path(root, "metadata"), pattern = "\\.tsv$", full.names = TRUE)
if (length(metadata) != 7L) stop("Expected seven reviewed metadata/reference tables")
deps <- c("Seurat", "SeuratObject", "Matrix", "GSVA", "ggplot2", "patchwork", "scales", "digest", "jsonlite", "ggtext", "png", "grid")
missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install dependencies separately: ", paste(missing, collapse = ", "))
scripts <- c("01_plot_normal_breast_candidates.R", "02_plot_er_positive_candidates.R",
             "03_plot_her2_candidates.R", "04_plot_triple_negative_candidates.R",
             "05_finalize_cohort_candidate_panels.R", "06_compute_patient_pam50_resemblance.R",
             "07_score_reviewed_cd4_states.R", "08_plot_patient_stromal_tcell_associations.R",
             "09_assemble_cohort_panels.R")
if (!all(file.exists(file.path(root, "stages", scripts)))) stop("Missing workflow stage")
if ("--preflight" %in% args) {
  message("Preflight PASS: explicit objects, reviewed metadata, dependencies and stage references. No single-cell computation executed.")
  quit(status = 0L)
}
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
if (!all(file.copy(metadata, file.path(out, "tables"), overwrite = FALSE))) stop("Failed to seed reviewed metadata")
Sys.setenv(ASC_BREAST_WORK_ROOT = out, ASC_BREAST_DATA_ROOT = data_dir,
           ASC_BREAST_NORMAL_OBJECT = normal, ASC_BREAST_FORMAL_OBJECT = formal,
           GSE161529_ASC_P_PROJECT_ROOT = out)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (script in scripts) {
  cmd <- c("--vanilla", shQuote(file.path(root, "stages", script)))
  if (script == "05_finalize_cohort_candidate_panels.R") cmd <- c(cmd, "--cohort=all")
  result <- system2(rscript, cmd)
  if (result != 0L) stop("Stage failed: ", script, "; preserve output for diagnosis and start a new output on retry")
}
message("Human breast evidence workflow completed in: ", out)

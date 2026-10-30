#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || !length(args)) {
  cat("Usage: Rscript analysis/bulk_rna/run_bulk_supplement.R --input EXPRESSION.csv --did-table PRIMARY_DID.csv --out NEW_OUTPUT_DIRECTORY [--preflight]\n")
  quit(status = if ("--help" %in% args) 0L else 2L)
}
value <- function(k) {
  i <- which(args == k)
  if (length(i) != 1L || i == length(args) || startsWith(args[i + 1L], "--")) stop("Required: ", k)
  args[i + 1L]
}
input <- normalizePath(value("--input"), mustWork = TRUE)
did <- normalizePath(value("--did-table"), mustWork = TRUE)
out <- normalizePath(value("--out"), mustWork = FALSE)
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty output directory")
if (tolower(out) %in% tolower(c(dirname(input), dirname(did)))) stop("Keep output separate from input directories")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- dirname(normalizePath(script, mustWork = TRUE))
stages <- c("01_fit_within_subset_responses.R", "02_summarize_response_overlap.R",
            "03_plot_response_volcano_and_camera.R")
preflight <- "--preflight" %in% args
if (preflight) stages <- stages[1L]
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (stage in stages) {
  cmd <- c("--vanilla", shQuote(file.path(root, "supplement", stage)),
           "--input", shQuote(input), "--out", shQuote(out), "--did-table", shQuote(did))
  if (preflight) cmd <- c(cmd, "--dry-run")
  status <- system2(rscript, cmd)
  if (status != 0L) stop("Supplement stage failed: ", stage)
}
if (preflight) message("First-stage input/model preflight passed; downstream checks require generated contrast tables")

# Explicit execution interface; statistical calculations below are unchanged.
if ("--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat("Usage: Rscript SCRIPT --input ANNOTATED_OBJECT.rds --out OUTPUT_DIRECTORY [--threads 4] [--preflight] [--resume]\n")
  quit(status = 0L)
}
if (!all(c("--input", "--out") %in% commandArgs(trailingOnly = TRUE)))
  stop("Both --input and --out must be supplied explicitly")
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}

input <- normalizePath(value_after("--input", NULL), winslash = "/", mustWork = TRUE)
out <- normalizePath(value_after("--out", "output"), winslash = "/", mustWork = FALSE)
threads <- as.integer(value_after("--threads", "4"))
preflight <- "--preflight" %in% args
resume <- "--resume" %in% args

if (!is.finite(threads) || threads < 1L) stop("--threads must be a positive integer")
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)) && !resume && !preflight) {
  stop("Output directory is not empty. Use a new --out directory or add --resume: ", out)
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cmd_all <- commandArgs(FALSE)
script_arg <- grep("^--file=", cmd_all, value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "run_mouse_analysis.R"
root <- normalizePath(dirname(script_file), winslash = "/", mustWork = TRUE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")

run_r <- function(script, extra = character()) {
  status <- system2(rscript, c("--vanilla", shQuote(file.path(root, "stages", script)), extra))
  if (!identical(status, 0L)) stop("Stage failed: ", script, " (exit status ", status, ")")
}

common <- c("--input", shQuote(input), "--out", shQuote(out), "--threads", as.character(threads))
if (preflight) {
  run_r("01_compute_signatures_and_trajectory.R", c(common, "--preflight"))
  message("Submission preflight PASS")
  quit(save = "no", status = 0L)
}

core_done <- all(file.exists(file.path(out, c(
  "seu_cd4_nonTreg.rds", "tf_by_cell.rds", "tf_act_long.rds", "path_scores.rds",
  "seu_cd4_nonTreg_final.rds", "seu_LR_final.rds"
))))
strict_done <- file.exists(file.path(out, "panel_validation.csv"))
liana_done <- file.exists(file.path(out, "liana_fixed", "LIANA_FIX_REPORT.md"))

if (resume && core_done) {
  message("[resume] Core analysis checkpoint found; skipping stage 02")
} else {
  run_r("01_compute_signatures_and_trajectory.R", common)
}

if (resume && strict_done) {
  message("[resume] Strict figure checkpoint found; skipping stage 03")
} else {
  run_r("02_finalize_figure_panels.R", common)
}

if (resume && liana_done) {
  message("[resume] LIANA completion report found; skipping stage 04")
} else {
  run_r("03_predict_asc_cd4_communication.R", c("--input", shQuote(file.path(out, "seu_LR_final.rds")),
                                    "--out", shQuote(file.path(out, "liana_fixed"))))
}
message("Submission workflow completed: ", out)

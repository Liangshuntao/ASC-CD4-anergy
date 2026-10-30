#!/usr/bin/env Rscript
# Export only the published-state inputs, not the complete historical workspace.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || !identical(args[c(1L, 3L)], c("--workspace", "--out")))
  stop("Usage: Rscript export_proteomics_checkpoint.R --workspace ORIGINAL.RData --out NEW_STATE.rds")
workspace <- normalizePath(args[2L], mustWork = TRUE)
out <- normalizePath(args[4L], mustWork = FALSE)
if (file.exists(out)) stop("Refusing to overwrite checkpoint")
env <- new.env(parent = baseenv())
load(workspace, envir = env)
required <- c("expr_imp_asc", "expr_imp_cd4", "sample_info_asc", "sample_info_cd4",
              "ann_asc", "ann_cd4", "asc_gsva", "asc_gsva_Hallmark", "asc_gsva_kegg",
              "cd4_gsva", "expr_asc_gene_median", "expr_cd4_for_gsva", "annot_col")
if (!all(required %in% ls(env))) stop("Workspace lacks required published-state objects")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(mget(required, envir = env, inherits = FALSE), out)
message("Whitelisted checkpoint exported; original workspace unchanged")

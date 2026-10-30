#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || "--help" %in% args) {
  cat("Usage: Rscript apply_reviewed_annotation.R --input CLUSTERED.rds --map REVIEWED.csv --reviewed-map --out NEW_ANNOTATED.rds\n")
  quit(status = if ("--help" %in% args) 0L else 2L)
}
value <- function(k) {
  i <- which(args == k)
  if (length(i) != 1L || i == length(args) || startsWith(args[i + 1L], "--")) stop("Required: ", k)
  args[i + 1L]
}
if (!"--reviewed-map" %in% args) stop("Reviewed manual annotation is required; no labels are inferred")
input <- normalizePath(value("--input"), mustWork = TRUE)
mapping <- normalizePath(value("--map"), mustWork = TRUE)
out <- normalizePath(value("--out"), mustWork = FALSE)
if (file.exists(out)) stop("Refusing to overwrite an object")
if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Install SeuratObject separately")
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE))
source(file.path(root, "lib", "reviewed_annotation.R"))
obj <- readRDS(input)
meta <- obj[[]]
labels <- map_reviewed_labels(meta, read.csv(mapping, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM"))
obj$celltype_v2 <- labels$celltype_v2
obj$major_celltype_v2 <- labels$major_celltype_v2
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, out)
message("Reviewed labels imported; assays, embeddings, clusters and input file unchanged")

#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
source_rda <- file.path(project_root, "provenance", "pam50", "genefu", "data", "pam50.scale.rda")
output_tsv <- file.path(project_root, "tables", "PAM50_genefu_2.34.0_scale_centroids.tsv")

stopifnot(file.exists(source_rda))
env <- new.env(parent = emptyenv())
load(source_rda, envir = env)
model <- env[["pam50.scale"]]
stopifnot(identical(model$std, "scale"), identical(model$method.cor, "spearman"))
stopifnot(identical(dim(model$centroids), c(50L, 5L)))

centroids <- data.frame(
  gene = rownames(model$centroids),
  LumA = model$centroids[, "LumA"],
  LumB = model$centroids[, "LumB"],
  HER2 = model$centroids[, "Her2"],
  Basal = model$centroids[, "Basal"],
  Normal = model$centroids[, "Normal"],
  check.names = FALSE
)
write.table(centroids, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

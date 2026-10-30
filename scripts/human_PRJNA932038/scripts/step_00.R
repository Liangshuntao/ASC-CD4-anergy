cfg <- load_config(); root <- project_root(); attempt <- Sys.getenv("ASC_P_ATTEMPT_DIR")
assert_packages(c("digest", "yaml", "Seurat", "Matrix", "dplyr"))

samples <- cfg$samples
expected_cells <- c(S53_N = 11719L, S53_T = 10252L, S76_N = 8467L, S76_T = 11931L, S81_N = 10597L, S81_T = 10484L)
audit_rows <- list(); file_rows <- list()
for (i in seq_len(nrow(samples))) {
  d <- file.path(cfg$input_root, samples$folder[[i]], "filtered_feature_bc_matrix")
  paths <- file.path(d, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"))
  if (!all(file.exists(paths))) stop("Missing 10x files for ", samples$sample_id[[i]])
  features <- read.delim(gzfile(paths[[2]]), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- readLines(gzfile(paths[[3]]))
  con <- gzfile(paths[[1]], open = "rt")
  repeat { line <- readLines(con, n = 1L); if (!length(line)) stop("MTX header missing"); if (!startsWith(line, "%")) break }
  close(con)
  dims <- as.integer(strsplit(line, "[[:space:]]+")[[1]])
  if (length(dims) != 3L || dims[[1]] != nrow(features) || dims[[2]] != length(barcodes)) stop("Dimension mismatch: ", samples$sample_id[[i]])
  if (length(barcodes) != expected_cells[[samples$sample_id[[i]]]]) stop("Unexpected cell count: ", samples$sample_id[[i]])
  audit_rows[[i]] <- cbind(samples[i, ], data.frame(matrix_path = d, genes = dims[[1]], cells = dims[[2]], nnz = dims[[3]], stringsAsFactors = FALSE))
  file_rows[[i]] <- data.frame(sample_id = samples$sample_id[[i]], file = paths,
                               bytes = file.info(paths)$size,
                               sha256 = vapply(paths, sha256_file, character(1)), stringsAsFactors = FALSE)
}
audit <- dplyr::bind_rows(audit_rows); manifest <- dplyr::bind_rows(file_rows)
stopifnot(nrow(audit) == 6L, sum(audit$cells) == 63450L, all(audit$genes == 38606L))
write.csv(audit, file.path(attempt, "tables", "input_audit.csv"), row.names = FALSE)
write.csv(manifest, file.path(attempt, "tables", "input_file_manifest.csv"), row.names = FALSE)
write.csv(samples, file.path(attempt, "tables", "sample_mapping.csv"), row.names = FALSE)
locked <- file.path(root, "config", "locked_package_versions.csv")
if (!file.exists(locked)) stop("Run install_dependencies.R before Step 00")
checkpoint <- file.path(attempt, "objects", "step00_input_audit.rds")
saveRDS(list(audit = audit, manifest = manifest, samples = samples, package_versions = read.csv(locked)), checkpoint)
write_waiting_review(attempt, 0, c(
  "## Required checks", "",
  "- Confirm six paired samples and 63,450 nuclei.",
  "- Confirm all matrices contain 38,606 features.",
  "- Confirm S81_N is flagged for 42% reads-in-cells and later sensitivity analysis.",
  "- Confirm the project-local package lock is acceptable.",
  "- Confirm filtered-only matrices prevent independent empty-droplet/ambient correction."
))
write_step_manifest(0, attempt, checkpoint, file.path(root, "scripts", "step_00.R"),
                    extra = list(input_manifest_sha256 = sha256_file(file.path(attempt, "tables", "input_file_manifest.csv"))))

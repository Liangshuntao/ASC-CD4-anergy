source(file.path(normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE), "scripts", "utils_GSE161529_ASC_P.R"))
cfg <- source_project_config()
ensure_packages(c("digest", "jsonlite", "patchwork"))
set.seed(cfg$seed)

author_code_paths <- file.path(cfg$author_code_dir, cfg$author_code_files)
source_object_paths <- c(file.path(cfg$input_dir, cfg$local_objects$filename),
                         file.path(cfg$author_asset_dir, cfg$required_author_assets$filename))
state <- begin_step(cfg, "01_map_author_annotations", c(author_code_paths[file.exists(author_code_paths)],
                                                         source_object_paths[file.exists(source_object_paths)]))

derive_missing_sub_object <- function(row) {
  destination <- file.path(cfg$updated_object_dir,
                           sub("\\.rds$", "_derived_from_author_total.rds", row$filename, ignore.case = TRUE))
  if (file.exists(destination)) return(list(path = destination, excluded = character(), retained = NA_integer_, status = "already_derived"))
  total_row <- cfg$local_objects[cfg$local_objects$object_id == row$total_object_id, ]
  if (nrow(total_row) != 1L) stop("Cannot resolve author Total object for ", row$object_id, call. = FALSE)
  total_path <- file.path(cfg$input_dir, total_row$filename)
  code_path <- file.path(cfg$author_code_dir, row$code_file)
  excluded <- extract_author_microenvironment_exclusions(code_path)
  if (!length(excluded)) return(list(path = NA_character_, excluded = character(), retained = NA_integer_, status = "author_filter_unresolved"))
  obj <- load_seurat_object(total_path)
  meta <- safe_meta_data(obj)
  cluster_col <- resolve_meta_column(meta, cfg$cluster_candidates, TRUE, "author cluster")
  author_cluster_index <- as.integer(meta[[cluster_col]])
  keep <- rownames(meta)[!author_cluster_index %in% as.integer(excluded)]
  if (!length(keep) || length(keep) == nrow(meta)) {
    stop("Parsed author exclusion did not produce a proper non-epithelial subset for ", row$object_id, call. = FALSE)
  }
  note <- list(
    status = "derived_from_author_total_object",
    source_total_rds = normalizePath(total_path, winslash = "/"),
    source_total_sha256 = sha256_file(total_path),
    author_code = normalizePath(code_path, winslash = "/"),
    author_code_sha256 = sha256_file(code_path),
    author_filter = paste0("Sub <- !(Cluster %in% c(", paste(excluded, collapse = ","), "))"),
    retained_original_cluster_identity = TRUE,
    no_qc_normalization_integration_reduction_or_clustering_rerun = TRUE
  )
  derived <- subset_legacy_seurat_rna(obj, keep, note)
  derived$derivation_status <- "derived_from_author_total_object"
  derived$derivation_author_code <- basename(code_path)
  save_rds_safe(derived, destination, cfg$overwrite, compress = TRUE)
  retained <- length(keep)
  rm(obj, derived); invisible(gc())
  list(path = destination, excluded = excluded, retained = retained, status = "derived_from_author_total_object")
}

derivation_records <- list()
for (i in seq_len(nrow(cfg$required_author_assets))) {
  row <- cfg$required_author_assets[i, ]
  released <- file.path(cfg$author_asset_dir, row$filename)
  status <- if (is_valid_file(released, 1024^2)) "author_released_object" else "fallback_not_attempted"
  path <- if (status == "author_released_object") released else NA_character_
  excluded <- character(); retained <- NA_integer_; note <- NA_character_
  if (status != "author_released_object") {
    result <- tryCatch(derive_missing_sub_object(row), error = identity)
    if (inherits(result, "error")) {
      status <- "fallback_failed"; note <- conditionMessage(result)
    } else {
      status <- result$status; path <- result$path; excluded <- result$excluded; retained <- result$retained
      if (is.na(path)) note <- "Could not parse the author's explicit Total-cluster exclusion used to define the microenvironment subset."
    }
  }
  derivation_records[[length(derivation_records) + 1L]] <- data.frame(
    object_id = row$object_id, status = status, path = path, author_code = row$code_file,
    author_excluded_cluster_indices = paste(excluded, collapse = ";"), retained_cell_number = retained,
    note = note, stringsAsFactors = FALSE
  )
}
derivation_table <- do.call(rbind, derivation_records)
derivation_out <- file.path(cfg$table_dir, "01_sub_object_resolution.tsv")
write_tsv_safe(derivation_table, derivation_out, cfg$overwrite)

registry <- rbind(localization_registry(cfg), analysis_registry(cfg)[, c("analysis_id", "subtype", "path", "source_status")])
registry <- registry[!duplicated(registry$analysis_id), ]

suggest_marker_compartment <- function(means) {
  scores <- c(
    `CAF/fibroblast` = mean(means[intersect(names(means), cfg$fibroblast_markers)], na.rm = TRUE),
    `pericyte/mural` = mean(means[intersect(names(means), cfg$mural_markers)], na.rm = TRUE),
    endothelial = mean(means[intersect(names(means), cfg$endothelial_markers)], na.rm = TRUE),
    contamination = mean(means[intersect(names(means), cfg$contamination_markers)], na.rm = TRUE)
  )
  scores[!is.finite(scores)] <- -Inf
  if (all(scores == -Inf)) return(NA_character_)
  names(which.max(scores))
}

mapping_rows <- list(); cluster_evidence <- list(); plot_cluster <- list(); plot_celltype <- list(); blocking <- character()
for (i in seq_len(nrow(registry))) {
  item <- registry[i, ]
  if (!file.exists(item$path)) {
    blocking <- c(blocking, paste0(item$analysis_id, ": required object is unavailable."))
    next
  }
  obj <- load_seurat_object(item$path)
  meta <- safe_meta_data(obj)
  cluster_col <- resolve_meta_column(meta, cfg$cluster_candidates, TRUE, "author cluster")
  patient_col <- resolve_meta_column(meta, cfg$patient_candidates, TRUE, "patient")
  celltype_col <- resolve_meta_column(meta, cfg$celltype_candidates, FALSE)
  clusters <- sort(unique(as.character(meta[[cluster_col]])))
  marker_mat <- extract_gene_matrix(obj, cfg$all_audit_genes, "data")
  # Test HER2 before ER because the string "HER2" contains "ER".
  code_name <- if (grepl("Normal", item$analysis_id)) "NormTotal.R" else if (grepl("HER2", item$analysis_id)) "HER2.R" else if (grepl("ER", item$analysis_id)) "ER.R" else "TNBC.R"
  code_path <- file.path(cfg$author_code_dir, code_name)
  excluded_total_indices <- extract_author_microenvironment_exclusions(code_path)

  for (cl in clusters) {
    cells <- rownames(meta)[as.character(meta[[cluster_col]]) == cl]
    author_index <- unique(as.integer(meta[cells, cluster_col]))
    author_excluded <- grepl("_Total$|^NormalTotal$", item$analysis_id) &&
      length(author_index) == 1L && as.character(author_index) %in% excluded_total_indices
    label <- if (!is.na(celltype_col)) {
      tab <- sort(table(as.character(meta[cells, celltype_col])), decreasing = TRUE); names(tab)[1]
    } else if (identical(item$analysis_id, "NormalFib")) {
      "author fibroblast-only object"
    } else if (author_excluded) {
      "author Total cluster excluded from microenvironment Sub"
    } else "UNRESOLVED"
    compartment <- if (identical(item$analysis_id, "NormalFib")) "CAF/fibroblast" else if (author_excluded) "epithelial/tumor" else classify_compartment(label)
    means <- Matrix::rowMeans(marker_mat[, cells, drop = FALSE])
    supporting <- paste(names(sort(means, decreasing = TRUE))[seq_len(min(6L, length(means)))], collapse = ";")
    patient_tab <- sort(table(as.character(meta[cells, patient_col])), decreasing = TRUE)
    dominance <- if (sum(patient_tab)) max(patient_tab) / sum(patient_tab) else NA_real_
    mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
      analysis_id = item$analysis_id, subtype = item$subtype, source_rds = item$path,
      source_status = item$source_status, author_cluster = cl, author_cluster_one_based_index = paste(author_index, collapse = ";"),
      author_label = label, compartment = compartment,
      evidence_source = if (!is.na(celltype_col)) paste0("metadata:", celltype_col) else if (item$analysis_id == "NormalFib") "author NormTotalFib fibroblast-only object" else if (author_excluded) paste0("author_code:", code_name, ":Sub exclusion") else "manual_review_required",
      marker_only_suggestion_not_annotation = suggest_marker_compartment(means), supporting_markers = supporting,
      cell_number = length(cells), patient_number = length(patient_tab), largest_patient_fraction = dominance,
      single_patient_dominant = is.finite(dominance) && dominance > cfg$dominance_threshold,
      mixed_flag = compartment == "unknown/mixed", approved = FALSE, stringsAsFactors = FALSE
    )
  }

  if (validate_r_code(code_path)) {
    code_lines <- readLines(code_path, warn = FALSE, encoding = "UTF-8")
    evidence_lines <- grep("Sub <-|RenameIdents|celltype|cell.type|fibro|CAF|pericy|mural|endothel|subset|idents", code_lines, ignore.case = TRUE, value = TRUE)
    cluster_evidence[[item$analysis_id]] <- utils::head(trimws(evidence_lines), 80)
  }
  red <- reduction_name(obj)
  if (!is.na(red)) {
    plot_cluster[[item$analysis_id]] <- embedding_group_plot(obj, red, meta[[cluster_col]], paste(item$analysis_id, "author clusters"))
    rows_now <- do.call(rbind, mapping_rows)
    rows_now <- rows_now[rows_now$analysis_id == item$analysis_id, ]
    label_map <- setNames(rows_now$compartment, rows_now$author_cluster)
    plot_celltype[[item$analysis_id]] <- embedding_group_plot(obj, red, unname(label_map[as.character(meta[[cluster_col]])]), paste(item$analysis_id, "mapped compartment"))
  }
  rm(obj, meta, marker_mat); invisible(gc())
}

mapping <- if (length(mapping_rows)) do.call(rbind, mapping_rows) else data.frame()
if (!nrow(mapping)) blocking <- c(blocking, "No mapping rows could be produced.")
if (nrow(mapping) && any(mapping$compartment == "unknown/mixed")) blocking <- c(blocking, "One or more author clusters remain unknown/mixed and require manual evidence-based mapping.")
generated_mapping_file <- cfg$generated_mapping_file %||% file.path(cfg$table_dir, "01_author_cluster_annotation_map_generated.tsv")
write_tsv_safe(mapping, generated_mapping_file, cfg$overwrite)

outputs <- c(derivation_out, generated_mapping_file)
if (length(plot_cluster)) outputs <- c(outputs, save_plot_dual(patchwork::wrap_plots(plot_cluster, ncol = 2), file.path(cfg$figure_dir, "audit", "01_original_embedding_by_cluster"), cfg, 14, 10))
if (length(plot_celltype)) outputs <- c(outputs, save_plot_dual(patchwork::wrap_plots(plot_celltype, ncol = 2), file.path(cfg$figure_dir, "audit", "01_original_embedding_by_author_celltype"), cfg, 14, 10))

evidence_text <- unlist(lapply(names(cluster_evidence), function(nm) c(paste0("### ", nm), "", paste0("    ", cluster_evidence[[nm]]), "")))
mapping_report <- file.path(cfg$report_dir, "01_AUTHOR_ANNOTATION_MAPPING.md")
write_lines_safe(c("# Author annotation mapping", "",
                   "Resolved annotations use only released metadata, the author fibroblast-only object, or the author's explicit Total-cluster exclusion.",
                   "Marker-only suggestions are review aids and are not promoted to annotations automatically.",
                   "Unknown/mixed clusters remain blocked pending human review.", "", evidence_text), mapping_report, cfg$overwrite)
outputs <- c(outputs, mapping_report)

checkpoint <- file.path(cfg$report_dir, "REVIEW_CLUSTER_BOUNDARY_CHECKPOINT.md")
checkpoint_lines <- c("# Review cluster boundary checkpoint", "",
                      paste0("- Active mapping SHA256: `", sha256_file(cfg$mapping_file), "`"),
                      "- Candidate analysis must not start until every row has a resolved compartment and approval is granted.",
                      "- CAF/fibroblast is the only candidate denominator; mural and endothelial remain specificity references.", "",
                      "## Blocking issues", "", if (length(blocking)) paste0("- ", unique(blocking)) else "- None detected; human review is still required.", "",
                       "## Review table", "", paste0("Inspect `", basename(cfg$mapping_file), "`, including evidence source, marker-only suggestions, patient composition, dominance, and mixed flags."))
write_lines_safe(checkpoint_lines, checkpoint, cfg$overwrite)
outputs <- c(outputs, checkpoint)

if (!file.exists(cfg$approval_file)) {
  write_lines_safe(c("approved: false", "approved_by: \"\"", "approved_at: \"\"",
                     paste0("mapping_sha256: \"", sha256_file(cfg$mapping_file), "\""),
                     "rationale: \"Pending human review of author cluster boundaries\""), cfg$approval_file, FALSE)
  outputs <- c(outputs, cfg$approval_file)
}
if (length(blocking)) outputs <- c(outputs, write_waiting_report(cfg, unique(blocking)))
end_step(cfg, state, outputs, status = if (length(blocking)) "BLOCKED" else "SUCCESS", notes = unique(blocking))
if (length(blocking)) stop("Precheck is biologically blocked; review the waiting/checkpoint reports.", call. = FALSE)

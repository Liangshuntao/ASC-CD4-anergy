source(file.path(normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE), "scripts", "utils_GSE161529_ASC_P.R"))
cfg <- source_project_config()
ensure_packages(c("digest", "jsonlite"))

state <- begin_step(cfg, "01b_refine_author_cluster_mapping", c(cfg$mapping_file, file.path(cfg$author_code_dir, cfg$author_code_files)))
mapping <- read_mapping(cfg, require_resolved = FALSE)
if (nrow(mapping) != 77L) stop("Expected 77 current mapping rows; found ", nrow(mapping), call. = FALSE)

code_for_analysis <- function(analysis_id) {
  if (grepl("Normal", analysis_id)) return("NormTotal.R")
  # HER2 must be tested before ER: HER2 contains the substring ER.
  if (grepl("HER2", analysis_id)) return("HER2.R")
  if (grepl("ER", analysis_id)) return("ER.R")
  "TNBC.R"
}

is_total <- grepl("_Total$|^NormalTotal$", mapping$analysis_id)
mapping$revised_compartment <- as.character(mapping$compartment)
mapping$revised_evidence_source <- as.character(mapping$evidence_source)
mapping$revised_mixed_flag <- as.character(mapping$mixed_flag)
mapping$revised_approval_status <- "not_approved"
mapping$approved <- FALSE
mapping$revision_flags <- ""
mapping$corrected_author_code <- ""
mapping$corrected_author_excluded <- FALSE

for (i in seq_len(nrow(mapping))) {
  aid <- mapping$analysis_id[[i]]
  code_name <- code_for_analysis(aid)
  code_path <- file.path(cfg$author_code_dir, code_name)
  mapping$corrected_author_code[[i]] <- code_name
  exclusions <- extract_author_microenvironment_exclusions(code_path)
  one_based <- trimws(unlist(strsplit(as.character(mapping$author_cluster_one_based_index[[i]]), ";", fixed = TRUE)))[1]
  excluded <- is_total[[i]] && one_based %in% exclusions
  mapping$corrected_author_excluded[[i]] <- excluded

  flags <- character()
  if (excluded) {
    mapping$revised_compartment[[i]] <- "epithelial/tumor"
    mapping$revised_evidence_source[[i]] <- paste0("author_code:", code_name, ":Sub exclusion (corrected route)")
    mapping$revised_mixed_flag[[i]] <- "FALSE"
    if (aid == "HER2_Total" && grepl("author_code:ER\\.R", mapping$evidence_source[[i]])) {
      flags <- c(flags, "HER2_Total_wrong_ER_route_corrected")
    }
  } else {
    # Remove stale ER.R exclusion labels from HER2_Total retained clusters.
    if (aid == "HER2_Total" && grepl("author_code:ER\\.R", mapping$evidence_source[[i]])) {
      mapping$revised_compartment[[i]] <- "unknown/mixed"
      mapping$revised_evidence_source[[i]] <- "mapping_review:stale_ER_R_exclusion_removed"
      mapping$revised_mixed_flag[[i]] <- "TRUE"
      flags <- c(flags, "HER2_Total_stale_ER_R_exclusion_removed")
    }
    # The released Sub objects are author-defined non-epithelial microenvironments.
    # Epithelial signal inside them is a boundary/doublet flag, not a final label.
    if (!is_total[[i]] && mapping$compartment[[i]] == "epithelial/tumor") {
      mapping$revised_compartment[[i]] <- "unknown/mixed"
      mapping$revised_evidence_source[[i]] <- "mapping_review:epithelial_signal_inside_non_epithelial_Sub"
      mapping$revised_mixed_flag[[i]] <- "TRUE"
      flags <- c(flags, "epithelial_signal_inside_non_epithelial_Sub")
    }
    # contamination and immune/contamination are not final compartments.
    if (mapping$compartment[[i]] %in% c("contamination", "immune/contamination")) {
      mapping$revised_compartment[[i]] <- "unknown/mixed"
      mapping$revised_evidence_source[[i]] <- "mapping_review:noncanonical_contamination_label"
      mapping$revised_mixed_flag[[i]] <- "TRUE"
      flags <- c(flags, "noncanonical_contamination_label_requires_lineage_review")
    }
  }

  # The old HER2_Total rows 5-7 were incorrectly treated as ER exclusions.
  # Retained-cluster marker suggestions are carried forward only as provisional review aids.
  if (aid == "HER2_Total" && !excluded && one_based == "5" &&
      identical(mapping$marker_only_suggestion_not_annotation[[i]], "CAF/fibroblast")) {
    mapping$revised_compartment[[i]] <- "CAF/fibroblast"
    mapping$revised_evidence_source[[i]] <- "mapping_review:HER2.R_retained_cluster_plus_CAF_marker_support"
    mapping$revised_mixed_flag[[i]] <- "FALSE"
    flags <- c(flags, "HER2_Total_retained_CAF_candidate_requires_human_confirmation")
  }
  if (aid == "HER2_Total" && !excluded && one_based == "7" &&
      identical(mapping$marker_only_suggestion_not_annotation[[i]], "endothelial")) {
    mapping$revised_compartment[[i]] <- "endothelial"
    mapping$revised_evidence_source[[i]] <- "mapping_review:HER2.R_retained_cluster_plus_endothelial_marker_support"
    mapping$revised_mixed_flag[[i]] <- "FALSE"
    flags <- c(flags, "HER2_Total_retained_endothelial_candidate_requires_human_confirmation")
  }
  mapping$revision_flags[[i]] <- paste(unique(flags), collapse = ";")
}

mapping$mapping_status <- "revised_draft_not_approved"
out_table <- file.path(cfg$table_dir, "01_author_cluster_annotation_map_revised_v4.tsv")
write_tsv_safe(mapping, out_table, overwrite = cfg$overwrite)
queue <- mapping[mapping$revised_compartment == "unknown/mixed" | nzchar(mapping$revision_flags),
                 c("analysis_id", "author_cluster", "author_cluster_one_based_index", "revised_compartment",
                   "revised_evidence_source", "marker_only_suggestion_not_annotation", "supporting_markers",
                   "cell_number", "patient_number", "largest_patient_fraction", "single_patient_dominant",
                   "revision_flags", "mapping_status"), drop = FALSE]
queue_out <- file.path(cfg$table_dir, "01_mapping_human_review_queue_v2.tsv")
write_tsv_safe(queue, queue_out, overwrite = cfg$overwrite)

changes <- mapping[nzchar(mapping$revision_flags) | mapping$revised_compartment != mapping$compartment, , drop = FALSE]
report <- c(
  "# Revised author-cluster mapping draft (not approved)", "",
  "This file preserves the original CodeBuddy mapping and applies only deterministic route correction plus conservative boundary handling.", "",
  "## Applied rules", "",
  "- HER2 is tested before ER when selecting the author script.",
  "- Author Total-object exclusion clusters are hard evidence for epithelial/tumor removal.",
  "- `contamination` and `immune/contamination` are converted to `unknown/mixed`; they are not accepted final compartments.",
  "- Epithelial signal in released non-epithelial Sub objects is retained as a boundary flag, not promoted to a final epithelial annotation.",
  "- All revised rows remain `not_approved`; no candidate analysis may use this draft until a human signs the approval file.", "",
  paste0("Current rows: ", nrow(mapping), "; rows changed or flagged: ", nrow(changes)), "",
  "## HER2 correction", "",
  "HER2_Total must use HER2.R exclusions (one-based clusters 1, 4, 8). The prior mapping used ER.R because `grepl(\"ER\", \"HER2_Total\")` was evaluated first.", "",
  "## Changed/flagged rows", "",
  if (nrow(changes)) paste0("- ", changes$analysis_id, " cluster ", changes$author_cluster, ": ", changes$revision_flags, "; ", changes$compartment, " -> ", changes$revised_compartment) else "- None"
)
out_report <- file.path(cfg$report_dir, "01_MAPPING_REFINEMENT_DRAFT_v4.md")
write_lines_safe(report, out_report, overwrite = cfg$overwrite)
end_step(cfg, state, c(out_table, queue_out, out_report), status = "REVIEW_BLOCKED", notes = c("Revised draft only; original mapping and approval were not overwritten.", "Human cluster boundary approval remains required."))

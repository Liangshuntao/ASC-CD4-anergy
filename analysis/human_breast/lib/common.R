suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[1])) y else x

script_project_root <- function() .work_root

source_project_config <- function() {
  root <- script_project_root()
  old_root <- Sys.getenv("GSE161529_ASC_P_PROJECT_ROOT", unset = NA_character_)
  Sys.setenv(GSE161529_ASC_P_PROJECT_ROOT = root)
  on.exit({
    if (is.na(old_root)) Sys.unsetenv("GSE161529_ASC_P_PROJECT_ROOT")
    else Sys.setenv(GSE161529_ASC_P_PROJECT_ROOT = old_root)
  }, add = TRUE)
  source(file.path(.code_root, "config", "analysis_parameters.R"), local = .GlobalEnv)
  invisible(get("cfg", envir = .GlobalEnv))
}

ensure_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

ensure_dirs <- function(cfg) {
  dirs <- c(cfg$author_asset_dir, cfg$author_code_dir, cfg$updated_object_dir,
            cfg$table_dir, cfg$figure_dir, cfg$report_dir, cfg$log_dir,
            cfg$manifest_dir, dirname(cfg$candidate_cache))
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

assert_new_path <- function(path, overwrite = FALSE) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_tsv_safe <- function(x, path, overwrite = FALSE) {
  assert_new_path(path, overwrite)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  invisible(path)
}

write_lines_safe <- function(x, path, overwrite = FALSE) {
  assert_new_path(path, overwrite)
  writeLines(enc2utf8(x), path, useBytes = TRUE)
  invisible(path)
}

save_rds_safe <- function(x, path, overwrite = FALSE, compress = TRUE) {
  assert_new_path(path, overwrite)
  saveRDS(x, path, compress = compress)
  invisible(path)
}

sha256_file <- function(path) {
  ensure_packages("digest")
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

md5_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

begin_step <- function(cfg, step_id, inputs = character()) {
  ensure_dirs(cfg)
  list(step_id = step_id, started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
       inputs = inputs, input_sha256 = setNames(vapply(inputs, sha256_file, character(1)), inputs))
}

end_step <- function(cfg, state, outputs = character(), status = "SUCCESS", notes = character()) {
  ensure_packages("jsonlite")
  state$finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  state$status <- status
  state$outputs <- outputs
  state$output_sha256 <- setNames(vapply(outputs, sha256_file, character(1)), outputs)
  state$notes <- notes
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
  path <- file.path(cfg$manifest_dir, sprintf("%s_%s.json", state$step_id, stamp))
  assert_new_path(path, FALSE)
  jsonlite::write_json(state, path, pretty = TRUE, auto_unbox = TRUE, na = "null")
  invisible(path)
}

save_session_info <- function(path, overwrite = FALSE) {
  write_lines_safe(capture.output(sessionInfo()), path, overwrite)
}

save_plot_dual <- function(plot, stem, cfg, width = cfg$plot_width, height = cfg$plot_height) {
  pdf <- paste0(stem, ".pdf")
  png <- paste0(stem, ".png")
  assert_new_path(pdf, cfg$overwrite)
  assert_new_path(png, cfg$overwrite)
  ggplot2::ggsave(pdf, plot, width = width, height = height, device = cairo_pdf, limitsize = FALSE)
  ggplot2::ggsave(png, plot, width = width, height = height, dpi = cfg$plot_dpi, limitsize = FALSE)
  c(pdf, png)
}

resolve_meta_column <- function(meta, candidates, required = TRUE, label = "metadata") {
  hit <- candidates[candidates %in% colnames(meta)]
  if (!length(hit)) {
    if (required) stop("Cannot resolve ", label, " column; checked: ", paste(candidates, collapse = ", "), call. = FALSE)
    return(NA_character_)
  }
  hit[[1]]
}

get_rna_matrix <- function(object, layer = c("data", "counts")) {
  layer <- match.arg(layer)
  assays <- safe_assay_names(object)
  if (!"RNA" %in% assays) stop("RNA assay is absent", call. = FALSE)
  rna <- tryCatch(slot(object, "assays")[["RNA"]], error = function(e) NULL)
  out <- NULL
  access_mode <- NA_character_
  if (!is.null(rna) && methods::is(rna, "Assay") && layer %in% methods::slotNames(rna)) {
    out <- methods::slot(rna, layer)
    access_mode <- paste0("legacy_Assay_slot:", layer)
  }
  if (is.null(out)) {
    out <- tryCatch(SeuratObject::LayerData(object, assay = "RNA", layer = layer), error = function(e) NULL)
    if (!is.null(out)) access_mode <- paste0("LayerData:", layer)
  }
  if (is.null(out)) {
    out <- tryCatch(SeuratObject::GetAssayData(object, assay = "RNA", layer = layer), error = function(e) NULL)
    if (!is.null(out)) access_mode <- paste0("GetAssayData_layer:", layer)
  }
  if (is.null(out)) {
    out <- tryCatch(SeuratObject::GetAssayData(object, assay = "RNA", slot = layer), error = function(e) NULL)
    if (!is.null(out)) access_mode <- paste0("GetAssayData_slot:", layer)
  }
  if (is.null(out)) stop("Unable to read RNA ", layer, " using legacy slots, layers, or GetAssayData", call. = FALSE)
  if (!inherits(out, c("matrix", "Matrix"))) stop("RNA ", layer, " is not matrix-like", call. = FALSE)
  attr(out, "gse161529_access_mode") <- access_mode
  out
}

extract_gene_matrix <- function(object, genes, layer = c("data", "counts"), cells = NULL) {
  layer <- match.arg(layer)
  mat <- get_rna_matrix(object, layer)
  present <- intersect(genes, rownames(mat))
  if (!length(present)) stop("None of the requested genes are present in RNA ", layer, call. = FALSE)
  if (is.null(cells)) cells <- colnames(mat)
  missing_cells <- setdiff(cells, colnames(mat))
  if (length(missing_cells)) stop("Requested cells are missing from RNA matrix", call. = FALSE)
  mat[present, cells, drop = FALSE]
}

load_seurat_object <- function(path) {
  if (!file.exists(path)) stop("Input RDS is missing: ", path, call. = FALSE)
  x <- readRDS(path)
  if (!inherits(x, "Seurat")) stop("RDS is not a Seurat object: ", path, call. = FALSE)
  # Fix missing images slot (Seurat v3 objects loaded under Seurat v5 runtime)
  img <- tryCatch(methods::slot(x, "images"), error = function(e) NULL)
  if (is.null(img)) methods::slot(x, "images") <- list()
  x
}

safe_meta_data <- function(object) {
  meta <- tryCatch(methods::slot(object, "meta.data"), error = function(e) NULL)
  if (is.null(meta)) meta <- tryCatch(object[[]], error = function(e) NULL)
  if (!is.data.frame(meta)) stop("Seurat metadata are unavailable or not a data.frame", call. = FALSE)
  meta
}

safe_assay_names <- function(object) {
  x <- tryCatch(names(methods::slot(object, "assays")), error = function(e) character())
  if (!length(x)) x <- tryCatch(SeuratObject::Assays(object), error = function(e) character())
  as.character(x)
}

safe_active_assay <- function(object) {
  x <- tryCatch(as.character(methods::slot(object, "active.assay")), error = function(e) character())
  if (!length(x) || !nzchar(x[[1]])) x <- tryCatch(SeuratObject::DefaultAssay(object), error = function(e) character())
  if (length(x)) x[[1]] else NA_character_
}

safe_reduction_names <- function(object) {
  x <- tryCatch(names(methods::slot(object, "reductions")), error = function(e) character())
  if (!length(x)) x <- tryCatch(SeuratObject::Reductions(object), error = function(e) character())
  as.character(x)
}

safe_cell_names <- function(object) {
  meta <- safe_meta_data(object)
  cells <- rownames(meta)
  if (!length(cells)) {
    counts <- get_rna_matrix(object, "counts")
    cells <- colnames(counts)
  }
  as.character(cells)
}

safe_embedding <- function(object, reduction) {
  reds <- tryCatch(methods::slot(object, "reductions"), error = function(e) NULL)
  emb <- NULL
  if (!is.null(reds) && reduction %in% names(reds)) {
    emb <- tryCatch(methods::slot(reds[[reduction]], "cell.embeddings"), error = function(e) NULL)
  }
  if (is.null(emb)) emb <- tryCatch(SeuratObject::Embeddings(object, reduction = reduction), error = function(e) NULL)
  if (!inherits(emb, c("matrix", "Matrix")) || ncol(emb) < 2L) {
    stop("Reduction ", reduction, " has no accessible >=2D cell embedding", call. = FALSE)
  }
  emb
}

object_version <- function(object) {
  as.character(utils::packageVersion("SeuratObject"))
}

reduction_name <- function(object) {
  red <- safe_reduction_names(object)
  preferred <- c("umap", "tsne", "pca")
  hit <- preferred[preferred %in% red]
  if (length(hit)) hit[[1]] else if (length(red)) red[[1]] else NA_character_
}

embedding_group_plot <- function(object, reduction, group, title = NULL, point_size = 0.25) {
  emb <- safe_embedding(object, reduction)
  meta <- safe_meta_data(object)
  if (is.character(group) && length(group) == 1L && group %in% colnames(meta)) group <- meta[[group]]
  if (length(group) != nrow(meta)) stop("Embedding group length does not match metadata rows", call. = FALSE)
  idx <- match(rownames(emb), rownames(meta))
  if (anyNA(idx)) stop("Embedding cells do not match metadata cell barcodes", call. = FALSE)
  tab <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], group = as.factor(group[idx]), check.names = FALSE)
  ggplot(tab, aes(dim1, dim2, colour = group)) +
    geom_point(size = point_size, alpha = 0.8) + coord_equal() +
    labs(x = colnames(emb)[1], y = colnames(emb)[2], colour = NULL, title = title) +
    theme_classic(base_size = 9) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
    guides(colour = guide_legend(ncol = 4, override.aes = list(size = 2)))
}

embedding_feature_plot <- function(object, reduction, genes, title = NULL, ncol = 3L) {
  ensure_packages("patchwork")
  emb <- safe_embedding(object, reduction)
  dat <- extract_gene_matrix(object, genes, "data", cells = rownames(emb))
  plots <- lapply(rownames(dat), function(gene) {
    tab <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2], expression = as.numeric(dat[gene, ]))
    tab <- tab[order(tab$expression), , drop = FALSE]
    ggplot(tab, aes(dim1, dim2, colour = expression)) +
      geom_point(size = 0.2) + coord_equal() +
      scale_colour_viridis_c(option = "C") +
      labs(x = NULL, y = NULL, colour = "RNA data", title = gene) +
      theme_classic(base_size = 8) + theme(legend.position = "right")
  })
  patchwork::wrap_plots(plots, ncol = ncol) + patchwork::plot_annotation(title = title)
}

subset_legacy_seurat_rna <- function(object, cells, derivation_note = list()) {
  meta <- safe_meta_data(object)
  missing <- setdiff(cells, rownames(meta))
  if (length(missing)) stop("Subset cells are missing from metadata", call. = FALSE)
  counts <- get_rna_matrix(object, "counts")[, cells, drop = FALSE]
  data_mat <- get_rna_matrix(object, "data")[, cells, drop = FALSE]
  project <- tryCatch(methods::slot(object, "project.name"), error = function(e) "GSE161529_derived")
  derived <- SeuratObject::CreateSeuratObject(counts = counts, project = project,
                                               meta.data = meta[cells, , drop = FALSE])
  derived <- SeuratObject::SetAssayData(derived, assay = "RNA", layer = "data", new.data = data_mat)
  for (red in safe_reduction_names(object)) {
    emb <- tryCatch(safe_embedding(object, red), error = function(e) NULL)
    if (is.null(emb) || !all(cells %in% rownames(emb))) next
    old_reds <- tryCatch(methods::slot(object, "reductions"), error = function(e) NULL)
    key <- if (!is.null(old_reds)) tryCatch(methods::slot(old_reds[[red]], "key"), error = function(e) paste0(red, "_")) else paste0(red, "_")
    derived[[red]] <- SeuratObject::CreateDimReducObject(embeddings = emb[cells, , drop = FALSE], key = key, assay = "RNA")
  }
  derived@misc$gse161529_derivation <- derivation_note
  derived
}

classify_compartment <- function(label) {
  z <- tolower(as.character(label))
  out <- rep("unknown/mixed", length(z))
  out[grepl("fibro|caf|stromal", z)] <- "CAF/fibroblast"
  out[grepl("pericy|mural|smooth muscle|smc", z)] <- "pericyte/mural"
  out[grepl("endothel|vascular|lymphatic", z)] <- "endothelial"
  out[grepl("immune|myeloid|lymph|t cell|b cell|macroph|monocyte|mast", z)] <- "immune"
  out[grepl("epithel|tumou?r|malignant|basal|luminal", z)] <- "epithelial/tumor"
  out
}

analysis_registry <- function(cfg) {
  local <- function(name) file.path(cfg$input_dir, name)
  asset <- function(name) file.path(cfg$author_asset_dir, name)
  derived <- function(name) file.path(cfg$updated_object_dir, name)
  rows <- data.frame(
    analysis_id = c("NormalFib", "ER", "HER2", "TNBC"),
    subtype = c("Normal", "ER+", "HER2+", "TNBC"),
    preferred = c(local("SeuratObject_NormTotalFib.rds"), asset("SeuratObject_ERTotalSub.rds"),
                  asset("SeuratObject_HER2Sub.rds"), asset("SeuratObject_TNBCSub.rds")),
    fallback = c(local("SeuratObject_NormTotalFib.rds"), derived("SeuratObject_ERTotalSub_derived_from_author_total.rds"),
                 derived("SeuratObject_HER2Sub_derived_from_author_total.rds"),
                 derived("SeuratObject_TNBCSub_derived_from_author_total.rds")),
    stringsAsFactors = FALSE
  )
  rows$path <- ifelse(file.exists(rows$preferred), rows$preferred, rows$fallback)
  rows$source_status <- ifelse(rows$path == rows$preferred, "author_released_object", "derived_from_author_total_object")
  rows
}

localization_registry <- function(cfg) {
  data.frame(
    analysis_id = c("NormalTotal", "ER_Total", "HER2_Total", "TNBC_Total"),
    subtype = c("Normal", "ER+", "HER2+", "TNBC"),
    path = file.path(cfg$input_dir, c("SeuratObject_NormTotal.rds", "SeuratObject_ERTotal.rds",
                                     "SeuratObject_HER2.rds", "SeuratObject_TNBC.rds")),
    source_status = "author_released_total_object",
    stringsAsFactors = FALSE
  )
}

read_mapping <- function(cfg, require_resolved = TRUE) {
  if (!file.exists(cfg$mapping_file)) stop("Cluster mapping is missing: ", cfg$mapping_file, call. = FALSE)
  x <- utils::read.delim(cfg$mapping_file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("analysis_id", "author_cluster", "compartment")
  if (!all(required %in% colnames(x))) stop("Mapping lacks required columns", call. = FALSE)
  if (require_resolved && any(x$compartment %in% c("", "unknown/mixed", "UNRESOLVED") | is.na(x$compartment))) {
    stop("Cluster mapping contains unresolved compartments", call. = FALSE)
  }
  x
}

attach_mapping <- function(meta, mapping, analysis_id, cluster_col) {
  m <- mapping[mapping$analysis_id == analysis_id, c("author_cluster", "compartment", "author_label"), drop = FALSE]
  if (!nrow(m)) stop("No mapping rows for ", analysis_id, call. = FALSE)
  idx <- match(as.character(meta[[cluster_col]]), as.character(m$author_cluster))
  if (anyNA(idx)) stop("At least one cell cluster is absent from the approved mapping for ", analysis_id, call. = FALSE)
  meta$author_compartment <- m$compartment[idx]
  meta$author_celltype <- m$author_label[idx]
  meta
}

safe_divide <- function(x, y) ifelse(y > 0, x / y, NA_real_)

read_simple_yaml <- function(path) {
  if (!file.exists(path)) stop("Approval file is missing: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[grepl(":", lines, fixed = TRUE)]
  keys <- trimws(sub(":.*$", "", lines))
  vals <- trimws(sub("^[^:]+:", "", lines))
  vals <- gsub("^[\"']|[\"']$", "", vals)
  setNames(as.list(vals), keys)
}

require_boundary_approval <- function(cfg) {
  approval <- read_simple_yaml(cfg$approval_file)
  approved <- tolower(approval$approved %||% "false") %in% c("true", "yes", "1")
  if (!approved) stop("Cluster boundary approval is not granted", call. = FALSE)
  current_hash <- sha256_file(cfg$mapping_file)
  if (!identical(tolower(approval$mapping_sha256 %||% ""), tolower(current_hash))) {
    stop("Approval mapping_sha256 does not match the current mapping", call. = FALSE)
  }
  if (!nzchar(approval$approved_by %||% "")) stop("approved_by is empty", call. = FALSE)
  invisible(approval)
}

candidate_status <- function(eng, thy1) {
  ifelse(eng > 0 & thy1 > 0, "ENGpos_THY1pos_candidate",
         ifelse(eng <= 0 & thy1 > 0, "ENGneg_THY1pos",
                ifelse(eng > 0 & thy1 <= 0, "ENGpos_THY1neg", "ENGneg_THY1neg")))
}

summarise_expression <- function(data_mat, counts_mat, groups, patient) {
  stopifnot(identical(colnames(data_mat), colnames(counts_mat)))
  genes <- intersect(rownames(data_mat), rownames(counts_mat))
  idx <- split(seq_along(groups), groups, drop = TRUE)
  out <- lapply(names(idx), function(g) {
    ii <- idx[[g]]
    do.call(rbind, lapply(genes, function(gene) {
      d <- as.numeric(data_mat[gene, ii])
      cts <- as.numeric(counts_mat[gene, ii])
      data.frame(group = g, gene = gene, pct_expressed = mean(cts > 0) * 100,
                 average_normalized_RNA = mean(d), median_raw_UMI = stats::median(cts),
                 mean_raw_UMI = mean(cts), cell_number = length(ii),
                 patient_support_number = length(unique(patient[ii][cts > 0])),
                 stringsAsFactors = FALSE)
    }))
  })
  do.call(rbind, out)
}

dotplot_from_summary <- function(tab, title = NULL) {
  ggplot(tab, aes(x = gene, y = group, size = pct_expressed, color = average_normalized_RNA)) +
    geom_point() +
    scale_color_viridis_c(option = "C") +
    scale_size_continuous(range = c(0.5, 8)) +
    labs(x = NULL, y = NULL, size = "% detected", color = "Mean RNA", title = title) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

is_valid_file <- function(path, min_bytes = 1, validator = NULL) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < min_bytes) return(FALSE)
  if (is.null(validator)) return(TRUE)
  isTRUE(tryCatch(validator(path), error = function(e) FALSE))
}

validate_r_code <- function(path) {
  if (!is_valid_file(path, min_bytes = 100)) return(FALSE)
  !inherits(try(parse(path, keep.source = FALSE), silent = TRUE), "try-error")
}

extract_author_microenvironment_exclusions <- function(code_path) {
  if (!validate_r_code(code_path)) return(character())
  text <- paste(readLines(code_path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  match <- regexec("(?is)Sub\\s*<-\\s*!\\s*\\(\\s*Cluster\\s*%in%\\s*c\\(([^)]*)\\)\\s*\\)", text, perl = TRUE)
  hit <- regmatches(text, match)[[1]]
  if (length(hit) < 2L) return(character())
  ids <- trimws(unlist(strsplit(gsub("['\"]", "", hit[[2]]), ",")))
  ids[grepl("^[0-9]+$", ids)]
}

copy_or_download <- function(urls, destination, min_bytes = 1, validator = NULL) {
  urls <- unique(urls[nzchar(urls)])
  if (!length(urls)) stop("No download URL was supplied for ", basename(destination), call. = FALSE)
  if (is_valid_file(destination, min_bytes, validator)) {
    return(list(path = destination, status = "already_present_validated", url = urls[[1]], quarantined = NA_character_))
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  quarantined <- NA_character_
  if (file.exists(destination)) {
    quarantined <- paste0(destination, ".invalid_", format(Sys.time(), "%Y%m%dT%H%M%S"))
    if (!file.rename(destination, quarantined)) stop("Could not quarantine invalid existing file: ", destination, call. = FALSE)
  }
  errors <- character()
  for (url in urls) {
    temporary <- paste0(destination, ".part")
    if (file.exists(temporary)) unlink(temporary)
    ok <- tryCatch({
      utils::download.file(url, temporary, mode = "wb", quiet = FALSE,
                           headers = c("User-Agent" = "GSE161529-ASC-P-validation/1.0"))
      is_valid_file(temporary, min_bytes, validator)
    }, error = function(e) { errors <<- c(errors, paste0(url, ": ", conditionMessage(e))); FALSE })
    if (isTRUE(ok)) {
      if (!file.rename(temporary, destination)) stop("Validated download could not be moved into place: ", destination, call. = FALSE)
      return(list(path = destination, status = "downloaded_validated", url = url, quarantined = quarantined))
    }
    if (file.exists(temporary)) unlink(temporary)
  }
  stop("All downloads failed validation for ", basename(destination), ": ", paste(errors, collapse = " | "),
       if (!is.na(quarantined)) paste0("; invalid prior file preserved at ", quarantined) else "", call. = FALSE)
}

write_waiting_report <- function(cfg, reasons) {
  path <- file.path(cfg$report_dir, "WAITING_FOR_AUTHOR_CLUSTER_MAPPING_APPROVAL.md")
  lines <- c("# Waiting for author cluster mapping approval", "",
             "Candidate detection is blocked. No biological identity was guessed.", "",
             "## Reasons", "", paste0("- ", reasons))
  if (!file.exists(path)) write_lines_safe(lines, path, FALSE)
  invisible(path)
}

archive_precheck_outputs <- function(cfg) {
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
  archive_root <- file.path(cfg$project_root, "provenance", "reruns", paste0("precheck_", stamp))
  if (dir.exists(archive_root)) stop("Precheck archive already exists: ", archive_root, call. = FALSE)
  fixed <- c(
    file.path(cfg$manifest_dir, "00b_author_asset_manifest.tsv"),
    file.path(cfg$table_dir, "00_GSE161529_object_inventory.tsv"),
    file.path(cfg$table_dir, "00_GSE161529_metadata_columns_by_object.tsv"),
    file.path(cfg$table_dir, "01_sub_object_resolution.tsv"), cfg$mapping_file,
    file.path(cfg$report_dir, "00_GSE161529_INPUT_AUDIT.md"),
    file.path(cfg$report_dir, "01_AUTHOR_ANNOTATION_MAPPING.md"),
    file.path(cfg$report_dir, "REVIEW_CLUSTER_BOUNDARY_CHECKPOINT.md"),
    file.path(cfg$report_dir, "WAITING_FOR_AUTHOR_CLUSTER_MAPPING_APPROVAL.md"),
    file.path(cfg$report_dir, "MANUAL_DOWNLOAD_REQUIRED.md"),
    file.path(cfg$log_dir, "00_audit_sessionInfo.txt"), cfg$approval_file
  )
  patterns <- c(
    list.files(file.path(cfg$figure_dir, "audit"), pattern = "^01_", full.names = TRUE),
    list.files(cfg$updated_object_dir, pattern = "_derived_from_author_total\\.rds$", full.names = TRUE)
  )
  paths <- unique(c(fixed, patterns))
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return(invisible(NA_character_))
  project_norm <- normalizePath(cfg$project_root, winslash = "/")
  relative <- substring(normalizePath(paths, winslash = "/"), nchar(project_norm) + 2L)
  archived <- file.path(archive_root, relative)
  for (i in seq_along(paths)) {
    dir.create(dirname(archived[[i]]), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(paths[[i]], archived[[i]])) stop("Failed to archive precheck output: ", paths[[i]], call. = FALSE)
  }
  manifest <- data.frame(original_path = paths, archived_path = archived,
                         archived_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), stringsAsFactors = FALSE)
  write_tsv_safe(manifest, file.path(archive_root, "ARCHIVE_MANIFEST.tsv"), FALSE)
  invisible(archive_root)
}

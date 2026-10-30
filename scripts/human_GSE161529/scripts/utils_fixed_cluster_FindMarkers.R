# Fixed-cluster, targeted FindMarkers runner for author-released GSE161529 objects.
# This helper never modifies an author RDS, and its cell-level tests are descriptive only.

fixed_marker_object_registry <- function(cfg) {
  rbind(
    data.frame(object_id = cfg$local_objects$object_id[cfg$local_objects$object_id != "NormTotal"],
               input_dir = cfg$input_dir,
               filename = cfg$local_objects$filename[cfg$local_objects$object_id != "NormTotal"],
               stringsAsFactors = FALSE),
    data.frame(object_id = cfg$required_author_assets$object_id,
               input_dir = cfg$author_asset_dir,
               filename = cfg$required_author_assets$filename,
               stringsAsFactors = FALSE)
  )
}

run_fixed_cluster_findmarkers <- function(object_id, cfg = source_project_config()) {
  ensure_packages(c("Seurat", "SeuratObject", "Matrix", "digest", "jsonlite", "pheatmap", "future"))
  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)
  set.seed(cfg$seed)
  if (identical(Sys.getenv("GSE161529_OVERWRITE"), "1")) cfg$overwrite <- TRUE

  registry <- fixed_marker_object_registry(cfg)
  row <- registry[registry$object_id == object_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Unsupported object_id: ", object_id, call. = FALSE)
  input_path <- file.path(row$input_dir[[1]], row$filename[[1]])
  if (!file.exists(input_path)) stop("Missing RDS for ", object_id, ": ", input_path, call. = FALSE)
  state <- begin_step(cfg, paste0("01b_", object_id, "_FindMarkers"), input_path)

  obj <- load_seurat_object(input_path)
  meta <- safe_meta_data(obj)
  cluster_col <- resolve_meta_column(meta, cfg$cluster_candidates, TRUE, paste0(object_id, " author cluster"))
  patient_col <- resolve_meta_column(meta, cfg$patient_candidates, TRUE, paste0(object_id, " patient"))
  if (!"RNA" %in% safe_assay_names(obj)) stop(object_id, " lacks an RNA assay", call. = FALSE)
  SeuratObject::DefaultAssay(obj) <- "RNA"
  cluster_values <- as.character(meta[[cluster_col]])
  clusters <- sort(unique(cluster_values))
  SeuratObject::Idents(obj) <- factor(cluster_values, levels = clusters)
  annotation_features <- intersect(cfg$all_audit_genes, rownames(obj[["RNA"]]))
  if (length(annotation_features) < 10L) stop("Too few audit marker genes present in ", object_id, call. = FALSE)

  marker_rows <- list()
  for (cl in clusters) {
    message("FindMarkers ", object_id, " cluster ", cl, " vs all remaining clusters")
    res <- Seurat::FindMarkers(
      object = obj, ident.1 = cl, assay = "RNA", slot = "data", features = annotation_features,
      only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25,
      test.use = "wilcox", max.cells.per.ident = 3000L, verbose = FALSE
    )
    if (!nrow(res)) next
    cells <- rownames(meta)[cluster_values == cl]
    res$gene <- rownames(res)
    res$analysis_id <- object_id
    res$author_cluster <- cl
    numeric_cluster <- suppressWarnings(as.integer(cl))
    res$author_cluster_one_based_index <- if (is.na(numeric_cluster)) NA_integer_ else numeric_cluster + 1L
    res$cell_number <- length(cells)
    res$patient_number <- length(unique(as.character(meta[cells, patient_col])))
    res <- res[, c("analysis_id", "author_cluster", "author_cluster_one_based_index", "gene",
                   setdiff(colnames(res), c("analysis_id", "author_cluster", "author_cluster_one_based_index", "gene"))), drop = FALSE]
    marker_rows[[length(marker_rows) + 1L]] <- res
  }
  markers <- if (length(marker_rows)) do.call(rbind, marker_rows) else data.frame()
  if (!nrow(markers)) stop("FindMarkers returned no markers for ", object_id, call. = FALSE)
  markers <- markers[order(markers$author_cluster, markers$p_val_adj, -markers$avg_log2FC), , drop = FALSE]
  marker_out <- file.path(cfg$table_dir, paste0(object_id, "_cluster_FindMarkers.tsv"))
  write_tsv_safe(markers, marker_out, cfg$overwrite)

  marker_sets <- list(
    CAF_fibroblast = cfg$fibroblast_markers,
    mural_pericyte = cfg$mural_markers,
    endothelial = cfg$endothelial_markers,
    immune = c("PTPRC", "CD3D", "CD3E", "TRAC", "CD8A", "CD4", "MS4A1", "CD79A", "CD68", "CD14", "TYROBP", "LST1", "FCER1G", "C1QC", "NKG7", "GNLY", "MZB1", "JCHAIN"),
    epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT14", "KRT5", "KRT17", "KRT7", "MUC1", "KRT1")
  )
  score_data <- extract_gene_matrix(obj, unique(unlist(marker_sets, use.names = FALSE)), "data")
  score_rows <- list()
  for (cl in clusters) {
    cells <- rownames(meta)[cluster_values == cl]
    numeric_cluster <- suppressWarnings(as.integer(cl))
    for (nm in names(marker_sets)) {
      genes <- intersect(marker_sets[[nm]], rownames(score_data))
      if (!length(genes)) next
      block <- score_data[genes, cells, drop = FALSE]
      score_rows[[length(score_rows) + 1L]] <- data.frame(
        analysis_id = object_id, author_cluster = cl,
        author_cluster_one_based_index = if (is.na(numeric_cluster)) NA_integer_ else numeric_cluster + 1L,
        lineage = nm, marker_gene_number = length(genes),
        percent_expressed = mean(as.matrix(block) > 0) * 100,
        mean_normalized_RNA = mean(as.matrix(block)), cell_number = length(cells),
        patient_number = length(unique(as.character(meta[cells, patient_col]))), stringsAsFactors = FALSE
      )
    }
  }
  scores <- do.call(rbind, score_rows)
  score_out <- file.path(cfg$table_dir, paste0(object_id, "_cluster_marker_scores.tsv"))
  write_tsv_safe(scores, score_out, cfg$overwrite)

  top <- do.call(rbind, lapply(split(markers, markers$author_cluster), function(x) utils::head(x, 15L)))
  top_out <- file.path(cfg$table_dir, paste0(object_id, "_cluster_top15_markers.tsv"))
  write_tsv_safe(top, top_out, cfg$overwrite)

  top_genes <- unique(top$gene)
  heat_data <- extract_gene_matrix(obj, top_genes, "data")
  cluster_means <- do.call(cbind, lapply(clusters, function(cl) {
    cells <- rownames(meta)[cluster_values == cl]
    Matrix::rowMeans(heat_data[, cells, drop = FALSE])
  }))
  colnames(cluster_means) <- paste0("Cluster_", clusters)
  heat_outputs <- character()
  if (nrow(cluster_means) > 1L) {
    heat_z <- t(scale(t(as.matrix(cluster_means))))
    heat_z[!is.finite(heat_z)] <- 0
    stem <- file.path(cfg$figure_dir, "identity_boundary", paste0(object_id, "_cluster_top15_marker_heatmap"))
    pdf_path <- paste0(stem, ".pdf"); png_path <- paste0(stem, ".png")
    assert_new_path(pdf_path, cfg$overwrite); assert_new_path(png_path, cfg$overwrite)
    grDevices::pdf(pdf_path, width = 10, height = 12)
    pheatmap::pheatmap(heat_z, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = TRUE,
                       fontsize_row = 5, main = paste0(object_id, " top marker genes by author cluster"))
    grDevices::dev.off()
    grDevices::png(png_path, width = 2400, height = 2800, res = 300)
    pheatmap::pheatmap(heat_z, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = TRUE,
                       fontsize_row = 5, main = paste0(object_id, " top marker genes by author cluster"))
    grDevices::dev.off()
    heat_outputs <- c(pdf_path, png_path)
  }

  report_out <- file.path(cfg$report_dir, paste0("01b_", object_id, "_FINDMARKERS_REPORT.md"))
  report <- c(
    paste0("# ", object_id, " fixed-cluster marker analysis"), "",
    "Read-only descriptive analysis of the author-released Seurat object.",
    "No QC, normalization, integration, dimensional reduction, clustering, or RDS write-back was performed.", "",
    paste0("Input SHA256: `", sha256_file(input_path), "`"),
    paste0("Cluster column: `", cluster_col, "`; clusters: ", paste(clusters, collapse = ", ")),
    paste0("Patient column: `", patient_col, "`; patients: ", length(unique(as.character(meta[[patient_col]])))),
    paste0("FindMarkers parameters: assay=RNA; slot=data; targeted audit marker panel (n=", length(annotation_features), "); test=wilcox; only.pos=TRUE; min.pct=0.10; logfc.threshold=0.25; max.cells.per.ident=3000; seed=", cfg$seed, "."), "",
    paste0("Marker rows: ", nrow(markers), "; top-marker rows: ", nrow(top), "; score rows: ", nrow(scores)), "",
    "Cell-level Wilcoxon output is descriptive evidence for author-cluster review, not patient-level inferential differential expression.",
    "Review lineages with patient support, mixed-lineage flags, and author Total/Sub context before approving a mapping."
  )
  write_lines_safe(report, report_out, cfg$overwrite)
  end_step(cfg, state, c(marker_out, score_out, top_out, heat_outputs, report_out), status = "SUCCESS",
           notes = c("Fixed-cluster targeted FindMarkers completed in memory.", "Author RDS was not modified."))
  rm(obj); invisible(gc())
}

# Panel C colored cluster labels + Panel E legend fix.
# Visual-only refinement: no scientific data changes.
# Creates new versioned outputs; never overwrites existing files.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so --file is available", call. = FALSE)
project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "scripts", "utils_GSE161529_ASC_P.R"))
cfg <- source_project_config()
if (identical(Sys.getenv("GSE161529_OVERWRITE"), "1")) cfg$overwrite <- TRUE
ensure_packages(c("SeuratObject", "Matrix", "digest", "jsonlite", "ggplot2", "patchwork", "scales", "ggtext", "png"))
set.seed(cfg$seed)

cohort_arg <- grep("^--cohort=", commandArgs(trailingOnly = TRUE), value = TRUE)
if (!length(cohort_arg)) stop("Specify --cohort=NormalTotal|ERTotal|HER2|TNBC|all", call. = FALSE)
cohort_val <- sub("^--cohort=", "", cohort_arg[[1]])
cohorts <- if (tolower(cohort_val) == "all") c("NormalTotal", "ERTotal", "HER2", "TNBC") else cohort_val
if (!all(cohorts %in% c("NormalTotal", "ERTotal", "HER2", "TNBC"))) stop("Unknown cohort(s): ", paste(setdiff(cohorts, c("NormalTotal", "ERTotal", "HER2", "TNBC")), collapse = ", "), call. = FALSE)

# ---- Cohort specifications ----
cohort_specs <- list(
  NormalTotal = list(
    rds = "SeuratObject_NormTotal.rds",
    map = "NormTotal_final_cluster_annotation.tsv",
    map_mode = "normtotal_annotation",
    denominator_mode = "fixed_cluster",
    denominator_cluster = "1",
    denominator_label = "Cluster 1 fibroblasts",
    tissue_label = "normal human breast tissue",
    panel_prefix = "NormalTotal",
    prior_panels = c(
      "Figure_A_NormalTotal_fixed_annotation_FINAL.png",
      "Figure_C_NormalTotal_ASC_P_like_candidate_map_FINAL.png",
      "Figure_D_NormalTotal_ASC_P_four_group_DotPlot_FINAL.png"
    )
  ),
  ERTotal = list(
    rds = "SeuratObject_ERTotal.rds",
    map = "ERTotal_cluster_annotation_map_approved.tsv",
    map_mode = "approved_map",
    denominator_mode = "approved_map",
    denominator_label = "approved CAF denominator",
    tissue_label = "ER-positive human breast tumors",
    panel_prefix = "ERTotal",
    prior_panels = c(
      "Figure_A_ERTotal_fixed_annotation_FINAL.png",
      "Figure_B_ERTotal_ASC_P_like_candidate_map_FINAL.png",
      "Figure_D_ERTotal_ASC_P_four_group_DotPlot_FINAL.png"
    )
  ),
  HER2 = list(
    rds = "SeuratObject_HER2.rds",
    map = "HER2_cluster_annotation_map_approved.tsv",
    map_mode = "approved_map",
    denominator_mode = "approved_map",
    denominator_label = "approved Cluster 4 CAF compartment",
    tissue_label = "HER2-positive human breast tumors",
    panel_prefix = "HER2",
    prior_panels = c(
      "Figure_A_HER2_fixed_annotation_FINAL.png",
      "Figure_B_HER2_ASC_P_like_candidate_map_FINAL.png",
      "Figure_D_HER2_ASC_P_three_group_DotPlot_FINAL.png"
    )
  ),
  TNBC = list(
    rds = "SeuratObject_TNBC.rds",
    map = "TNBC_cluster_annotation_map_approved.tsv",
    map_mode = "approved_map",
    denominator_mode = "approved_map",
    denominator_label = "approved Cluster 5 CAF compartment",
    tissue_label = "TNBC human breast tumors",
    panel_prefix = "TNBC",
    prior_panels = c(
      "Figure_A_TNBC_fixed_annotation_FINAL.png",
      "Figure_B_TNBC_ASC_P_like_candidate_map_FINAL.png",
      "Figure_D_TNBC_ASC_P_three_group_DotPlot_FINAL.png"
    )
  )
)

# ---- Expected audit values ----
expected_values <- list(
  NormalTotal = list(candidates = 545, core = 196, pct = 35.96),
  ERTotal = list(candidates = 447, core = 201, pct = 44.97),
  HER2 = list(candidates = 258, core = 99, pct = 38.37),
  TNBC = list(candidates = 674, core = 317, pct = 47.03)
)

asc_genes <- c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM")

lineage_sets <- list(
  CAF = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "DPT", "SFRP2", "PDGFRA", "COL10A1", "COMP"),
  Mural = c("PDGFRB", "RGS5", "CSPG4", "ACTA2", "TAGLN", "MYH11", "NOTCH3", "MCAM", "CAVIN3"),
  Endothelial = c("PECAM1", "VWF", "CLDN5", "EMCN", "KDR", "PLVAP", "ESAM", "RAMP2"),
  Immune = c("PTPRC", "CD3D", "CD3E", "TRAC", "CD79A", "MS4A1", "NKG7", "GNLY", "LST1", "TYROBP"),
  Epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT14", "KRT5", "KRT17", "MUC1")
)

# ---- Helper functions ----
luminance <- function(hex) {
  rgb_val <- col2rgb(hex)
  0.299 * rgb_val[1] + 0.587 * rgb_val[2] + 0.114 * rgb_val[3]
}

text_color_for <- function(bg) {
  ifelse(luminance(bg) > 128, "black", "white")
}

save_plot <- function(plot, stem, width, height) {
  paths <- paste0(stem, c(".pdf", ".png", ".tiff"))
  for (path in paths) assert_new_path(path, cfg$overwrite)
  # Use the base PDF device for this raster-backed panel export; Cairo can fail
  # on long Windows paths while writing an otherwise valid figure.
  pdf_device <- "pdf"
  ggplot2::ggsave(paths[[1]], plot, width = width, height = height, units = "in", device = pdf_device, bg = "white", limitsize = FALSE)
  # Use cairo_png for ggtext::element_markdown support
  grDevices::png(paths[[2]], width = width, height = height, units = "in", res = 300, type = if (isTRUE(capabilities("cairo"))) "cairo" else "windows")
  print(plot); grDevices::dev.off()
  # TIFF with cairo
  grDevices::tiff(paths[[3]], width = width, height = height, units = "in", res = 600, compression = "lzw", type = if (isTRUE(capabilities("cairo"))) "cairo" else "windows")
  print(plot); grDevices::dev.off()
  paths
}

save_tsv <- function(x, path) {
  assert_new_path(path, cfg$overwrite)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  path
}

read_panel_png <- function(path) {
  if (!is_valid_file(path, min_bytes = 1000)) stop("Missing or empty reviewed panel PNG: ", path, call. = FALSE)
  patchwork::wrap_elements(full = grid::rasterGrob(png::readPNG(path), interpolate = TRUE), clip = FALSE)
}

parse_denominator <- function(map, expected_clusters) {
  required <- c("cluster", "celltype_final", "compartment", "candidate_denominator")
  if (!all(required %in% colnames(map))) stop("Approved map must contain: ", paste(required, collapse = ", "), call. = FALSE)
  map$cluster <- trimws(as.character(map$cluster))
  map$compartment <- tolower(trimws(as.character(map$compartment)))
  if (anyDuplicated(map$cluster) || !setequal(map$cluster, expected_clusters) || nrow(map) != length(expected_clusters)) {
    stop("Approved map must cover every author cluster exactly once", call. = FALSE)
  }
  denom_text <- tolower(trimws(as.character(map$candidate_denominator)))
  true_values <- c("true", "1", "yes", "y")
  false_values <- c("false", "0", "no", "n")
  if (anyNA(denom_text) || any(!denom_text %in% c(true_values, false_values))) {
    stop("candidate_denominator must be TRUE/FALSE", call. = FALSE)
  }
  denom <- denom_text %in% true_values
  if (!any(denom)) stop("No candidate_denominator=TRUE cluster found", call. = FALSE)
  list(denominator_clusters = map$cluster[denom], map = map)
}

out_dir <- file.path(cfg$figure_dir, "reviewer_ASC_P")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Process each cohort ----
for (cohort in cohorts) {
  spec <- cohort_specs[[cohort]]
  message("\n=== [03d ", cohort, "] Starting Panel C color labels + Panel E legend fix ===")

  input_path <- file.path(cfg$input_dir, spec$rds)
  map_path <- file.path(cfg$table_dir, spec$map)
  prior_panel_paths <- file.path(out_dir, spec$prior_panels)

  if (!file.exists(input_path)) stop("Missing RDS for ", cohort, ": ", input_path, call. = FALSE)
  if (!file.exists(map_path)) stop("Missing map for ", cohort, ": ", map_path, call. = FALSE)
  if (!all(file.exists(prior_panel_paths))) stop("Missing prior A/B/D panel PNGs for ", cohort, call. = FALSE)

  state <- begin_step(cfg, paste0("03d_", cohort, "_panelC_colorlabels_panelE_legend_fix"), c(input_path, map_path, prior_panel_paths))

  # ---- Load object ----
  message("[03d ", cohort, "] Stage 1/5: loading object and annotations")
  obj <- load_seurat_object(input_path)
  meta <- safe_meta_data(obj)
  cluster_col <- resolve_meta_column(meta, cfg$cluster_candidates, TRUE, paste0(cohort, " author cluster"))
  if (!"RNA" %in% safe_assay_names(obj)) stop(cohort, " lacks RNA assay", call. = FALSE)
  counts <- get_rna_matrix(obj, "counts")
  data_mat <- get_rna_matrix(obj, "data")
  cells <- Reduce(intersect, list(rownames(meta), colnames(counts), colnames(data_mat)))
  meta <- meta[cells, , drop = FALSE]; counts <- counts[, cells, drop = FALSE]; data_mat <- data_mat[, cells, drop = FALSE]
  cluster_values <- trimws(as.character(meta[[cluster_col]]))
  sort_cluster_ids <- function(x) {
    numeric_x <- suppressWarnings(as.numeric(x))
    if (length(x) && all(is.finite(numeric_x))) x[order(numeric_x)] else sort(x)
  }
  clusters <- sort_cluster_ids(unique(cluster_values))

  # ---- Parse map / denominator ----
  if (identical(spec$map_mode, "normtotal_annotation")) {
    fixed_annotation <- utils::read.delim(map_path, stringsAsFactors = FALSE, check.names = FALSE)
    celltype_final <- fixed_annotation$celltype_final[match(clusters, as.character(fixed_annotation$cluster))]
    compartment <- rep("unknown", length(clusters))
    compartment[celltype_final %in% c("Fibroblast")] <- "caf"
    compartment[grepl("endothel", celltype_final, ignore.case = TRUE)] <- "endothelial"
    compartment[grepl("mural|pericyte", celltype_final, ignore.case = TRUE)] <- "mural"
    compartment[grepl("immune", celltype_final, ignore.case = TRUE)] <- "immune"
    compartment[grepl("luminal|basal|epithel", celltype_final, ignore.case = TRUE)] <- "epithelial"
    denom_clusters <- spec$denominator_cluster
  } else {
    parsed <- parse_denominator(utils::read.delim(map_path, stringsAsFactors = FALSE, check.names = FALSE), clusters)
    map <- parsed$map
    map <- map[match(clusters, map$cluster), , drop = FALSE]
    celltype_final <- map$celltype_final
    compartment <- map$compartment
    denom_clusters <- parsed$denominator_clusters
  }

  # ---- Candidate definition ----
  message("[03d ", cohort, "] Stage 2/5: computing candidates and ENG+THY1 core")
  missing_asc <- setdiff(asc_genes, rownames(counts))
  if (length(missing_asc)) stop("Missing ASC_P genes: ", paste(missing_asc, collapse = ", "), call. = FALSE)
  detected <- as.matrix(counts[asc_genes, cells, drop = FALSE]) > 0
  n_detected <- as.integer(colSums(detected))
  cell_is_denom <- cluster_values %in% denom_clusters
  candidate <- cell_is_denom & n_detected >= 3L
  eng_thy1_core <- candidate & detected["ENG", ] & detected["THY1", ]
  n_denom <- sum(cell_is_denom)
  n_candidate <- sum(candidate)
  n_core <- sum(eng_thy1_core)
  pct_candidate <- 100 * n_candidate / n_denom
  pct_core <- 100 * n_core / n_candidate

  exp <- expected_values[[cohort]]
  if (n_candidate != exp$candidates || n_core != exp$core || abs(pct_core - exp$pct) > 0.1) {
    stop(sprintf("Audit mismatch for %s: expected %d/%d/%.2f%%, got %d/%d/%.2f%%",
                 cohort, exp$candidates, exp$core, exp$pct, n_candidate, n_core, pct_core), call. = FALSE)
  }
  message(sprintf("[03d %s] Audit OK: candidates=%d, core=%d, core%%=%.2f%%", cohort, n_candidate, n_core, pct_core))

  # ---- Get embedding ----
  reds <- safe_reduction_names(obj)
  reduction <- if ("tsne" %in% reds) "tsne" else if ("umap" %in% reds) "umap" else stop("No tsne/umap")
  reduction_label <- if (reduction == "tsne") "t-SNE" else "UMAP"
  emb <- safe_embedding(obj, reduction)
  emb <- emb[rownames(emb) %in% cells, , drop = FALSE]
  candidate <- candidate[match(rownames(emb), cells)]
  eng_thy1_core <- eng_thy1_core[match(rownames(emb), cells)]

  # ---- Build Panel A palette (same as 03b) ----
  labels <- setNames(paste0("Cluster ", clusters, ": ", celltype_final), clusters)
  if (identical(spec$map_mode, "normtotal_annotation")) {
    # NormalTotal uses a specific palette
    norm_palette <- c(
      "Luminal progenitor" = "#009E73", "Fibroblast" = "#D55E00",
      "Basal/myoepithelial" = "#CC79A7", "Mature luminal" = "#56B4E9",
      "Vascular endothelial" = "#0072B2", "Mural/pericyte-VSMC" = "#E69F00",
      "Immune" = "#999999", "Unresolved PTPRC+ENG cluster" = "#7B61A8"
    )
    cluster_palette <- setNames(norm_palette[celltype_final], labels)
  } else {
    cluster_palette <- setNames(grDevices::hcl.colors(length(clusters), "Dark 3"), labels)
  }
  # Map cluster ID to its Panel A color
  cluster_to_color <- setNames(unname(cluster_palette), clusters)

  # ---- Panel A: fixed annotation with the legend on the right ----
  # Rebuild from the same embedding, author clusters and approved annotation
  # used by the reviewed Panel A.  Only legend placement is changed.
  emb_cluster <- cluster_values[match(rownames(emb), cells)]
  annotation_df <- data.frame(
    dim1 = emb[, 1], dim2 = emb[, 2],
    cluster = emb_cluster,
    annotation = factor(unname(labels[emb_cluster]), levels = names(cluster_palette)),
    stringsAsFactors = FALSE
  )
  annotation_centres <- stats::aggregate(
    annotation_df[, c("dim1", "dim2")],
    by = list(cluster = annotation_df$cluster), FUN = stats::median
  )
  pA <- ggplot2::ggplot(annotation_df, ggplot2::aes(dim1, dim2, colour = annotation)) +
    ggplot2::geom_point(size = 0.15, alpha = 0.7) +
    ggplot2::geom_label(
      data = annotation_centres,
      ggplot2::aes(dim1, dim2, label = cluster),
      inherit.aes = FALSE, colour = "black", fill = "white", linewidth = 0.22,
      size = 5, fontface = "bold", alpha = 0.8
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_manual(values = cluster_palette, drop = FALSE) +
    ggplot2::labs(
      title = paste0("Broad cell-type annotation of ", spec$tissue_label),
      x = paste0(reduction_label, "-1"), y = paste0(reduction_label, "-2"), colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 6.5),
      legend.key.height = grid::unit(0.38, "cm"),
      legend.key.width = grid::unit(0.32, "cm")
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1, override.aes = list(size = 2)))

  # ---- Panel C: lineage DotPlot with Panel-A-coloured lineage strips ----
  message("[03d ", cohort, "] Stage 3/5: creating Panel C with coloured lineage strips")
  lineage_order <- unique(unlist(lineage_sets, use.names = FALSE))
  strip_lineages <- c("CAF", "Endothelial", "Epithelial_tumor", "Immune", "Mural")
  lineage_rows <- list()
  for (cl in clusters) {
    cl_cells <- cells[cluster_values == cl]
    if (!length(cl_cells)) next
    for (gene in lineage_order) {
      if (!gene %in% rownames(counts) || !gene %in% rownames(data_mat)) next
      raw <- as.numeric(counts[gene, cl_cells])
      norm <- as.numeric(data_mat[gene, cl_cells])
      lineage_rows[[length(lineage_rows) + 1L]] <- data.frame(
        group = cl, gene = gene,
        percent_expressed_counts = 100 * mean(raw > 0),
        average_normalized_RNA = mean(norm),
        stringsAsFactors = FALSE
      )
    }
  }
  lineage_summary <- do.call(rbind, lineage_rows)
  gene_to_lineage <- setNames(rep(names(lineage_sets), lengths(lineage_sets)), unlist(lineage_sets, use.names = FALSE))
  lineage_summary$lineage <- factor(
    unname(gene_to_lineage[as.character(lineage_summary$gene)]),
    levels = strip_lineages
  )
  lineage_summary$gene <- factor(lineage_summary$gene, levels = lineage_order)
  lineage_summary$group <- factor(lineage_summary$group, levels = rev(clusters))

  # A lineage can contain more than one biological cluster.  For each strip,
  # use the Panel A colour of the cluster with the highest mean expression
  # across that lineage's marker set. This also handles cohorts without a
  # separately annotated mural cluster without inventing an identity.
  strip_cluster <- vapply(strip_lineages, function(lineage_name) {
    x <- lineage_summary[as.character(lineage_summary$lineage) == lineage_name, , drop = FALSE]
    by_cluster <- tapply(x$average_normalized_RNA, as.character(x$group), mean, na.rm = TRUE)
    names(by_cluster)[which.max(by_cluster)][[1]]
  }, character(1))
  strip_fill <- unname(cluster_to_color[strip_cluster])
  if (anyNA(strip_fill)) stop("Unable to resolve Panel A colour for one or more Panel C lineage strips", call. = FALSE)

  set_strip_grob_style <- function(grob, fill, text_colour) {
    if (inherits(grob, "rect")) {
      grob$gp$fill <- fill
      grob$gp$col <- NA
    }
    if (inherits(grob, "text")) {
      grob$gp$col <- text_colour
    }
    if (!is.null(grob$children)) grob$children <- lapply(grob$children, set_strip_grob_style, fill = fill, text_colour = text_colour)
    if (!is.null(grob$grobs)) grob$grobs <- lapply(grob$grobs, set_strip_grob_style, fill = fill, text_colour = text_colour)
    grob
  }

  cluster_label_lookup <- stats::setNames(as.character(celltype_final), as.character(clusters))
  if (anyNA(cluster_label_lookup) || any(!nzchar(cluster_label_lookup))) {
    stop("Missing cluster annotation label for Panel C y-axis", call. = FALSE)
  }

  pC <- ggplot2::ggplot(lineage_summary, ggplot2::aes(gene, group, size = percent_expressed_counts, colour = average_normalized_RNA)) +
    ggplot2::geom_point() +
    ggplot2::scale_size(range = c(0.4, 7), name = "Cells detected (%)") +
    ggplot2::scale_colour_viridis_c(option = "C", name = "Mean normalized expression") +
    ggplot2::facet_grid(. ~ lineage, scales = "free_x", space = "free_x") +
    ggplot2::labs(title = "Lineage-marker expression across cell clusters", x = NULL, y = NULL) +
    ggplot2::scale_y_discrete(labels = function(x) {
      x_chr <- as.character(x)
      ann <- unname(cluster_label_lookup[x_chr])
      paste0(x_chr, ": ", ann)
    }) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.ticks.y = ggplot2::element_line(),
      strip.background = ggplot2::element_rect(fill = "grey94", colour = NA),
      strip.text.x = ggplot2::element_text(face = "bold", size = 9),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
  pC_grob <- ggplot2::ggplotGrob(pC)
  strip_indices <- grep("^strip-t", pC_grob$layout$name)
  if (length(strip_indices) != length(strip_lineages)) stop("Unexpected Panel C strip layout", call. = FALSE)
  for (i in seq_along(strip_indices)) {
    pC_grob$grobs[[strip_indices[[i]]]] <- set_strip_grob_style(
      pC_grob$grobs[[strip_indices[[i]]]],
      fill = strip_fill[[i]],
      text_colour = text_color_for(strip_fill[[i]])
    )
  }
  pC_render <- patchwork::wrap_elements(full = pC_grob, clip = FALSE)

  # ---- Panel E: ENG+THY1 core with fixed legend ----
  message("[03d ", cohort, "] Stage 4/5: creating Panel E with fixed legend")
  core_group <- rep("Other cells", nrow(emb))
  core_group[eng_thy1_core] <- "ENG+THY1 core\nwithin >=3/6"
  core_levels <- c("Other cells", "ENG+THY1 core\nwithin >=3/6")
  core_group <- factor(core_group, levels = core_levels)

  panel_e_data <- data.frame(
    cell = rownames(emb), dim1 = emb[, 1], dim2 = emb[, 2],
    ASC_P_n_detected = n_detected[match(rownames(emb), cells)],
    ASC_P_like_candidate_ge3 = candidate,
    ENG_THY1_core_within_ge3 = eng_thy1_core,
    panelE_group = as.character(core_group),
    stringsAsFactors = FALSE
  )

  pE <- ggplot2::ggplot(panel_e_data, ggplot2::aes(dim1, dim2, colour = panelE_group)) +
    ggplot2::geom_point(size = 0.45, alpha = 0.5) +
    ggplot2::scale_colour_manual(values = c(
      "Other cells" = "grey85",
      "ENG+THY1 core\nwithin >=3/6" = "#E31A1C"
    ), drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "ENG+THY1 within >=3/6 ASC_P-like candidates",
      subtitle = sprintf("n = %d of %d candidates (%.2f%%)", n_core, n_candidate, pct_core),
      x = paste0(reduction_label, "-1"), y = paste0(reduction_label, "-2"), colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      aspect.ratio = 1, legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(5, 5, 15, 5)
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.2, alpha = 1)))

  # ---- Save Panel C and Panel E ----
  message("[03d ", cohort, "] Stage 5/5: exporting outputs")
  panel_a_stem <- file.path(out_dir, paste0("Figure_A_", spec$panel_prefix, "_fixed_annotation_LEGEND_RIGHT_FINAL"))
  panel_c_stem <- file.path(out_dir, paste0("Figure_C_", spec$panel_prefix, "_gene_level_lineage_DotPlot_COLORED_LINEAGE_STRIPS_REVISED_FINAL"))
  panel_e_suffix <- "_CORE_ONLY_RED_FINAL"
  panel_e_stem <- file.path(out_dir, paste0("Figure_E_", spec$panel_prefix, "_ASC_P_ENG_THY1_within_3of6_LEGEND_FIXED_REVISED_FINAL", panel_e_suffix))
  panel_a_outputs <- save_plot(pA, panel_a_stem, 12, 9)
  panel_c_outputs <- save_plot(pC_render, panel_c_stem, 16, 7)
  panel_e_outputs <- save_plot(pE, panel_e_stem, 10, 9)

  # Save data TSVs
  data_path <- file.path(out_dir, paste0("Figure_E_", spec$panel_prefix, "_ENG_THY1_within_3of6_embedding_data_LEGEND_FIXED_REVISED.tsv"))
  summary_path <- file.path(out_dir, paste0(spec$panel_prefix, "_ASC_P_ENG_THY1_within_3of6_summary_LEGEND_FIXED_REVISED.tsv"))
  summary_table <- data.frame(
    cohort = cohort, reduction = reduction_label, denominator_label = spec$denominator_label,
    denominator_cells = n_denom, candidate_ge3_cells = n_candidate,
    candidate_ge3_percent_of_denominator = pct_candidate,
    ENG_THY1_core_within_ge3_cells = n_core,
    ENG_THY1_core_percent_of_candidate_ge3 = pct_core,
    stringsAsFactors = FALSE
  )
  table_outputs <- c(save_tsv(panel_e_data, data_path), save_tsv(summary_table, summary_path))

  # ---- A-E composite using the new legend-right A and existing B/D PNGs ----
  pB <- read_panel_png(prior_panel_paths[[2]])
  pD <- read_panel_png(prior_panel_paths[[3]])
  composite <- ((pA | pB) / pC_render / (pD | pE)) +
    patchwork::plot_annotation(tag_levels = "A", theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))) +
    patchwork::plot_layout(heights = c(1.15, 0.85, 0.85))
  composite_stem <- file.path(out_dir, paste0("Reviewer_Figure_", spec$panel_prefix, "_ASC_P_like_with_panelE_COLORED_LINEAGE_STRIPS_FINAL"))
  composite_outputs <- save_plot(composite, composite_stem, 18, 17)

  # ---- Legend ----
  legend_path <- file.path(out_dir, paste0("Reviewer_Figure_", spec$panel_prefix, "_ASC_P_like_with_panelE_COLORED_LINEAGE_STRIPS_FINAL_legend.md"))
  legend_lines <- c(
    "# Figure legend", "",
    paste0("**Partial ASC_P-like transcriptional candidates in ", spec$tissue_label, ".**"), "",
    "(A) Broad cell-type annotation on the author-provided t-SNE/UMAP embedding. Numeric labels identify author clusters; the fixed annotation legend is shown at the right of the embedding.", "",
    paste0("(B) Localization of cells meeting the partial ASC_P-like candidate definition. Candidates are defined only within the ", spec$denominator_label, " as cells with RNA counts >0 for at least three of six genes (ENG, THY1, NT5E, PDGFRA, VCAM1, MCAM): n = ", n_candidate, " (", sprintf("%.2f", pct_candidate), "% of the denominator)."), "",
    "(C) Lineage-marker context across author clusters. Rows are author-cluster numbers; dot size denotes the percentage of cells with detectable RNA and colour denotes mean normalized RNA. The CAF, Endothelial, Epithelial_tumor, Immune, and Mural headers are visual guides: each header background uses the Panel A colour of the cluster with the highest average expression of that marker set. These header colours do not constitute additional cell-type calls.", "",
    "(D) Six-gene expression context across the prespecified candidate, non-candidate CAF, vascular, and/or mural comparison groups used for this cohort.", "",
    paste0("(E) Nested ENG+THY1 core within the same >=3/6 candidates: RNA-count detection of both ENG and THY1, n = ", n_core, " of ", n_candidate, " (", sprintf("%.2f", pct_core), "%). This is a descriptive transcriptomic summary, not independent validation or evidence of protein-defined CD105+/CD90+ cells. These data support partial ASC_P-like transcriptomic candidates, not an identical human ASC_P niche.")
  )
  assert_new_path(legend_path, cfg$overwrite)
  write_lines_safe(legend_lines, legend_path, cfg$overwrite)

  all_outputs <- c(table_outputs, panel_a_outputs, panel_c_outputs, panel_e_outputs, composite_outputs, legend_path)
  for (f in all_outputs) {
    ok <- is_valid_file(f, min_bytes = 50)
    sz <- if (file.exists(f)) file.info(f)$size else NA
    message(sprintf("  CHECK: %s -> exists=%s size=%s ok=%s", basename(f), file.exists(f), sz, ok))
  }
  if (!all(vapply(setdiff(all_outputs, legend_path), is_valid_file, logical(1), min_bytes = 50)) || !is_valid_file(legend_path, min_bytes = 30)) {
    stop("One or more 03d outputs are missing or empty for ", cohort, call. = FALSE)
  }

  end_step(cfg, state, all_outputs, status = "SUCCESS", notes = c(
    "Plotting-only visual refinement: Panel C numeric cluster axis and Panel-A-coloured lineage strips + Panel E legend fix.",
    sprintf(">=3/6 candidates: %d (%.2f%% of denominator); ENG+THY1 core: %d (%.2f%% of candidates).", n_candidate, pct_candidate, n_core, pct_core),
    "Scientific data unchanged. A-E composite uses a legend-right redraw of Panel A from the same object/map and reuses reviewed B/D PNG panels; C and E are vector sources.",
    "Original RDS and approved maps were not modified."
  ))
  rm(obj, counts, data_mat, meta)
  gc()
  message(sprintf("[03d %s] SUCCESS: all outputs generated.", cohort))
}

message("\n=== All cohorts completed ===\n")

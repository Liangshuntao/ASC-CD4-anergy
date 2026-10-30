# Explicit execution context supplied by run_breast_analysis.R.
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.code_root <- normalizePath(file.path(dirname(.script), ".."), winslash = "/", mustWork = TRUE)
.work_root <- Sys.getenv("ASC_BREAST_WORK_ROOT")
.data_root <- Sys.getenv("ASC_BREAST_DATA_ROOT")
if (!nzchar(.work_root) || !nzchar(.data_root)) stop("Use run_breast_analysis.R with explicit input/output arguments")
#!/usr/bin/env Rscript
# Assemble approved final figure panels only; no scRNA analysis is rerun.

suppressPackageStartupMessages({
  library(grid)
  library(png)
  library(digest)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
script_path <- normalizePath(file.path(.code_root,"stages","09_assemble_cohort_panels.R"), winslash = "/", mustWork = TRUE)
project_dir <- .work_root
figure_dir <- file.path(project_dir, "figures", "reviewer_ASC_P")
hash_file <- function(path) digest(file = path, algo = "sha256")
read_panel <- function(path) {
  if (!file.exists(path)) stop("Missing required panel: ", path, call. = FALSE)
  png::readPNG(path)
}
draw_panel <- function(img) {
  grid.rect(gp = gpar(fill = "white", col = NA))
  ar <- dim(img)[2] / dim(img)[1]
  panel_ar <- convertWidth(unit(1, "npc"), "in", valueOnly = TRUE) / convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  if (ar >= panel_ar) {
    grid.raster(img, width = unit(1, "npc"), height = unit(panel_ar / ar, "npc"), interpolate = TRUE)
  } else {
    grid.raster(img, width = unit(ar / panel_ar, "npc"), height = unit(1, "npc"), interpolate = TRUE)
  }
}
draw_panel_compact_x <- function(img, width_fraction = 0.78) {
  # Preserve the original panel aspect ratio; this narrows the D-row footprint
  # without non-uniformly deforming DotPlot circles or source-panel text.
  pushViewport(viewport(x = 0.5, y = 0.5, width = unit(width_fraction, "npc"), height = unit(1, "npc")))
  draw_panel(img)
  popViewport()
}
draw_abde <- function(output, device = c("png", "pdf", "tiff")) {
  device <- match.arg(device)
  if (device == "png") ragg::agg_png(output, width = 10800, height = 7560, units = "px", res = 600)
  else if (device == "tiff") ragg::agg_tiff(output, width = 10800, height = 7560, units = "px", res = 600, compression = "lzw")
  else pdf(output, width = 18, height = 12.6, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  cohorts <- c("Normal", "ER-positive", "HER2-positive", "TNBC")
  cohort_keys <- c("NormalTotal", "ERTotal", "HER2", "TNBC")
  row_labels <- c("A", "B", "C", "D")
  row_keys <- c("A", "B", "D", "E")
  files <- list(
    A = c(
      NormalTotal = "Figure_A_NormalTotal_fixed_annotation_LEGEND_RIGHT_FINAL.png",
      ERTotal = "Figure_A_ERTotal_fixed_annotation_LEGEND_RIGHT_FINAL.png",
      HER2 = "Figure_A_HER2_fixed_annotation_LEGEND_RIGHT_FINAL.png",
      TNBC = "Figure_A_TNBC_fixed_annotation_LEGEND_RIGHT_FINAL.png"),
    B = c(
      NormalTotal = "Figure_C_NormalTotal_ASC_P_like_candidate_map_FINAL.png",
      ERTotal = "Figure_B_ERTotal_ASC_P_like_candidate_map_FINAL.png",
      HER2 = "Figure_B_HER2_ASC_P_like_candidate_map_FINAL.png",
      TNBC = "Figure_B_TNBC_ASC_P_like_candidate_map_FINAL.png"),
    D = c(
      NormalTotal = "Figure_D_NormalTotal_ASC_P_four_group_DotPlot_FINAL.png",
      ERTotal = "Figure_D_ERTotal_ASC_P_four_group_DotPlot_FINAL.png",
      HER2 = "Figure_D_HER2_ASC_P_three_group_DotPlot_FINAL.png",
      TNBC = "Figure_D_TNBC_ASC_P_three_group_DotPlot_FINAL.png"),
    E = c(
      NormalTotal = "Figure_E_NormalTotal_ASC_P_ENG_THY1_within_3of6_LEGEND_FIXED_REVISED_FINAL_CORE_ONLY_RED_FINAL.png",
      ERTotal = "Figure_E_ERTotal_ASC_P_ENG_THY1_within_3of6_LEGEND_FIXED_REVISED_FINAL_CORE_ONLY_RED_FINAL.png",
      HER2 = "Figure_E_HER2_ASC_P_ENG_THY1_within_3of6_LEGEND_FIXED_REVISED_FINAL_CORE_ONLY_RED_FINAL.png",
      TNBC = "Figure_E_TNBC_ASC_P_ENG_THY1_within_3of6_LEGEND_FIXED_REVISED_FINAL_CORE_ONLY_RED_FINAL.png"))
  layout <- grid.layout(5, 5, widths = unit.c(unit(0.6, "in"), rep(unit(4.35, "in"), 4)), heights = unit.c(unit(0.6, "in"), rep(unit(3, "in"), 4)))
  pushViewport(viewport(layout = layout))
  for (j in seq_along(cohorts)) {
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = j + 1)); grid.text(cohorts[j], gp = gpar(fontface = "bold", fontsize = 28)); popViewport()
  }
  for (i in seq_along(row_keys)) {
    pushViewport(viewport(layout.pos.row = i + 1, layout.pos.col = 1)); grid.text(row_labels[i], gp = gpar(fontface = "bold", fontsize = 36)); popViewport()
    for (j in seq_along(cohort_keys)) {
      path <- file.path(figure_dir, files[[row_keys[i]]][[cohort_keys[j]]])
      pushViewport(viewport(layout.pos.row = i + 1, layout.pos.col = j + 1))
      if (row_keys[i] == "D") draw_panel_compact_x(read_panel(path)) else draw_panel(read_panel(path))
      popViewport()
    }
  }
  popViewport()
  invisible(unlist(files, use.names = FALSE))
}
draw_c <- function(output, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") ragg::agg_png(output, width = 7200, height = 12300, units = "px", res = 600)
  else pdf(output, width = 12, height = 20.5, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  cohorts <- c("NormalTotal", "ERTotal", "HER2", "TNBC")
  files <- c(
    NormalTotal = "Figure_C_NormalTotal_gene_level_lineage_DotPlot_COLORED_LINEAGE_STRIPS_REVISED_FINAL.png",
    ERTotal = "Figure_C_ERTotal_gene_level_lineage_DotPlot_COLORED_LINEAGE_STRIPS_REVISED_FINAL.png",
    HER2 = "Figure_C_HER2_gene_level_lineage_DotPlot_COLORED_LINEAGE_STRIPS_REVISED_FINAL.png",
    TNBC = "Figure_C_TNBC_gene_level_lineage_DotPlot_COLORED_LINEAGE_STRIPS_REVISED_FINAL.png")
  layout <- grid.layout(5, 2, widths = unit.c(unit(0.45, "in"), unit(11.55, "in")), heights = unit.c(unit(0.5, "in"), rep(unit(5, "in"), 4)))
  pushViewport(viewport(layout = layout))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2)); grid.text("C. Lineage-marker context across author clusters", x = unit(0, "npc"), just = "left", gp = gpar(fontface = "bold", fontsize = 14)); popViewport()
  for (i in seq_along(cohorts)) {
    pushViewport(viewport(layout.pos.row = i + 1, layout.pos.col = 1)); grid.text(cohorts[i], rot = 90, gp = gpar(fontface = "bold", fontsize = 11)); popViewport()
    pushViewport(viewport(layout.pos.row = i + 1, layout.pos.col = 2)); draw_panel(read_panel(file.path(figure_dir, files[[cohorts[i]]]))); popViewport()
  }
  popViewport()
  invisible(unname(files))
}

abde_png <- file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_ABDE_summary.png")
abde_pdf <- file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_ABDE_summary.pdf")
abde_tiff <- file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_ABDE_summary.tiff")
c_png <- file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_C_lineage_summary.png")
c_pdf <- file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_C_lineage_summary.pdf")
abde_inputs <- draw_abde(abde_png, "png"); draw_abde(abde_pdf, "pdf"); draw_abde(abde_tiff, "tiff")
c_inputs <- draw_c(c_png, "png"); draw_c(c_pdf, "pdf")

abde_legend <- c(
  "# Figure legend",
  "",
  "**Four-cohort summary of descriptive ASC_P-like transcriptomic localization in GSE161529.**",
  "",
  "Columns show NormalTotal, ERTotal, HER2 and TNBC, respectively. All cohorts use the same descriptive candidate rule within their preapproved fibroblast/CAF denominator: RNA-count detection of at least three of ENG, THY1, NT5E, PDGFRA, VCAM1 and MCAM.",
  "",
  "(A) Author-provided t-SNE embeddings with fixed broad cell-type annotation. (B) Localization of partial ASC_P-like transcriptomic candidates. (C) Six-gene expression across prespecified candidate, non-candidate fibroblast/CAF, vascular and, where applicable, mural comparison groups. (D) ENG+THY1 RNA-count co-detection nested within the same >=3/6 candidate set.",
  "",
  "These figures support partial ASC_P-like transcriptomic candidates only. They do not demonstrate protein-defined CD105+/CD90+ ASC cells, an identical human ASC_P niche, or conserved mechanism." )
c_legend <- c(
  "# Figure legend",
  "",
  "**Lineage-marker context across author clusters in four GSE161529 cohorts.**",
  "",
  "Rows show NormalTotal, ERTotal, HER2 and TNBC, respectively. Dot size denotes the percentage of cells with detectable RNA and colour denotes mean normalized RNA. CAF, Endothelial, Epithelial_tumor, Immune and Mural header strips are visual guides inherited from the final cohort-specific panels; strip colours do not constitute additional cell-type calls.",
  "",
  "The figure provides descriptive lineage context for the candidate analysis and does not establish a protein-defined ASC population or an identical human ASC_P niche." )
writeLines(abde_legend, file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_ABDE_summary_legend.md"))
writeLines(c_legend, file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_C_lineage_summary_legend.md"))

manifest <- data.frame(
  record_type = c(rep("input", length(c(abde_inputs, c_inputs))), rep("output", 7)),
  path = c(file.path(figure_dir, c(abde_inputs, c_inputs)), abde_png, abde_pdf, abde_tiff, c_png, c_pdf,
           file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_ABDE_summary_legend.md"), file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_C_lineage_summary_legend.md")),
  sha256 = NA_character_, stringsAsFactors = FALSE)
manifest$sha256 <- vapply(manifest$path, hash_file, character(1))
manifest$generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
manifest$script_sha256 <- hash_file(script_path)
write.table(manifest, file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_summary_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(figure_dir, "Figure_GSE161529_ASC_P_like_four_cohort_summary_sessionInfo.txt"))

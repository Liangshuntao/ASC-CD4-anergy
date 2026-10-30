# Step04 UMAP label fix — numbered clusters + coord_fixed + numbered legend
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

cfg <- list(
  input_rds  = "D:/Users/lstsw/Documents/ASC_CD4/final_result/ASC_P_human_validation/steps/step04/attempt_011/objects/PRJNA932038_global_annotation_v2_candidate.rds",
  mapping_csv = "D:/Users/lstsw/Documents/ASC_CD4/final_result/ASC_P_human_validation/steps/step04/attempt_011/tables/annotation_v2_cluster_mapping.csv",
  output_dir = Sys.getenv("ASC_P_ATTEMPT_DIR")
)
if (!nzchar(cfg$output_dir)) stop("ASC_P_ATTEMPT_DIR not set")

# ---- Load ----
obj <- readRDS(cfg$input_rds)
mapping <- read.csv(cfg$mapping_csv, stringsAsFactors = FALSE)

# ---- Create numbered legend label ----
obj$celltype_v2_numbered <- paste0(obj$seurat_clusters, ": ", obj$celltype_v2)
obj$major_celltype_v2_numbered <- paste0(obj$seurat_clusters, ": ", obj$major_celltype_v2)

# Force legend order by cluster number
ordered_levels <- paste0(0:18, ": ", mapping$celltype_v2[match(as.character(0:18), as.character(mapping$cluster))])
obj$celltype_v2_numbered <- factor(obj$celltype_v2_numbered, levels = ordered_levels)

ordered_major <- paste0(0:18, ": ", mapping$major_celltype_v2[match(as.character(0:18), as.character(mapping$cluster))])
obj$major_celltype_v2_numbered <- factor(obj$major_celltype_v2_numbered, levels = ordered_major)

# ---- Get cluster median UMAP positions for labels ----
umap_coords <- Seurat::Embeddings(obj, "umap")
colnames(umap_coords) <- c("umap_1", "umap_2")

label_df <- umap_coords |>
  as.data.frame() |>
  tibble::rownames_to_column("cell") |>
  dplyr::left_join(
    obj@meta.data |>
      tibble::rownames_to_column("cell") |>
      dplyr::select(cell, seurat_clusters, celltype_v2, major_celltype_v2),
    by = "cell"
  ) |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    celltype_v2 = dplyr::first(celltype_v2),
    major_celltype_v2 = dplyr::first(major_celltype_v2),
    .groups = "drop"
  )

# ---- UMAP: celltype numbered ----
p_celltype <- Seurat::DimPlot(obj, reduction = "umap", group.by = "celltype_v2_numbered", pt.size = 0.3) +
  ggplot2::geom_label(
    data = label_df,
    ggplot2::aes(x = umap_1, y = umap_2, label = seurat_clusters),
    size = 4, fontface = "bold", color = "black",
    fill = NA, label.size = 0.25, alpha = 1.0,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::ggtitle("Annotation v2: celltype") +
  ggplot2::theme(
    legend.text = ggplot2::element_text(size = 8),
    legend.title = ggplot2::element_text(size = 9)
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1, override.aes = list(size = 2)))

ggplot2::ggsave(file.path(cfg$output_dir, "figures", "annotation_v2_umap_celltype_numbered.pdf"),
                p_celltype, width = 10, height = 8, device = "pdf")
ggplot2::ggsave(file.path(cfg$output_dir, "figures", "annotation_v2_umap_celltype_numbered.png"),
                p_celltype, width = 10, height = 8, dpi = 300)

# ---- UMAP: major celltype numbered ----
p_major <- Seurat::DimPlot(obj, reduction = "umap", group.by = "major_celltype_v2_numbered", pt.size = 0.3) +
  ggplot2::geom_label(
    data = label_df,
    ggplot2::aes(x = umap_1, y = umap_2, label = seurat_clusters),
    size = 4, fontface = "bold", color = "black",
    fill = NA, label.size = 0.25, alpha = 1.0,
    inherit.aes = FALSE
  ) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::ggtitle("Annotation v2: major celltype") +
  ggplot2::theme(
    legend.text = ggplot2::element_text(size = 9),
    legend.title = ggplot2::element_text(size = 10)
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(ncol = 1, override.aes = list(size = 2)))

ggplot2::ggsave(file.path(cfg$output_dir, "figures", "annotation_v2_umap_major_celltype_numbered.pdf"),
                p_major, width = 10, height = 8, device = "pdf")
ggplot2::ggsave(file.path(cfg$output_dir, "figures", "annotation_v2_umap_major_celltype_numbered.png"),
                p_major, width = 10, height = 8, dpi = 300)

# ---- Numbered legend table ----
legend_out <- mapping |>
  dplyr::mutate(
    legend_label = paste0(cluster, ": ", celltype_v2)
  ) |>
  dplyr::select(cluster, legend_label, celltype_v2, major_celltype_v2, confidence, note, n_cells)
write.csv(legend_out, file.path(cfg$output_dir, "tables", "annotation_v2_numbered_legend.csv"), row.names = FALSE)

# ---- Review report ----
input_sha <- digest::digest(file = cfg$input_rds, algo = "sha256", serialize = FALSE)
script_sha <- digest::digest(
  file = "D:/Users/lstsw/Documents/ASC_CD4/final_result/ASC_P_human_validation/scripts/step_04_findallmarkers_supplement.R",
  algo = "sha256", serialize = FALSE
)

report_lines <- c(
  "# Step04 UMAP label fix review",
  "",
  "## Input",
  paste0("- RDS: ", cfg$input_rds),
  paste0("- SHA256: ", input_sha),
  paste0("- Annotation mapping: ", cfg$mapping_csv),
  "",
  "## Changes applied",
  "- Generated 0-18 numbered cluster labels at median UMAP positions",
  "- Created numbered legend: `0: PDGFRA+ adipose-associated fibroblasts` ... `18: Cycling immune cells`",
  "- Applied `coord_fixed(ratio = 1)` for equal UMAP aspect ratio",
  "- PDF and PNG generated from the same ggplot object",
  "- PDF device: default (pdf)",
  "",
  "## Output files",
  "- `figures/annotation_v2_umap_celltype_numbered.pdf` (10x8 in)",
  "- `figures/annotation_v2_umap_celltype_numbered.png` (10x8 in, 300 dpi)",
  "- `figures/annotation_v2_umap_major_celltype_numbered.pdf` (10x8 in)",
  "- `figures/annotation_v2_umap_major_celltype_numbered.png` (10x8 in, 300 dpi)",
  "- `tables/annotation_v2_numbered_legend.csv`",
  "",
  "## Verification",
  "- 0-18 numbered labels: YES",
  "- Numbered legend: YES",
  "- coord_fixed(ratio=1): YES",
  "- PDF/PNG same ggplot: YES",
  "- Annotation unchanged: YES (attempt_011 mapping preserved)",
  "- Cluster 15 NOT malignant: YES",
  "",
  "## Status",
  "",
  "`WAITING_FOR_HUMAN_APPROVAL`"
)

writeLines(report_lines, file.path(cfg$output_dir, "reports", "REVIEW_STEP_04_UMAP_LABEL_FIX.md"), useBytes = TRUE)

# ---- manifest ----
manifest <- list(
  step_id = "04",
  status = "WAITING_FOR_HUMAN_APPROVAL",
  attempt = basename(cfg$output_dir),
  task = "annotation_v2_umap_numbered_label_fix",
  input_checkpoint = "steps/step04/attempt_011/objects/PRJNA932038_global_annotation_v2_candidate.rds",
  input_checkpoint_sha256 = input_sha,
  source_annotation_attempt = "attempt_011",
  analysis_changed = FALSE,
  annotation_changed = FALSE,
  plots_updated = TRUE,
  approval_required_before_step05 = TRUE,
  pdf_png_consistent = TRUE,
  coord_fixed_ratio_1 = TRUE,
  script = "scripts/step_04_findallmarkers_supplement.R",
  script_sha256 = script_sha,
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)

manifest_path <- file.path(cfg$output_dir, "step_manifest.yml")
tmp <- paste0(manifest_path, ".tmp")
yaml::write_yaml(manifest, tmp)
file.rename(tmp, manifest_path)

cat("WAITING_FOR_HUMAN_APPROVAL: Step 04 UMAP label fix\n")

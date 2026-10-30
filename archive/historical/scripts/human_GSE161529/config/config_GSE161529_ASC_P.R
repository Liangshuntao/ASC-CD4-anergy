# GSE161529 ASC_P-like validation: centralized configuration
# This project is descriptive. Cells are never treated as biological replicates.

.root_override <- Sys.getenv("GSE161529_ASC_P_PROJECT_ROOT", unset = "")
if (nzchar(.root_override)) {
  .project_root <- normalizePath(.root_override, winslash = "/", mustWork = TRUE)
  .config_file <- file.path(.project_root, "config", "config_GSE161529_ASC_P.R")
} else {
  .config_file <- normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE)
  .project_root <- normalizePath(file.path(dirname(.config_file), ".."), winslash = "/", mustWork = TRUE)
}
.public_analysis_root <- normalizePath(file.path(.project_root, ".."), winslash = "/", mustWork = TRUE)

cfg <- list(
  project_root = .project_root,
  input_dir = normalizePath(
    file.path(.public_analysis_root, "data", "processed", "GSE161529_figshare"),
    winslash = "/", mustWork = TRUE
  ),
  author_asset_dir = file.path(.project_root, "inputs", "author_assets"),
  author_code_dir = file.path(.project_root, "provenance", "author_code"),
  updated_object_dir = file.path(.project_root, "objects", "updated_copies"),
  table_dir = file.path(.project_root, "tables"),
  figure_dir = file.path(.project_root, "figures"),
  report_dir = file.path(.project_root, "reports"),
  log_dir = file.path(.project_root, "logs"),
  manifest_dir = file.path(.project_root, "manifests"),
  config_dir = file.path(.project_root, "config"),
  figshare_article_id = "17058077",
  figshare_api_url = "https://api.figshare.com/v2/articles/17058077",
  github_repo = "https://github.com/yunshun/HumanBreast10X",
  github_raw_roots = c(
    "https://raw.githubusercontent.com/yunshun/HumanBreast10X/main/RCode",
    "https://raw.githubusercontent.com/yunshun/HumanBreast10X/master/RCode"
  ),
  figshare_article_url = "https://figshare.com/articles/dataset/Data_R_code_and_output_Seurat_Objects_for_single_cell_RNA-seq_analysis_of_human_breast_tissues/17058077",
  rscript = "E:/R/R-4.3.3/bin/Rscript.exe",
  species = "Homo sapiens",
  genome_build = "GRCh38 (Cell Ranger human reference v3.0.0)",
  statistical_unit = "patient/sample",
  inference = "descriptive_not_inferential",
  patient_candidates = c("patient_id", "patient", "donor", "group"),
  sample_candidates = c("sample_id", "sample", "orig.ident", "group"),
  cluster_candidates = c("seurat_clusters", "cluster", "integrated_snn_res.0.1", "integrated_snn_res.0.07", "integrated_snn_res.0.05"),
  celltype_candidates = c("author_celltype", "celltype", "cell_type", "CellType", "annotation", "cell_annotation"),
  approval_file = file.path(.project_root, "config", "cluster_boundary_approval.yml"),
  # Active mapping is the reviewed draft; the original CodeBuddy output remains preserved.
  mapping_file = file.path(.project_root, "tables", "01_author_cluster_annotation_map_revised_v5.tsv"),
  generated_mapping_file = file.path(.project_root, "tables", "01_author_cluster_annotation_map_generated.tsv"),
  candidate_cache = file.path(.project_root, "objects", "03_candidate_cell_metadata.rds"),
  prjna932038_dir = "D:/Users/lstsw/Documents/ASC_CD4/final_result/ASC_P_human_validation",
  seed = 1234L,
  overwrite = FALSE,
  boundary_pct_threshold = 10,
  boundary_min_markers = 2L,
  dominance_threshold = 0.70,
  plot_width = 9,
  plot_height = 6,
  plot_dpi = 300
)

cfg$signature_genes <- c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM")
cfg$fibroblast_markers <- c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "DPT", "SFRP2", "PDGFRA", "COL10A1", "COMP")
cfg$mural_markers <- c("PDGFRB", "RGS5", "CSPG4", "ACTA2", "TAGLN", "MYH11", "NOTCH3", "MCAM", "CAVIN3")
cfg$endothelial_markers <- c("PECAM1", "VWF", "CLDN5", "EMCN", "KDR", "PLVAP", "ESAM", "RAMP2")
cfg$contamination_markers <- c("PTPRC", "EPCAM", "KRT8", "KRT18", "KRT19")
cfg$all_audit_genes <- unique(c(cfg$signature_genes, cfg$fibroblast_markers, cfg$mural_markers,
                                cfg$endothelial_markers, cfg$contamination_markers))

cfg$local_objects <- data.frame(
  object_id = c("NormTotal", "NormTotalFib", "ERTotal", "HER2", "TNBC"),
  subtype = c("Normal", "Normal", "ER+", "HER2+", "TNBC"),
  role = c("total", "fibroblast", "total", "total", "total"),
  filename = c("SeuratObject_NormTotal.rds", "SeuratObject_NormTotalFib.rds",
               "SeuratObject_ERTotal.rds", "SeuratObject_HER2.rds", "SeuratObject_TNBC.rds"),
  stringsAsFactors = FALSE
)

cfg$required_author_assets <- data.frame(
  object_id = c("ERTotalSub", "HER2Sub", "TNBCSub"),
  subtype = c("ER+", "HER2+", "TNBC"),
  filename = c("SeuratObject_ERTotalSub.rds", "SeuratObject_HER2Sub.rds", "SeuratObject_TNBCSub.rds"),
  total_object_id = c("ERTotal", "HER2", "TNBC"),
  code_file = c("ER.R", "HER2.R", "TNBC.R"),
  stringsAsFactors = FALSE
)

cfg$author_code_files <- c("NormTotal.R", "ER.R", "HER2.R", "TNBC.R")
cfg$analysis_order <- c("NormalFib", "ER", "HER2", "TNBC")

rm(.config_file, .project_root, .public_analysis_root, .root_override)

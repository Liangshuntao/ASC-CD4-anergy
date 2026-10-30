# Materialize and document the selected Phase 1 Harmony solution.
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(harmony)
  library(ggplot2); library(patchwork); library(digest)
})

project_root <- "D:/Users/lstsw/Documents/Adipose_derived stromal_metabolic niches/adipose_stromal_metabolic_niche_lit/public_data_analysis/GSE161529_ASC_P_validation"
input_rds <- file.path(project_root, "outputs", "ERTotal_cluster1_TNK_reclustering_harmony", "ERTotal_cluster1_phase1_harmony.rds")
output_dir <- file.path(project_root, "outputs", "ERTotal_cluster1_TNK_Harmony_optimized")
grid_dir <- output_dir
script_file <- file.path(project_root, "scripts", "05b7_ERTotal_cluster1_Harmony_optimized_finalize.R")
if (!file.exists(input_rds)) stop("Input RDS missing: ", input_rds)
selected_file <- file.path(grid_dir, "ERTotal_cluster1_optimized_selected_parameters.tsv")
if (!file.exists(selected_file)) stop("Run 05b6 first; selected-parameter file is missing.")

write_tsv <- function(x, name) write.table(x, file.path(output_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
selected <- read.delim(selected_file, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(selected) != 1L) stop("Selected-parameter table must contain exactly one row.")
dims <- as.integer(selected$dims); k <- as.integer(selected$k_param); resolution <- as.numeric(selected$resolution)
seed <- 20260716L
set.seed(seed)

obj0 <- readRDS(input_rds)
counts <- LayerData(obj0, assay = "RNA", layer = "counts")
meta0 <- obj0[[]]
patient <- as.character(meta0[colnames(counts), "patient_id"])
names(patient) <- colnames(counts)
if (anyNA(patient)) stop("patient_id is missing for one or more cells.")

sets <- list(
  pan_t = c("CD3D", "CD3E", "TRAC"), helper = c("IL7R", "CCR7", "TCF7", "MAL", "CD40LG"),
  treg = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF18", "IKZF2"), cd8 = c("CD8A", "CD8B"),
  cytotoxic = c("CCL5", "GZMK", "GZMB", "PRF1", "NKG7")
)
if (!all(unlist(sets) %in% rownames(counts))) stop("One or more anchor genes are absent from the RNA assay.")
det <- counts[unique(unlist(sets)), , drop = FALSE] > 0
n_detect <- function(genes) Matrix::colSums(det[genes, , drop = FALSE])
pan_n <- n_detect(sets$pan_t); helper_n <- n_detect(sets$helper); treg_n <- n_detect(sets$treg)
cd8_n <- n_detect(sets$cd8); cyto_n <- n_detect(sets$cytotoxic)
anchor3 <- rep("Other", ncol(counts)); names(anchor3) <- colnames(counts)
anchor3[pan_n >= 1L & cd8_n >= 1L & cyto_n >= 2L & helper_n <= 1L & treg_n < 3L] <- "CD8_cytotoxic"
anchor3[pan_n >= 1L & helper_n >= 2L & cyto_n <= 1L & cd8_n == 0L & treg_n < 3L] <- "CD4_helper"
anchor3[pan_n >= 1L & treg_n >= 3L] <- "Treg"

obj <- CreateSeuratObject(counts = counts, meta.data = data.frame(patient_id = patient, row.names = colnames(counts)), project = "ERTotal_cluster1_optimized")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000L, verbose = FALSE)
hvg <- VariableFeatures(obj); upper <- toupper(hvg)
cc <- toupper(unique(c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)))
exclude <- grepl("^MT-", upper) | grepl("^RPL", upper) | grepl("^RPS", upper) | grepl("^IG[HKL]", upper) | grepl("^TR[ABDG][VJ]", upper) | upper %in% cc
pca_features <- hvg[!exclude]
if (length(pca_features) < 1000L) stop("Too few retained PCA features.")
obj <- ScaleData(obj, features = pca_features, verbose = FALSE)
obj <- RunPCA(obj, features = pca_features, npcs = 50L, seed.use = seed, verbose = FALSE)
set.seed(seed + dims * 100L)
obj <- RunHarmony(obj, group.by.vars = "patient_id", reduction.use = "pca", dims.use = seq_len(dims), reduction.save = "harmony_optimized", max_iter = 20L, plot_convergence = FALSE, verbose = FALSE)
set.seed(seed + dims * 100L + k)
obj <- FindNeighbors(obj, reduction = "harmony_optimized", dims = seq_len(dims), k.param = k, graph.name = c("optimized_nn", "optimized_snn"), verbose = FALSE)
obj <- FindClusters(obj, graph.name = "optimized_snn", resolution = resolution, algorithm = 1, random.seed = seed, verbose = FALSE)
cluster_col <- paste0("optimized_snn_res.", format(resolution, nsmall = 1L))
if (!cluster_col %in% colnames(obj[[]])) stop("Selected clustering column is missing: ", cluster_col)
obj$optimized_cluster <- factor(obj[[cluster_col, drop = TRUE]])
obj$anchor_state_for_evaluation <- factor(anchor3, levels = c("CD4_helper", "Treg", "CD8_cytotoxic", "Other"))
obj <- RunUMAP(obj, reduction = "harmony_optimized", dims = seq_len(dims), reduction.name = "umap_optimized", reduction.key = "UMAPOPT_", seed.use = seed, verbose = FALSE)

cluster_sizes <- as.data.frame(table(cluster = obj$optimized_cluster), stringsAsFactors = FALSE); names(cluster_sizes)[2] <- "n_cells"
write_tsv(cluster_sizes, "ERTotal_cluster1_optimized_selected_cluster_sizes.tsv")
composition <- as.data.frame.matrix(table(cluster = obj$optimized_cluster, anchor_state = obj$anchor_state_for_evaluation), stringsAsFactors = FALSE)
composition$cluster <- rownames(composition); rownames(composition) <- NULL
composition <- composition[, c("cluster", setdiff(names(composition), "cluster")), drop = FALSE]
write_tsv(composition, "ERTotal_cluster1_optimized_selected_anchor_composition.tsv")
patient_qc <- do.call(rbind, lapply(levels(obj$optimized_cluster), function(cl) {
  ix <- which(as.character(obj$optimized_cluster) == cl); x <- sort(table(patient[ix]), decreasing = TRUE)
  data.frame(cluster = cl, total_cells = length(ix), contributing_patients = length(x), dominant_patient = names(x)[1], dominant_patient_fraction = as.numeric(x[1]) / length(ix), stringsAsFactors = FALSE)
}))
write_tsv(patient_qc, "ERTotal_cluster1_optimized_selected_patient_QC.tsv")

Idents(obj) <- "optimized_cluster"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25, verbose = FALSE)
write_tsv(markers, "ERTotal_cluster1_optimized_selected_cluster_markers.tsv")
top_markers <- do.call(rbind, lapply(split(markers, markers$cluster), function(x) head(x[order(-x$avg_log2FC, x$p_val_adj), , drop = FALSE], 25L)))
write_tsv(top_markers, "ERTotal_cluster1_optimized_selected_cluster_top25_markers.tsv")

p_cluster <- DimPlot(obj, reduction = "umap_optimized", group.by = "optimized_cluster", label = TRUE, repel = TRUE, label.size = 3.3, pt.size = 0.25) + NoLegend() + coord_fixed() + ggtitle("A  Optimized Harmony clusters") + theme(plot.title = element_text(face = "bold"))
p_anchor <- DimPlot(obj, reduction = "umap_optimized", group.by = "anchor_state_for_evaluation", cols = c("#2C7FB8", "#7A0177", "#D95F0E", "grey82"), pt.size = 0.25) + coord_fixed() + ggtitle("B  Predefined evaluation anchors") + theme(plot.title = element_text(face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8))
dot_genes <- intersect(c("CD3D", "TRAC", "IL7R", "CCR7", "TCF7", "MAL", "CD40LG", "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "CD8A", "CD8B", "CCL5", "GZMK", "NKG7", "PRF1"), rownames(obj))
p_dot <- DotPlot(obj, features = dot_genes, group.by = "optimized_cluster", cols = c("grey92", "#B2182B"), dot.scale = 5) + RotatedAxis() + ggtitle("C  Canonical marker distribution") + xlab("Genes") + ylab("Optimized cluster") + theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 55, hjust = 1, size = 8), legend.position = "right")
main_fig <- (p_cluster | p_anchor | p_dot) + plot_layout(widths = c(1.1, 1.25, 1.55))
ggsave(file.path(output_dir, "Reviewer_Figure_ERTotal_cluster1_Harmony_optimized.png"), main_fig, width = 16, height = 5.2, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "Reviewer_Figure_ERTotal_cluster1_Harmony_optimized.pdf"), main_fig, width = 16, height = 5.2, bg = "white")

grid <- read.delim(file.path(grid_dir, "ERTotal_cluster1_optimized_parameter_grid_metrics.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
grid <- grid[grid$candidate_id != "baseline_HVG2000_CCincluded_PC20_k20", , drop = FALSE]
selected_plot <- selected[, c("dims", "k_param", "resolution"), drop = FALSE]
heat <- function(metric, title) {
  ggplot(grid, aes_string(x = "factor(k_param)", y = "factor(dims)", fill = metric)) + geom_tile(color = "white", linewidth = 0.5) + geom_tile(data = selected_plot, aes(x = factor(k_param), y = factor(dims)), fill = NA, color = "black", linewidth = 1.15, inherit.aes = FALSE) + facet_wrap(~ resolution, nrow = 1) + scale_fill_viridis_c(option = "C", name = metric) + labs(title = title, x = "k-nearest neighbors", y = "Harmony PCs") + theme_bw(base_size = 9) + theme(strip.background = element_rect(fill = "grey94"), panel.grid = element_blank(), plot.title = element_text(face = "bold"))
}
grid_fig <- heat("broad_anchor_ARI", "A  Broad CD4/CD8 anchor concordance") / heat("anchor3_ARI", "B  Three-state anchor concordance")
ggsave(file.path(output_dir, "Reviewer_Figure_ERTotal_cluster1_Harmony_optimized_grid.png"), grid_fig, width = 15, height = 6, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "Reviewer_Figure_ERTotal_cluster1_Harmony_optimized_grid.pdf"), grid_fig, width = 15, height = 6, bg = "white")

saveRDS(obj, file.path(output_dir, "ERTotal_cluster1_Harmony_optimized.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "ERTotal_cluster1_optimized_finalize_sessionInfo.txt"))
writeLines(c("# Optimized Phase 1 Harmony solution", "", paste0("Input SHA256: ", digest(file = input_rds, algo = "sha256")), paste0("Grid script SHA256: ", digest(file = file.path(project_root, "scripts", "05b6_ERTotal_cluster1_Harmony_parameter_optimization.R"), algo = "sha256")), paste0("Finalize script SHA256: ", digest(file = script_file, algo = "sha256")), paste0("Selected configuration: ", selected$candidate_id, "; resolution=", resolution, "; clusters=", nrow(cluster_sizes)), "Representation: 3,000 VST HVGs; mitochondrial, ribosomal, immunoglobulin, TCR V/J, and cell-cycle genes excluded from PCA; Harmony corrected patient_id.", "Evaluation anchors are predefined marker-detection rules used only to compare configurations. They are not final cell-type labels.", "This Phase 1 result is a lineage-resolution optimization. Formal doublet assessment and final biological annotation remain separate checkpoints."), file.path(output_dir, "ERTotal_cluster1_optimized_FINAL_RUN_REPORT.md"))
cat("Finalized optimized object: ", file.path(output_dir, "ERTotal_cluster1_Harmony_optimized.rds"), "\n", sep = "")

# Explicit execution context supplied by run_breast_analysis.R.
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.code_root <- normalizePath(file.path(dirname(.script), ".."), winslash = "/", mustWork = TRUE)
.work_root <- Sys.getenv("ASC_BREAST_WORK_ROOT")
.data_root <- Sys.getenv("ASC_BREAST_DATA_ROOT")
if (!nzchar(.work_root) || !nzchar(.data_root)) stop("Use run_breast_analysis.R with explicit input/output arguments")
#!/usr/bin/env Rscript
# Read-only, patient/sample-level malignant-epithelial pseudobulk resemblance analysis.
# It does not assign clinical PAM50 subtypes and never overwrites any source RDS.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(script_file) || !nzchar(script_file)) stop("Run with Rscript so --file is available.", call. = FALSE)
project_root <- .work_root
helper_path <- file.path(.code_root, "lib", "common.R")
source(helper_path)

# ============================= USER CONFIGURATION =============================
# Replace UNRESOLVED only with values documented in the required-configuration TSV.
# `unit_type` must be "patient" only when the field is demonstrated to be a patient ID;
# otherwise use "sample" with an explicit justification. The script blocks UNRESOLVED.
cfg_pam50 <- list(
  output_dir = file.path(project_root, "outputs", "ThreeTumor_PAM50_LumB_similarity"),
  centroid_tsv = file.path(project_root, "tables", "PAM50_genefu_2.34.0_scale_centroids.tsv"),
  centroid_source = "Bioconductor 3.18 genefu 2.34.0 pam50.scale; source package SHA256 BB55478C0AF19FC6EC9CD750423E94E63A9A6EC56C0CE838EF9784632CE2D663",
  centroid_preprocessing = "genefu pam50.scale: log2(CPM+1), then gene-wise mean centering and sample-SD scaling jointly across eligible pseudobulks; official method.cor=spearman",
  min_cells = c(20L, 50L, 100L), primary_min_cells = 50L,
  uncertain_margin = 0.05, primary_correlation = "spearman",
  correlation_methods = c("pearson", "spearman"), seed = 1234L, write_figures = TRUE,
  datasets = list(
    ERTotal = list(dataset = "ERTotal", rds = file.path(.data_root, "SeuratObject_ERTotal.rds"), annotation_map = file.path(project_root, "tables", "ERTotal_cluster_annotation_map_approved.tsv"), unit_map = file.path(project_root, "tables", "ThreeTumor_sample_patient_mapping.tsv"), cluster_col = "seurat_clusters", patient_or_sample_col = "group", unit_type = "patient", unit_evidence = "Scientific Data 2022 Supplementary Table 1 maps each author ERTotal group/sample to an explicit Patient value.", approved_malignant_labels = c("ER-positive luminal malignant epithelial cells", "Cycling ER-positive malignant epithelial cells")),
    HER2 = list(dataset = "HER2", rds = file.path(.data_root, "SeuratObject_HER2.rds"), annotation_map = file.path(project_root, "tables", "HER2_cluster_annotation_map_approved.tsv"), unit_map = file.path(project_root, "tables", "ThreeTumor_sample_patient_mapping.tsv"), cluster_col = "seurat_clusters", patient_or_sample_col = "group", unit_type = "patient", unit_evidence = "Scientific Data 2022 Supplementary Table 1 maps each author HER2 group/sample to an explicit Patient value.", approved_malignant_labels = c("HER2-positive luminal-like malignant epithelial cells", "Cycling HER2-positive malignant epithelial cells")),
    TNBC = list(dataset = "TNBC", rds = file.path(.data_root, "SeuratObject_TNBC.rds"), annotation_map = file.path(project_root, "tables", "TNBC_cluster_annotation_map_approved.tsv"), unit_map = file.path(project_root, "tables", "ThreeTumor_sample_patient_mapping.tsv"), cluster_col = "seurat_clusters", patient_or_sample_col = "group", unit_type = "patient", unit_evidence = "Scientific Data 2022 Supplementary Table 1 maps each author TNBC group/sample to an explicit Patient value.", approved_malignant_labels = c("Epithelial/tumor cells", "Cycling epithelial/tumor cells"))
  )
)
# ==============================================================================

ensure_packages(c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "scales", "digest"))
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
out <- function(x) file.path(cfg_pam50$output_dir, x)
write_tsv <- function(x, name) utils::write.table(x, out(name), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
programs <- list(Luminal_ER = c("ESR1","PGR","GATA3","FOXA1","KRT8","KRT18","KRT19","MUC1"), HER2 = c("ERBB2","GRB7","STARD3","PGAP3"), Basal_TNBC = c("KRT5","KRT14","KRT17","EGFR","KRT6A","KRT6B"), Proliferation = c("MKI67","TOP2A","UBE2C","BIRC5","CCNB1","CDC20","AURKA","TYMS"))
centroid_cols <- c("LumA", "LumB", "HER2", "Basal", "Normal")

collapse_duplicate_genes <- function(mat) {
  gene <- toupper(trimws(rownames(mat)))
  assert(!anyNA(gene) && all(nzchar(gene)), "RNA matrix contains missing gene symbols.")
  if (!anyDuplicated(gene)) { rownames(mat) <- gene; return(mat) }
  unique_gene <- unique(gene)
  collapse_map <- Matrix::sparseMatrix(i = match(gene, unique_gene), j = seq_along(gene), x = 1,
                                        dims = c(length(unique_gene), length(gene)), dimnames = list(unique_gene, gene))
  collapsed <- collapse_map %*% mat
  colnames(collapsed) <- colnames(mat)
  collapsed
}

normalise_logcpm <- function(pb) {
  lib <- Matrix::colSums(pb)
  assert(all(is.finite(lib) & lib > 0), "Pseudobulk library size is zero or non-finite.")
  # Right multiplication applies each pseudobulk column's own library size.
  norm <- pb %*% Matrix::Diagonal(x = 1e6 / lib)
  dimnames(norm) <- dimnames(pb)
  result <- log2(norm + 1)
  dimnames(result) <- dimnames(pb)
  result
}
check_logcpm <- function() {
  x <- Matrix::Matrix(c(10, 10, 0, 20), nrow = 2, sparse = TRUE)
  y <- normalise_logcpm(x)
  assert(isTRUE(all.equal(as.numeric(y[1, 1]), log2(500001), tolerance = 1e-10)), "Internal CPM column-normalisation check failed.")
  assert(isTRUE(all.equal(as.numeric(y[2, 2]), log2(1000001), tolerance = 1e-10)), "Internal CPM column-normalisation check failed.")
}

read_centroids <- function() {
  assert(!is.na(cfg_pam50$centroid_tsv) && file.exists(cfg_pam50$centroid_tsv), "BLOCKED: reviewed local PAM50 centroid TSV is required.")
  assert(!is.na(cfg_pam50$centroid_source) && nzchar(cfg_pam50$centroid_source), "BLOCKED: centroid source/version is required.")
  assert(!is.na(cfg_pam50$centroid_preprocessing) && nzchar(cfg_pam50$centroid_preprocessing), "BLOCKED: centroid preprocessing specification is required.")
  z <- utils::read.delim(cfg_pam50$centroid_tsv, check.names = FALSE, stringsAsFactors = FALSE)
  assert(identical(names(z), c("gene", centroid_cols)), "Centroid TSV columns must exactly be: gene, LumA, LumB, HER2, Basal, Normal.")
  z$gene <- toupper(trimws(z$gene)); assert(!anyNA(z$gene) && all(nzchar(z$gene)) && !anyDuplicated(z$gene), "Centroid genes must be non-empty and unique.")
  z[, centroid_cols] <- lapply(z[, centroid_cols, drop = FALSE], as.numeric)
  assert(all(vapply(z[, centroid_cols, drop = FALSE], function(x) all(is.finite(x)) && stats::var(x) > 0, logical(1))), "Every centroid must be numeric, finite, and non-constant.")
  assert(nrow(z) == 50L, paste0("Centroid contains ", nrow(z), " unique genes, not the standard 50; PAM50 analysis is blocked."))
  rownames(z) <- z$gene; as.matrix(z[, centroid_cols, drop = FALSE])
}

read_map <- function(d) {
  assert(file.exists(d$annotation_map), paste0(d$dataset, ": annotation map missing."))
  z <- utils::read.delim(d$annotation_map, check.names = FALSE, stringsAsFactors = FALSE)
  assert(all(c("cluster", "celltype_final", "malignant_status") %in% names(z)), paste0(d$dataset, ": map must contain cluster, celltype_final, malignant_status."))
  z$cluster <- as.character(z$cluster); z$celltype_final <- as.character(z$celltype_final); z$malignant_status <- as.character(z$malignant_status)
  assert(!anyDuplicated(z$cluster), paste0(d$dataset, ": duplicated cluster in annotation map.")); z
}

validate_dataset_config <- function(d) {
  assert(file.exists(d$rds), paste0(d$dataset, ": RDS missing."))
  assert(file.exists(d$unit_map), paste0(d$dataset, ": sample-to-patient mapping missing."))
  for (nm in c("cluster_col", "patient_or_sample_col", "unit_type", "unit_evidence")) assert(!is.na(d[[nm]]) && nzchar(d[[nm]]) && d[[nm]] != "UNRESOLVED", paste0("BLOCKED: ", d$dataset, " configuration ", nm, " is UNRESOLVED."))
  assert(d$unit_type %in% c("patient", "sample"), paste0(d$dataset, ": unit_type must be patient or sample."))
  labels <- unique(as.character(d$approved_malignant_labels)); labels <- labels[!is.na(labels) & nzchar(labels)]
  assert(length(labels) > 0L, paste0("BLOCKED: ", d$dataset, " approved_malignant_labels are unresolved.")); labels
}

extract_pseudobulk <- function(d) {
  labels <- validate_dataset_config(d); map <- read_map(d)
  selected <- map[map$celltype_final %in% labels, , drop = FALSE]
  assert(nrow(selected) > 0L && all(selected$malignant_status == "Malignant"), paste0("BLOCKED: every selected ", d$dataset, " map row must have malignant_status == Malignant."))
  obj <- load_seurat_object(d$rds) # helper repairs legacy Seurat v3 only in memory.
  meta <- safe_meta_data(obj); assert(all(c(d$cluster_col, d$patient_or_sample_col) %in% names(meta)), paste0(d$dataset, ": configured metadata columns absent."))
  source_unit <- trimws(as.character(meta[[d$patient_or_sample_col]]))
  assert(!anyNA(source_unit) && all(nzchar(source_unit)) && length(unique(source_unit)) > 1L, paste0(d$dataset, ": source unit field is missing, blank, or constant."))
  unit_map <- utils::read.delim(d$unit_map, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  assert(all(c("dataset", "group_id", "patient_id") %in% names(unit_map)), paste0(d$dataset, ": unit map lacks dataset/group_id/patient_id."))
  unit_map <- unit_map[unit_map$dataset == d$dataset, , drop = FALSE]
  assert(!anyDuplicated(unit_map$group_id) && !anyNA(unit_map$patient_id) && all(nzchar(unit_map$patient_id)), paste0(d$dataset, ": invalid unit map."))
  unit <- unit_map$patient_id[match(source_unit, unit_map$group_id)]
  assert(!anyNA(unit) && all(nzchar(unit)) && length(unique(unit)) > 1L, paste0(d$dataset, ": not every source group maps to a patient."))
  counts <- get_rna_matrix(obj, "counts"); assert(identical(colnames(counts), rownames(meta)), paste0(d$dataset, ": RNA columns and metadata rows differ."))
  assert(inherits(counts, "Matrix") && all(is.finite(counts@x)) && all(counts@x >= 0) && all(counts@x == round(counts@x)), paste0(d$dataset, ": RNA counts must be sparse nonnegative integers."))
  cluster <- as.character(meta[[d$cluster_col]]); assert(all(unique(cluster) %in% map$cluster), paste0(d$dataset, ": annotation map does not cover observed clusters."))
  mapped <- map[match(cluster, map$cluster), , drop = FALSE]; keep <- mapped$celltype_final %in% labels
  assert(any(keep), paste0(d$dataset, ": no cells in approved malignant labels."))
  counts <- collapse_duplicate_genes(counts[, keep, drop = FALSE]); unit <- unit[keep]
  membership <- Matrix::sparseMatrix(i = seq_along(unit), j = match(unit, unique(unit)), x = 1, dims = c(length(unit), length(unique(unit))), dimnames = list(NULL, unique(unit)))
  pb <- counts %*% membership; raw_patient <- colnames(pb); colnames(pb) <- paste0(d$dataset, "__", raw_patient)
  qc <- data.frame(dataset = d$dataset, unit_type = d$unit_type, unit_id = raw_patient, patient_id = colnames(pb), malignant_cells = as.integer(Matrix::colSums(membership)), total_UMI = as.numeric(Matrix::colSums(pb)), detected_genes = as.integer(Matrix::colSums(pb > 0)), stringsAsFactors = FALSE)
  for (n in cfg_pam50$min_cells) qc[[paste0("included_ge", n)]] <- qc$malignant_cells >= n
  inclusion <- data.frame(dataset = d$dataset, cluster = map$cluster, celltype_final = map$celltype_final, malignant_status = map$malignant_status, included_as_malignant = map$celltype_final %in% labels, inclusion_basis = ifelse(map$celltype_final %in% labels, "explicit approved_malignant_labels", "not selected"), stringsAsFactors = FALSE)
  inclusion$cell_count <- as.integer(table(factor(cluster, levels = map$cluster)))
  list(pb = pb, qc = qc, inclusion = inclusion, rds_sha256 = sha256_file(d$rds), map_sha256 = sha256_file(d$annotation_map))
}

usable_expression <- function(pbs, centroids) {
  common_features <- Reduce(intersect, lapply(pbs, rownames))
  # Reorder every cohort to the shared feature vector before column-binding.
  common_pbs <- lapply(pbs, function(x) x[common_features, , drop = FALSE])
  logcpm <- lapply(common_pbs, normalise_logcpm); combined <- do.call(cbind, logcpm)
  expressed <- rownames(combined)[Matrix::rowSums(do.call(cbind, common_pbs)) > 0]
  row_sd <- apply(combined, 1L, stats::sd); scaled <- sweep(sweep(combined[row_sd > 0, , drop = FALSE], 1L, rowMeans(combined)[row_sd > 0], "-"), 1L, row_sd[row_sd > 0], "/")
  usable <- Reduce(intersect, list(rownames(scaled), rownames(centroids), common_features, expressed))
  coverage <- data.frame(gene = rownames(centroids), centroid_present = TRUE, common_feature = rownames(centroids) %in% common_features, expressed_in_pseudobulk = rownames(centroids) %in% expressed, nonzero_variance = rownames(centroids) %in% rownames(scaled), usable = rownames(centroids) %in% usable, stringsAsFactors = FALSE)
  list(scaled = scaled, usable = usable, coverage = coverage)
}

save_final_figure <- function(plot) {
  stem <- out("Reviewer_Figure_ThreeTumor_LumB_similarity_FINAL")
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = 10, height = 12, device = cairo_pdf)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = 10, height = 12, dpi = 300, device = "png")
  ggplot2::ggsave(paste0(stem, ".tiff"), plot, width = 10, height = 12, dpi = 300, device = "tiff", compression = "lzw")
  files <- paste0(stem, c(".pdf", ".png", ".tiff"))
  assert(all(file.info(files)$size > 0), "Figure export produced an empty file.")
}

score_programs <- function(scaled) {
  availability <- do.call(rbind, lapply(names(programs), function(nm) data.frame(program = nm, gene = programs[[nm]], common_feature = programs[[nm]] %in% rownames(scaled), zero_variance = !(programs[[nm]] %in% rownames(scaled)), used = programs[[nm]] %in% rownames(scaled), stringsAsFactors = FALSE)))
  score <- vapply(names(programs), function(nm) { genes <- intersect(programs[[nm]], rownames(scaled)); assert(length(genes) > 0L, paste0("No usable genes for fixed ", nm, " program.")); colMeans(scaled[genes, , drop = FALSE]) }, numeric(ncol(scaled)))
  colnames(score) <- names(programs); program_tab <- data.frame(patient_id = rownames(score), score, check.names = FALSE, stringsAsFactors = FALSE)
  assert(nrow(program_tab) == ncol(scaled) && identical(program_tab$patient_id, colnames(scaled)) && all(vapply(program_tab[names(programs)], function(x) all(is.finite(x)), logical(1))), "Program-score matrix validation failed.")
  program_tab$dataset <- sub("__.*$", "", program_tab$patient_id); program_tab$unit_id <- sub("^[^_]+__", "", program_tab$patient_id)
  list(scores = program_tab, coverage = availability)
}

correlate_centroids <- function(scaled, usable, centroids, method) {
  expression_matrix <- as.matrix(scaled[usable, , drop = FALSE])
  storage.mode(expression_matrix) <- "double"
  centroid_matrix <- as.matrix(centroids[usable, , drop = FALSE])
  storage.mode(centroid_matrix) <- "double"
  r <- stats::cor(expression_matrix, centroid_matrix, method = method)
  rownames(r) <- colnames(expression_matrix)
  assert(nrow(r) == length(colnames(expression_matrix)), "Centroid correlation/unit identifier dimension mismatch.")
  z <- data.frame(patient_id = rownames(r), r, check.names = FALSE, stringsAsFactors = FALSE); z$dataset <- sub("__.*$", "", z$patient_id); z$unit_id <- sub("^[^_]+__", "", z$patient_id); z$method <- method
  z
}

annotate_centroids <- function(x, usable_n) {
  values <- as.matrix(x[, centroid_cols, drop = FALSE]); ord <- t(apply(values, 1L, order, decreasing = TRUE)); x$best_centroid <- centroid_cols[ord[,1]]; x$best_correlation <- values[cbind(seq_len(nrow(x)), ord[,1])]; x$second_best_correlation <- values[cbind(seq_len(nrow(x)), ord[,2])]
  x$classification_margin <- x$best_correlation - x$second_best_correlation; x$LumB_correlation <- x$LumB; x$LumB_margin_vs_others <- x$LumB - apply(values[, setdiff(centroid_cols,"LumB"), drop = FALSE], 1L, max)
  x$uncertain <- !is.finite(x$classification_margin) | x$classification_margin <= cfg_pam50$uncertain_margin | usable_n < 40L
  x$nearest_centroid_reported <- ifelse(x$uncertain, "uncertain", x$best_centroid); x
}

build_correlation_profile_pca <- function(primary) {
  x <- as.matrix(primary[, centroid_cols, drop = FALSE])
  storage.mode(x) <- "double"
  assert(all(is.finite(x)), "PCA input contains non-finite centroid correlations.")
  assert(all(apply(x, 2L, stats::sd) > 0), "PCA input contains a constant centroid-correlation column.")
  fit <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  variance <- fit$sdev^2 / sum(fit$sdev^2)
  scores <- data.frame(
    patient_id = primary$patient_id,
    dataset = primary$dataset,
    unit_id = primary$unit_id,
    nearest_centroid_reported = primary$nearest_centroid_reported,
    PC1 = fit$x[, "PC1"],
    PC2 = fit$x[, "PC2"],
    stringsAsFactors = FALSE
  )
  loadings <- data.frame(centroid = rownames(fit$rotation), fit$rotation, check.names = FALSE, row.names = NULL)
  variance_table <- data.frame(
    component = paste0("PC", seq_along(variance)),
    variance_explained = variance,
    cumulative_variance = cumsum(variance),
    stringsAsFactors = FALSE
  )
  assert(nrow(scores) == nrow(primary) && all(is.finite(scores$PC1)) && all(is.finite(scores$PC2)), "PCA score validation failed.")
  list(scores = scores, loadings = loadings, variance = variance_table)
}

build_nearest_centroid_dot_matrix <- function(primary) {
  centroid_levels <- c("LumA", "LumB", "HER2", "Basal", "Normal", "uncertain")
  z <- primary[, c("patient_id", "dataset", "unit_id", "nearest_centroid_reported"), drop = FALSE]
  z$dataset <- factor(z$dataset, levels = c("ERTotal", "HER2", "TNBC"))
  z$nearest_centroid_reported <- factor(z$nearest_centroid_reported, levels = centroid_levels)
  z <- z[order(z$dataset, z$nearest_centroid_reported, z$unit_id), , drop = FALSE]
  z$slot <- ave(seq_len(nrow(z)), z$dataset, FUN = seq_along)
  z$dataset_plot <- factor(as.character(z$dataset), levels = c("TNBC", "HER2", "ERTotal"))
  summary_rows <- do.call(rbind, lapply(levels(z$dataset), function(ds) {
    d <- z[z$dataset == ds, , drop = FALSE]
    counts <- table(factor(d$nearest_centroid_reported, levels = centroid_levels))
    counts <- counts[counts > 0]
    data.frame(
      dataset = ds,
      dataset_plot = factor(ds, levels = c("TNBC", "HER2", "ERTotal")),
      total_n = nrow(d),
      composition_label = paste(paste0(names(counts), " ", as.integer(counts), "/", nrow(d)), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
  assert(nrow(z) == nrow(primary) && all(table(z$dataset) == table(primary$dataset)), "Nearest-centroid dot-matrix validation failed.")
  list(points = z, summary = summary_rows)
}

group_statistics <- function(primary, formal) {
  contrasts <- list(c("ERTotal","HER2"), c("ERTotal","TNBC")); do.call(rbind, lapply(contrasts, function(g) {
    a <- primary$LumB_correlation[primary$dataset == g[1]]; b <- primary$LumB_correlation[primary$dataset == g[2]]
    if (!formal) return(data.frame(contrast=paste(g,collapse=" vs "), n_group1=length(a), n_group2=length(b), analysis_status="descriptive_only", stringsAsFactors=FALSE))
    w <- stats::wilcox.test(a,b, exact=FALSE, conf.int=TRUE); u <- unname(w$statistic)
    data.frame(contrast=paste(g,collapse=" vs "), n_group1=length(a), n_group2=length(b), wilcoxon_W=u, hodges_lehmann=unname(w$estimate), ci_low=w$conf.int[1], ci_high=w$conf.int[2], rank_biserial=2*u/(length(a)*length(b))-1, median_difference=stats::median(a)-stats::median(b), p_value=w$p.value, analysis_status="exploratory_patient_level", stringsAsFactors=FALSE)
  }))
}

main <- function() {
  assert(identical(names(cfg_pam50$datasets), c("ERTotal","HER2","TNBC")), "Dataset order must be ERTotal, HER2, TNBC."); assert(cfg_pam50$primary_min_cells %in% cfg_pam50$min_cells && cfg_pam50$primary_correlation %in% cfg_pam50$correlation_methods, "Invalid primary analysis configuration.")
  check_logcpm(); centroids <- read_centroids(); pieces <- lapply(cfg_pam50$datasets, extract_pseudobulk); qc <- do.call(rbind,lapply(pieces,`[[`,"qc")); inclusion <- do.call(rbind,lapply(pieces,`[[`,"inclusion")); write_tsv(qc,"ThreeTumor_patient_pseudobulk_QC.tsv"); write_tsv(inclusion,"ThreeTumor_malignant_epithelial_cell_inclusion.tsv")
  all_methods <- list(); sensitivity <- list(); leave_one_out <- list(); primary <- NULL; program_out <- NULL
  primary_ue <- NULL; primary_coverage <- NULL; coverage_by_threshold <- list()
  for (threshold in cfg_pam50$min_cells) {
    eligible <- table(factor(qc$dataset[qc[[paste0("included_ge",threshold)]]], levels=names(cfg_pam50$datasets))); eligible_ok <- all(eligible >= 3L)
    if (!eligible_ok) { sensitivity[[as.character(threshold)]] <- data.frame(min_cells=threshold,method=cfg_pam50$primary_correlation,dataset=names(eligible),eligible_units=as.integer(eligible),median_LumB=NA_real_,IQR_LumB=NA_real_,usable_pam50_n=NA_integer_,ERTotal_minus_HER2=NA_real_,ERTotal_minus_TNBC=NA_real_,analysis_status="not_interpretable_fewer_than_3_units",stringsAsFactors=FALSE); next }
    pbs <- lapply(seq_along(pieces), function(i) { q <- pieces[[i]]$qc; ids <- q$patient_id[q[[paste0("included_ge",threshold)]]]; pieces[[i]]$pb[,ids,drop=FALSE] })
    ue <- usable_expression(pbs,centroids); usable_n <- length(ue$usable)
    coverage_threshold <- ue$coverage
    coverage_threshold$min_cells <- threshold
    coverage_by_threshold[[as.character(threshold)]] <- coverage_threshold
    if (usable_n < 30L) { sensitivity[[as.character(threshold)]] <- data.frame(min_cells=threshold,method=cfg_pam50$primary_correlation,dataset=names(eligible),eligible_units=as.integer(eligible),median_LumB=NA_real_,IQR_LumB=NA_real_,usable_pam50_n=usable_n,ERTotal_minus_HER2=NA_real_,ERTotal_minus_TNBC=NA_real_,analysis_status="not_interpretable_PAM50_coverage_under_30",stringsAsFactors=FALSE); next }
    by_method <- lapply(cfg_pam50$correlation_methods,function(m) correlate_centroids(ue$scaled,ue$usable,centroids,m)); names(by_method)<-cfg_pam50$correlation_methods; tab <- do.call(rbind,by_method); tab$min_cells <- threshold; all_methods[[as.character(threshold)]] <- tab
    summaries <- do.call(rbind,lapply(cfg_pam50$correlation_methods,function(m) {
      z <- tab[tab$method == m, , drop = FALSE]
      med <- vapply(names(cfg_pam50$datasets), function(ds) median(z$LumB[z$dataset == ds]), numeric(1))
      do.call(rbind,lapply(names(cfg_pam50$datasets),function(ds) {
        v <- z$LumB[z$dataset == ds]
        data.frame(
          min_cells=threshold, method=m, dataset=ds, eligible_units=length(v),
          median_LumB=median(v), IQR_LumB=IQR(v), usable_pam50_n=usable_n,
          ERTotal_minus_HER2=unname(med["ERTotal"]-med["HER2"]),
          ERTotal_minus_TNBC=unname(med["ERTotal"]-med["TNBC"]),
          analysis_status=ifelse(usable_n>=40,"exploratory","exploratory_low_coverage")
        )
      }))
    })); sensitivity[[as.character(threshold)]] <- summaries
    if (threshold == cfg_pam50$primary_min_cells) {
      primary <- annotate_centroids(tab[tab$method==cfg_pam50$primary_correlation,,drop=FALSE],usable_n)
      primary_ue <- ue
      primary_coverage <- coverage_threshold
      program_out <- score_programs(ue$scaled)
    }
    for (pid in tab$patient_id[tab$method==cfg_pam50$primary_correlation]) {
      z <- tab[tab$method==cfg_pam50$primary_correlation & tab$patient_id!=pid,]
      leave_one_out[[length(leave_one_out)+1L]] <- data.frame(
        min_cells=threshold,
        method=cfg_pam50$primary_correlation,
        removed_patient=pid,
        ERTotal_minus_HER2=median(z$LumB[z$dataset=="ERTotal"])-median(z$LumB[z$dataset=="HER2"]),
        ERTotal_minus_TNBC=median(z$LumB[z$dataset=="ERTotal"])-median(z$LumB[z$dataset=="TNBC"]),
        stringsAsFactors=FALSE
      )
    }
  }
  assert(!is.null(primary) && !is.null(primary_ue) && !is.null(primary_coverage), "BLOCKED: primary threshold did not yield >=3 units per cohort and >=30 usable PAM50 genes.")
  pca_out <- build_correlation_profile_pca(primary); dot_out <- build_nearest_centroid_dot_matrix(primary)
  write_tsv(primary_coverage,"ThreeTumor_PAM50_gene_coverage.tsv"); write_tsv(do.call(rbind,coverage_by_threshold),"ThreeTumor_PAM50_gene_coverage_all_thresholds.tsv"); write_tsv(program_out$coverage,"ThreeTumor_subtype_program_gene_coverage.tsv"); write_tsv(program_out$scores,"ThreeTumor_subtype_program_scores_by_patient.tsv"); write_tsv(primary,"ThreeTumor_PAM50_centroid_correlations_by_patient.tsv"); write_tsv(do.call(rbind,all_methods),"ThreeTumor_PAM50_centroid_correlations_all_methods.tsv"); sensitivity_tab <- do.call(rbind,sensitivity); write_tsv(sensitivity_tab,"ThreeTumor_PAM50_sensitivity_summary.tsv"); leave_one_out_tab <- do.call(rbind,leave_one_out); write_tsv(leave_one_out_tab,"ThreeTumor_PAM50_leave_one_out_sensitivity.tsv")
  primary_pam50_file <- out("ThreeTumor_PAM50_centroid_correlations_by_patient.tsv")
  anchor_manifest <- data.frame(
    manifest_version = "1",
    artifact = basename(primary_pam50_file),
    artifact_sha256 = sha256_file(primary_pam50_file),
    n_patients = nrow(primary),
    datasets = paste(names(cfg_pam50$datasets), collapse = ";"),
    primary_min_cells = cfg_pam50$primary_min_cells,
    primary_correlation = cfg_pam50$primary_correlation,
    usable_pam50_genes = length(primary_ue$usable),
    generating_script = basename(script_file),
    generating_script_sha256 = sha256_file(script_file),
    stringsAsFactors = FALSE
  )
  write_tsv(anchor_manifest, "ThreeTumor_PAM50_primary_anchor_manifest.tsv")
  write_tsv(pca_out$scores,"ThreeTumor_PAM50_correlation_profile_PCA_scores.tsv"); write_tsv(pca_out$loadings,"ThreeTumor_PAM50_correlation_profile_PCA_loadings.tsv"); write_tsv(pca_out$variance,"ThreeTumor_PAM50_correlation_profile_PCA_variance.tsv"); write_tsv(dot_out$points,"ThreeTumor_nearest_centroid_patient_dot_matrix.tsv"); write_tsv(dot_out$summary,"ThreeTumor_nearest_centroid_composition_summary.tsv")
  summary <- do.call(rbind,lapply(split(primary$LumB_correlation,primary$dataset),function(x)data.frame(n=length(x),median=median(x),IQR=IQR(x),min=min(x),max=max(x)))); summary$dataset<-rownames(summary);rownames(summary)<-NULL; write_tsv(summary,"ThreeTumor_LumB_group_summary.tsv")
  eligible_primary <- table(factor(qc$dataset[qc[[paste0("included_ge",cfg_pam50$primary_min_cells)]]],levels=names(cfg_pam50$datasets))); stats <- group_statistics(primary,all(eligible_primary>=5L)); if("p_value"%in%names(stats))stats$p_adjust_holm<-p.adjust(stats$p_value,"holm");write_tsv(stats,"ThreeTumor_LumB_patient_level_statistics.tsv")
  if (isTRUE(cfg_pam50$write_figures)) {
    dataset_palette <- c(ERTotal = "#D55E00", HER2 = "#009E73", TNBC = "#0072B2")
    centroid_palette <- c(LumA = "#009E73", LumB = "#00A6D6", HER2 = "#E69F00", Basal = "#D55E00", Normal = "#56B4E9", uncertain = "#CC79A7")
    centroid_shapes <- c(LumA = 16, LumB = 17, HER2 = 15, Basal = 18, Normal = 8, uncertain = 4)
    set.seed(cfg_pam50$seed)
    long <- do.call(rbind,lapply(names(programs),function(nm)data.frame(dataset=program_out$scores$dataset,program=nm,score=program_out$scores[[nm]])))
    pA <- ggplot2::ggplot(long,ggplot2::aes(dataset,score,colour=dataset))+ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::geom_jitter(width=.12)+ggplot2::facet_wrap(~program,scales="free_y")+ggplot2::scale_colour_manual(values=dataset_palette)+ggplot2::theme_classic(base_size=11)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30,hjust=1))+ggplot2::labs(title="A. Fixed program resemblance",x=NULL,y="Joint gene-z mean",colour="Cohort")
    heat <- stats::reshape(primary[,c("patient_id","dataset",centroid_cols)],varying=centroid_cols,v.names="correlation",timevar="centroid",times=centroid_cols,direction="long")
    pB <- ggplot2::ggplot(heat,ggplot2::aes(centroid,patient_id,fill=correlation))+ggplot2::geom_tile()+ggplot2::scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0)+ggplot2::facet_grid(dataset~.,scales="free_y",space="free_y")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(axis.title=ggplot2::element_blank(),panel.grid=ggplot2::element_blank())+ggplot2::labs(title="B. PAM50 centroid correlations",fill="Spearman r")
    pC <- ggplot2::ggplot(primary,ggplot2::aes(dataset,LumB_correlation,colour=dataset))+ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::geom_jitter(width=.12)+ggplot2::scale_colour_manual(values=dataset_palette)+ggplot2::theme_classic(base_size=11)+ggplot2::theme(legend.position="none")+ggplot2::labs(title="C. Luminal B transcriptional resemblance",x=NULL,y="Spearman correlation")
    pca_scores <- pca_out$scores
    pca_scores$dataset <- factor(pca_scores$dataset,levels=c("ERTotal","HER2","TNBC"))
    pca_scores$nearest_centroid_reported <- factor(pca_scores$nearest_centroid_reported,levels=names(centroid_shapes))
    hulls <- do.call(rbind,lapply(split(pca_scores,pca_scores$dataset),function(d)d[chull(d$PC1,d$PC2),,drop=FALSE]))
    pc1_pct <- 100*pca_out$variance$variance_explained[pca_out$variance$component=="PC1"]
    pc2_pct <- 100*pca_out$variance$variance_explained[pca_out$variance$component=="PC2"]
    pD <- ggplot2::ggplot(pca_scores,ggplot2::aes(PC1,PC2,colour=dataset,shape=nearest_centroid_reported))+ggplot2::geom_hline(yintercept=0,linewidth=.3,linetype=3,colour="grey75")+ggplot2::geom_vline(xintercept=0,linewidth=.3,linetype=3,colour="grey75")+ggplot2::geom_polygon(data=hulls,ggplot2::aes(PC1,PC2,fill=dataset,group=dataset),inherit.aes=FALSE,alpha=.10,linewidth=.45)+ggplot2::geom_point(size=3,stroke=1)+ggplot2::scale_colour_manual(values=dataset_palette)+ggplot2::scale_fill_manual(values=dataset_palette,guide="none")+ggplot2::scale_shape_manual(values=centroid_shapes)+ggplot2::coord_equal()+ggplot2::theme_classic(base_size=11)+ggplot2::labs(title="D. Patient PAM50 correlation-profile PCA",x=sprintf("PC1 (%.1f%%)",pc1_pct),y=sprintf("PC2 (%.1f%%)",pc2_pct),colour="Cohort",shape="Nearest centroid")
    dot_points <- dot_out$points; dot_summary <- dot_out$summary
    pE <- ggplot2::ggplot(dot_points,ggplot2::aes(slot,dataset_plot,fill=nearest_centroid_reported))+ggplot2::geom_point(shape=21,size=4.2,stroke=.45,colour="grey20")+ggplot2::geom_text(data=dot_summary,ggplot2::aes(x=14.2,y=dataset_plot,label=composition_label),inherit.aes=FALSE,hjust=0,size=3.4)+ggplot2::scale_fill_manual(values=centroid_palette,drop=FALSE)+ggplot2::scale_x_continuous(limits=c(.5,23),breaks=NULL,expand=c(0,0))+ggplot2::theme_classic(base_size=11)+ggplot2::theme(axis.line.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),legend.position="bottom",legend.box="horizontal")+ggplot2::guides(fill=ggplot2::guide_legend(nrow=1,byrow=TRUE))+ggplot2::labs(title="E. Nearest-centroid composition (one dot per patient)",x=NULL,y=NULL,fill="Nearest centroid")
    final <- (pA | pB) / (pC | pD) / pE + patchwork::plot_layout(heights=c(1.1,1,0.55))
    save_final_figure(final)
  }
  primary_summary <- sensitivity_tab[
    sensitivity_tab$min_cells == cfg_pam50$primary_min_cells &
      sensitivity_tab$method == cfg_pam50$primary_correlation,
    , drop = FALSE
  ]
  ert_minus_her2 <- unique(primary_summary$ERTotal_minus_HER2)
  ert_minus_tnbc <- unique(primary_summary$ERTotal_minus_TNBC)
  evidence <- if(all(eligible_primary>=5L)&&length(primary_ue$usable)>=40L) "Suggestive" else "Insufficient"
  report<-c(
    "# Three-tumor PAM50/Luminal B resemblance run report","",
    paste0("Evidence level: **",evidence,"**. This evaluates transcriptional resemblance, not clinical PAM50 diagnosis."),"",
    paste0("Helper: ",helper_path," SHA256=",sha256_file(helper_path)),
    paste0("Centroid source: ",cfg_pam50$centroid_source),
    paste0("Centroid preprocessing: ",cfg_pam50$centroid_preprocessing),
    paste0("Usable PAM50 genes at primary threshold: ",length(primary_ue$usable),"/",nrow(centroids)),
    paste0("Primary median LumB difference, ERTotal minus HER2: ",format(ert_minus_her2,digits=4)),
    paste0("Primary median LumB difference, ERTotal minus TNBC: ",format(ert_minus_tnbc,digits=4)),"",
    "## Required interpretation limits",
    "- ERTotal is the prespecified ER-positive, luminal-enriched single-cell context for the primary downstream analysis; these results do not show that ERTotal has the highest continuous LumB correlation among all three cohorts.",
    "- Cross-cohort processing/batch effects may confound similarity.",
    "- Any sample-level configuration cannot support patient-level claims.",
    "- No result clinically labels ERTotal or any tumour as PAM50 Luminal B."
  );writeLines(report,out("ThreeTumor_PAM50_LumB_similarity_RUN_REPORT.md"));writeLines(capture.output(sessionInfo()),out("ThreeTumor_PAM50_LumB_similarity_sessionInfo.txt"))
}

dir.create(cfg_pam50$output_dir, recursive=TRUE, showWarnings=FALSE)
tryCatch(main(),error=function(e){writeLines(c("# Three-tumor PAM50/Luminal B resemblance status","","## BLOCKED",conditionMessage(e),"","No fallback to cells, unapproved epithelial labels, downloaded centroids, or altered source RDS was performed."),out("ThreeTumor_PAM50_LumB_similarity_RUN_REPORT.md"));stop(e)})

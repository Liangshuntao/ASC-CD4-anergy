# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 457:508
clin_patient <- read.delim(file.path(metabric_dir, "data_clinical_patient.txt"), skip = 4, header = TRUE, stringsAsFactors = FALSE)
clin_sample <- read.delim(file.path(metabric_dir, "data_clinical_sample.txt"), skip = 4, header = TRUE, stringsAsFactors = FALSE)

metabric_clin <- merge(clin_patient, clin_sample, by = "PATIENT_ID")
rownames(metabric_clin) <- metabric_clin$SAMPLE_ID
metabric_clin$PAM50_Subtype <- metabric_clin$CLAUDIN_SUBTYPE 

# ------------------------------------------------------------------------------
# Read the local expression matrix while preserving sample names
# ------------------------------------------------------------------------------
message("Reading the microarray expression matrix with sample names preserved...")
# Preserve original names with check.names = FALSE
expr_raw <- read.delim(file.path(metabric_dir, "data_mrna_illumina_microarray.txt"), 
                       header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Remove duplicate and empty identifiers, then set row names
expr_raw <- expr_raw[!duplicated(expr_raw$Hugo_Symbol) & expr_raw$Hugo_Symbol != "", ]
rownames(expr_raw) <- expr_raw$Hugo_Symbol
metabric_expr <- as.matrix(expr_raw[, -(1:2)])

# ------------------------------------------------------------------------------
# Select Luminal B samples and prepare survival variables
# ------------------------------------------------------------------------------
message("Selecting Luminal B samples and aligning the expression matrix...")
metabric_lumb_clin <- metabric_clin %>%
  mutate(
    OS.time = as.numeric(OS_MONTHS),
    OS.status = ifelse(grepl("DECEASED", OS_STATUS, ignore.case = TRUE), 1, 0)
  ) %>%
  filter(grepl("LumB", PAM50_Subtype, ignore.case = TRUE) | grepl("Luminal B", PAM50_Subtype, ignore.case = TRUE)) %>%
  filter(!is.na(OS.time) & OS.time > 0)

# Check the number of aligned samples
valid_samples <- intersect(rownames(metabric_lumb_clin), colnames(metabric_expr))
message(paste("Aligned samples:", length(valid_samples)))

metabric_lumb_clin <- metabric_lumb_clin[valid_samples, ]
expr_lumb_meta <- as.matrix(metabric_expr[, valid_samples])
# ------------------------------------------------------------------------------
# Calculate ssGSEA scores for the microarray data
# ------------------------------------------------------------------------------
message("Calculating ssGSEA scores for the microarray data...")
meta_genes <- rownames(expr_lumb_meta)
sig_asc_p_meta <- intersect(c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM"), meta_genes)
sig_cd4t_meta <- intersect(c("CD4", "CD3D", "CD3E"), meta_genes)
metabric_sigs <- list(ASC_P = sig_asc_p_meta, CD4T = sig_cd4t_meta)

if("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  meta_gsva <- gsva(ssgseaParam(expr_lumb_meta, metabric_sigs))
} else {
  meta_gsva <- gsva(expr_lumb_meta, metabric_sigs, method="ssgsea")
}

# SOURCE END 457:508

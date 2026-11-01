# Original cell-cycle-adjusted model and CAMERA sensitivity
# Executed by run_bulk_analysis.R; not a standalone entry point.

# SOURCE BEGIN 589:708
# Step 4. Robustness: adjust DID for a cell-cycle score (no new experiments)

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(msigdbr)
  library(tibble)
})

## -------------------------
# Required objects from the primary model
# expr_mat: genes by samples, containing log2(normalized expression + 1)
# sample_info: sample, subtype (N/P) and timepoint (Before/After)
# The interaction design is reconstructed below
## -------------------------

# Retrieve cell-cycle Hallmark sets mapped to mouse genes
m_h <- msigdbr(species = "Mus musculus", category = "H")

gs_g2m <- m_h %>% filter(gs_name == "HALLMARK_G2M_CHECKPOINT") %>% pull(gene_symbol) %>% unique()
gs_e2f <- m_h %>% filter(gs_name == "HALLMARK_E2F_TARGETS") %>% pull(gene_symbol) %>% unique()

cc_genes <- unique(c(gs_g2m, gs_e2f))
cc_genes <- intersect(cc_genes, rownames(expr_mat))

cat("Cell-cycle genes in expr_mat:", length(cc_genes), "\n")
stopifnot(length(cc_genes) >= 10)  # Require at least 10 matched cell-cycle genes

# Compute the mean of gene-wise standardized expression values
zmat <- t(scale(t(expr_mat[cc_genes, , drop = FALSE])))   # row-wise z-score
cc_score <- colMeans(zmat, na.rm = TRUE)

df_cc <- data.frame(sample = names(cc_score), cc_score = as.numeric(cc_score)) %>%
  left_join(sample_info, by = "sample")

print(df_cc %>% arrange(desc(cc_score)))

# Add cc_score as a covariate while retaining the interaction term
sample_info2 <- sample_info %>%
  left_join(df_cc %>% select(sample, cc_score), by = "sample") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("Before","After")),
    subtype   = factor(subtype,   levels = c("N","P"))
  )

design_cc <- model.matrix(~ 0 + subtype * timepoint + cc_score, data = sample_info2)
colnames(design_cc) <- make.names(colnames(design_cc))
print(colnames(design_cc))

# Use the same contrasts with the adjusted design
cont_cc <- makeContrasts(
  Pre_P_vs_N  = subtypeP - subtypeN,
  Post_P_vs_N = (subtypeP + subtypeP.timepointAfter) - subtypeN,
  DID_P_specific_change = subtypeP.timepointAfter,
  P_After_vs_Before = timepointAfter + subtypeP.timepointAfter,
  N_After_vs_Before = timepointAfter,
  levels = design_cc
)

fit_cc  <- lmFit(expr_mat, design_cc)
fit_cc2 <- eBayes(contrasts.fit(fit_cc, cont_cc))

get_res <- function(fit_obj, coef_name){
  topTable(fit_obj, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("Gene") %>%
    mutate(FDR = p.adjust(P.Value, "BH"))
}

res_did_cc  <- get_res(fit_cc2, "DID_P_specific_change")
res_pre_cc  <- get_res(fit_cc2, "Pre_P_vs_N")
res_post_cc <- get_res(fit_cc2, "Post_P_vs_N")

write.csv(res_did_cc,  "ASC_RNA_DID_P_specific_change_ccAdjusted.csv", row.names = FALSE)

# Compare primary and cell-cycle-adjusted DID estimates
# Read the primary DID table if it is not already loaded
if(!exists("res_did")){
  res_did <- read.csv("ASC_RNA_DID_P_specific_change_all.csv", check.names = FALSE)
}

cmp <- res_did %>%
  select(Gene, logFC_raw = logFC, FDR_raw = FDR) %>%
  inner_join(res_did_cc %>% select(Gene, logFC_cc = logFC, FDR_cc = FDR),
             by = "Gene")

cat("Correlation of DID logFC (raw vs cc-adjusted): ",
    cor(cmp$logFC_raw, cmp$logFC_cc, use = "pairwise.complete.obs"), "\n")

# Selected genes for comparison
genes_mech <- c(
  "Prelp","Thbs1","Thbs2","Timp3","Has2","Cemip","Npnt","Itga11",
  "Alcam","Lrrc32","Jag1","Bmp2","Ltbp2","Igfbp3","Angptl2","Slc7a5","Id1","Klf4","Serpina3n","Cd34","Ptprb"
)

cmp_mech <- cmp %>%
  filter(Gene %in% genes_mech) %>%
  arrange(FDR_cc)

print(cmp_mech)

write.csv(cmp_mech, "ASC_DID_keyGenes_raw_vs_ccAdjusted.csv", row.names = FALSE)

# Run Hallmark CAMERA with the cell-cycle-adjusted design
hallmark_list <- split(m_h$gene_symbol, m_h$gs_name)
idx_h <- ids2indices(hallmark_list, rownames(expr_mat))

did_contrast_cc <- cont_cc[, "DID_P_specific_change"]

cam_asc_did_cc <- camera(
  y = expr_mat,
  index = idx_h,
  design = design_cc,
  contrast = did_contrast_cc
) %>%
  rownames_to_column("Pathway") %>%
  mutate(FDR = p.adjust(PValue, "BH")) %>%
  arrange(FDR)

write.csv(cam_asc_did_cc, "ASC_DID_CAMERA_Hallmark_ccAdjusted.csv", row.names = FALSE)
print(head(cam_asc_did_cc[, c("Pathway","Direction","PValue","FDR")], 20))

# SOURCE END 589:708

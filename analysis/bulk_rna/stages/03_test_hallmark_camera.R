# Original primary Hallmark CAMERA
# Executed by run_bulk_analysis.R; not a standalone entry point.

# SOURCE BEGIN 381:408
# Hallmark CAMERA analysis of the DID contrast

library(msigdbr)
library(limma)
library(dplyr)
library(tibble)

# CAMERA requires the expression matrix, design matrix and contrast
# Reuse design and cont from the primary model
# Extract the DID contrast vector
did_contrast <- cont[, "DID_P_specific_change"]

m_h <- msigdbr(species="Mus musculus", category="H")
hallmark_list <- split(m_h$gene_symbol, m_h$gs_name)
idx_h <- ids2indices(hallmark_list, rownames(expr_mat))

cam_asc_did <- camera(
  y = expr_mat,
  index = idx_h,
  design = design,
  contrast = did_contrast
) %>%
  rownames_to_column("Pathway") %>%
  mutate(FDR = p.adjust(PValue, "BH")) %>%
  arrange(FDR)

head(cam_asc_did, 20)
write.csv(cam_asc_did, "ASC_DID_CAMERA_Hallmark.csv", row.names = FALSE)

# SOURCE END 381:408

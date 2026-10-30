# Original primary Hallmark CAMERA
# Executed by run_bulk_analysis.R; not a standalone entry point.

# SOURCE BEGIN 381:408
#### 第二步：对 res_did 做 Hallmark 富集（更系统地把模块补齐）

library(msigdbr)
library(limma)
library(dplyr)
library(tibble)

# gene-level fit 已经在 fit2 里了，但 camera 需要 design 与 contrast
# 复用你建模时的 design / cont（确保仍在环境中）
# DID contrast 向量：
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

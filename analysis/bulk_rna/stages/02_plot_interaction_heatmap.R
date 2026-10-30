# Figure 8C source block; final visual edits not reproduced
# Executed by run_bulk_analysis.R; not a standalone entry point.

# SOURCE BEGIN 320:379
###第一步：把 DID 的 top genes 做热图，并用 annotation_row 标模块（投稿风格）
library(dplyr)
library(pheatmap)

# 1) 选 top genes（比如前 60 个，或 FDR<0.01 且 |logFC|>0.5）
top_genes <- res_did %>%
  filter(FDR < 0.01) %>%
  arrange(FDR) %>%
  slice(1:min(60, n())) %>%
  pull(Gene)

mat_top <- expr_mat[top_genes, , drop=FALSE]

write.csv(mat_top,  "ASC_DID_genes.csv",  row.names = T)

# 2) 列注释（和之前一致）
annot_col <- data.frame(
  subtype   = sample_info$subtype,
  timepoint = sample_info$timepoint,
  row.names = sample_info$sample
)

# 3) 行注释（模块标注：先对你这批强信号做）
module_map <- tibble(
  Gene = top_genes,
  module = case_when(
    Gene %in% c("Prelp","Npnt","Itga11","Timp3","Thbs1","Thbs2","Cemip") ~ "ECM_remodeling",
    Gene %in% c("Has2","Cemip") ~ "HA_GAG",
    Gene %in% c("Alcam") ~ "Adhesion_contact",
    Gene %in% c("Lrrc32") ~ "GARP_TGFb_axis",
    Gene %in% c("Gpx3") ~ "Redox",
    Gene %in% c("Igfbp3","Angptl2") ~ "Secreted_growth_inflammation",
    TRUE ~ "Other"
  )
)

annot_row <- data.frame(
  module = module_map$module,
  row.names = module_map$Gene
)

# 4) 输出 PDF（投稿风格）
pdf("ASC_DID_topgenes_heatmap.pdf", width = 11, height = 13)
pheatmap(
  mat_top,
  scale = "row",
  annotation_col = annot_col,
  annotation_row = annot_row,
  show_colnames = T,
  show_rownames = T,
  cellwidth = 14,
  cellheight = 10,
  fontsize_row = 7,
  fontsize_col = 7,
  border_color = "white",
  treeheight_row = 0,
  treeheight_col = 0,
  main = "ASC DID (P-specific coculture response): top genes"
)
dev.off()

# SOURCE END 320:379

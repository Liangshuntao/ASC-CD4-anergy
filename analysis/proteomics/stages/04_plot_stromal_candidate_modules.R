pdf("ASC_candidate_modules.pdf", width = 8, height = 12)
# SOURCE BEGIN 681:707
asc_suppressor_genes <- c(
  ## 1) TGF-β / BMP / growth factors
  "Tgfb1","Tgfb2","Tgfb3","Thbs1","Ltbp1","Ltbp3","Eng","Bmp4","Bmp2","Bmp7",
  "Bmpr1a","Bmpr2","Smad2","Smad3","Smad7",
  ## 2) Checkpoints / death ligands
  "Cd274","Pdcd1lg2","Fasl","Tnfsf10","Lgals9","Lgals1","Cd47",
  ## 3) Metabolic suppressors
  "Ido1","Arg1","Nos2","Entpd1","Nt5e","Ptgs2","Ptges","Aldh1a1","Aldh1a2",
  ## 4) ROS / redox
  "Nox4","Cyba","Sod2","Gpx1","Prdx1","Prdx2","Cat",
  ## 5) Notch / Wnt ligands
  "Jag1","Jag2","Dll1","Dll4","Wnt5a","Wnt2","Wnt3a","Dkk1","Rspo3",
  ## 6) Complement / coagulation
  "C3","C4b","Cfb","Cfd","Vwf","Serpine1","Serpine2","Thbd","Fga","Fgb","Fgg",
  ## 7) ECM / adhesion
  "Fn1","Col1a1","Col1a2","Col3a1","Col4a1","Col4a2",
  "Lama4","Lamb1","Lamb2","Vcam1","Icam1","Sele","Selp"
)

## 和 ASC 表的基因列取交集
asc_genes_in_data <- intersect(asc_suppressor_genes, rownames(expr_asc_gene_median))

rownames(expr_asc_gene_median)
expr_asc_gene_median

asc_mat <- expr_asc_gene_median[asc_genes_in_data, ]  # 或 expr_imp_asc


# SOURCE END 681:707

# SOURCE BEGIN 729:886
## 1. 定义每个基因所属模块（按你现在这张图里出现的基因来列）
gene_module_map <- c(
  ## TGF-β / BMP 轴
  "Tgfb1"  = "TGFb/BMP",
  "Tgfb2"  = "TGFb/BMP",
  "Tgfb3"  = "TGFb/BMP",
  "Thbs1"  = "TGFb/BMP",
  "Ltbp1"  = "TGFb/BMP",
  "Ltbp3"  = "TGFb/BMP",
  "Eng"    = "TGFb/BMP",
  "Bmp2"   = "TGFb/BMP",
  "Bmp4"   = "TGFb/BMP",
  "Bmp7"   = "TGFb/BMP",
  "Smad2"  = "TGFb/BMP",
  "Smad3"  = "TGFb/BMP",
  "Smad7"  = "TGFb/BMP",
  
  ## Checkpoint / death ligand
  "Cd274"    = "Checkpoint",
  "Pdcd1lg2" = "Checkpoint",
  "Fasl"     = "Checkpoint",
  "Tnfsf10"  = "Checkpoint",
  "Lgals9"   = "Checkpoint",
  "Lgals1"   = "Checkpoint",
  "Cd47"     = "Checkpoint",
  
  ## Metabolic immune suppression (IDO/ARG/腺苷/RA/前列腺素)
  "Ido1"    = "Metabolic",
  "Arg1"    = "Metabolic",
  "Nos2"    = "Metabolic",
  "Entpd1"  = "Metabolic",
  "Nt5e"    = "Metabolic",
  "Ptgs2"   = "Metabolic",
  "Ptges"   = "Metabolic",
  "Aldh1a1" = "Metabolic",
  "Aldh1a2" = "Metabolic",
  
  ## ROS / redox
  "Nox4"  = "ROS/redox",
  "Cyba"  = "ROS/redox",
  "Sod2"  = "ROS/redox",
  "Gpx1"  = "ROS/redox",
  "Prdx1" = "ROS/redox",
  "Prdx2" = "ROS/redox",
  "Cat"   = "ROS/redox",
  
  ## Notch / Wnt ligands
  "Jag1"  = "Notch/Wnt",
  "Jag2"  = "Notch/Wnt",
  "Dll1"  = "Notch/Wnt",
  "Dll4"  = "Notch/Wnt",
  "Wnt5a" = "Notch/Wnt",
  "Wnt2"  = "Notch/Wnt",
  "Wnt3a" = "Notch/Wnt",
  "Dkk1"  = "Notch/Wnt",
  "Rspo3" = "Notch/Wnt",
  
  ## Complement / coagulation
  "C3"      = "Complement/coag",
  "C4b"     = "Complement/coag",
  "Cfb"     = "Complement/coag",
  "Cfd"     = "Complement/coag",
  "Vwf"     = "Complement/coag",
  "Serpine1"= "Complement/coag",
  "Serpine2"= "Complement/coag",
  "Thbd"    = "Complement/coag",
  "Fga"     = "Complement/coag",
  "Fgb"     = "Complement/coag",
  "Fgg"     = "Complement/coag",
  
  ## ECM / adhesion
  "Fn1"    = "ECM/adhesion",
  "Col1a1" = "ECM/adhesion",
  "Col1a2" = "ECM/adhesion",
  "Col3a1" = "ECM/adhesion",
  "Col4a1" = "ECM/adhesion",
  "Col4a2" = "ECM/adhesion",
  "Lama4"  = "ECM/adhesion",
  "Lamb1"  = "ECM/adhesion",
  "Lamb2"  = "ECM/adhesion",
  "Vcam1"  = "ECM/adhesion",
  "Icam1"  = "ECM/adhesion",
  "Sele"   = "ECM/adhesion",
  "Selp"   = "ECM/adhesion"
)

## 2. 和表达矩阵取交集，只保留在 asc_mat 里其实有的基因
common_genes <- intersect(names(gene_module_map), rownames(asc_mat))

asc_mat_sub <- asc_mat[common_genes, , drop = FALSE]

## 3. 构建行注释：module
row_annot <- data.frame(
  module = gene_module_map[common_genes]
)
rownames(row_annot) <- common_genes

## 4. 按模块顺序 + 模块内平均表达差排序（ASC_P vs ASC_N）
# 简单做法：先按 module，再按 ASC_P 平均值排序
# 假设 group_asc 是 factor(c("ASC_N","ASC_P",...))，长度 = ncol(asc_mat_sub)
asc_mean <- rowMeans(asc_mat_sub[, group_asc == "ASC_P", drop = FALSE], na.rm = TRUE)

row_annot$module <- factor(
  row_annot$module,
  levels = c("TGFb/BMP","Checkpoint","Metabolic","ROS/redox",
             "Notch/Wnt","Complement/coag","ECM/adhesion")
)

ord <- order(row_annot$module, asc_mean, decreasing = TRUE)
asc_mat_sub <- asc_mat_sub[ord, , drop = FALSE]
row_annot   <- row_annot[ord, , drop = FALSE]

## 5. 列注释：condition（和之前 CD4/ASC 热图保持一致配色）
ann_col <- data.frame(
  condition = group_asc
)
rownames(ann_col) <- colnames(asc_mat_sub)

## 6. 注释颜色：module + condition
ann_colors <- list(
  module = c(
    "TGFb/BMP"        = "#E41A1C",
    "Checkpoint"      = "#377EB8",
    "Metabolic"       = "#4DAF4A",
    "ROS/redox"       = "#984EA3",
    "Notch/Wnt"       = "#FF7F00",
    "Complement/coag" = "#A65628",
    "ECM/adhesion"    = "#999999"
  ),
  condition = c(
    "ASC_N" = "#00BFC4",
    "ASC_P" = "#F8766D"
  )
)

## 7. 统一“投稿风格”的 pheatmap
pheatmap(
  asc_mat_sub,
  scale = "row",
#  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annot,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,      # 行很多，字体略小
  fontsize_col = 8,
  border_color = "white",
  cellwidth = 20,
  cellheight = 14,
  treeheight_row = 0,
  treeheight_col = 0,
  main = "ASC-derived suppressive/abortive pathways (ASC_P vs ASC_N)"
)



# SOURCE END 729:886

dev.off()

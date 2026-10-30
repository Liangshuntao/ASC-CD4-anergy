#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

getwd()
####设置工作路径
setwd("D:/R/20251126 ASC protein/")
#######1.加载所需要的包########

## 安装（第一次用时运行）
# install.packages(c("tidyverse", "data.table"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("limma", "imputeLCMD", "clusterProfiler", "org.Mm.eg.db"))
BiocManager::install("imputeLCMD")

library(tidyverse)
library(data.table)
library(limma)
library(imputeLCMD)
library(clusterProfiler)
library(org.Mm.eg.db)  # 小鼠；如果将来做人，就改 org.Hs.eg.db
library(dplyr)
library(pheatmap)

## ========== 1. 读入原始 CSV + 构建基础对象 ==========

dat <- fread("ASC CD4 protein.csv", header = TRUE, sep = ",", data.table = FALSE)

# 基础注释（后面要用）
ann_all <- dat[, c("Accession", "Description")]
rownames(ann_all) <- ann_all$Accession

# 1.2 提取表达矩阵（去掉注释列）
feature_cols <- c("Accession", "Description")
expr_mat <- dat[, !(names(dat) %in% feature_cols)]
expr_mat <- as.data.frame(expr_mat)
expr_mat  <- dat$Accession

# 确保为数值
expr_mat <- mutate_all(expr_mat, as.numeric)

# 样本名
sample_names <- colnames(expr_mat)
print(sample_names)

## 按你之前的说明，这里手动指定 condition
## 3,4,5: ASC_CD105⁻CD90⁺（记 ASC_N）
## 6,7,8: ASC_CD105⁺CD90⁺（记 ASC_P）
## 9,10,11: CD4_N (与 ASC_N 共培养)
## 12,13,14: CD4_P (与 ASC_P 共培养)

condition <- c(
  rep("ASC_N", 3),
  rep("ASC_P", 3),
  rep("CD4_N", 3),
  rep("CD4_P", 3)
)

stopifnot(length(condition) == length(sample_names))

sample_info <- data.frame(
  sample    = sample_names,
  condition = factor(condition,
                     levels = c("ASC_N", "ASC_P", "CD4_N", "CD4_P")),
  stringsAsFactors = FALSE
)

print(sample_info)


# 1.3 log2 转换
expr_log2 <- log2(expr_mat + 1)
rownames(expr_log2) <- dat$Accession


## ========== 2. 拆分 ASC / CD4 两个子矩阵 ==========

asc_idx  <- which(sample_info$condition %in% c("ASC_N", "ASC_P"))
cd4_idx  <- which(sample_info$condition %in% c("CD4_N", "CD4_P"))

expr_log2_asc   <- expr_log2[, asc_idx, drop = FALSE]
expr_log2_cd4   <- expr_log2[, cd4_idx, drop = FALSE]

sample_info_asc <- sample_info[asc_idx, , drop = FALSE]
sample_info_cd4 <- sample_info[cd4_idx, , drop = FALSE]


## ========== 3. 封装：按组过滤 + MinProb 插补 的函数 ==========

filter_and_impute <- function(expr_log2_sub, group_vector,
                              min_nonNA_per_group = 2,
                              q_minprob = 0.01) {
  groups <- factor(group_vector)
  
  # 统计每个蛋白在每组的非 NA 数
  count_non_na_by_group <- function(mat, groups) {
    res <- lapply(levels(groups), function(g){
      idx <- which(groups == g)
      rowSums(!is.na(mat[, idx, drop = FALSE]))
    })
    res <- do.call(cbind, res)
    colnames(res) <- levels(groups)
    res
  }
  
  non_na_counts <- count_non_na_by_group(expr_log2_sub, groups)
  
  # 保留：至少在任一组中，非 NA 数 >= min_nonNA_per_group
  keep <- apply(non_na_counts, 1, function(x) any(x >= min_nonNA_per_group))
  expr_filt <- expr_log2_sub[keep, , drop = FALSE]
  
  # MinProb 插补（左截断）
  mat_for_imp <- as.matrix(expr_filt)
  mat_imp <- impute.MinProb(mat_for_imp, q = q_minprob)
  expr_imp  <- as.data.frame(mat_imp)
  
  return(list(
    expr_imp       = expr_imp,       # 插补后的表达矩阵
    keep_idx_logical = keep,         # 相对于原始行的逻辑向量
    non_na_counts  = non_na_counts   # 每组非 NA 统计
  ))
}

## ========== 4. ASC 分支：ASC_P vs ASC_N ==========

cat("==== ASC 分支开始 ====\n")

# 4.1 过滤 + 插补
asc_proc <- filter_and_impute(
  expr_log2_sub       = expr_log2_asc,
  group_vector        = sample_info_asc$condition,
  min_nonNA_per_group = 2,
  q_minprob           = 0.01
)

expr_imp_asc <- asc_proc$expr_imp
keep_asc     <- asc_proc$keep_idx_logical

# 匹配注释
ann_asc <- ann_all[keep_asc, ]
stopifnot(nrow(expr_imp_asc) == nrow(ann_asc))
rownames(expr_imp_asc) <- ann_asc$Accession

# 4.2 PCA 看一下（可选）
pca_asc <- prcomp(t(expr_imp_asc), scale. = FALSE)
pca_asc_df <- data.frame(
  PC1 = pca_asc$x[, 1],
  PC2 = pca_asc$x[, 2],
  sample_info_asc
)

p.asc.pca=ggplot(pca_asc_df,
       aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  theme_bw() +
  ggtitle("PCA - ASC only")

p.asc.pca

ggsave(filename = "ASC protein PCA.plot.pdf",
       plot =p.asc.pca,
       width = 8,
       height = 6,
       dpi = 300)


# 4.3 limma 差异分析：ASC_P vs ASC_N
group_asc <- droplevels(sample_info_asc$condition)  # ASC_N, ASC_P

design_asc <- model.matrix(~ 0 + group_asc)
colnames(design_asc) <- levels(group_asc)
print(design_asc)

contrast_asc <- makeContrasts(
  ASC_P_vs_ASC_N = ASC_P - ASC_N,
  levels = design_asc
)

fit_asc <- lmFit(expr_imp_asc, design_asc)
fit_asc <- contrasts.fit(fit_asc, contrast_asc)
fit_asc <- eBayes(fit_asc)

res_ASC <- topTable(fit_asc, coef = "ASC_P_vs_ASC_N",
                    number = Inf, sort.by = "P")

res_ASC_anno <- res_ASC %>%
  rownames_to_column("Accession") %>%
  left_join(ann_asc, by = "Accession") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.585, "sig", "ns")
  )

write.csv(res_ASC_anno,
          "DE_ASC_P_vs_ASC_N_asc_only_all.csv",
          row.names = FALSE)

filter(res_ASC_anno, sig == "sig")

write.csv(filter(res_ASC_anno, sig == "sig"),
          "DE_ASC_P_vs_ASC_N_asc_only_sig.csv",
          row.names = FALSE)

# 4.4 火山图 + 热图（ASC）
colnames(res_ASC_anno)

p.asc.volcano=ggplot(res_ASC_anno,
       aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.7, size = 0.5) +
  scale_color_manual(values = c("ns" = "grey70", "sig" = "red")) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  ggtitle("Volcano - ASC_P vs ASC_N")

p.asc.volcano

ggsave(filename = "ASC protein Volcano.plot.pdf",,
       plot = p.asc.volcano,
       width = 7,
       height =6 )



# 取 top 50 显著差异蛋白画热图
top50_asc <- res_ASC_anno %>%
  filter(sig == "sig") %>%
  arrange(FDR) %>%
  slice(1:min(50, n())) %>%
  pull(Accession)

mat_top50_asc <- expr_imp_asc[top50_asc, ]

library(pheatmap)

pdf("ASC_P vs ASC_N top50 DEPs.pdf", width = 4, height = 5)
pheatmap(mat_top50_asc,
         scale = "row",
         annotation_col = data.frame(condition = group_asc,
                                     row.names = colnames(mat_top50_asc)),
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = FALSE,
         main = "ASC_P vs ASC_N - top50 DEPs")


dev.off()

pdf("CD105+CD90+ ASC vs CD105-CD90+ ASC DEPs.pdf", width = 4, height = 5)
p.all.asc=pheatmap(expr_imp_asc,
         scale = "row",
         annotation_col = data.frame(condition = group_asc,
                                     row.names = colnames(mat_top50_asc)),
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = FALSE,
         main = "CD105+CD90+ ASC vs CD105-CD90+ ASC",
         border_color="white")
dev.off()

##4.5 基因名整理（从 Description 里 GN= 提取）
library(dplyr)
library(tibble)

## 4.5 从 Description 抽基因名 GN=
extract_gene_from_desc <- function(desc) {
  ifelse(
    grepl(" GN=", desc),
    sub(".* GN=([^ ]+).*", "\\1", desc),
    NA
  )
}
res_ASC_anno

res_ASC_anno_gene <- res_ASC_anno %>%
  mutate(Gene = extract_gene_from_desc(Description)) %>%
  filter(!is.na(Gene) & Gene != "")






## 如果一个基因对应多条蛋白记录，这里简单保留 P 值最小的那条
res_ASC_gene <- res_ASC_anno_gene %>%
  arrange(P.Value) %>%
  group_by(Gene) %>%
  slice(1) %>%
  ungroup()

write.csv(res_ASC_gene,
          "DE_ASC_P_vs_ASC_N_asc_only_geneLevel_all.csv",
          row.names = FALSE)


## 4.6 ASC GO/KEGG富集分析

library(clusterProfiler)
library(org.Mm.eg.db)

## 4.6.1 SYMOBL -> ENTREZ 映射
id_map_asc <- bitr(
  res_ASC_gene$Gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

## 背景基因 universe：所有有 ENTREZID 的 ASC 差异分析基因
gene_universe_asc <- id_map_asc$ENTREZID

## 4.6.2 选显著差异基因（ASC_P vs ASC_N）
deg_sig_asc <- res_ASC_gene %>%
  filter(FDR < 0.05 & abs(logFC) >= 0.585) #0.585 变化倍数1.5

deg_sig_asc_id <- deg_sig_asc %>%
  inner_join(id_map_asc, by = c("Gene" = "SYMBOL")) %>%
  pull(ENTREZID) %>%
  unique()

length(deg_sig_asc_id)  # 看看显著 DEG 有多少用于富集

## 4.6.3 GO BP 富集
ego_all_asc <- enrichGO(
  gene          = deg_sig_asc_id,
  universe      = gene_universe_asc,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

ego_all_asc_df <- as.data.frame(ego_all_asc)
write.csv(ego_all_asc_df,
          "ASC_GO_All_enrich_ASC_P_vs_ASC_N.csv",
          row.names = FALSE)

## 4.6.4 KEGG 富集（mmu）
ekegg_asc <- enrichKEGG(
  gene          = deg_sig_asc_id,
  universe      = gene_universe_asc,
  organism      = "mmu",
  keyType       = "kegg",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

ekegg_asc <- setReadable(
  ekegg_asc,             # 需要转换的富集结果对象
  OrgDb = org.Mm.eg.db,    # 指定物种注释数据库（这里是小鼠）
  keyType = "ENTREZID"     # 指定输入 ID 的类型（默认为 ENTREZID，通常不需要修改）
)


ekegg_asc_df <- as.data.frame(ekegg_asc@result)
write.csv(ekegg_asc_df,
          "ASC_KEGG_enrich_ASC_P_vs_ASC_N.csv",
          row.names = FALSE)

dotplot(ego_all_asc, showCategory = 20)
dotplot(ekegg_asc, showCategory = 20)


##4.7 ASC GSVA （Hallmark，ASC_P vs ASC_N）

library(GSVA)
library(msigdbr)

## 4.7.1 从 expr_imp_asc + ann_asc 得到 ASC 基因层表达矩阵
library(dplyr)
library(tibble)

expr_imp_asc_with_gene2 <- expr_imp_asc %>%
  as.data.frame() %>%
  rownames_to_column("Accession") %>%   # 把行名变成一列
  inner_join(
    res_ASC_anno_gene[, c("Accession", "Gene")] %>% distinct(),
    by = "Accession"
  ) %>%
  relocate(Gene, .after = Accession)    # 把 Gene 列挪到 Accession 后面

head(expr_imp_asc_with_gene2)

library(dplyr)


## 按 Gene 聚合，相同基因的多个 Accession 取中位数
expr_asc_gene_median <- expr_imp_asc_with_gene2 %>%
  group_by(Gene) %>%
  summarise(
    across(
      .cols = where(is.numeric),      # 所有数值列：ASC_N,1~ASC_P,3
      .fns  = ~ median(.x, na.rm = TRUE)
    )
  ) %>%
  ungroup()

## 变成基因×样本矩阵
expr_asc_gene_median <- as.data.frame(expr_asc_gene_median)
rownames(expr_asc_gene_median) <- expr_asc_gene_median$Gene
expr_asc_gene_median$Gene <- NULL

dim(expr_asc_gene_median)
head(expr_asc_gene_median[, 1:3])

write.csv(expr_asc_gene_median,
          "expr_asc_gene_median.csv",
          row.names = T)

expr_asc_for_gsva <- as.matrix(expr_asc_gene_median)  # 行 = 基因，列 = ASC 样本


## 4.7.2 Hallmark 基因集（小鼠）
m_hallmark_asc <- msigdbr(species = "Mus musculus", 
                          category = "H")

hallmark_list_asc <- split(m_hallmark_asc$gene_symbol, m_hallmark_asc$gs_name)

## 确认交集
length(intersect(rownames(expr_asc_for_gsva),
                 unique(m_hallmark_asc$gene_symbol)))

## 4.7.3 GSVA：ASC 样本
set.seed(123)
asc_gsva <- gsva(
  expr          = expr_asc_for_gsva,
  gset.idx.list = hallmark_list_asc,
  method        = "gsva",
  kcdf          = "Gaussian",
  mx.diff       = TRUE,
  parallel.sz   = 1
)

dim(asc_gsva)  # 行 = hallmark 通路，列 = ASC 样本

write.csv(asc_gsva,
          "asc_gsva_Hallmark.csv",
          row.names = T)
asc_gsva_Hallmark= as.data.frame(asc_gsva)

## 4.7.4 GSVA limma分析：ASC_P vs ASC_N

## 保证列顺序与 sample_info_asc 一致
asc_gsva <- asc_gsva[, sample_info_asc$sample, drop = FALSE]
group_asc <- droplevels(sample_info_asc$condition)  # ASC_N / ASC_P

design_asc_gsva <- model.matrix(~ 0 + group_asc)
colnames(design_asc_gsva) <- levels(group_asc)

contrast_asc_gsva <- makeContrasts(
  ASC_P_vs_ASC_N = ASC_P - ASC_N,
  levels = design_asc_gsva
)

fit_asc_gsva <- lmFit(asc_gsva, design_asc_gsva)
fit_asc_gsva <- contrasts.fit(fit_asc_gsva, contrast_asc_gsva)
fit_asc_gsva <- eBayes(fit_asc_gsva)

res_ASC_gsva <- topTable(
  fit_asc_gsva,
  coef    = "ASC_P_vs_ASC_N",
  number  = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("Pathway") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.2, "sig", "ns")
  )

write.csv(res_ASC_gsva,
          "GSVA_ASC_P_vs_ASC_N_Hallmark_all.csv",
          row.names = FALSE)

write.csv(filter(res_ASC_gsva, sig == "sig"),
          "GSVA_ASC_P_vs_ASC_N_Hallmark_sig.csv",
          row.names = FALSE)

## KEGG GSVA 的基因集
m_c2_mouse <- msigdbr(
  species  = "Mus musculus",
  category = "C2"
)

unique(m_c2_mouse$gs_subcat)


m_kegg_asc <- msigdbr(
  species     = "Mus musculus",
  category    = "C2",
)


kegg_list_asc <- split(m_kegg_asc$gene_symbol,
                       m_kegg_asc$gs_name)

asc_gsva_kegg <- gsva(
  expr          = expr_asc_for_gsva,
  gset.idx.list = kegg_list_asc,
  method        = "gsva",
  kcdf          = "Gaussian",
  mx.diff       = TRUE,
  parallel.sz   = 1
)

dim(asc_gsva_kegg)  # 行 = hallmark 通路，列 = ASC 样本

write.csv(asc_gsva_kegg,
          "asc_gsva_kegg.csv",
          row.names = T)

asc_gsva_kegg= as.data.frame(asc_gsva_kegg)




## 保证列顺序与 sample_info_asc 一致
asc_gsva_kegg <- asc_gsva_kegg[, sample_info_asc$sample, drop = FALSE]
group_asc <- droplevels(sample_info_asc$condition)  # ASC_N / ASC_P


fit_asc_gsva_kegg <- lmFit(asc_gsva_kegg, design_asc_gsva)
fit_asc_gsva_kegg <- contrasts.fit(fit_asc_gsva_kegg, contrast_asc_gsva)
fit_asc_gsva_kegg <- eBayes(fit_asc_gsva_kegg)

res_ASC_gsva_kegg <- topTable(
  fit_asc_gsva_kegg,
  coef    = "ASC_P_vs_ASC_N",
  number  = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("Pathway") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.2, "sig", "ns")
  )

write.csv(res_ASC_gsva_kegg,
          "GSVA_ASC_P_vs_ASC_N__kegg_all.csv",
          row.names = FALSE)

write.csv(filter(res_ASC_gsva_kegg, sig == "sig"),
          "GSVA_ASC_P_vs_ASC_N_kegg_sig.csv",
          row.names = FALSE)


###提取筛选后的GSVA
df_asc_gsva_sel <- read.csv("GSVA_ASC_P_vs_ASC_N__sig_slected.csv")


library(dplyr)
library(ggplot2)

df_plot <- df_asc_gsva_sel %>%
  mutate(
    direction = ifelse(logFC > 0, "Up_in_ASC_P", "Down_in_ASC_P"),
    Pathway   = factor(Pathway, levels = Pathway[order(logFC)]),
    FDR_label = sprintf("FDR=%.3g", FDR)
  )

p.ASC.gava.plot=ggplot(df_plot, aes(x = logFC, y = Pathway, fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(
    aes(label = FDR_label,
        hjust  = ifelse(logFC >= 0, -0.05, 1.05)),
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "Up_in_ASC_P"   = "#D95F02",
      "Down_in_ASC_P" = "#1B9E77"
    ),
    name = NULL,
    labels = c(
      "Up_in_ASC_P"   = "Higher in CD105⁺CD90⁺ ASC",
      "Down_in_ASC_P" = "Lower in CD105⁺CD90⁺ ASC"
    )
  ) +
  xlab("GSVA logFC (ASC_P vs ASC_N)") +
  ylab("") +
  ggtitle("Key ASC pathways associated with CD4 T-cell anergy") +
  theme_bw() +
  theme(
    plot.title  = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  ) +
  coord_cartesian(
    xlim = c(min(df_plot$logFC) - 0.3,
             max(df_plot$logFC) + 0.3)
  )


p.ASC.gava.plot

ggsave(filename = "ASC protein GSVA.plot.pdf",plot =p.ASC.gava.plot,width = 16,height = 6,dpi = 600)

###热图展示

##按照列合并
asc_gsva_Hallmark_kegg <- rbind(asc_gsva_Hallmark, asc_gsva_kegg)

df_asc_gsva_sel[,1]


ASC_GSVA_selected =asc_gsva_Hallmark_kegg[df_asc_gsva_sel[,1],]


# 列注释：ASC_N / ASC_P
## 1. 设定 condition 的顺序：先 ASC_N，再 ASC_P
sample_info_asc$condition <- factor(sample_info_asc$condition,
                                    levels = c("ASC_N", "ASC_P"))

## 2. 按 condition 排序样本
sample_info_asc <- sample_info_asc[order(sample_info_asc$condition), ]

## 3. 按样本顺序重排热图矩阵的列
ASC_GSVA_selected <- ASC_GSVA_selected[, sample_info_asc$sample]

## 4. 构建列注释
annot_col_asc <- data.frame(
  condition = sample_info_asc$condition
)
rownames(annot_col_asc) <- sample_info_asc$sample


## 先定义你想要的行顺序
row_order <- c(
  "WP_PLATELETMEDIATED_INTERACTIONS_WITH_VASCULAR_AND_CIRCULATING_CELLS",
  "REACTOME_ARACHIDONATE_PRODUCTION_FROM_DAG",
  "MIZUKAMI_HYPOXIA_UP",
  "RAMJAUN_APOPTOSIS_BY_TGFB1_VIA_SMAD4_UP",
  "KEGG_MEDICUS_REFERENCE_LECTIN_PATHWAY_OF_COMPLEMENT_CASCADE_C4_C2_TO_C3_CONVERTASE_FORMATION",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

## 按这个顺序重排矩阵行（顺便防一下有缺的）
row_order_use <- intersect(row_order, rownames(ASC_GSVA_selected))
ASC_GSVA_selected <- ASC_GSVA_selected[row_order_use, , drop = FALSE]



## 5. 画图：关键是这里 cluster_cols = FALSE
## 画图：关键是 cluster_rows = FALSE，保持你手动设定的行顺序
pdf("ASC_GSVA_heatmap.pdf", width = 10, height = 4)
pheatmap(
  ASC_GSVA_selected,
  scale = "row",
  cluster_rows = FALSE,              # ← 不对行聚类，行顺序照 row_order 来
  cluster_cols = FALSE,              # ← 列顺序前面已经处理过
  annotation_col = annot_col_asc,
  main = "GSVA in ASC_N vs ASC_P",
  fontsize_row = 6,
  fontsize_col = 4,
  cellwidth = 20,
  cellheight = 14,
  show_rownames = TRUE,
  show_colnames = FALSE,
  border_color = "white",
  treeheight_row = 0,
  treeheight_col = 0
)
dev.off()

graphics.off()



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

# 画热图（投稿风格参数你之前已经在用）
library(pheatmap)
pheatmap(
  asc_mat,
  scale = "row",
  annotation_col = data.frame(
    condition = group_asc,
    row.names = colnames(asc_mat)
  ),
  treeheight_row = 0,
  treeheight_col = 0,
  show_rownames  = TRUE,
  main = "ASC-derived suppressive/abortive pathways (ASC_P vs ASC_N)",
  border_color   = "white",
  fontsize_row   = 6,
  fontsize_col   = 8
)

library(dplyr)
library(pheatmap)

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














## ========== 5. CD4 分支：CD4_P vs CD4_N ==========

cat("==== CD4 分支开始 ====\n")

# 5.1 过滤 + 插补
cd4_proc <- filter_and_impute(
  expr_log2_sub       = expr_log2_cd4,
  group_vector        = sample_info_cd4$condition,
  min_nonNA_per_group = 2,
  q_minprob           = 0.01
)

expr_imp_cd4 <- cd4_proc$expr_imp
keep_cd4     <- cd4_proc$keep_idx_logical

ann_cd4 <- ann_all[keep_cd4, ]
stopifnot(nrow(expr_imp_cd4) == nrow(ann_cd4))
rownames(expr_imp_cd4) <- ann_cd4$Accession

# 5.2 PCA（CD4）
pca_cd4 <- prcomp(t(expr_imp_cd4), scale. = FALSE)
pca_cd4_df <- data.frame(
  PC1 = pca_cd4$x[, 1],
  PC2 = pca_cd4$x[, 2],
  sample_info_cd4
)

p.cd4.pca=ggplot(pca_cd4_df,
       aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  theme_bw() +
  ggtitle("PCA - CD4 only")

ggsave(filename = "CD4 protein PCA.plot.pdf",
       plot =p.cd4.pca,
       width = 8,
       height = 6,
       dpi = 300)




# 5.3 limma 差异分析：CD4_P vs CD4_N
group_cd4 <- droplevels(sample_info_cd4$condition)  # CD4_N, CD4_P

design_cd4 <- model.matrix(~ 0 + group_cd4)
colnames(design_cd4) <- levels(group_cd4)
print(design_cd4)

contrast_cd4 <- makeContrasts(
  CD4_P_vs_CD4_N = CD4_P - CD4_N,
  levels = design_cd4
)

fit_cd4 <- lmFit(expr_imp_cd4, design_cd4)
fit_cd4 <- contrasts.fit(fit_cd4, contrast_cd4)
fit_cd4 <- eBayes(fit_cd4)

res_CD4 <- topTable(fit_cd4, coef = "CD4_P_vs_CD4_N",
                    number = Inf, sort.by = "P")

res_CD4_anno <- res_CD4 %>%
  rownames_to_column("Accession") %>%
  left_join(ann_cd4, by = "Accession") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.585, "sig", "ns")
  )

write.csv(res_CD4_anno,
          "DE_CD4_P_vs_CD4_N_cd4_only_all.csv",
          row.names = FALSE)

write.csv(filter(res_CD4_anno, sig == "sig"),
          "DE_CD4_P_vs_CD4_N_cd4_only_sig.csv",
          row.names = FALSE)

# 5.4 火山图 + 热图（CD4）
p.cd4.volcano=ggplot(res_CD4_anno,
       aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.7, size = 0.5) +
  scale_color_manual(values = c("ns" = "grey70", "sig" = "red")) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed") +
  geom_hline(yintercept = -1*log10(0.05), linetype = "dashed") +
  theme_bw() +
  ggtitle("Volcano - CD4_P vs CD4_N")


p.cd4.volcano

ggsave(filename = "CD4 protein Volcano.plot.pdf",,
       plot = p.cd4.volcano,
       width = 7,
       height =6 )





top50_cd4 <- res_CD4_anno %>%
  filter(sig == "sig") %>%
  arrange(FDR) %>%
  slice(1:min(50, n())) %>%
  pull(Accession)

mat_top50_cd4 <- expr_imp_cd4[top50_cd4, ]

pdf("CD4_P vs CD4_N top50 DEPs.pdf", width = 4, height = 5)
pheatmap(mat_top50_cd4,
         scale = "row",
         annotation_col = data.frame(condition = group_cd4,
                                     row.names = colnames(mat_top50_cd4)),
         treeheight_row = 0,
         treeheight_col = 0,
        
         show_rownames = FALSE,
         main = "CD4_P vs CD4_N - top50 DEPs")
dev.off()


pdf("CD4_P vs CD4_N DEPs.pdf", width = 4, height = 5)
p.all.asc=pheatmap(expr_imp_cd4,
                   scale = "row",
                   annotation_col = data.frame(condition = group_cd4,
                                               row.names = colnames(mat_top50_cd4)),
                   treeheight_row = 0,
                   treeheight_col = 0,
                   show_rownames = FALSE,
                   main = "CD4_P vs CD4_N",
                   border_color="white")
dev.off()





## ========== 6. （可选）UniProt -> 基因名 & 富集分析骨架 ==========

## 一次性映射所有 Accession 到 SYMBOL/ENTREZID
uniprot_ids <- ann_all$Accession |> unique() |> na.omit()
id_map <- bitr(uniprot_ids,
               fromType = "UNIPROT",
               toType   = c("SYMBOL", "ENTREZID"),
               OrgDb    = org.Mm.eg.db)

## 把 SYMBOL 合并到 ASC / CD4 DE 结果中
res_ASC_anno <- res_ASC_anno %>%
  left_join(id_map, by = c("Accession" = "UNIPROT"))

res_CD4_anno <- res_CD4_anno %>%
  left_join(id_map, by = c("Accession" = "UNIPROT"))

write.csv(res_ASC_anno,
          "DE_ASC_P_vs_ASC_N_asc_only_all_withSymbol.csv",
          row.names = FALSE)
write.csv(res_CD4_anno,
          "DE_CD4_P_vs_CD4_N_cd4_only_all_withSymbol.csv",
          row.names = FALSE)

## 以 CD4 上调基因为例做 GO/KEGG
genes_cd4_up <- res_CD4_anno %>%
  filter(sig == "sig", logFC > 0, !is.na(SYMBOL.x)) %>%
  pull(SYMBOL.x) %>%
  unique()

ent_cd4_up <- bitr(genes_cd4_up,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Mm.eg.db)

ego_cd4_up <- enrichGO(
  gene         = ent_cd4_up$ENTREZID,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "ALL", # 选择GO本体类型 ("BP"生物过程, "MF"分子功能, "CC"细胞组分, "ALL"全部
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable     = TRUE
)

write.csv(as.data.frame(ego_cd4_up),
          "GO_CD4_P_vs_CD4_N_up.csv",
          row.names = FALSE)

dotplot(ego_cd4_up, showCategory = 20)

kegg_CD4_up <- enrichKEGG(gene         = ent_cd4_up$ENTREZID,         # 您的基因列表
                           organism     = 'mmu',            # 物种代码："mmu"鼠'hsa' 代表智人 (Homo sapiens)
                           pvalueCutoff = 0.05,             # P值阈值
                           qvalueCutoff = 0.2,              # Q值阈值
                           pAdjustMethod = "BH",            # P值校正方法
)

# 查看富集结果的前几行
head(kegg_CD4_up)


kegg_results_readable <- setReadable(
  kegg_CD4_up,             # 需要转换的富集结果对象
  OrgDb = org.Mm.eg.db,    # 指定物种注释数据库（这里是小鼠）
  keyType = "ENTREZID"     # 指定输入 ID 的类型（默认为 ENTREZID，通常不需要修改）
)

write.csv(as.data.frame(kegg_results_readable),
          "kegg_CD4_P_vs_CD4_N_up.csv",
          row.names = FALSE)


## ================== CD4 GSVA（基于 gene-level expr_imp_cd4） ==================

library(GSVA)
library(msigdbr)

library(clusterProfiler)
library(org.Mm.eg.db)

## 保证 ann_cd4 与 expr_imp_cd4 对齐
stopifnot(nrow(expr_imp_cd4) == nrow(ann_cd4))

## 1) 从 Description 里提取 GN= 后面的基因符号
extract_gene_from_desc <- function(desc) {
  # 有 GN= 才抽，否则 NA
  ifelse(
    grepl(" GN=", desc),
    sub(".* GN=([^ ]+).*", "\\1", desc),
    NA
  )
}

gene_cd4 <- extract_gene_from_desc(ann_cd4$Description)

## 2) 把基因名挂到表达矩阵上
cd4_gene_df <- cbind(Gene = gene_cd4, as.data.frame(expr_imp_cd4))

## 去掉没有 GN 的行
cd4_gene_df <- cd4_gene_df[!is.na(cd4_gene_df$Gene) & cd4_gene_df$Gene != "", ]

library(dplyr)

## 3) 每个基因可能有多条蛋白记录，使用中位数聚合
cd4_gene_df2 <- cd4_gene_df %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
  ungroup()

## 4) 转成 data.frame 再设行名（这样不会有 tibble 的 warning）
cd4_gene_df2 <- as.data.frame(cd4_gene_df2)

rownames(cd4_gene_df2) <- cd4_gene_df2$Gene
cd4_gene_df2$Gene <- NULL

## 最终的基因层表达矩阵：
expr_cd4_for_gsva <- as.matrix(cd4_gene_df2)

dim(expr_cd4_for_gsva)
head(rownames(expr_cd4_for_gsva))

library(msigdbr)
library(GSVA)

## 1) 准备基因集：小鼠 Hallmark
m_hallmark <- msigdbr(species = "Mus musculus", category = "H")
hallmark_list <- split(m_hallmark$gene_symbol, m_hallmark$gs_name)

## 2) 简单检查：交集基因数不为零才靠谱
length(intersect(rownames(expr_cd4_for_gsva),
                 unique(m_hallmark$gene_symbol)))
# 正常应该是几百/上千，>0 即可

## 3) GSVA
set.seed(123)
cd4_gsva <- gsva(
  expr          = expr_cd4_for_gsva,  # 基因 × 样本
  gset.idx.list = hallmark_list,
  method        = "gsva",
  kcdf          = "Gaussian",
  mx.diff       = TRUE,
  parallel.sz   = 1
)

dim(cd4_gsva)  # 行 = 通路，列 = 样本

group_cd4 <- droplevels(sample_info_cd4$condition)  # CD4_N / CD4_P

design_cd4_gsva <- model.matrix(~ 0 + group_cd4)
colnames(design_cd4_gsva) <- levels(group_cd4)

contrast_cd4_gsva <- makeContrasts(
  CD4_P_vs_CD4_N = CD4_P - CD4_N,
  levels = design_cd4_gsva
)

fit_cd4_gsva <- lmFit(cd4_gsva, design_cd4_gsva)
fit_cd4_gsva <- contrasts.fit(fit_cd4_gsva, contrast_cd4_gsva)
fit_cd4_gsva <- eBayes(fit_cd4_gsva)

res_cd4_gsva <- topTable(
  fit_cd4_gsva,
  coef    = "CD4_P_vs_CD4_N",
  number  = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("Pathway") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.2, "sig", "ns")
  )
res_cd4_gsva


write.csv(res_cd4_gsva,
          "GSVA_CD4_P_vs_CD4_N_Hallmark_all.csv",
          row.names = FALSE)
write.csv(filter(res_cd4_gsva, sig == "sig"),
          "GSVA_CD4_P_vs_CD4_N_Hallmark_sig.csv",
          row.names = FALSE)

colnames(res_cd4_gsva)
library(dplyr)
res_cd4_gsva_sig=filter(res_cd4_gsva, sig == "sig")

cd4_gsva_sig_selected <- res_cd4_gsva_sig %>%
  filter(Pathway %in% c("HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                        "HALLMARK_COAGULATION",
                        "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                        "HALLMARK_TGF_BETA_SIGNALING",
                        "HALLMARK_XENOBIOTIC_METABOLISM",
                        "HALLMARK_PROTEIN_SECRETION",
                        "HALLMARK_COMPLEMENT",
                        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                        "HALLMARK_E2F_TARGETS",
                        "HALLMARK_MYC_TARGETS_V1",
                        "HALLMARK_MYC_TARGETS_V2",
                        "HALLMARK_G2M_CHECKPOINT"))
cd4_gsva_sig_selected

## 打开一个 PDF 设备
pdf("CD4_proliferation_heatmap.pdf", width = 6, height = 4)

cd4_gsva=as.matrix.data.frame(cd4_gsva)

cd4_GSVA_selected =cd4_gsva[c("HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                              "HALLMARK_COAGULATION",
                              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                              "HALLMARK_TGF_BETA_SIGNALING",
                              "HALLMARK_XENOBIOTIC_METABOLISM",
                              "HALLMARK_PROTEIN_SECRETION",
                              "HALLMARK_COMPLEMENT",
                              "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                              "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                              "HALLMARK_E2F_TARGETS",
                              "HALLMARK_MYC_TARGETS_V1",
                              "HALLMARK_MYC_TARGETS_V2",
                              "HALLMARK_G2M_CHECKPOINT"),]
## 打开一个 PDF 设备
pdf("CD4_GSVA_heatmap.pdf", width = 6, height = 4)
## 在这个设备上画图
pheatmap(
  cd4_GSVA_selected,
  scale = "row",  # 行标准化，看相对高低
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  annotation_col = annot_col,
  main = "GSVA in CD4_N vs CD4_P",
  fontsize_row = 6,
  fontsize_col = 4,
  cellwidth = 20,
  cellheight =14,
  show_rownames = T,
  show_colnames = FALSE,
  border_color = "white",
  treeheight_row = 0,
  treeheight_col = 0,
)
## 关闭设备，文件写入完成
dev.off()

library(dplyr)
library(ggplot2)

## 1. 选取四个增殖相关 hallmark
paths_cycle <- c(
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_COAGULATION",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_G2M_CHECKPOINT"
)

df_cycle <- res_cd4_gsva %>%
  filter(Pathway %in% paths_cycle) %>%
  select(Pathway, logFC, FDR)

## 2. 调整顺序：按 logFC 从高到低排，方便视觉
df_cycle <- df_cycle %>%
  arrange(logFC) %>%               # 这里是从低到高，方便做横图
  mutate(Pathway = factor(Pathway, levels = Pathway))

colnames(df_cycle)

write.csv(
  df_cycle,
  "GSVA_CD4_P_vs_CD4_N_selected_pathways.csv",
  row.names = FALSE
)

## 3. 画图：横向条形图（CD4_P vs CD4_N）
library(dplyr)
library(ggplot2)

## 3.1. 构造方向列 + 排序
df_cycle2 <- df_cycle %>%
  mutate(
    direction = ifelse(logFC > 0, "Up_in_CD4_P", "Down_in_CD4_P"),
    Pathway   = factor(Pathway, levels = Pathway[order(logFC)]),
    label     = sprintf("FDR=%.3g", FDR)   # 只标注 FDR
  )

## 3.2. 画横向条形图
p.CD4.protein=ggplot(df_cycle2, aes(x = logFC, y = Pathway, fill = direction)) +
                     geom_col(width = 0.7) +
                      geom_vline(xintercept = 0, linetype = "dashed") +
                    geom_text(
                              aes(label = label,
                               hjust  = ifelse(logFC >= 0, -0.1, 1.1)),
                          size = 2
                              ) +
                     scale_fill_manual(
              values = c(
             "Up_in_CD4_P"   = "#D95F02",
           "Down_in_CD4_P" = "#1B9E77"
                         ),
                name = NULL,
                labels = c("Up_in_CD4_P" = "Higher in CD4_P",
               "Down_in_CD4_P" = "Lower in CD4_P")
               ) +
          xlab("logFC (CD4_P vs CD4_N)") +
          ylab("") +
         ggtitle("Proliferation-related Hallmark GSVA (CD4_P vs CD4_N)") +
          theme_bw() +
             theme(
           plot.title = element_text(hjust = 0.5),
              axis.text.y = element_text(size = 9),
               legend.position = "bottom"
              ) +
             coord_cartesian(
                 xlim = c(min(df_cycle2$logFC) - 0.3,
                        max(df_cycle2$logFC) + 0.3)
                     ) 

p.CD4.protein

ggsave(filename = "CD4 protein GSVA.pdf",plot =p.CD4.protein,width = 8,height = 6,dpi = 600)

library(dplyr)
library(pheatmap)

## 一组典型的增殖/细胞周期基因（人符号，小鼠看你数据是人还是鼠）
genes_prolif <- c(
  "MKI67",   # Ki-67
  "PCNA",
  "CDK1",
  "CCNB1",
  "CCNB2",
  "CDC20",
  "TOP2A",
  "MCM2",
  "MCM4",
  "MCM6"
)
##在表达矩阵里匹配这些基因（大小写不敏感）
# 假定 expr_cd4_for_gsva 的行名是 gene symbol
all_genes      <- rownames(expr_cd4_for_gsva)
all_genes_up   <- toupper(all_genes)
genes_prolif_up <- toupper(genes_prolif)

idx <- all_genes_up %in% genes_prolif_up

prolif_genes_in_mat <- all_genes[idx]
prolif_genes_in_mat
# 看看实际匹配到了哪些

mat_prolif <- expr_cd4_for_gsva[prolif_genes_in_mat, , drop = FALSE]

# 列注释：CD4_N / CD4_P
annot_col <- data.frame(
  condition = sample_info_cd4$condition
)
rownames(annot_col) <- sample_info_cd4$sample  # 要和列名一致

# 确保顺序一致
mat_prolif <- mat_prolif[, rownames(annot_col), drop = FALSE]

# 画热图（row scale 展示相对高低）
## 打开一个 PDF 设备
pdf("CD4_proliferation_heatmap.pdf", width = 6, height = 4)

## 在这个设备上画图
pheatmap(
  mat_prolif,
  scale = "row",  # 行标准化，看相对高低
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  annotation_col = annot_col,
  main = "Proliferation-related genes in CD4_N vs CD4_P",
  fontsize_row = 9,
  fontsize_col = 4,
  cellwidth = 20,
  cellheight =14,
  show_rownames = T,
  show_colnames = FALSE,
  border_color = "white",
  treeheight_row = 0,
  treeheight_col = 0,
)
## 关闭设备，文件写入完成
dev.off()


## KEGG GSVA for CD4的基因集
m_c2_mouse <- msigdbr(
  species  = "Mus musculus",
  category = "C2"
)

unique(m_c2_mouse$gs_subcat)


m_kegg_asc <- msigdbr(
  species     = "Mus musculus",
  category    = "C2",
)


kegg_list_asc <- split(m_kegg_asc$gene_symbol,
                       m_kegg_asc$gs_name)

CD4_gsva_kegg <- gsva(
  expr          = expr_cd4_for_gsva,
  gset.idx.list = kegg_list_asc,
  method        = "gsva",
  kcdf          = "Gaussian",
  mx.diff       = TRUE,
  parallel.sz   = 1
)

dim(CD4_gsva_kegg)  # 行 = 通路，列 = 样本

group_cd4 <- droplevels(sample_info_cd4$condition)  # CD4_N / CD4_P

design_cd4_gsva <- model.matrix(~ 0 + group_cd4)
colnames(design_cd4_gsva) <- levels(group_cd4)

contrast_cd4_gsva <- makeContrasts(
  CD4_P_vs_CD4_N = CD4_P - CD4_N,
  levels = design_cd4_gsva
)

fit_cd4_gsva_kegg <- lmFit(CD4_gsva_kegg, design_cd4_gsva)
fit_cd4_gsva_kegg <- contrasts.fit(fit_cd4_gsva_kegg, contrast_cd4_gsva)
fit_cd4_gsva_kegg <- eBayes(fit_cd4_gsva_kegg)

library(tibble)

res_cd4_gsva_kegg <- topTable(
  fit_cd4_gsva_kegg,
  coef    = "CD4_P_vs_CD4_N",
  number  = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("Pathway") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.2, "sig", "ns")
  )
res_cd4_gsva_kegg

write.csv(res_cd4_gsva_kegg,
          "GSVA_CD4_P_vs_CD4_N_kegg_all.csv",
          row.names = FALSE)

write.csv(filter(res_cd4_gsva_kegg, sig == "sig"),
          "GSVA_CD4_P_vs_CD4_N_kegg_sig.csv",
          row.names = FALSE)

write.csv(CD4_gsva_kegg,
          "GSVA_CD4_P_vs_CD4_N_kegg.csv",
          row.names = T)


####提取模块可视化###
## 1. 细胞周期 / E2F / 有丝分裂相关（proliferative arrest）

cd4_cycle_paths <- c(
  "REACTOME_MITOTIC_TELOPHASE_CYTOKINESIS",
  "REACTOME_TRANSCRIPTION_OF_E2F_TARGETS_UNDER_NEGATIVE_CONTROL_BY_P107_RBL1_AND_P130_RBL2_IN_COMPLEX_WITH_HDAC1",
  "REACTOME_TRANSCRIPTION_OF_E2F_TARGETS_UNDER_NEGATIVE_CONTROL_BY_DREAM_COMPLEX",
  "REN_BOUND_BY_E2F",
  "KONG_E2F3_TARGETS"
)


## 2. IFN / JAK–STAT / cytokine signaling 下调

cd4_ifn_jak_stat_paths <- c(
  "KEGG_MEDICUS_REFERENCE_TYPE_II_INTERFERON_TO_JAK_STAT_SIGNALING_PATHWAY",
  "WP_TYPE_II_INTERFERON_SIGNALING",
  "ZHANG_INTERFERON_RESPONSE",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "GRANDVAUX_IFN_RESPONSE_NOT_VIA_IRF3",
  "DER_IFN_ALPHA_RESPONSE_UP",
  "DER_IFN_BETA_RESPONSE_UP",
  "DER_IFN_GAMMA_RESPONSE_UP",
  "KEGG_MEDICUS_REFERENCE_CYTOKINE_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_MEDICUS_REFERENCE_IFN_RIPK1_3_SIGNALING_PATHWAY",
  "KEGG_MEDICUS_REFERENCE_EGF_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_MEDICUS_VARIANT_ERBB2_OVEREXPRESSION_TO_EGF_JAK_STAT_SIGNALING_PATHWAY",
  "REACTOME_STAT3_NUCLEAR_EVENTS_DOWNSTREAM_OF_ALK_SIGNALING",
  "REACTOME_STAT5_ACTIVATION",
  "REACTOME_ERYTHROPOIETIN_ACTIVATES_STAT5"
)


## 3. TGF-β / SMAD 相关通路上调

cd4_tgfb_paths <- c(
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
  "REACTOME_SIGNALING_BY_TGFBR3",
  "REACTOME_TGFBR3_PTM_REGULATION",
  "VERRECCHIA_EARLY_RESPONSE_TO_TGFB1",
  "VERRECCHIA_RESPONSE_TO_TGFB1_C1",
  "VERRECCHIA_RESPONSE_TO_TGFB1_C4",
  "PLASARI_TGFB1_SIGNALING_VIA_NFIC_10HR_UP",
  "PLASARI_TGFB1_SIGNALING_VIA_NFIC_10HR_DN"
)


## 4. ROS / 缺氧 / OXPHOS 失衡（hypoxia–ROS↑ + OXPHOS↓）

cd4_ros_hypoxia_oxphos_paths <- c(
  "WP_HIF1A_AND_PPARG_IN_HYPERTROPHIC_CARDIOMYOPATHY",
  "FARDIN_HYPOXIA_11",
  "WINTER_HYPOXIA_DN",
  "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
  "GROSS_HYPOXIA_VIA_HIF1A_ONLY",
  "GROSS_HIF1A_TARGETS_UP",
  "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES",
  "WP_MITOCHONDRIAL_COMPLEX_I_ASSEMBLY_MODEL_OXPHOS_SYSTEM",
  "WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA",
  "KEGG_MEDICUS_REFERENCE_MITOCHONDRIAL_COMPLEX_UCP1_IN_THERMOGENESIS",
  "WP_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "MOOTHA_VOXPHOS"
)


## 5. 脂质 / 胆固醇 / 脂介质相关（lipid metabolic rewiring）

cd4_lipid_paths <- c(
  "GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_LIGHTGREEN_DN",
  "GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_MAGENTA_DN",
  "GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_GREEN_UP",
  "WUNDER_INFLAMMATORY_RESPONSE_AND_CHOLESTEROL_DN",
  "REACTOME_NR1H3_NR1H2_REGULATE_GENE_EXPRESSION_LINKED_TO_CHOLESTEROL_TRANSPORT_AND_EFFLUX",
  "WP_FATTY_ACID_BIOSYNTHESIS",
  "WP_SREBF_AND_MIR33_IN_CHOLESTEROL_AND_LIPID_HOMEOSTASIS",
  "WP_LIPID_METABOLISM_IN_SENESCENT_CELLS",
  "WP_LIPID_PARTICLES_COMPOSITION",
  "REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS_VIA_7ALPHA_HYDROXYCHOLESTEROL",
  "REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS_VIA_24_HYDROXYCHOLESTEROL",
  "REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS_VIA_27_HYDROXYCHOLESTEROL",
  "WP_FAMILIAL_HYPERLIPIDEMIA_TYPE_1",
  "WP_FAMILIAL_HYPERLIPIDEMIA_TYPE_3",
  "WP_FAMILIAL_HYPERLIPIDEMIA_TYPE_4",
  "WP_FAMILIAL_HYPERLIPIDEMIA_TYPE_5",
  "REACTOME_LYSOSPHINGOLIPID_AND_LPA_RECEPTORS"
)


## 6. 补体 / 凝血 / 血小板互作（complement–platelet niche）

cd4_complement_paths <- c(
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "BIOCARTA_LECTIN_PATHWAY",
  "KEGG_MEDICUS_REFERENCE_LECTIN_PATHWAY_OF_COMPLEMENT_CASCADE_C4_C2_TO_C3_CONVERTASE_FORMATION",
  "REACTOME_ALTERNATIVE_COMPLEMENT_ACTIVATION",
  "WP_COMPLEMENTMEDIATED_INFLAMMATION_OF_PULMONARY_ALVEOLUS_IN_COVID19_HYPOTHETICAL_PATHWAY",
  "REACTOME_SYNTHESIS_OF_PROSTAGLANDINS_PG_AND_THROMBOXANES_TX",
  "WP_PLATELETMEDIATED_INTERACTIONS_WITH_VASCULAR_AND_CIRCULATING_CELLS",
  "WP_COVID19_THROMBOSIS_AND_ANTICOAGULATION"
)


## 7. Notch / FOXP3/Treg-like 程序

cd4_notch_foxp3_paths <- c(
  "WP_APOPTOSISRELATED_NETWORK_DUE_TO_ALTERED_NOTCH3_IN_OVARIAN_CANCER",
  "REACTOME_SIGNALING_BY_NOTCH4",
  "REACTOME_NEGATIVE_REGULATION_OF_NOTCH4_SIGNALING",
  "WP_NOTCH1_REGULATION_OF_ENDOTHELIAL_CELL_CALCIFICATION",
  "BRAUNE_GEIST_20_GENE_NOTCH_SIG_BREAST_CANCER",
  "KEGG_MEDICUS_REFERENCE_NOTCH_MESP2_SIGNALING",
  "KEGG_MEDICUS_VARIANT_NOTCH_OVEREXPRESSION_TO_NOTCH_SIGNALING_PATHWAY",
  "ZHENG_FOXP3_TARGETS_IN_T_LYMPHOCYTE_DN",
  "GAVIN_FOXP3_TARGETS_CLUSTER_P3",
  "ZHENG_FOXP3_TARGETS_DN"
)


## 8. Tryptophan / IDO 轴（用来说明“不是这条路”）

cd4_trp_catabolism_paths <- c(
  "REACTOME_TRYPTOPHAN_CATABOLISM"
)


## 9. 合并为一个总的“CD4 anergy 机制通路列表”

cd4_mech_paths <- unique(c(
  cd4_cycle_paths,
  cd4_ifn_jak_stat_paths,
  cd4_tgfb_paths,
  cd4_ros_hypoxia_oxphos_paths,
  cd4_lipid_paths,
  cd4_complement_paths,
  cd4_notch_foxp3_paths,
  cd4_trp_catabolism_paths
))
length(cd4_mech_paths)


###1. 在 df_kegg 里标记每个通路属于哪个机制模块
###现在给CD4_gsva_kegg_sig(df_kegg)打标签：
library(dplyr)
library(ggplot2)
library(forcats)

##CD4_gsva_kegg_sig(df_kegg)是你 GSVA + limma 后的结果表
## 必须至少包含: Pathway, logFC, FDR

df_kegg=filter(res_cd4_gsva_kegg, sig == "sig")


df_cd4_mech <- df_kegg %>%
  filter(Pathway %in% cd4_mech_paths) %>%     # 只保留机制相关通路
  mutate(
    module = case_when(
      Pathway %in% cd4_cycle_paths              ~ "Cell cycle / E2F",
      Pathway %in% cd4_ifn_jak_stat_paths       ~ "IFN / JAK-STAT",
      Pathway %in% cd4_tgfb_paths               ~ "TGF-β / SMAD",
      Pathway %in% cd4_ros_hypoxia_oxphos_paths ~ "Hypoxia / ROS / OXPHOS",
      Pathway %in% cd4_lipid_paths              ~ "Lipid / cholesterol",
      Pathway %in% cd4_complement_paths         ~ "Complement / coagulation",
      Pathway %in% cd4_notch_foxp3_paths        ~ "Notch / FOXP3-like",
      Pathway %in% cd4_trp_catabolism_paths     ~ "Tryptophan catabolism",
      TRUE                                      ~ "Other"
    )
  ) %>%
  # 按模块和 logFC 排序 path，方便画图
  group_by(module) %>%
  arrange(logFC, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    Pathway = fct_inorder(Pathway)  # 按当前顺序固定 factor
  )

## 导出一份机制子集表到补充数据
write.csv(df_cd4_mech,
          "Supp_Table_CD4_GSVA_mechanism_modules.csv",
          row.names = FALSE)

#这一步做完，你就有一张干净的表：Pathway + logFC + FDR + module，可以作为一个补充表（Supp Table）



##2. 补充图：按模块分面展示 logFC（CD4_P vs CD4_N）


ggplot(df_cd4_mech, 
       aes(x = Pathway, y = logFC, fill = logFC)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(
    aes(
      label = ifelse(FDR < 0.01,
                     sprintf("FDR=%.2g", FDR),
                     ""),
      hjust = ifelse(logFC >= 0, -0.05, 1.05)
    ),
    size  = 2.5
  ) +
  coord_flip() +
  facet_wrap(~ module, scales = "free_y", ncol = 2) +
  scale_fill_gradient2(
    low  = "#1B9E77",
    mid  = "white",
    high = "#D95F02",
    midpoint = 0,
    name = "logFC\n(CD4_P vs CD4_N)"
  ) +
  xlab("") +
  ylab("GSVA logFC (CD4_P vs CD4_N)") +
  ggtitle("Mechanism-related pathways in CD4 GSVA (C2, FDR<0.05, |logFC|≥0.2)") +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 6),
    strip.text   = element_text(size = 9),
    legend.position = "right"
  )


library(dplyr)
library(forcats)
library(ggplot2)

## 1. 预处理：按 module 内 logFC 排序，固定 Pathway 顺序
df_cd4_mech_plot <- df_cd4_mech %>%
  group_by(module) %>%
  arrange(logFC, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    Pathway = fct_inorder(Pathway),
    FDR_label = sprintf("FDR=%.2g", FDR)   # 提前算好标签
  )

## 2. 画图：按模块分面、横向条形图
p.cd4.mech.plot=ggplot(df_cd4_mech_plot,
       aes(x = logFC, y = Pathway, fill = logFC)) +
  geom_col(width = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_text(
    aes(
      label = FDR_label,
      hjust = ifelse(logFC >= 0, -0.1, 1.1)
    ),
    size = 1.8
  ) +
  facet_wrap(~ module, scales = "free_y", ncol = 2) +
  scale_fill_gradient2(
    low  = "#1B9E77",
    mid  = "white",
    high = "#D95F02",
    midpoint = 0,
    name = "logFC\n(CD4_P vs CD4_N)"
  ) +
  scale_x_continuous(
    name   = "GSVA logFC (CD4_P vs CD4_N)",
    expand = expansion(mult = 0.15)
  ) +
  ylab("") +
  ggtitle("Mechanism-related pathways in CD4 GSVA (C2, FDR<0.05, |logFC|≥0.2)") +
  theme_bw(base_size = 9) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    axis.text.x       = element_text(size = 6),
    axis.text.y       = element_text(size = 6),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text        = element_text(size = 9),
    legend.position   = "right"
  )

p.cd4.mech.plot

ggsave(filename = "CD4 protein GSVA.mech.plot.pdf",plot =p.cd4.mech.plot,width = 24,height = 12,dpi = 600)


############################################################
## 结束
###########################################################







getwd()
setwd("D:/R/20260126/")

# ==============================================================================
# 1. 环境准备与数据读取
# ==============================================================================
# 检查并安装必要的包
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("ggrepel")) install.packages("ggrepel")

library(tidyverse)
library(pheatmap)
library(ggrepel)

# 读取数据
raw_df <- read.csv("20260126 bulk ASC.csv", stringsAsFactors = FALSE)

# 数据清洗：去除基因名中的引号
raw_df$Gene.Symbol <- gsub("'", "", raw_df$Gene.Symbol)

# 定义样本列名 (根据您的最新定义)
cols_base_P <- c("Normalized.read.count.ASC_P_1", "Normalized.read.count.ASC_P_2", "Normalized.read.count.ASC_P_3")
cols_base_N <- c("Normalized.read.count.ASC_N_1", "Normalized.read.count.ASC_N_2", "Normalized.read.count.ASC_N_3")
cols_after_P <- c("Normalized.read.count.after_ASC_P_1", "Normalized.read.count.after_ASC_P_2", "Normalized.read.count.after_ASC_P_3")
cols_after_N <- c("Normalized.read.count.after_ASC_N_1", "Normalized.read.count.after_ASC_N_2", "Normalized.read.count.after_ASC_N_3")

# 过滤低表达基因 (所有样本平均 count > 10)
all_cols <- c(cols_base_P, cols_base_N, cols_after_P, cols_after_N)
raw_df$mean_expr <- rowMeans(raw_df[, all_cols])
df <- raw_df %>% filter(mean_expr > 10)

# ==============================================================================
# 2. 差异分析 (Differential Expression Analysis)
# ==============================================================================

# 辅助函数：计算 Log2FC 和 P-value
calc_deg <- function(data, group1_cols, group2_cols) {
  # Group1 vs Group2 (LogFC > 0 mean Group1 is higher)
  mean1 <- rowMeans(data[, group1_cols])
  mean2 <- rowMeans(data[, group2_cols])
  
  # Log2 Fold Change (添加 1 防止除零)
  logfc <- log2((mean1 + 1) / (mean2 + 1))
  
  # T-test P-value
  pvals <- apply(data, 1, function(row) {
    tryCatch(t.test(log1p(as.numeric(row[group1_cols])), 
                    log1p(as.numeric(row[group2_cols])))$p.value, 
             error = function(e) 1)
  })
  
  return(data.frame(Gene = data$Gene.Symbol, LogFC = logfc, PValue = pvals))
}

# --- 分析 1: 身份验证 (Baseline: ASC_P vs ASC_N) ---
res_identity <- calc_deg(df, cols_base_P, cols_base_N)
res_identity$Type <- "Identity (Base P vs N)"

# --- 分析 2: 功能机制 (Co-culture: After_P vs After_N) ---
# 这是解释抑制功能差异的核心
res_mechanism <- calc_deg(df, cols_after_P, cols_after_N)
res_mechanism$Type <- "Mechanism (After P vs N)"

# --- 分析 3: 动态响应 (Induction: After_P vs Base_P) ---
# 解释基质细胞如何响应 T 细胞
res_induction <- calc_deg(df, cols_after_P, cols_base_P)
res_induction$Type <- "Induction (After P vs Base P)"

# ==============================================================================
# 3. 结果可视化：火山图 (针对机制分析)
# ==============================================================================

# 标记显著性
res_mechanism$Sig <- ifelse(res_mechanism$PValue < 0.05 & abs(res_mechanism$LogFC) > 1,
                            ifelse(res_mechanism$LogFC > 1, "Up in After_P", "Down in After_P"), "NS")

# 挑选需要标注的基因 (基于您的课题假设)
# 身份基因: Eng (CD105)
# 耐受基因: Jag1 (Notch), Nt5e (CD73), Aldh1a2 (RA), Lgals9 (Galectin)
# 血管基因: Emcn, Pecam1
genes_to_label <- c("Eng", "Jag1", "Nt5e", "Aldh1a2", "Il33", "Lgals9", "Emcn", "Tgfb1")

ggplot(res_mechanism, aes(x = LogFC, y = -log10(PValue))) +
  geom_point(aes(color = Sig), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Down in After_P" = "gray", "NS" = "gray90", "Up in After_P" = "#E41A1C")) +
  geom_text_repel(data = subset(res_mechanism, Gene %in% genes_to_label),
                  aes(label = Gene), box.padding = 0.5, max.overlaps = Inf, fontface = "bold") +
  theme_classic() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Mechanism: After_ASC_P vs After_ASC_N",
       subtitle = "Genes enriched in the suppressive niche after co-culture",
       x = "Log2 Fold Change (P/N)", y = "-Log10 P-value")

ggsave("Fig_Mechanism_Volcano.pdf", width = 6, height = 5)

# ==============================================================================
# 4. 结果可视化：重点基因热图 (Heatmap)
# ==============================================================================

# 提取重点基因的表达量
heatmap_genes <- c("Eng", "Thy1",           # Markers
                   "Jag1", "Dll4",          # Notch Ligands
                   "Nt5e", "Aldh1a2",       # Enzymes (CD73, RA)
                   "Tgfb1", "Il33",         # Cytokines
                   "Emcn", "Pecam1")        # Vascular

plot_data <- df %>% 
  filter(Gene.Symbol %in% heatmap_genes) %>%
  select(Gene.Symbol, all_of(cols_base_N), all_of(cols_base_P), 
         all_of(cols_after_N), all_of(cols_after_P))

# 整理数据用于绘图
rownames(plot_data) <- plot_data$Gene.Symbol
plot_data <- plot_data[, -1] # 去除 Gene 列
plot_mat <- log1p(plot_data) # Log 转换

# 列注释 (分组信息)
annotation_col <- data.frame(
  Group = factor(c(rep("Base_N", 3), rep("Base_P", 3), 
                   rep("After_N", 3), rep("After_P", 3)),
                 levels = c("Base_N", "Base_P", "After_N", "After_P"))
)
rownames(annotation_col) <- colnames(plot_mat)

# 绘图
pheatmap(plot_mat,
         annotation_col = annotation_col,
         cluster_cols = FALSE, # 保持分组顺序
         cluster_rows = TRUE,
         scale = "row",        # 按行归一化，突出表达趋势
         show_colnames = FALSE,
         main = "Expression Dynamics of Key Regulatory Genes",
         color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
         filename = "Fig_Heatmap_KeyGenes.pdf")

# ==============================================================================
# 5. 导出关键数据表
# ==============================================================================
# 将分析结果合并导出，方便您在论文正文中引用数值
final_table <- res_mechanism %>% 
  select(Gene, LogFC_After_PvsN = LogFC, Pval_After_PvsN = PValue) %>%
  left_join(res_induction %>% select(Gene, LogFC_Induction = LogFC, Pval_Induction = PValue), by = "Gene") %>%
  left_join(res_identity %>% select(Gene, LogFC_Identity = LogFC), by = "Gene")

write.csv(final_table, "Table_S1_Detailed_Stats.csv", row.names = FALSE)




###chatgpt###

#第一步：读取数据，识别样本列与分组

library(tidyverse)

asc <- read.csv("20260126 bulk ASC.csv", check.names = FALSE)

# 1) 指定列
id_col   <- "Gene ID"
sym_col  <- "Gene Symbol"

stopifnot(all(c(id_col, sym_col) %in% colnames(asc)))

# 2) 清理 Gene Symbol 的单引号与空白
asc <- asc %>%
  mutate(
    !!sym_col := str_replace_all(.data[[sym_col]], "^'+|'+$", ""),  # 去掉首尾单引号（可能多个）
    !!sym_col := str_replace_all(.data[[sym_col]], "\"", ""),      # 防御性：去掉双引号
    !!sym_col := str_trim(.data[[sym_col]])
  )

# 3) 处理 Gene Symbol 缺失：用 Gene ID 兜底
asc <- asc %>%
  mutate(
    GeneKey = ifelse(is.na(.data[[sym_col]]) | .data[[sym_col]] == "",
                     .data[[id_col]],
                     .data[[sym_col]])
  )

# 4) 样本列：排除 Gene ID / Gene Symbol / GeneKey
sample_cols <- setdiff(colnames(asc), c(id_col, sym_col, "GeneKey"))

# 5) 构建表达矩阵（normalized read count）
mat <- as.matrix(asc[, sample_cols])
rownames(mat) <- asc$GeneKey

# 6) Gene Symbol 去重：同名基因取中位数（或均值）
#    推荐中位数，抗异常值；你也可改成 mean
expr <- log2(mat + 1)

expr_df <- as.data.frame(expr) %>%
  rownames_to_column("GeneKey") %>%
  left_join(asc %>% select(all_of(id_col), all_of(sym_col), GeneKey) %>% distinct(),
            by = "GeneKey") %>%
  mutate(GeneSymbol_clean = ifelse(!is.na(.data[[sym_col]]) & .data[[sym_col]] != "",
                                   .data[[sym_col]], GeneKey))

expr_gene <- expr_df %>%
  select(GeneSymbol_clean, all_of(sample_cols)) %>%
  group_by(GeneSymbol_clean) %>%
  summarise(across(all_of(sample_cols), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
  rename(Gene = GeneSymbol_clean)

expr_mat <- as.matrix(expr_gene[, sample_cols])
rownames(expr_mat) <- expr_gene$Gene

dim(expr_mat)
expr_mat[1:3, 1:6]


#第二步：清理列名 + 构建 sample_info（适配你真实列名）

library(tidyverse)

# 1) 清理列名：去掉统一前缀
colnames(expr_mat) <- sub("^Normalized read count\\s+", "", colnames(expr_mat))

head(colnames(expr_mat), 10)

# 2) 生成 sample_info
sample_info <- data.frame(sample = colnames(expr_mat), stringsAsFactors = FALSE) %>%
  mutate(
    timepoint = ifelse(str_detect(sample, "^after_"), "After", "Before"),
    subtype   = case_when(
      str_detect(sample, "ASC_P_") ~ "P",
      str_detect(sample, "ASC_N_") ~ "N",
      TRUE ~ NA_character_
    ),
    group = paste0(timepoint, "_ASC_", subtype)  # 方便做颜色/注释
  )

print(sample_info)
stopifnot(!any(is.na(sample_info$subtype)))

sample_info$timepoint <- factor(sample_info$timepoint, levels = c("Before","After"))
sample_info$subtype   <- factor(sample_info$subtype,   levels = c("N","P"))
sample_info$group     <- factor(sample_info$group)

table(sample_info$subtype, sample_info$timepoint)


# 第三步 PCA 快速看 batch/timepoint 是否主导（可选但建议）

library(ggplot2)

pca <- prcomp(t(expr_mat), scale. = FALSE)
pca_df <- data.frame(pca$x[, 1:2], sample_info)

p=ggplot(pca_df, aes(PC1, PC2, color = subtype, shape = timepoint, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  theme_bw() +
  ggtitle("ASC bulk RNA-seq PCA (log2 normalized counts)")


print(p)

ggsave("ASC_PCA.pdf", plot = p, width = 10, height = 8, dpi = 300)



# 第4步：limma 交互模型（规避 before/after 批次效应的核心）

library(limma)
library(tibble)
library(dplyr)

design <- model.matrix(~ 0 + subtype * timepoint, data = sample_info)
colnames(design) <- make.names(colnames(design))
design
colnames(design)

fit <- lmFit(expr_mat, design)

# 交互项模型下的核心 contrasts
cont <- makeContrasts(
  # 共培养前：P vs N（Before 的固有差异）
  Pre_P_vs_N  = subtypeP - subtypeN,
  
  # 共培养后：P vs N（After 的差异）
  # After 时：P = subtypeP + subtypeP.timepointAfter；N = subtypeN + 0
  Post_P_vs_N = (subtypeP + subtypeP.timepointAfter) - subtypeN,
  
  # DID：共培养诱导变化是否 P 特异（最关键机制池）
  DID_P_specific_change = subtypeP.timepointAfter,
  
  # 各自 after-before（可选）
  P_After_vs_Before = timepointAfter + subtypeP.timepointAfter,
  N_After_vs_Before = timepointAfter,
  
  levels = design
)

fit2 <- eBayes(contrasts.fit(fit, cont))

get_res <- function(coef_name){
  topTable(fit2, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("Gene") %>%
    mutate(FDR = p.adjust(P.Value, "BH"))
}

res_pre  <- get_res("Pre_P_vs_N")
res_post <- get_res("Post_P_vs_N")
res_did  <- get_res("DID_P_specific_change")

# 保存
write.csv(res_pre,  "ASC_RNA_Pre_P_vs_N.csv",  row.names = FALSE)
write.csv(res_post, "ASC_RNA_Post_P_vs_N.csv", row.names = FALSE)
write.csv(res_did,  "ASC_RNA_DID_P_specific_change.csv", row.names = FALSE)

# 先看 DID top hits（机制候选库）
nrow(res_did)
head(res_did[, c("Gene","logFC","FDR")], 20)


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


# 第一步：把 ASC 的三个对比都跑齐并落表（Pre、Post、DID）
sig_cut_fdr  <- 0.05
sig_cut_lfc  <- 0.5

write_sig <- function(res, prefix){
  res_sig <- res %>%
    mutate(sig = ifelse(FDR < sig_cut_fdr & abs(logFC) >= sig_cut_lfc, "sig", "ns"))
  write.csv(res_sig, paste0(prefix, "_all.csv"), row.names = FALSE)
  write.csv(filter(res_sig, sig=="sig"), paste0(prefix, "_sig.csv"), row.names = FALSE)
  invisible(res_sig)
}

res_pre2  <- write_sig(res_pre,  "ASC_RNA_Pre_P_vs_N")
res_post2 <- write_sig(res_post, "ASC_RNA_Post_P_vs_N")
res_did2  <- write_sig(res_did,  "ASC_RNA_DID_P_specific_change")


sig_cut_fdr  <- 0.05
sig_cut_lfc  <- 0.5

write_sig <- function(res, prefix){
  res_sig <- res %>%
    mutate(sig = ifelse(FDR < sig_cut_fdr & abs(logFC) >= sig_cut_lfc, "sig", "ns"))
  write.csv(res_sig, paste0(prefix, "_all.csv"), row.names = FALSE)
  write.csv(filter(res_sig, sig=="sig"), paste0(prefix, "_sig.csv"), row.names = FALSE)
  invisible(res_sig)
}

res_pre2  <- write_sig(res_pre,  "ASC_RNA_Pre_P_vs_N")
res_post2 <- write_sig(res_post, "ASC_RNA_Post_P_vs_N")
res_did2  <- write_sig(res_did,  "ASC_RNA_DID_P_specific_change")

# 第二步：主文级别的 3 张图（ASC-only）
# 1. PCA（你已完成 PDF：ASC_PCA.pdf）
# 2. DID 火山图（最关键，强调“P-specific response”）
# 3. “机制候选基因模块热图”（你已生成 ASC_DID_topgenes_heatmap.pdf，但建议再补一张更聚焦的“候选机制基因热图”，用于主文或补图）


library(ggrepel)
library(ggplot2)
library(dplyr)

top_lab <- res_did %>%
  arrange(FDR) %>%
  slice(1:15) %>%
  pull(Gene)

dfv <- res_did %>%
  mutate(sig = case_when(
    FDR < 0.05 & logFC >=  0.5 ~ "Up",
    FDR < 0.05 & logFC <= -0.5 ~ "Down",
    TRUE ~ "ns"
  ))

p_vol <- ggplot(dfv, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = sig), alpha = 0.75, size = 1.6) +
  scale_color_manual(values = c("Up"="#D95F02","Down"="#1B9E77","ns"="grey70")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = dfv %>% filter(Gene %in% top_lab),
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_bw() +
  ggtitle("ASC DID (P-specific coculture response): differential genes") +
  xlab("log2 fold change (interaction term)") +
  ylab(expression(-log[10](FDR))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

p_vol

ggsave("ASC_DID_volcano.pdf", p_vol, width = 15, height = 10)





library(ggrepel)
library(ggplot2)
library(dplyr)

# 1. 准备 Top 15 基因 (保持原逻辑)
top_lab <- res_did %>%
  arrange(FDR) %>%
  slice(1:15) %>%
  pull(Gene)

# 2. 【关键修改】定义您想强制标记的基因 (手动添加 Nt5e 等)
my_interest_genes <- c("Nt5e", "Jag1", "Eng", "Has2", "Alcam","Emcn","Pecam1")

# 3. 合并标签列表 (Top 15 + 您感兴趣的基因，并去重)
final_labels <- unique(c(top_lab, my_interest_genes))

# 4. 准备绘图数据 (保持原逻辑)
dfv <- res_did %>%
  mutate(sig = case_when(
    FDR < 0.05 & logFC >=  0.5 ~ "Up",
    FDR < 0.05 & logFC <= -0.5 ~ "Down",
    TRUE ~ "ns"
  ))

# 5. 绘图
p_vol <- ggplot(dfv, aes(x = logFC, y = -log10(FDR))) +
  # (1) 基础点图层
  geom_point(aes(color = sig), alpha = 0.6, size = 1.6) +
  scale_color_manual(values = c("Up"="#D95F02", "Down"="#1B9E77", "ns"="grey85")) +
  
  # (2) 【新增】高亮重点基因的点 (加黑边、填黄色，醒目)
  geom_point(data = dfv %>% filter(Gene %in% my_interest_genes),
             color = "black",   # 边框黑色
             fill = "yellow",   # 填充黄色
             shape = 21,        # 空心圆模式(支持fill)
             size = 3.5,        # 稍微大一点
             stroke = 1) +      # 边框加粗
  
  # (3) 辅助线
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  
  # (4) 【修改】标签层 (使用合并后的 final_labels)
  ggrepel::geom_text_repel(
    data = dfv %>% filter(Gene %in% final_labels), # 这里改成了 final_labels
    aes(label = Gene),
    size = 4,                # 字稍微大一点
    fontface = "bold",       # 加粗
    box.padding = 0.5,       # 文本框边距
    point.padding = 0.5,     # 点边距
    max.overlaps = Inf,      # 允许重叠无限大(强制显示)
    min.segment.length = 0   # 总是显示连线
  ) +
  
  # (5) 主题修饰
  theme_bw() +
  ggtitle("ASC DID (P-specific coculture response): differential genes") +
  xlab("log2 fold change (interaction term)") +
  ylab(expression(-log[10](FDR))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        panel.grid.major = element_blank(), # 去除网格线更清爽
        panel.grid.minor = element_blank())

# 6. 显示并保存
print(p_vol)

ggsave("ASC_DID_volcano_highlight.3.pdf", p_vol, width = 10, height = 8)







# 第三步：把 Hallmark DID 机制子集做成 Supplement 级别图 + 表

df_top15 <- cam_asc_did %>%
  arrange(FDR) %>%
  slice(1:15) %>%
  mutate(score = ifelse(Direction=="Up", 1, -1) * (-log10(FDR)),
         Pathway = factor(Pathway, levels = rev(Pathway)))

p_top15 <- ggplot(df_top15, aes(x = score, y = Pathway, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  xlab(expression(paste("Signed ", -log[10], "(FDR)"))) +
  ylab("") +
  ggtitle("ASC DID: Top Hallmark programs") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom")
p_top15

ggsave("ASC_DID_Hallmark_top15.pdf", p_top15, width = 6.0, height = 4.6)


# Step 4. Robustness: adjust DID for a cell-cycle score (no new experiments)

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(msigdbr)
  library(tibble)
})

## -------------------------
## 0) 前置：你应当已经有这些对象
## expr_mat: 行=GeneSymbol，列=样本（ASC_N_1... after_ASC_P_3），数值=log2(normalized+1) 或类似
## sample_info: 包含 sample / subtype(N,P) / timepoint(Before,After)
## design: 你之前的 model.matrix(~0 + subtype*timepoint)（如果没有也没关系，下面会重建）
## -------------------------

## 1) 取 Hallmark 里的 cell-cycle gene sets（Mus musculus）
m_h <- msigdbr(species = "Mus musculus", category = "H")

gs_g2m <- m_h %>% filter(gs_name == "HALLMARK_G2M_CHECKPOINT") %>% pull(gene_symbol) %>% unique()
gs_e2f <- m_h %>% filter(gs_name == "HALLMARK_E2F_TARGETS") %>% pull(gene_symbol) %>% unique()

cc_genes <- unique(c(gs_g2m, gs_e2f))
cc_genes <- intersect(cc_genes, rownames(expr_mat))

cat("Cell-cycle genes in expr_mat:", length(cc_genes), "\n")
stopifnot(length(cc_genes) >= 10)  # 少于10说明映射/基因名有问题

## 2) 计算 cell-cycle score（建议用：先按基因做 row-zscore，再对基因取均值）
zmat <- t(scale(t(expr_mat[cc_genes, , drop = FALSE])))   # row-wise z-score
cc_score <- colMeans(zmat, na.rm = TRUE)

df_cc <- data.frame(sample = names(cc_score), cc_score = as.numeric(cc_score)) %>%
  left_join(sample_info, by = "sample")

print(df_cc %>% arrange(desc(cc_score)))

## 3) 重建设计矩阵：加入 cc_score 作为协变量（主效应），仍保留 DID 交互项
sample_info2 <- sample_info %>%
  left_join(df_cc %>% select(sample, cc_score), by = "sample") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("Before","After")),
    subtype   = factor(subtype,   levels = c("N","P"))
  )

design_cc <- model.matrix(~ 0 + subtype * timepoint + cc_score, data = sample_info2)
colnames(design_cc) <- make.names(colnames(design_cc))
print(colnames(design_cc))

## 4) contrasts：与原来一致（但 DID 项现在仍然是 subtypeP.timepointAfter）
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

## 5) 对比：原始 DID vs cc-adjusted DID（整体相关 + 重点候选基因）
# 如果你原始 res_did 不在环境里，就从你之前导出的 CSV 读入：
if(!exists("res_did")){
  res_did <- read.csv("ASC_RNA_DID_P_specific_change_all.csv", check.names = FALSE)
}

cmp <- res_did %>%
  select(Gene, logFC_raw = logFC, FDR_raw = FDR) %>%
  inner_join(res_did_cc %>% select(Gene, logFC_cc = logFC, FDR_cc = FDR),
             by = "Gene")

cat("Correlation of DID logFC (raw vs cc-adjusted): ",
    cor(cmp$logFC_raw, cmp$logFC_cc, use = "pairwise.complete.obs"), "\n")

# 你关心的“机制候选基因”清单（可按需增删）
genes_mech <- c(
  "Prelp","Thbs1","Thbs2","Timp3","Has2","Cemip","Npnt","Itga11",
  "Alcam","Lrrc32","Jag1","Bmp2","Ltbp2","Igfbp3","Angptl2","Slc7a5","Id1","Klf4","Serpina3n","Cd34","Ptprb"
)

cmp_mech <- cmp %>%
  filter(Gene %in% genes_mech) %>%
  arrange(FDR_cc)

print(cmp_mech)

write.csv(cmp_mech, "ASC_DID_keyGenes_raw_vs_ccAdjusted.csv", row.names = FALSE)

## 6) 重新做 Hallmark CAMERA（cc-adjusted 版本）
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

library(tidyverse)
library(ggplot2)

# 读入 ccAdjusted CAMERA 结果
cc_path <- "ASC_DID_CAMERA_Hallmark_ccAdjusted.csv"
cc <- read.csv(cc_path, check.names = FALSE)

stopifnot(all(c("Pathway","Direction","PValue","FDR") %in% colnames(cc)))

cc2 <- cc %>%
  mutate(
    Direction = factor(Direction, levels = c("Down","Up")),
    score = ifelse(Direction == "Up", 1, -1) * (-log10(FDR + 1e-300))
  )

# -------------------------
# S5D 主图：Top 15 Hallmarks
# -------------------------
top_n <- 15
df_top <- cc2 %>%
  arrange(FDR) %>%
  slice(1:min(top_n, n())) %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))

pdf("Figure_S5D_Hallmark_CAMERA_ccAdjusted_top15.pdf", width = 9, height = 4.8)
print(
  ggplot(df_top, aes(x = score, y = Pathway, fill = Direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text(
      aes(label = sprintf("FDR=%.2g", FDR),
          hjust = ifelse(score >= 0, -0.05, 1.05)),
      size = 2.7
    ) +
    scale_fill_manual(
      values = c("Up" = "#D95F02", "Down" = "#1B9E77"),
      name = NULL,
      labels = c("Down" = "Down in DID", "Up" = "Up in DID")
    ) +
    xlab("Signed -log10(FDR)") +
    ylab("") +
    ggtitle("ASC DID Hallmark CAMERA (cell-cycle adjusted) - Top 15") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    ) +
    coord_cartesian(
      xlim = c(min(df_top$score) - 0.6, max(df_top$score) + 0.6)
    )
)
dev.off()

# -------------------------
# S5D 主图：ALL Hallmarks
# -------------------------
top_n <- 50
df_top <- cc2 %>%
  arrange(FDR) %>%
  slice(1:min(top_n, n())) %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))

pdf("Figure_S5D_Hallmark_CAMERA_ccAdjusted_top50.pdf", width = 9, height = 16)
print(
  ggplot(df_top, aes(x = score, y = Pathway, fill = Direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text(
      aes(label = sprintf("FDR=%.2g", FDR),
          hjust = ifelse(score >= 0, -0.05, 1.05)),
      size = 2.7
    ) +
    scale_fill_manual(
      values = c("Up" = "#D95F02", "Down" = "#1B9E77"),
      name = NULL,
      labels = c("Down" = "Down in DID", "Up" = "Up in DID")
    ) +
    xlab("Signed -log10(FDR)") +
    ylab("") +
    ggtitle("ASC DID Hallmark CAMERA (cell-cycle adjusted)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    ) +
    coord_cartesian(
      xlim = c(min(df_top$score) - 0.6, max(df_top$score) + 0.6)
    )
)
dev.off()



# ---------------------------------------
# 可选：机制聚焦（IFNα/IFNγ/ROS）小图
# ---------------------------------------
mech_paths <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")

df_mech <- cc2 %>%
  filter(Pathway %in% mech_paths) %>%
  mutate(Pathway = factor(Pathway, levels = rev(mech_paths)))

pdf("Figure_S5D_Hallmark_CAMERA_ccAdjusted_IFN_ROS.pdf", width = 6.6, height = 3.2)
print(
  ggplot(df_mech, aes(x = score, y = Pathway, fill = Direction)) +
    geom_col(width = 0.70) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text(
      aes(label = sprintf("FDR=%.2g", FDR),
          hjust = ifelse(score >= 0, -0.06, 1.06)),
      size = 3.0
    ) +
    scale_fill_manual(values = c("Up" = "#D95F02", "Down" = "#1B9E77"), name = NULL) +
    xlab("Signed -log10(FDR)") +
    ylab("") +
    ggtitle("Mechanism-focused Hallmarks (cell-cycle adjusted)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 9),
      legend.position = "none"
    )
)
dev.off()


## 7) 生成一个“raw vs cc-adjusted logFC”散点图（审稿人友好）
library(ggplot2)

p_scatter <- ggplot(cmp, aes(x = logFC_raw, y = logFC_cc)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  ggtitle("ASC DID logFC: raw vs cell-cycle-adjusted") +
  xlab("DID logFC (raw)") +
  ylab("DID logFC (cc-adjusted)") +
  theme(plot.title = element_text(hjust = 0.5))
p_scatter
ggsave("ASC_DID_raw_vs_ccAdjusted_scatter.pdf", p_scatter, width = 4.8, height = 4.2)

## 8) 可选：只标注关键基因（在散点图上写字）
library(ggrepel)

p_scatter_lab <- p_scatter +
  ggrepel::geom_text_repel(
    data = cmp_mech,
    aes(label = Gene),
    color = "#FF6600",          # ✅ 标记文字红色
    size = 3,
    max.overlaps = Inf
  )
p_scatter_lab
ggsave("ASC_DID_raw_vs_ccAdjusted_scatter_labeled.pdf", p_scatter_lab, width = 5.6, height = 4.6)



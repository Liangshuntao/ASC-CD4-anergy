# SOURCE BEGIN 1241:1277
cd4_gsva=as.matrix(cd4_gsva)

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

# SOURCE END 1241:1277

# SOURCE BEGIN 1279:1426
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

# SOURCE END 1279:1426

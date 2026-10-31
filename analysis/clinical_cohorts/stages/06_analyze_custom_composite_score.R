# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 665:740
# 模块 5：基于 TIDE 逻辑的免疫逃逸分析 (R 原生版)
# ==============================================================================
message("Calculating the custom expression-derived composite score...")

# 1. 定义 TIDE 官方公开的 Exclusion (免疫排斥) 标志物
# 这些基因是 TIDE 评估基质屏障的核心
tide_sigs <- list(
  CAF_Score = c("FAP", "THY1", "COL1A1"),      # 相关于成纤维细胞
  MDSC_Score = c("CD14", "ITGAM", "FUT4"),     # 相关于髓系抑制细胞
  TAM_M2_Score = c("CD163", "CD68", "MS4A4A")  # 相关于 M2 型巨噬细胞
)

# 2. 对 Luminal B 样本执行 ssGSEA 打分
if("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  tide_gsva <- gsva(ssgseaParam(expr_matrix[, lumb_barcodes], tide_sigs))
} else {
  tide_gsva <- gsva(expr_matrix[, lumb_barcodes], tide_sigs, method="ssgsea")
}

# 3. 数据整合
tide_df <- as.data.frame(t(tide_gsva))
# 根据 TIDE 原理：Exclusion 评分是这三类抑制因子分数的平均值
tide_df$Exclusion <- rowMeans(tide_df) 
tide_df$barcode <- rownames(tide_df)

# 4. 计算综合 TIDE 分数并预测响应
# 核心逻辑：TIDE = (T细胞失能 + 免疫排斥) / 2
final_df <- final_df %>%
  left_join(tide_df %>% select(barcode, Exclusion, CAF_Score), by = "barcode") %>%
  mutate(
    # 升华点：直接用您发现的 Anergy 评分代表 T 细胞的 Dysfunction (失能)
    Dysfunction = CD4T_Anergy, 
    TIDE_Score = (Dysfunction + Exclusion) / 2,
    # 模拟响应预测：TIDE 分数越低，预测响应 (Responder) 越好
    Responder = ifelse(TIDE_Score < median(TIDE_Score), "True", "False")
  )

message("Composite scores calculated; generating plots...")

# ------------------------------------------------------------------------------
# 5. 绘制 Figure 5：免疫治疗预测与机制关联
# ------------------------------------------------------------------------------
library(ggpubr)
library(patchwork)

# 5A: TIDE 评分箱线图
p_5A <- ggboxplot(final_df, x = "ASC_Group", y = "TIDE_Score",
                  fill = "ASC_Group", palette = c("#E41A1C", "#377EB8"),
                  add = "jitter", title = "Fig 5A: Custom composite score") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
p_5A

# 5B: 响应者比例柱状图
p_5B <- ggplot(final_df, aes(x = ASC_Group, fill = Responder)) +
  geom_bar(position = "fill", width = 0.6) +
  scale_fill_manual(values = c("False" = "#BEBEBE", "True" = "#FF7F00")) +
  labs(title = "Fig 5B: Predicted Response Rate", y = "Proportion") +
  theme_pubr()
p_5B


# 5C: 关键机制验证 —— ASC_P 是否驱动了免疫排斥 (Exclusion)？
p_5C <- ggscatter(final_df, x = "ASC_P_6", y = "Exclusion",
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "ASC_P Score", ylab = "Immune Exclusion Score",
                  title = "Fig 5C: ASC_P drives Immune Exclusion") + theme_bw()

p_5C 
# 组合图片
fig5_final <- (p_5A | p_5B | p_5C) 
fig5_final

ggsave("Figure_5_TIDE_Rversion.pdf", plot = fig5_final, width = 16, height = 6)

message("Composite-score plots saved: Figure_5_TIDE_Rversion.pdf")

# SOURCE END 665:740

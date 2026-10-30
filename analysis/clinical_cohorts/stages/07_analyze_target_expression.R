# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 744:794
# 模块 6：核心药物靶点与耐药特征分析 (Figure 6)
# ==============================================================================
message(">>> 正在分析 ASC_P 相关亚群的药物靶点图谱...")

library(ggpubr)
library(patchwork)
library(dplyr)
library(tidyr)

# 1. 定义乳腺癌关键药物及其对应的靶点基因
# 内分泌治疗(ESR1/PGR), CDK4/6抑制剂(CDK4/6), PI3K抑制剂(PIK3CA), 化疗(TOP2A/MKI67)
target_genes <- c("ESR1", "PGR", "CDK4", "CDK6", "PIK3CA", "TOP2A", "MKI67")

# 2. 提取并准备数据
# 确保只选取 Luminal B 样本
existing_genes <- intersect(target_genes, rownames(expr_matrix))
target_expr <- as.data.frame(t(expr_matrix[existing_genes, final_df$barcode]))
target_expr$barcode <- rownames(target_expr)

# 合并临床分组信息
plot_df <- left_join(final_df %>% select(barcode, ASC_Group, ASC_P_6), target_expr, by = "barcode")

# 3. 绘制 Figure 6A：核心靶点表达差异
# 这将展示哪些靶点在 ASC_High 组里是活跃的或失活的
p_6A_list <- list()
for(gene in existing_genes){
  p <- ggboxplot(plot_df, x = "ASC_Group", y = gene,
                 fill = "ASC_Group", palette = c("#E41A1C", "#377EB8"),
                 add = "jitter", title = paste0(gene, " Expression")) +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    theme(legend.position = "none")
  p_6A_list[[gene]] <- p
}

# 4. 绘制 Figure 6B：机制关联验证
# 验证 ASC_P 评分是否与 CDK4/6 这一乳腺癌核心靶点相关
p_6B <- ggscatter(plot_df, x = "ASC_P_6", y = "ESR1",
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  color = "#377EB8", alpha = 0.5,
                  xlab = "ASC_P Score", ylab = "ESR1 Expression",
                  title = "Fig 6D: ASC_P score and ESR1 expression") + 
  theme_bw()

# 5. 组合展示并保存
fig6_final <- (p_6A_list[["ESR1"]] | p_6A_list[["CDK4"]] | p_6A_list[["PIK3CA"]] | p_6B) 
fig6_final

ggsave("Figure_6_Target_Vulnerability_Final.pdf", plot = fig6_final, width = 21.33, height = 6)

message("🎉 Figure 6 已生成：Figure_6_Target_Vulnerability_Final.pdf")

# SOURCE END 744:794

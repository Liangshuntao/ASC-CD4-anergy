# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 744:794
# Selected target-expression profiles
# ==============================================================================
message("Generating target-expression profiles for ASC_P-associated groups...")

library(ggpubr)
library(patchwork)
library(dplyr)
library(tidyr)

# Define selected target and proliferation-marker genes
# ESR1/PGR, CDK4/CDK6, PIK3CA and TOP2A/MKI67
target_genes <- c("ESR1", "PGR", "CDK4", "CDK6", "PIK3CA", "TOP2A", "MKI67")

# Extract and prepare expression data
# Select Luminal B samples
existing_genes <- intersect(target_genes, rownames(expr_matrix))
target_expr <- as.data.frame(t(expr_matrix[existing_genes, final_df$barcode]))
target_expr$barcode <- rownames(target_expr)

# Join clinical group labels
plot_df <- left_join(final_df %>% select(barcode, ASC_Group, ASC_P_6), target_expr, by = "barcode")

# Compare selected target-gene expression between groups
# Display expression levels in the ASC score groups
p_6A_list <- list()
for(gene in existing_genes){
  p <- ggboxplot(plot_df, x = "ASC_Group", y = gene,
                 fill = "ASC_Group", palette = c("#E41A1C", "#377EB8"),
                 add = "jitter", title = paste0(gene, " Expression")) +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    theme(legend.position = "none")
  p_6A_list[[gene]] <- p
}

# Plot score-expression correlations
# Examine the association between ASC_P score and target expression
p_6B <- ggscatter(plot_df, x = "ASC_P_6", y = "ESR1",
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  color = "#377EB8", alpha = 0.5,
                  xlab = "ASC_P Score", ylab = "ESR1 Expression",
                  title = "Fig 6D: ASC_P score and ESR1 expression") + 
  theme_bw()

# Assemble and save the panels
fig6_final <- (p_6A_list[["ESR1"]] | p_6A_list[["CDK4"]] | p_6A_list[["PIK3CA"]] | p_6B) 
fig6_final

ggsave("Figure_6_Target_Vulnerability_Final.pdf", plot = fig6_final, width = 21.33, height = 6)

message("Target-expression plots saved: Figure_6_Target_Vulnerability_Final.pdf")

# SOURCE END 744:794

# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 365:427
# Transcriptomic and pathway analyses
# ==============================================================================
message("Generating pathway analysis plots...")

# GSEA comparing high and low score groups
group_info <- factor(final_df$ASC_Group, levels = c("ASC_Low", "ASC_High"))
design <- model.matrix(~ group_info)
fit <- lmFit(expr_matrix[, final_df$barcode], design) %>% eBayes()
diff_res <- topTable(fit, coef = 2, number = Inf)
diff_res$gene <- rownames(diff_res)
gene_list <- setNames(sign(diff_res$logFC) * -log10(diff_res$P.Value), diff_res$gene) %>% sort(decreasing = TRUE)

h_df <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name, gene_symbol)
gsea_res <- GSEA(gene_list, TERM2GENE = h_df, pvalueCutoff = 0.05, eps = 0)

p_3A <- dotplot(gsea_res, showCategory = 20,# Show the top 20 categories
                split = ".sign") + facet_grid(.~.sign) +
  scale_fill_gradientn(colours = c("#377EB8", "white", "#E41A1C")) +
  theme_bw() + labs(title = "Fig 3A: Hallmark GSEA (High vs Low)")

p_3A

ggsave("Fig 3A GSEA Dotplot.pdf", plot = p_3A, width = 10, height = 10)

# Pathway heatmap ordered by anergy score
h_list <- h_df %>% split(x = .$gene_symbol, f = .$gs_name)
pathway_ssgsea <- gsva(ssgseaParam(expr_matrix[, final_df$barcode], h_list))

selected_pathways <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                       "HALLMARK_TGF_BETA_SIGNALING",
                       "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
                       "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                       "HALLMARK_ANGIOGENESIS", 
                       "HALLMARK_HYPOXIA")

# Order the matrix and annotation rows
annotation_final <- final_df %>%
  select(barcode, ASC_Group, CD4T_Anergy, Treg_Infiltration, CD4T_Infiltration, OS.status) %>%
  arrange(ASC_Group, desc(CD4T_Anergy)) %>% as.data.frame()
rownames(annotation_final) <- annotation_final$barcode
annotation_final$barcode <- NULL
annotation_final$OS.status <- factor(annotation_final$OS.status, labels = c("Alive", "Dead"))

plot_mat <- pathway_ssgsea[selected_pathways, rownames(annotation_final)] %>% t() %>% scale() %>% t()
rownames(plot_mat) <- gsub("HALLMARK_", "", rownames(plot_mat))

ann_colors = list(
  ASC_Group = c(ASC_Low = "#377EB8", ASC_High = "#E41A1C"),
  OS.status = c(Alive = "#BEBEBE", Dead = "black"),
  CD4T_Anergy = colorRampPalette(c("blue",  "white","red"))(100), # Blue-white-red score scale
  Treg_Infiltration = colorRampPalette(c("white", "#984EA3"))(100)
)
ann_colors
pdf("Fig 3B Ultimate Heatmap.pdf", width = 8, height = 6)
pheatmap(plot_mat, annotation_col = annotation_final, annotation_colors = ann_colors,
         show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = TRUE,
         treeheight_row = 0,    # Hide the row dendrogram
         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
         main = "Fig 3B: Integrated Immune Landscape", gaps_col = table(annotation_final$ASC_Group)[1])
dev.off()
graphics.off()
message("Pathway analysis plots saved.")

# SOURCE END 365:427

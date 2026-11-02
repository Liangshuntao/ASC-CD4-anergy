# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 665:740
# Custom expression-derived immune-state composite; not the official TIDE implementation
# ==============================================================================
message("Calculating the custom expression-derived composite score...")

# Define the selected exclusion-related signatures
# Signature groups represent fibroblast, myeloid and macrophage programs
tide_sigs <- list(
  CAF_Score = c("FAP", "THY1", "COL1A1"),      # Fibroblast-related genes
  MDSC_Score = c("CD14", "ITGAM", "FUT4"),     # Myeloid suppressor-related genes
  TAM_M2_Score = c("CD163", "CD68", "MS4A4A")  # M2 macrophage-related genes
)

# Calculate ssGSEA scores in Luminal B samples
if("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  tide_gsva <- gsva(ssgseaParam(expr_matrix[, lumb_barcodes], tide_sigs))
} else {
  tide_gsva <- gsva(expr_matrix[, lumb_barcodes], tide_sigs, method="ssgsea")
}

# Combine the scores
tide_df <- as.data.frame(t(tide_gsva))
# Custom exclusion score: mean of the three signature scores
tide_df$Exclusion <- rowMeans(tide_df) 
tide_df$barcode <- rownames(tide_df)

# Calculate the custom composite and median-based labels
# Custom composite = (anergy score + exclusion score) / 2
final_df <- final_df %>%
  left_join(tide_df %>% select(barcode, Exclusion, CAF_Score), by = "barcode") %>%
  mutate(
    # Use CD4T_Anergy as the dysfunction-related component
    Dysfunction = CD4T_Anergy, 
    TIDE_Score = (Dysfunction + Exclusion) / 2,
    # Historical Responder labels indicate median-based score groups, not observed treatment responses
    Responder = ifelse(TIDE_Score < median(TIDE_Score), "True", "False")
  )

message("Composite scores calculated; generating plots...")

# ------------------------------------------------------------------------------
# Plot the custom composite and its associations
# ------------------------------------------------------------------------------
library(ggpubr)
library(patchwork)

# Composite-score boxplot
p_5A <- ggboxplot(final_df, x = "ASC_Group", y = "TIDE_Score",
                  fill = "ASC_Group", palette = c("#E41A1C", "#377EB8"),
                  add = "jitter", ylab = "Composite score", title = "Fig 5A: Custom composite score") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
p_5A

# Proportions of median-based score labels
p_5B <- ggplot(final_df, aes(x = ASC_Group, fill = Responder)) +
  geom_bar(position = "fill", width = 0.6) +
  scale_fill_manual(values = c("False" = "#BEBEBE", "True" = "#FF7F00"),
                    name = "Group", labels = c("False" = "High", "True" = "Low")) +
  labs(title = "Fig 5B: Median-based composite groups", y = "Proportion") +
  theme_pubr()
p_5B


# Association between ASC_P and the exclusion score
p_5C <- ggscatter(final_df, x = "ASC_P_6", y = "Exclusion",
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "ASC_P Score", ylab = "Custom exclusion score",
                  title = "Fig 5C: ASC_P and custom exclusion score") + theme_bw()

p_5C 
# Assemble the panels
fig5_final <- (p_5A | p_5B | p_5C) 
fig5_final

ggsave("Figure_5_TIDE_Rversion.pdf", plot = fig5_final, width = 16, height = 6)

message("Composite-score plots saved: Figure_5_TIDE_Rversion.pdf")

# SOURCE END 665:740

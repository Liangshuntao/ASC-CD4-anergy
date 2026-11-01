# SOURCE BEGIN 554:675
# Extract the selected GSVA results
df_asc_gsva_sel <- read.csv(asc_pathways_file, fileEncoding = "UTF-8-BOM")


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

# Plot the pathway heatmap

# Combine pathway matrices by column
asc_gsva_Hallmark_kegg <- rbind(asc_gsva_Hallmark, asc_gsva_kegg)

df_asc_gsva_sel[,1]


ASC_GSVA_selected =asc_gsva_Hallmark_kegg[df_asc_gsva_sel[,1],]


# Sample annotations: ASC_N and ASC_P
# Set condition order to ASC_N followed by ASC_P
sample_info_asc$condition <- factor(sample_info_asc$condition,
                                    levels = c("ASC_N", "ASC_P"))

# Order samples by condition
sample_info_asc <- sample_info_asc[order(sample_info_asc$condition), ]

# Align matrix columns to the sample order
ASC_GSVA_selected <- ASC_GSVA_selected[, sample_info_asc$sample]

# Construct column annotations
annot_col_asc <- data.frame(
  condition = sample_info_asc$condition
)
rownames(annot_col_asc) <- sample_info_asc$sample


# Specify the pathway row order
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

# Reorder rows using available pathway labels
row_order_use <- intersect(row_order, rownames(ASC_GSVA_selected))
ASC_GSVA_selected <- ASC_GSVA_selected[row_order_use, , drop = FALSE]



# Plot with column clustering disabled
# Disable row clustering to preserve the selected pathway order
pdf("ASC_GSVA_heatmap.pdf", width = 10, height = 4)
pheatmap(
  ASC_GSVA_selected,
  scale = "row",
  cluster_rows = FALSE,              # Preserve row_order
  cluster_cols = FALSE,              # Preserve the ordered sample columns
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

# SOURCE END 554:675

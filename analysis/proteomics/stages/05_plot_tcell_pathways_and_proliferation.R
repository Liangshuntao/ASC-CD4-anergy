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
# Open a PDF graphics device
pdf("CD4_GSVA_heatmap.pdf", width = 6, height = 4)
# Draw the plot on the active device
pheatmap(
  cd4_GSVA_selected,
  scale = "row",  # Standardize rows to display relative expression
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
# Close the device and finish writing the file
dev.off()

# SOURCE END 1241:1277

# SOURCE BEGIN 1279:1426
library(dplyr)
library(ggplot2)

# Select four proliferation-associated Hallmark pathways
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

# Order pathways by log fold change
df_cycle <- df_cycle %>%
  arrange(logFC) %>%               # Ascending order for the horizontal plot
  mutate(Pathway = factor(Pathway, levels = Pathway))

colnames(df_cycle)

write.csv(
  df_cycle,
  "GSVA_CD4_P_vs_CD4_N_selected_pathways.csv",
  row.names = FALSE
)

# Plot horizontal bars for CD4_P versus CD4_N
library(dplyr)
library(ggplot2)

# Define direction labels and ordering
df_cycle2 <- df_cycle %>%
  mutate(
    direction = ifelse(logFC > 0, "Up_in_CD4_P", "Down_in_CD4_P"),
    Pathway   = factor(Pathway, levels = Pathway[order(logFC)]),
    label     = sprintf("FDR=%.3g", FDR)   # Label FDR values
  )

# Draw the horizontal bar plot
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

# Selected proliferation and cell-cycle genes
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
# Match genes case-insensitively
# The matrix row names contain gene symbols
all_genes      <- rownames(expr_cd4_for_gsva)
all_genes_up   <- toupper(all_genes)
genes_prolif_up <- toupper(genes_prolif)

idx <- all_genes_up %in% genes_prolif_up

prolif_genes_in_mat <- all_genes[idx]
prolif_genes_in_mat
# Inspect the matched genes

mat_prolif <- expr_cd4_for_gsva[prolif_genes_in_mat, , drop = FALSE]

# Sample annotations: CD4_N and CD4_P
annot_col <- data.frame(
  condition = sample_info_cd4$condition
)
rownames(annot_col) <- sample_info_cd4$sample  # Match annotation row names to expression columns

# Verify sample ordering
mat_prolif <- mat_prolif[, rownames(annot_col), drop = FALSE]

# Draw the row-standardized heatmap
# Open a PDF graphics device
pdf("CD4_proliferation_heatmap.pdf", width = 6, height = 4)

# Draw on the active device
pheatmap(
  mat_prolif,
  scale = "row",  # Standardize rows to display relative expression
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
# Close the graphics device
dev.off()

# SOURCE END 1279:1426

# Figure 8D/E original plot sources; visual edits not reproduced
# Executed by run_bulk_analysis.R; not a standalone entry point.

# SOURCE BEGIN 710:762
library(tidyverse)
library(ggplot2)

# Read cell-cycle-adjusted CAMERA results
cc_path <- "ASC_DID_CAMERA_Hallmark_ccAdjusted.csv"
cc <- read.csv(cc_path, check.names = FALSE)

stopifnot(all(c("Pathway","Direction","PValue","FDR") %in% colnames(cc)))

cc2 <- cc %>%
  mutate(
    Direction = factor(Direction, levels = c("Down","Up")),
    score = ifelse(Direction == "Up", 1, -1) * (-log10(FDR + 1e-300))
  )

# -------------------------
# Plot the top 15 Hallmark pathways
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

# SOURCE END 710:762

# SOURCE BEGIN 840:867
# Compare primary and cell-cycle-adjusted log fold changes
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

# Label selected genes on the comparison plot
library(ggrepel)

p_scatter_lab <- p_scatter +
  ggrepel::geom_text_repel(
    data = cmp_mech,
    aes(label = Gene),
    color = "#FF6600",          # Highlight label color
    size = 3,
    max.overlaps = Inf
  )
p_scatter_lab
ggsave("ASC_DID_raw_vs_ccAdjusted_scatter_labeled.pdf", p_scatter_lab, width = 5.6, height = 4.6)

# SOURCE END 840:867

# Figure 8B highlighted volcano source; visual edits not reproduced
# Executed by run_bulk_analysis.R; not a standalone entry point.

# SOURCE BEGIN 491:558
library(ggrepel)
library(ggplot2)
library(dplyr)

# Select the top 15 genes
top_lab <- res_did %>%
  arrange(FDR) %>%
  slice(1:15) %>%
  pull(Gene)

# Specify additional genes to label
my_interest_genes <- c("Nt5e", "Jag1", "Eng", "Has2", "Alcam","Emcn","Pecam1")

# Combine and deduplicate the label lists
final_labels <- unique(c(top_lab, my_interest_genes))

# Prepare plotting data
dfv <- res_did %>%
  mutate(sig = case_when(
    FDR < 0.05 & logFC >=  0.5 ~ "Up",
    FDR < 0.05 & logFC <= -0.5 ~ "Down",
    TRUE ~ "ns"
  ))

# Construct the volcano plot
p_vol <- ggplot(dfv, aes(x = logFC, y = -log10(FDR))) +
  # Base point layer
  geom_point(aes(color = sig), alpha = 0.6, size = 1.6) +
  scale_color_manual(values = c("Up"="#D95F02", "Down"="#1B9E77", "ns"="grey85")) +
  
  # Highlight selected genes
  geom_point(data = dfv %>% filter(Gene %in% my_interest_genes),
             color = "black",   # Black border
             fill = "yellow",   # Yellow fill
             shape = 21,        # Point shape supporting separate fill and border
             size = 3.5,        # Point size
             stroke = 1) +      # Border width
  
  # Reference lines
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  
  # Label genes in final_labels
  ggrepel::geom_text_repel(
    data = dfv %>% filter(Gene %in% final_labels), # Use the combined label list
    aes(label = Gene),
    size = 4,                # Label size
    fontface = "bold",       # Bold labels
    box.padding = 0.5,       # Label-box padding
    point.padding = 0.5,     # Point padding
    max.overlaps = Inf,      # Retain all selected labels
    min.segment.length = 0   # Draw label connectors regardless of length
  ) +
  
  # Plot theme
  theme_bw() +
  ggtitle("ASC DID (P-specific coculture response): differential genes") +
  xlab("log2 fold change (interaction term)") +
  ylab(expression(-log[10](FDR))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank())

# Display and save the plot
print(p_vol)

ggsave("ASC_DID_volcano_highlight.3.pdf", p_vol, width = 10, height = 8)

# SOURCE END 491:558

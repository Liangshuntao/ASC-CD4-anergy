# PCA revision source; current input reproduces Figure 8A variance (89.26%, 4.42%).
# Figure 8A: PCA with variance explained on axis labels
# Matches original Figure 8A pipeline from 20260126 bulk ASC.R (lines 155-262)
# Adds PC1/PC2 variance percentage to axis labels per reviewer request

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# ---- Step 1: Load and prepare data (exactly as original script) ----
asc <- read.csv(input_file, check.names = FALSE)

id_col  <- "Gene ID"
sym_col <- "Gene Symbol"

stopifnot(all(c(id_col, sym_col) %in% colnames(asc)))

asc <- asc %>%
  mutate(
    !!sym_col := str_replace_all(.data[[sym_col]], "^'+|'+$", ""),
    !!sym_col := str_replace_all(.data[[sym_col]], "\"", ""),
    !!sym_col := str_trim(.data[[sym_col]])
  )

asc <- asc %>%
  mutate(
    GeneKey = ifelse(is.na(.data[[sym_col]]) | .data[[sym_col]] == "",
                     .data[[id_col]],
                     .data[[sym_col]])
  )

sample_cols <- setdiff(colnames(asc), c(id_col, sym_col, "GeneKey"))
mat <- as.matrix(asc[, sample_cols])
rownames(mat) <- asc$GeneKey

expr <- log2(mat + 1)

expr_df <- as.data.frame(expr) %>%
  rownames_to_column("GeneKey") %>%
  left_join(asc %>% select(all_of(id_col), all_of(sym_col), GeneKey) %>% distinct(),
            by = "GeneKey") %>%
  mutate(GeneSymbol_clean = ifelse(!is.na(.data[[sym_col]]) & .data[[sym_col]] != "",
                                   .data[[sym_col]], GeneKey))

expr_gene <- expr_df %>%
  select(GeneSymbol_clean, all_of(sample_cols)) %>%
  group_by(GeneSymbol_clean) %>%
  summarise(across(all_of(sample_cols), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
  rename(Gene = GeneSymbol_clean)

expr_mat <- as.matrix(expr_gene[, sample_cols])
rownames(expr_mat) <- expr_gene$Gene

# ---- Step 2: Clean column names and build sample_info (exactly as original) ----
colnames(expr_mat) <- sub("^Normalized read count\\s+", "", colnames(expr_mat))

sample_info <- data.frame(sample = colnames(expr_mat), stringsAsFactors = FALSE) %>%
  mutate(
    timepoint = ifelse(str_detect(sample, "^after_"), "After", "Before"),
    subtype   = case_when(
      str_detect(sample, "ASC_P_") ~ "P",
      str_detect(sample, "ASC_N_") ~ "N",
      TRUE ~ NA_character_
    ),
    group = paste0(timepoint, "_ASC_", subtype)
  )

sample_info$timepoint <- factor(sample_info$timepoint, levels = c("Before", "After"))
sample_info$subtype   <- factor(sample_info$subtype,   levels = c("N", "P"))
sample_info$group     <- factor(sample_info$group)

# ---- Step 3: PCA (exactly as original) ----
pca <- prcomp(t(expr_mat), scale. = FALSE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 2)

cat(sprintf("PC1: %.2f%% variance\n", var_pct[1]))
cat(sprintf("PC2: %.2f%% variance\n", var_pct[2]))

pca_df <- data.frame(pca$x[, 1:2], sample_info)

# ---- Step 4: Plot with variance on axis labels ----
p <- ggplot(pca_df, aes(PC1, PC2, color = subtype, shape = timepoint, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  theme_bw() +
  ggtitle("ASC bulk RNA-seq PCA (log2 normalized counts)") +
  xlab(sprintf("PC1 (%.2f%% variance)", var_pct[1])) +
  ylab(sprintf("PC2 (%.2f%% variance)", var_pct[2])) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p)

# ---- Step 5: Save ----
output_dir <- getwd()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(output_dir, "Figure8A_PCA_with_variance.pdf"),
       plot = p, width = 10, height = 8, dpi = 300, device = cairo_pdf)

cat("Saved: Figure8A_PCA_with_variance.pdf\n")

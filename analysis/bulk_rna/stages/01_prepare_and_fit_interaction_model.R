# Final unadjusted limma/DID model and result exports; exploratory t-test excluded
# Executed by run_bulk_analysis.R; not a standalone entry point.

# SOURCE BEGIN 158:243
library(tidyverse)

asc <- read.csv(input_file, check.names = FALSE)

# Specify identifier columns
id_col   <- "Gene ID"
sym_col  <- "Gene Symbol"

stopifnot(all(c(id_col, sym_col) %in% colnames(asc)))

# Remove quotation marks and surrounding whitespace from gene symbols
asc <- asc %>%
  mutate(
    !!sym_col := str_replace_all(.data[[sym_col]], "^'+|'+$", ""),  # Remove leading and trailing single quotes
    !!sym_col := str_replace_all(.data[[sym_col]], "\"", ""),      # Remove double quotes
    !!sym_col := str_trim(.data[[sym_col]])
  )

# Use Gene ID when Gene Symbol is missing
asc <- asc %>%
  mutate(
    GeneKey = ifelse(is.na(.data[[sym_col]]) | .data[[sym_col]] == "",
                     .data[[id_col]],
                     .data[[sym_col]])
  )

# Identify sample columns by excluding annotation columns
sample_cols <- setdiff(colnames(asc), c(id_col, sym_col, "GeneKey"))

# Construct the normalized-expression matrix
mat <- as.matrix(asc[, sample_cols])
rownames(mat) <- asc$GeneKey

# Aggregate duplicate gene symbols by median after log2 transformation

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

dim(expr_mat)
expr_mat[1:3, 1:6]


# Clean sample names and construct sample_info

library(tidyverse)

# Remove the common sample-column prefix
colnames(expr_mat) <- sub("^Normalized read count\\s+", "", colnames(expr_mat))

head(colnames(expr_mat), 10)

# Construct sample metadata
sample_info <- data.frame(sample = colnames(expr_mat), stringsAsFactors = FALSE) %>%
  mutate(
    timepoint = ifelse(str_detect(sample, "^after_"), "After", "Before"),
    subtype   = case_when(
      str_detect(sample, "ASC_P_") ~ "P",
      str_detect(sample, "ASC_N_") ~ "N",
      TRUE ~ NA_character_
    ),
    group = paste0(timepoint, "_ASC_", subtype)  # Group labels for plot annotations
  )

print(sample_info)
stopifnot(!any(is.na(sample_info$subtype)))

sample_info$timepoint <- factor(sample_info$timepoint, levels = c("Before","After"))
sample_info$subtype   <- factor(sample_info$subtype,   levels = c("N","P"))
sample_info$group     <- factor(sample_info$group)

table(sample_info$subtype, sample_info$timepoint)

# SOURCE END 158:243

# SOURCE BEGIN 268:317
library(limma)
library(tibble)
library(dplyr)

design <- model.matrix(~ 0 + subtype * timepoint, data = sample_info)
colnames(design) <- make.names(colnames(design))
design
colnames(design)

fit <- lmFit(expr_mat, design)

# Contrasts from the subtype-by-timepoint interaction model
cont <- makeContrasts(
  # ASC_P versus ASC_N before coculture
  Pre_P_vs_N  = subtypeP - subtypeN,
  
  # ASC_P versus ASC_N after coculture
  # After: P = subtypeP + subtypeP.timepointAfter; N = subtypeN
  Post_P_vs_N = (subtypeP + subtypeP.timepointAfter) - subtypeN,
  
  # Difference-in-differences: subtype-biased coculture-associated change
  DID_P_specific_change = subtypeP.timepointAfter,
  
  # Within-subtype after-versus-before contrasts
  P_After_vs_Before = timepointAfter + subtypeP.timepointAfter,
  N_After_vs_Before = timepointAfter,
  
  levels = design
)

fit2 <- eBayes(contrasts.fit(fit, cont))

get_res <- function(coef_name){
  topTable(fit2, coef = coef_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("Gene") %>%
    mutate(FDR = p.adjust(P.Value, "BH"))
}

res_pre  <- get_res("Pre_P_vs_N")
res_post <- get_res("Post_P_vs_N")
res_did  <- get_res("DID_P_specific_change")

# Export result tables
write.csv(res_pre,  "ASC_RNA_Pre_P_vs_N.csv",  row.names = FALSE)
write.csv(res_post, "ASC_RNA_Post_P_vs_N.csv", row.names = FALSE)
write.csv(res_did,  "ASC_RNA_DID_P_specific_change.csv", row.names = FALSE)

# Inspect the leading DID results
nrow(res_did)
head(res_did[, c("Gene","logFC","FDR")], 20)

# SOURCE END 268:317

# SOURCE BEGIN 411:425
# Export complete and significance-filtered tables for the three primary contrasts
sig_cut_fdr  <- 0.05
sig_cut_lfc  <- 0.5

write_sig <- function(res, prefix){
  res_sig <- res %>%
    mutate(sig = ifelse(FDR < sig_cut_fdr & abs(logFC) >= sig_cut_lfc, "sig", "ns"))
  write.csv(res_sig, paste0(prefix, "_all.csv"), row.names = FALSE)
  write.csv(filter(res_sig, sig=="sig"), paste0(prefix, "_sig.csv"), row.names = FALSE)
  invisible(res_sig)
}

res_pre2  <- write_sig(res_pre,  "ASC_RNA_Pre_P_vs_N")
res_post2 <- write_sig(res_post, "ASC_RNA_Post_P_vs_N")
res_did2  <- write_sig(res_did,  "ASC_RNA_DID_P_specific_change")

# SOURCE END 411:425

# Final unadjusted limma/DID model and result exports; exploratory t-test excluded
# Executed by run_bulk_analysis.R; not a standalone entry point.

# SOURCE BEGIN 158:243
library(tidyverse)

asc <- read.csv(input_file, check.names = FALSE)

# 1) 指定列
id_col   <- "Gene ID"
sym_col  <- "Gene Symbol"

stopifnot(all(c(id_col, sym_col) %in% colnames(asc)))

# 2) 清理 Gene Symbol 的单引号与空白
asc <- asc %>%
  mutate(
    !!sym_col := str_replace_all(.data[[sym_col]], "^'+|'+$", ""),  # 去掉首尾单引号（可能多个）
    !!sym_col := str_replace_all(.data[[sym_col]], "\"", ""),      # 防御性：去掉双引号
    !!sym_col := str_trim(.data[[sym_col]])
  )

# 3) 处理 Gene Symbol 缺失：用 Gene ID 兜底
asc <- asc %>%
  mutate(
    GeneKey = ifelse(is.na(.data[[sym_col]]) | .data[[sym_col]] == "",
                     .data[[id_col]],
                     .data[[sym_col]])
  )

# 4) 样本列：排除 Gene ID / Gene Symbol / GeneKey
sample_cols <- setdiff(colnames(asc), c(id_col, sym_col, "GeneKey"))

# 5) 构建表达矩阵（normalized read count）
mat <- as.matrix(asc[, sample_cols])
rownames(mat) <- asc$GeneKey

# 6) Gene Symbol 去重：同名基因取中位数（或均值）
#    推荐中位数，抗异常值；你也可改成 mean
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


#第二步：清理列名 + 构建 sample_info（适配你真实列名）

library(tidyverse)

# 1) 清理列名：去掉统一前缀
colnames(expr_mat) <- sub("^Normalized read count\\s+", "", colnames(expr_mat))

head(colnames(expr_mat), 10)

# 2) 生成 sample_info
sample_info <- data.frame(sample = colnames(expr_mat), stringsAsFactors = FALSE) %>%
  mutate(
    timepoint = ifelse(str_detect(sample, "^after_"), "After", "Before"),
    subtype   = case_when(
      str_detect(sample, "ASC_P_") ~ "P",
      str_detect(sample, "ASC_N_") ~ "N",
      TRUE ~ NA_character_
    ),
    group = paste0(timepoint, "_ASC_", subtype)  # 方便做颜色/注释
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

# 交互项模型下的核心 contrasts
cont <- makeContrasts(
  # 共培养前：P vs N（Before 的固有差异）
  Pre_P_vs_N  = subtypeP - subtypeN,
  
  # 共培养后：P vs N（After 的差异）
  # After 时：P = subtypeP + subtypeP.timepointAfter；N = subtypeN + 0
  Post_P_vs_N = (subtypeP + subtypeP.timepointAfter) - subtypeN,
  
  # DID：共培养诱导变化是否 P 特异（最关键机制池）
  DID_P_specific_change = subtypeP.timepointAfter,
  
  # 各自 after-before（可选）
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

# 保存
write.csv(res_pre,  "ASC_RNA_Pre_P_vs_N.csv",  row.names = FALSE)
write.csv(res_post, "ASC_RNA_Post_P_vs_N.csv", row.names = FALSE)
write.csv(res_did,  "ASC_RNA_DID_P_specific_change.csv", row.names = FALSE)

# 先看 DID top hits（机制候选库）
nrow(res_did)
head(res_did[, c("Gene","logFC","FDR")], 20)

# SOURCE END 268:317

# SOURCE BEGIN 411:425
# 第一步：把 ASC 的三个对比都跑齐并落表（Pre、Post、DID）
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

# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 457:508
clin_patient <- read.delim(file.path(metabric_dir, "data_clinical_patient.txt"), skip = 4, header = TRUE, stringsAsFactors = FALSE)
clin_sample <- read.delim(file.path(metabric_dir, "data_clinical_sample.txt"), skip = 4, header = TRUE, stringsAsFactors = FALSE)

metabric_clin <- merge(clin_patient, clin_sample, by = "PATIENT_ID")
rownames(metabric_clin) <- metabric_clin$SAMPLE_ID
metabric_clin$PAM50_Subtype <- metabric_clin$CLAUDIN_SUBTYPE 

# ------------------------------------------------------------------------------
# 2. 本地读取表达矩阵 (修复连字符变小数点的问题)
# ------------------------------------------------------------------------------
message(">>> 正在读取微阵列芯片表达矩阵 (加上 check.names = FALSE 保持样本名原样)...")
# 【核心修复】：末尾加上 check.names = FALSE
expr_raw <- read.delim(file.path(metabric_dir, "data_mrna_illumina_microarray.txt"), 
                       header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# 清洗矩阵：去重、去空、设行名
expr_raw <- expr_raw[!duplicated(expr_raw$Hugo_Symbol) & expr_raw$Hugo_Symbol != "", ]
rownames(expr_raw) <- expr_raw$Hugo_Symbol
metabric_expr <- as.matrix(expr_raw[, -(1:2)])

# ------------------------------------------------------------------------------
# 3. 提取 Luminal B 并清洗生存数据
# ------------------------------------------------------------------------------
message(">>> 正在提取 Luminal B 亚型并对齐表达矩阵...")
metabric_lumb_clin <- metabric_clin %>%
  mutate(
    OS.time = as.numeric(OS_MONTHS),
    OS.status = ifelse(grepl("DECEASED", OS_STATUS, ignore.case = TRUE), 1, 0)
  ) %>%
  filter(grepl("LumB", PAM50_Subtype, ignore.case = TRUE) | grepl("Luminal B", PAM50_Subtype, ignore.case = TRUE)) %>%
  filter(!is.na(OS.time) & OS.time > 0)

# 这一次，交集就不会是 0 了！
valid_samples <- intersect(rownames(metabric_lumb_clin), colnames(metabric_expr))
message(paste(">>> 成功对齐样本数：", length(valid_samples), "例"))

metabric_lumb_clin <- metabric_lumb_clin[valid_samples, ]
expr_lumb_meta <- as.matrix(metabric_expr[, valid_samples])
# ------------------------------------------------------------------------------
# 4. ssGSEA 跨平台打分
# ------------------------------------------------------------------------------
message(">>> 正在执行芯片数据的 ssGSEA 打分...")
meta_genes <- rownames(expr_lumb_meta)
sig_asc_p_meta <- intersect(c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM"), meta_genes)
sig_cd4t_meta <- intersect(c("CD4", "CD3D", "CD3E"), meta_genes)
metabric_sigs <- list(ASC_P = sig_asc_p_meta, CD4T = sig_cd4t_meta)

if("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  meta_gsva <- gsva(ssgseaParam(expr_lumb_meta, metabric_sigs))
} else {
  meta_gsva <- gsva(expr_lumb_meta, metabric_sigs, method="ssgsea")
}

# SOURCE END 457:508

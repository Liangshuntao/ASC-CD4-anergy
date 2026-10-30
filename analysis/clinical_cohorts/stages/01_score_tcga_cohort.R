# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 44:88
clinical_base <- as.data.frame(colData(brca_data)) %>%
  mutate(patient_id = substr(barcode, 1, 12)) %>%
  filter(sample_type == "Primary Tumor") %>%
  left_join(brca_subtype, by = c("patient_id" = "patient")) %>%
  mutate(
    OS.time = ifelse(vital_status == "Dead", as.numeric(days_to_death), as.numeric(days_to_last_follow_up)),
    OS.status = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(OS.time) & OS.time > 0 & !is.na(BRCA_Subtype_PAM50))

# 表达矩阵准备
counts_tpm <- assay(brca_data, "tpm_unstrand")
gene_info <- rowData(brca_data)
valid_rows <- which(!is.na(gene_info$gene_name) & !duplicated(gene_info$gene_name))
expr_matrix <- as.matrix(counts_tpm[valid_rows, clinical_base$barcode])
rownames(expr_matrix) <- gene_info$gene_name[valid_rows]

# 定义 Signature 集 (包含 ASC_P, CD4T, Anergy, Treg)
sig_list_all <- list(
  ASC_P_6 = c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM"),
  CD4T_Infiltration = c("CD4", "CD3D", "CD3E"),
  CD4T_Anergy = c("RNF128", "CBLB", "EGR2", "DGKA", "ITCH", "IZUMO1R"),
  Treg_Infiltration = c("FOXP3", "IL2RA", "IKZF2", "CCR4"),
  CD4T_Exhaustion = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4")
)

# 执行多维度 GSVA 打分
if("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  gsva_res <- gsva(ssgseaParam(expr_matrix, sig_list_all))
} else {
  gsva_res <- gsva(expr_matrix, sig_list_all, method="ssgsea")
}

# 整合进 Luminal B 核心表
gsva_df <- as.data.frame(t(gsva_res))
gsva_df$barcode <- rownames(gsva_df)

final_df <- clinical_base %>%
  left_join(gsva_df, by = "barcode") %>%
  filter(BRCA_Subtype_PAM50 == "LumB") %>%
  mutate(
    ASC_Group = ifelse(ASC_P_6 > median(ASC_P_6), "ASC_High", "ASC_Low"),
    CD4_Group = ifelse(CD4T_Infiltration > median(CD4T_Infiltration), "CD4_High", "CD4_Low"),
    Combined_Group = paste0(ASC_Group, " / ", CD4_Group)
  )

# SOURCE END 44:88

# SOURCE BEGIN 171:206
# limma differential-protein analysis: ASC_P versus ASC_N
group_asc <- droplevels(sample_info_asc$condition)  # ASC_N, ASC_P

design_asc <- model.matrix(~ 0 + group_asc)
colnames(design_asc) <- levels(group_asc)
print(design_asc)

contrast_asc <- makeContrasts(
  ASC_P_vs_ASC_N = ASC_P - ASC_N,
  levels = design_asc
)

fit_asc <- lmFit(expr_imp_asc, design_asc)
fit_asc <- contrasts.fit(fit_asc, contrast_asc)
fit_asc <- eBayes(fit_asc)

res_ASC <- topTable(fit_asc, coef = "ASC_P_vs_ASC_N",
                    number = Inf, sort.by = "P")

res_ASC_anno <- res_ASC %>%
  rownames_to_column("Accession") %>%
  left_join(ann_asc, by = "Accession") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.585, "sig", "ns")
  )

write.csv(res_ASC_anno,
          "DE_ASC_P_vs_ASC_N_asc_only_all.csv",
          row.names = FALSE)

filter(res_ASC_anno, sig == "sig")

write.csv(filter(res_ASC_anno, sig == "sig"),
          "DE_ASC_P_vs_ASC_N_asc_only_sig.csv",
          row.names = FALSE)

# SOURCE END 171:206

# SOURCE BEGIN 943:975
group_cd4 <- droplevels(sample_info_cd4$condition)  # CD4_N, CD4_P

design_cd4 <- model.matrix(~ 0 + group_cd4)
colnames(design_cd4) <- levels(group_cd4)
print(design_cd4)

contrast_cd4 <- makeContrasts(
  CD4_P_vs_CD4_N = CD4_P - CD4_N,
  levels = design_cd4
)

fit_cd4 <- lmFit(expr_imp_cd4, design_cd4)
fit_cd4 <- contrasts.fit(fit_cd4, contrast_cd4)
fit_cd4 <- eBayes(fit_cd4)

res_CD4 <- topTable(fit_cd4, coef = "CD4_P_vs_CD4_N",
                    number = Inf, sort.by = "P")

res_CD4_anno <- res_CD4 %>%
  rownames_to_column("Accession") %>%
  left_join(ann_cd4, by = "Accession") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.585, "sig", "ns")
  )

write.csv(res_CD4_anno,
          "DE_CD4_P_vs_CD4_N_cd4_only_all.csv",
          row.names = FALSE)

write.csv(filter(res_CD4_anno, sig == "sig"),
          "DE_CD4_P_vs_CD4_N_cd4_only_sig.csv",
          row.names = FALSE)

# SOURCE END 943:975

# SOURCE BEGIN 448:484
# limma pathway-score analysis: ASC_P versus ASC_N

# Align columns with sample_info_asc
asc_gsva <- asc_gsva[, sample_info_asc$sample, drop = FALSE]
group_asc <- droplevels(sample_info_asc$condition)  # ASC_N / ASC_P

design_asc_gsva <- model.matrix(~ 0 + group_asc)
colnames(design_asc_gsva) <- levels(group_asc)

contrast_asc_gsva <- makeContrasts(
  ASC_P_vs_ASC_N = ASC_P - ASC_N,
  levels = design_asc_gsva
)

fit_asc_gsva <- lmFit(asc_gsva, design_asc_gsva)
fit_asc_gsva <- contrasts.fit(fit_asc_gsva, contrast_asc_gsva)
fit_asc_gsva <- eBayes(fit_asc_gsva)

res_ASC_gsva <- topTable(
  fit_asc_gsva,
  coef    = "ASC_P_vs_ASC_N",
  number  = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("Pathway") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.2, "sig", "ns")
  )

write.csv(res_ASC_gsva,
          "GSVA_ASC_P_vs_ASC_N_Hallmark_all.csv",
          row.names = FALSE)

write.csv(filter(res_ASC_gsva, sig == "sig"),
          "GSVA_ASC_P_vs_ASC_N_Hallmark_sig.csv",
          row.names = FALSE)

# SOURCE END 448:484

# SOURCE BEGIN 1183:1216
group_cd4 <- droplevels(sample_info_cd4$condition)  # CD4_N / CD4_P

design_cd4_gsva <- model.matrix(~ 0 + group_cd4)
colnames(design_cd4_gsva) <- levels(group_cd4)

contrast_cd4_gsva <- makeContrasts(
  CD4_P_vs_CD4_N = CD4_P - CD4_N,
  levels = design_cd4_gsva
)

fit_cd4_gsva <- lmFit(cd4_gsva, design_cd4_gsva)
fit_cd4_gsva <- contrasts.fit(fit_cd4_gsva, contrast_cd4_gsva)
fit_cd4_gsva <- eBayes(fit_cd4_gsva)

res_cd4_gsva <- topTable(
  fit_cd4_gsva,
  coef    = "CD4_P_vs_CD4_N",
  number  = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("Pathway") %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    sig = ifelse(FDR < 0.05 & abs(logFC) >= 0.2, "sig", "ns")
  )
res_cd4_gsva


write.csv(res_cd4_gsva,
          "GSVA_CD4_P_vs_CD4_N_Hallmark_all.csv",
          row.names = FALSE)
write.csv(filter(res_cd4_gsva, sig == "sig"),
          "GSVA_CD4_P_vs_CD4_N_Hallmark_sig.csv",
          row.names = FALSE)

# SOURCE END 1183:1216

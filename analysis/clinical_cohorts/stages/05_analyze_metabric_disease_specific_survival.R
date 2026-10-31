# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 576:660
# METABRIC disease-specific survival analysis in joint quartile groups
# ==============================================================================
message("Starting METABRIC disease-specific survival analysis...")

# 1. Prepare disease-specific survival variables
metabric_dss_clin <- metabric_clin %>%
  mutate(
    # Follow-up time in months
    DSS.time = as.numeric(OS_MONTHS),
    # Disease deaths are events; other nonmissing statuses are censored
    DSS.status = ifelse(VITAL_STATUS == "Died of Disease", 1, 0) 
  ) %>%
  filter(grepl("LumB", PAM50_Subtype, ignore.case = TRUE) | grepl("Luminal B", PAM50_Subtype, ignore.case = TRUE)) %>%
  filter(!is.na(DSS.time) & DSS.time > 0)

# 2. Align clinical samples and expression scores
valid_dss_samples <- intersect(rownames(metabric_dss_clin), colnames(metabric_expr))
metabric_dss_clin <- metabric_dss_clin[valid_dss_samples, ]

# Attach the previously computed ssGSEA scores
metabric_dss_clin$ASC_Score <- as.numeric(meta_gsva["ASC_P", valid_dss_samples])
metabric_dss_clin$CD4_Score <- as.numeric(meta_gsva["CD4T", valid_dss_samples])

# 3. Define the upper and lower score quartiles
q75_asc <- quantile(metabric_dss_clin$ASC_Score, 0.75, na.rm=TRUE)
q25_asc <- quantile(metabric_dss_clin$ASC_Score, 0.25, na.rm=TRUE)
q75_cd4 <- quantile(metabric_dss_clin$CD4_Score, 0.75, na.rm=TRUE)
q25_cd4 <- quantile(metabric_dss_clin$CD4_Score, 0.25, na.rm=TRUE)

metabric_extreme <- metabric_dss_clin %>%
  # Retain samples in the upper or lower quartile of both scores
  filter((ASC_Score >= q75_asc | ASC_Score <= q25_asc) & 
           (CD4_Score >= q75_cd4 | CD4_Score <= q25_cd4)) %>%
  mutate(
    ASC_Group = ifelse(ASC_Score >= q75_asc, "ASC_High", "ASC_Low"),
    CD4_Group = ifelse(CD4_Score >= q75_cd4, "CD4_High", "CD4_Low"),
    Combined_Group = paste0(ASC_Group, " / ", CD4_Group)
  )

# 4. Plot disease-specific survival curves
graphics.off()
fit_extreme <- survfit(Surv(DSS.time, DSS.status) ~ Combined_Group, data = metabric_extreme)

p_4A_Rescue <- ggsurvplot(fit_extreme, data = metabric_extreme, pval = TRUE, palette = "jco",
                          xlim = c(0, 180), break.time.by = 60, 
                          title = "Fig 4A Rescue: METABRIC DSS (Top/Bottom 25%)",
                          xlab = "Disease-Specific Survival (Months)", risk.table = TRUE,
                          ggtheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold")))

pdf("Fig_4A_METABRIC_Rescue_DSS.pdf", width = 12, height = 8)
print(p_4A_Rescue)
dev.off()

message("Survival curves saved: Fig_4A_METABRIC_Rescue_DSS.pdf")


# ==============================================================================
# Cox model for disease-specific survival in the joint quartile subset
# ==============================================================================
message("Fitting the METABRIC Cox model and generating the forest plot...")

library(survival)
library(survminer)

# 1. Use the joint quartile subset (metabric_extreme)
# Set ASC_Low / CD4_High as the reference group
metabric_extreme$Combined_Group <- factor(metabric_extreme$Combined_Group)
metabric_extreme$Combined_Group <- relevel(metabric_extreme$Combined_Group, ref = "ASC_Low / CD4_High")

# 2. Fit the Cox model using DSS.time and DSS.status
cox_extreme <- coxph(Surv(DSS.time, DSS.status) ~ Combined_Group, data = metabric_extreme)

# 3. Generate the forest plot
p_4B_Rescue <- ggforest(cox_extreme, data = metabric_extreme, 
                        main = "Fig 4B Rescue: METABRIC DSS Hazard Ratio (Top/Bottom 25%)")

# 4. Display and save the plot
print(p_4B_Rescue)

ggsave("Fig_4B_METABRIC_Rescue_DSS_Forest.pdf", 
       plot = p_4B_Rescue, 
       width = 10, 
       height = 7)

message("Forest plot saved: Fig_4B_METABRIC_Rescue_DSS_Forest.pdf")

# SOURCE END 576:660

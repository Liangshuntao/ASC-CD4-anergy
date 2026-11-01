# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 90:361
# ==============================================================================
# Subtype-stratified and joint-group survival analyses
# ==============================================================================
message("Generating clinical analysis plots...")

# ------------------------------------------------------------------------------
# Prepare the table containing all subtypes
# ------------------------------------------------------------------------------
all_subtypes_df <- clinical_base %>%
  left_join(gsva_df, by = "barcode")

# ------------------------------------------------------------------------------
# Kaplan-Meier analysis using upper and lower score quartiles
# ------------------------------------------------------------------------------
run_subgroup_km <- function(data, subtype_name) {
  # Select the requested subtype
  sub_df <- data %>% dplyr::filter(BRCA_Subtype_PAM50 == subtype_name)
  
  if(nrow(sub_df) < 20) {
    message(paste("Skipping", subtype_name, ": insufficient samples"))
    return(NULL)
  }
  
  # Calculate the 25th and 75th percentiles of ASC_P_6
  q_h <- quantile(sub_df$ASC_P_6, 0.75, na.rm = TRUE)
  q_l <- quantile(sub_df$ASC_P_6, 0.25, na.rm = TRUE)
  
  # Retain and label the upper and lower quartile groups
  sub_extreme <- sub_df %>%
    dplyr::filter(ASC_P_6 >= q_h | ASC_P_6 <= q_l) %>%
    mutate(Group = ifelse(ASC_P_6 >= q_h, "High", "Low"))
  
  # Check the number of death events
  if(sum(sub_extreme$OS.status) < 3) {
    message(paste("Skipping", subtype_name, ": insufficient death events"))
    return(NULL)
  }
  
  # Fit and plot survival curves
  fit <- survfit(Surv(OS.time, OS.status) ~ Group, data = sub_extreme)
  p <- ggsurvplot(fit, data = sub_extreme, pval = TRUE, conf.int = TRUE,
                  palette = c("#E41A1C", "#377EB8"), # Red: high score; blue: low score
                  xlim = c(0, 2000), break.time.by = 500, # Display follow-up through 2,000 days
                  title = paste("Subtype:", subtype_name),
                  xlab = "Time (Days)", legend.labs = c("High", "Low"),
                  ggtheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  
  return(p)
}

# ------------------------------------------------------------------------------
# Generate and assemble subtype-specific survival plots
# ------------------------------------------------------------------------------
subtypes <- c("Basal", "LumA", "LumB", "Her2")
plot_list_final <- list()

for(st in subtypes) {
  p_obj <- run_subgroup_km(all_subtypes_df, st)
  if(!is.null(p_obj)) {
    # Place the legend above each plot
    p_ggplot <- p_obj$plot + theme(legend.position = "top")
    plot_list_final[[st]] <- p_ggplot
    message(paste("Survival curves generated for subtype:", st))
  }
}

# Assemble four panels in a 2-by-2 layout
figure_1_ABCD <- plot_grid(
  plotlist = plot_list_final,
  labels = c("A", "B", "C", "D"),
  ncol = 2, 
  nrow = 2,
  align = "hv"
)

print(figure_1_ABCD)

# Save the survival panels as PDF
ggsave("Fig 1A-D Subtypes KM.pdf", 
       plot = figure_1_ABCD, 
       width = 10, 
       height = 10, 
       device = "pdf")
message("Survival curve panels saved.")

# ------------------------------------------------------------------------------
# Joint-group survival analysis in Luminal B samples
# ------------------------------------------------------------------------------
message("Generating joint-group survival curves...")

fit_combined <- survfit(Surv(OS.time, OS.status) ~ Combined_Group, data = final_df)
p_1E <- ggsurvplot(fit_combined, data = final_df, pval = TRUE, palette = "jco",
                   title = "Fig 1E: LumB ASC_P & CD4T Interaction",
                   xlab = "Time (Days)", risk.table = TRUE)
print(p_1E)

pdf("Fig 1E Interaction KM.pdf", width = 12, height = 8)
print(p_1E)
dev.off()
message("Joint-group survival plot saved.")



# ------------------------------------------------------------------------------
# Cox regression forest plot
# ------------------------------------------------------------------------------
message("Generating the Cox model forest plot...")

final_df$Combined_Group <- factor(final_df$Combined_Group)
final_df$Combined_Group <- relevel(final_df$Combined_Group, ref = "ASC_Low / CD4_High")
cox_model <- coxph(Surv(OS.time, OS.status) ~ Combined_Group, data = final_df)

p_1F <- ggforest(cox_model, data = final_df, main = "Fig 1F: Hazard Ratio for Interaction")
print(p_1F)

ggsave("Fig 1F Cox Forest.pdf", plot = p_1F, width = 10, height = 7)
message("Cox model forest plot saved.")

# ==============================================================================
# Correlations between stromal and T-cell scores
# ==============================================================================
message("Generating cell-score correlation plots...")

# ASC_P versus anergy score
p_2A <- ggscatter(final_df, x = "ASC_P_6", y = "CD4T_Anergy", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", 
                  xlab = "ASC_P Signature Score", 
                  ylab = "CD4T Anergy Score", 
                  title = "Fig 2A: ASC_P vs. CD4T Anergy",
                  color = "#E41A1C", shape = 21, fill = "lightgray", size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_2A

ggsave("Fig_2A_ASC_P_vs_Anergy.pdf", 
       plot = p_2A, 
       width = 6, 
       height = 5, 
       device = "pdf")



# ASC_P versus Treg score
p_2B <- ggscatter(final_df, x = "ASC_P_6", y = "Treg_Infiltration", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", 
                  xlab = "ASC_P Signature Score", 
                  ylab = "Treg Infiltration Score (FOXP3+)", 
                  title = "Fig 2B: ASC_P vs. Treg Infiltration",
                  color = "#984EA3", shape = 21, fill = "lightgray", size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_2B

ggsave("Fig_2B_ASC_P_vs_Treg.pdf", 
       plot = p_2B, 
       width = 6, 
       height = 5, 
       device = "pdf")


# ASC_P versus CD4 T-cell infiltration score
p_2C <- ggscatter(final_df, x = "ASC_P_6", y = "CD4T_Infiltration", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", 
                  xlab = "ASC_P Signature Score", 
                  ylab = "CD4T Infiltration Score", 
                  title = "Fig 2C: ASC_P vs. CD4T Infiltration",
                  color = "#377EB8", shape = 21, fill = "lightgray", size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p_2C

ggsave("Fig_2C_ASC_P_vs_Infiltration.pdf", 
       plot = p_2C, 
       width = 6, 
       height = 5, 
       device = "pdf")



# CD4 T-cell infiltration versus anergy score
p_2D <- ggscatter(final_df, x = "CD4T_Infiltration", y = "CD4T_Anergy", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", 
                  xlab = "CD4T Infiltration Score", 
                  ylab = "CD4T Anergy Score", 
                  title = "Fig 2D: Infiltration vs. Anergy",
                  color = "#999999", shape = 21, fill = "white", size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p_2D

ggsave("Fig_2D_Infiltration_vs_Anergy.pdf", 
       plot = p_2D, 
       width = 6, 
       height = 5, 
       device = "pdf")


# ==============================================================================
# ASC_P and CD4 T-cell exhaustion score correlation
# ==============================================================================
message("Generating the CD4 T-cell exhaustion correlation plot...")

library(ggpubr)
library(cowplot)

# ASC_P versus exhaustion score
p_2E <- ggscatter(final_df, x = "ASC_P_6", y = "CD4T_Exhaustion", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", 
                  xlab = "ASC_P Signature Score", 
                  ylab = "CD4T Exhaustion Score", 
                  title = "Fig 2E: ASC_P vs. CD4T Exhaustion",
                  color = "#FF7F00", shape = 21, fill = "lightgray", size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_2E)

ggsave("Fig_2E_ASC_P_vs_Exhaustion.pdf", 
       plot = p_2E, 
       width = 6, 
       height = 5, 
       device = "pdf")

# ==============================================================================
# Correlation between anergy and exhaustion scores
# ==============================================================================
message("Generating the anergy versus exhaustion correlation plot...")

library(ggpubr)

# Anergy versus exhaustion score
p_2F <- ggscatter(final_df, x = "CD4T_Anergy", y = "CD4T_Exhaustion", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", 
                  xlab = "CD4T Anergy Score (Early Dysfunction)", 
                  ylab = "CD4T Exhaustion Score (Terminal State)", 
                  title = "Fig 2F: Anergy vs. Exhaustion",
                  color = "#4DAF4A", shape = 21, fill = "lightgray", size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Display the plot
print(p_2F)

# Save the individual correlation plot
ggsave("Fig_2F_Anergy_vs_Exhaustion.pdf", 
       plot = p_2F, 
       width = 6, 
       height = 5, 
       device = "pdf")

message("Anergy versus exhaustion correlation plot saved.")




# Assemble a 2-by-2 panel layout
figure_2_ABCD <- plot_grid(p_2A, p_2B, p_2C, p_2D, 
                           labels = c("A", "B", "C", "D"), 
                           ncol = 2, nrow = 2, align = "hv")

print(figure_2_ABCD)

# Save as PDF
ggsave("Fig 2ABCD Cellular Interaction_Premium.pdf", 
       plot = figure_2_ABCD, 
       width = 12, 
       height = 10, 
       device = "pdf")

message("Cell-score correlation panels saved.")

# SOURCE END 90:361

# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 90:361
# ==============================================================================
# Figure 1：临床预后与交互分析 (现象篇)
# ==============================================================================
message("Generating clinical analysis plots...")

# ------------------------------------------------------------------------------
# 准备步骤：生成包含所有亚型的数据框 (用于 1A-1D)
# ------------------------------------------------------------------------------
all_subtypes_df <- clinical_base %>%
  left_join(gsva_df, by = "barcode")

# ------------------------------------------------------------------------------
# 1A-1D：定义全亚型分层生存 KM 函数 (25% 极值划分法)
# ------------------------------------------------------------------------------
run_subgroup_km <- function(data, subtype_name) {
  # 1. 过滤对应亚型
  sub_df <- data %>% dplyr::filter(BRCA_Subtype_PAM50 == subtype_name)
  
  if(nrow(sub_df) < 20) {
    message(paste("Skipping", subtype_name, ": insufficient samples"))
    return(NULL)
  }
  
  # 2. 计算 ASC_P_6 得分的上下 25% 分位数
  q_h <- quantile(sub_df$ASC_P_6, 0.75, na.rm = TRUE)
  q_l <- quantile(sub_df$ASC_P_6, 0.25, na.rm = TRUE)
  
  # 3. 提取极值样本并分组
  sub_extreme <- sub_df %>%
    dplyr::filter(ASC_P_6 >= q_h | ASC_P_6 <= q_l) %>%
    mutate(Group = ifelse(ASC_P_6 >= q_h, "High", "Low"))
  
  # 4. 确保有足够的死亡事件进行统计
  if(sum(sub_extreme$OS.status) < 3) {
    message(paste("Skipping", subtype_name, ": insufficient death events"))
    return(NULL)
  }
  
  # 5. 生存分析模型拟合与绘图
  fit <- survfit(Surv(OS.time, OS.status) ~ Group, data = sub_extreme)
  p <- ggsurvplot(fit, data = sub_extreme, pval = TRUE, conf.int = TRUE,
                  palette = c("#E41A1C", "#377EB8"), # 红蓝配色，红高蓝低
                  xlim = c(0, 2000), break.time.by = 500, # 核心修改：截断 2000 天
                  title = paste("Subtype:", subtype_name),
                  xlab = "Time (Days)", legend.labs = c("High", "Low"),
                  ggtheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  
  return(p)
}

# ------------------------------------------------------------------------------
# 执行循环并合成 Fig 1A-D 拼图
# ------------------------------------------------------------------------------
subtypes <- c("Basal", "LumA", "LumB", "Her2")
plot_list_final <- list()

for(st in subtypes) {
  p_obj <- run_subgroup_km(all_subtypes_df, st)
  if(!is.null(p_obj)) {
    # 提取 ggplot 对象，统一将图例放顶部以方便拼图
    p_ggplot <- p_obj$plot + theme(legend.position = "top")
    plot_list_final[[st]] <- p_ggplot
    message(paste("Survival curves generated for subtype:", st))
  }
}

# 使用 cowplot 将四张图拼在一起 (2x2布局)
figure_1_ABCD <- plot_grid(
  plotlist = plot_list_final,
  labels = c("A", "B", "C", "D"),
  ncol = 2, 
  nrow = 2,
  align = "hv"
)

print(figure_1_ABCD)

# 保存 1A-1D 高清 PDF
ggsave("Fig 1A-D Subtypes KM.pdf", 
       plot = figure_1_ABCD, 
       width = 10, 
       height = 10, 
       device = "pdf")
message("Survival curve panels saved.")

# ------------------------------------------------------------------------------
# 1E：LumB 交互生存分析 (使用只包含 LumB 的 final_df)
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
# 1F：Cox 回归森林图
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
# Figure 2：细胞互作相关性分析 (互作篇)
# ==============================================================================
message("Generating cell-score correlation plots...")

# 2A：ASC_P vs Anergy (红色系，代表功能废黜风险)
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



# 2B：ASC_P vs Treg (紫色系，代表特异性抑制细胞的招募)
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


# 2C：ASC_P vs Infiltration (蓝色系，代表整体免疫招募)
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



# 2D：Infiltration vs Anergy (灰色系，特异性对照，证明失能非数量引起)
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
# 补充模块：Fig 2E - ASC_P 与 CD4T Exhaustion (耗竭) 的相关性分析
# ==============================================================================
message("Generating the CD4 T-cell exhaustion correlation plot...")

library(ggpubr)
library(cowplot)

# 2E：ASC_P vs Exhaustion (橙色系，代表慢性抗原刺激下的免疫耗竭)
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
# 补充模块：Fig 2F - Anergy (失能) 与 Exhaustion (耗竭) 的内在相关性
# ==============================================================================
message("Generating the anergy versus exhaustion correlation plot...")

library(ggpubr)

# 2F：Anergy vs Exhaustion (墨绿色系，代表 T 细胞命运的衰竭轨迹)
p_2F <- ggscatter(final_df, x = "CD4T_Anergy", y = "CD4T_Exhaustion", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", 
                  xlab = "CD4T Anergy Score (Early Dysfunction)", 
                  ylab = "CD4T Exhaustion Score (Terminal State)", 
                  title = "Fig 2F: Anergy vs. Exhaustion",
                  color = "#4DAF4A", shape = 21, fill = "lightgray", size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 打印显示在屏幕上
print(p_2F)

# 单独保存这张极具机制说明力的图
ggsave("Fig_2F_Anergy_vs_Exhaustion.pdf", 
       plot = p_2F, 
       width = 6, 
       height = 5, 
       device = "pdf")

message("Anergy versus exhaustion correlation plot saved.")




# 合成 2x2 拼图
figure_2_ABCD <- plot_grid(p_2A, p_2B, p_2C, p_2D, 
                           labels = c("A", "B", "C", "D"), 
                           ncol = 2, nrow = 2, align = "hv")

print(figure_2_ABCD)

# 保存高分辨率 PDF
ggsave("Fig 2ABCD Cellular Interaction_Premium.pdf", 
       plot = figure_2_ABCD, 
       width = 12, 
       height = 10, 
       device = "pdf")

message("Cell-score correlation panels saved.")

# SOURCE END 90:361

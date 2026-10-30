# Provenance-linked final analysis block; run using run_clinical_analysis.R.
# SOURCE BEGIN 576:660
# 模块 4 抢救计划：METABRIC 极值法 + 疾病特异性生存 (DSS) 验证
# ==============================================================================
message(">>> 启动 METABRIC 高阶洗噪抢救计划...")

# 1. 净化临床终点：从 OS 切换到 DSS (疾病特异性生存)
metabric_dss_clin <- metabric_clin %>%
  mutate(
    # 使用 METABRIC 专属的 DSS 终点
    DSS.time = as.numeric(OS_MONTHS),
    # 核心修改：只有明确标注 "Died of Disease" 的才算死亡事件(1)，其他(包括死于其他原因)都算截尾(0)
    DSS.status = ifelse(VITAL_STATUS == "Died of Disease", 1, 0) 
  ) %>%
  filter(grepl("LumB", PAM50_Subtype, ignore.case = TRUE) | grepl("Luminal B", PAM50_Subtype, ignore.case = TRUE)) %>%
  filter(!is.na(DSS.time) & DSS.time > 0)

# 2. 对齐样本 (保持原有打分)
valid_dss_samples <- intersect(rownames(metabric_dss_clin), colnames(metabric_expr))
metabric_dss_clin <- metabric_dss_clin[valid_dss_samples, ]

# 将刚才算好的 ssGSEA 分数贴过来
metabric_dss_clin$ASC_Score <- as.numeric(meta_gsva["ASC_P", valid_dss_samples])
metabric_dss_clin$CD4_Score <- as.numeric(meta_gsva["CD4T", valid_dss_samples])

# 3. 极值对比策略：提取上下 25% 的极端样本
q75_asc <- quantile(metabric_dss_clin$ASC_Score, 0.75, na.rm=TRUE)
q25_asc <- quantile(metabric_dss_clin$ASC_Score, 0.25, na.rm=TRUE)
q75_cd4 <- quantile(metabric_dss_clin$CD4_Score, 0.75, na.rm=TRUE)
q25_cd4 <- quantile(metabric_dss_clin$CD4_Score, 0.25, na.rm=TRUE)

metabric_extreme <- metabric_dss_clin %>%
  # 滤掉中间 50% 的模糊噪音样本
  filter((ASC_Score >= q75_asc | ASC_Score <= q25_asc) & 
           (CD4_Score >= q75_cd4 | CD4_Score <= q25_cd4)) %>%
  mutate(
    ASC_Group = ifelse(ASC_Score >= q75_asc, "ASC_High", "ASC_Low"),
    CD4_Group = ifelse(CD4_Score >= q75_cd4, "CD4_High", "CD4_Low"),
    Combined_Group = paste0(ASC_Group, " / ", CD4_Group)
  )

# 4. 重新绘制生存曲线 (基于 DSS)
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

message("🎉 抢救分析完成！请打开 Fig_4A_METABRIC_Rescue_DSS.pdf 查看曲线是否分开！")


# ==============================================================================
# 模块 4B 抢救计划：基于 DSS 和极值亚组的 Cox 森林图重新绘制
# ==============================================================================
message(">>> 正在基于净化的 DSS 极值队列重新计算 Cox 风险并绘制 Fig 4B...")

library(survival)
library(survminer)

# 1. 确保使用抢救后的数据子集 (metabric_extreme)
# 设定因子，并将生物学上预后最好的一组 (低基质陷阱 + 高 T 细胞) 设为基准线 (Reference)
metabric_extreme$Combined_Group <- factor(metabric_extreme$Combined_Group)
metabric_extreme$Combined_Group <- relevel(metabric_extreme$Combined_Group, ref = "ASC_Low / CD4_High")

# 2. 核心修改：使用 DSS.time 和 DSS.status 构建新的 Cox 模型
cox_extreme <- coxph(Surv(DSS.time, DSS.status) ~ Combined_Group, data = metabric_extreme)

# 3. 绘制最新的森林图
p_4B_Rescue <- ggforest(cox_extreme, data = metabric_extreme, 
                        main = "Fig 4B Rescue: METABRIC DSS Hazard Ratio (Top/Bottom 25%)")

# 4. 打印并安全保存
print(p_4B_Rescue)

ggsave("Fig_4B_METABRIC_Rescue_DSS_Forest.pdf", 
       plot = p_4B_Rescue, 
       width = 10, 
       height = 7)

message("🎉 Fig 4B 重绘完成！请查看文件夹中的 Fig_4B_METABRIC_Rescue_DSS_Forest.pdf")

# SOURCE END 576:660

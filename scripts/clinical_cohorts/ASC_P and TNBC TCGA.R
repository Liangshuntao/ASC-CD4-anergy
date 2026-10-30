##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)

getwd()
setwd("D:/R/20260317/")


# ==============================================================================
# 课题：全亚型分层分析 —— ASC_P 基质检查点 ssGSEA 预后评估
# 核心：全量计算打分 + 手动亚型拆分 + 5年聚焦生存分析
# ==============================================================================

# ------------------------------------------------------------------------------
# 模块 1：全量数据加载与环境准备
# ------------------------------------------------------------------------------
message(">>> 模块 1：正在加载官方 BRCA 数据...")
options(timeout = 3600)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(GSVA)
  library(survival)
  library(survminer)
  library(dplyr)
})

# 下载全量 BRCA 数据
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query)
brca_data <- GDCprepare(query)

# ------------------------------------------------------------------------------
# 模块 2：构建全量临床背景 (含亚型与生存)
# ------------------------------------------------------------------------------
message(">>> 模块 2：整合临床亚型与生存数据...")

# 获取官方亚型列表
brca_subtype <- TCGAquery_subtype(tumor = "brca") %>%
  select(patient, BRCA_Subtype_PAM50)

# 提取并清洗全量临床表
clinical_base <- as.data.frame(colData(brca_data)) %>%
  mutate(patient_id = substr(barcode, 1, 12)) %>%
  filter(sample_type == "Primary Tumor") %>%
  # 注入亚型
  left_join(brca_subtype, by = c("patient_id" = "patient")) %>%
  # 计算生存
  mutate(
    OS.time = ifelse(vital_status == "Dead", as.numeric(days_to_death), as.numeric(days_to_last_follow_up)),
    OS.status = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(OS.time) & OS.time > 0 & !is.na(BRCA_Subtype_PAM50))

# ------------------------------------------------------------------------------
# 模块 3 & 4：全背景 ssGSEA 算法打分
# ------------------------------------------------------------------------------
message(">>> 模块 3 & 4：全量样本 ssGSEA 计算中 (请稍候)...")

counts_tpm <- assay(brca_data, "tpm_unstrand")
gene_info <- rowData(brca_data)
valid_rows <- which(!is.na(gene_info$gene_name) & !duplicated(gene_info$gene_name))

# 构建全背景矩阵 (仅保留有亚型标注的病人)
expr_full <- as.matrix(counts_tpm[valid_rows, clinical_base$barcode])
rownames(expr_full) <- gene_info$gene_name[valid_rows]

# 您的专属 10 基因 Signature
asc_p_list <- list(ASC_P = intersect(c("ENG", "THY1",'NT5E',"PDGFRA","VCAM1","MCAM"), rownames(expr_full))) 


# 执行打分
if(packageVersion("GSVA") >= "1.50.0") {
  gsva_par <- ssgseaParam(expr_full, asc_p_list)
  ssgsea_res <- gsva(gsva_par)
} else {
  ssgsea_res <- gsva(expr_full, asc_p_list, method="ssgsea")
}

# 合并分数到临床总表
ssgsea_df <- as.data.frame(t(ssgsea_res))
ssgsea_df$barcode <- rownames(ssgsea_df)
clinical_final <- clinical_base %>% 
  left_join(ssgsea_df, by = "barcode") %>%
  rename(ASC_P_Score = ASC_P)

# ==============================================================================
# 模块 5 修复版：显式调用 dplyr::filter，规避函数名冲突
# ==============================================================================
message(">>> 模块 5：正在生成各亚型生存曲线图 (已修复函数名冲突)...")

# 定义分析函数
run_subgroup_km <- function(data, subtype_name) {
  # 【核心修复】使用 dplyr:: 确保调用正确的方法
  sub_df <- data %>% dplyr::filter(BRCA_Subtype_PAM50 == subtype_name)
  
  if(nrow(sub_df) < 20) {
    message(paste("跳过", subtype_name, "：样本量太少"))
    return(NULL)
  }
  
  # 25% 极值划分
  q_h <- quantile(sub_df$ASC_P_Score, 0.75, na.rm = TRUE)
  q_l <- quantile(sub_df$ASC_P_Score, 0.25, na.rm = TRUE)
  
  sub_extreme <- sub_df %>%
    dplyr::filter(ASC_P_Score >= q_h | ASC_P_Score <= q_l) %>%
    mutate(Group = ifelse(ASC_P_Score >= q_h, "High", "Low"))
  
  # 检查死亡事件
  if(sum(sub_extreme$OS.status) < 3) {
    message(paste("跳过", subtype_name, "：死亡事件不足"))
    return(NULL)
  }
  
  fit <- survfit(Surv(OS.time, OS.status) ~ Group, data = sub_extreme)
  
  p <- ggsurvplot(fit, data = sub_extreme, pval = TRUE, conf.int = TRUE,
                  xlim = c(0, 2000), break.time.by = 500, risk.table = TRUE,
                  palette = c("#E7B800", "#2E9FDF"),
                  title = paste("Subtype:", subtype_name, "(ASC_P 5-Year OS)"),
                  xlab = "Time (Days)", legend.labs = c("High", "Low"),
                  ggtheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5)))
  return(p)
}

# 依次生成结果
subtypes <- c("Basal", "LumA", "LumB", "Her2")
for(st in subtypes) {
  p <- run_subgroup_km(clinical_final, st)
  if(!is.null(p)) {
    message(paste("✅ 成功绘制亚型曲线：", st))
    print(p)
  }
}

# ==============================================================================
# 模块 6：Figure 1 合成 —— 全亚型生存分析拼图
# ==============================================================================
message(">>> 模块 6：正在合成 Figure 1 拼图...")

library(survminer)
library(cowplot) # 如果没装请先 install.packages("cowplot")

# 1. 预备一个列表来存储生成的图
plot_list_final <- list()

# 2. 循环生成并存储图表 (去掉风险表以便拼图美观，如果需要可以后期加回)
for(st in subtypes) {
  p_obj <- run_subgroup_km(clinical_final, st)
  
  if(!is.null(p_obj)) {
    # 提取 ggplot 对象进行统一调整
    # 我们暂时去掉风险表，保持拼图整洁；增加标签辨识度
    p_ggplot <- p_obj$plot + 
      theme(legend.position = "top") + 
      labs(title = paste0("Subtype: ", st))
    
    plot_list_final[[st]] <- p_ggplot
  }
}

# 3. 使用 cowplot 进场大合成
# 设置为 2x2 布局
figure_1 <- plot_grid(
  plotlist = plot_list_final,
  labels = c("A", "B", "C", "D"), # 自动添加大写字母编号
  ncol = 2, 
  nrow = 2,
  align = "hv"
)

# 4. 预览并保存高分辨率图片 (供论文投稿使用)
print(figure_1)

ggsave("Figure_1_ASC_P_Subtypes_KM_CD90_CD105_CD73_PDGFRA_VCAM1_CD.pdf", 
       plot = figure_1, 
       width = 12 ,
       height = 10, 
       device = "pdf")

message("🎉 Figure 1 已合成并保存为高分辨率 PDF！")

# ------------------------------------------------------------------------------
# 模块 7：计算 CD4+ T 细胞浸润并与 ASC_P 执行交互分析
# ==============================================================================
# 终极模块：ASC_P (6-Gene) 与 CD4T 免疫交互分析 (全流程自启修正版)
# ==============================================================================
message(">>> 正在初始化数据环境...")

# 1. 确保基础矩阵准备就绪
counts_tpm <- assay(brca_data, "tpm_unstrand")
gene_info <- rowData(brca_data)
valid_rows <- which(!is.na(gene_info$gene_name) & !duplicated(gene_info$gene_name))

# 提取全量样本（Primary Tumor）
all_barcodes <- colData(brca_data) %>% 
  as.data.frame() %>% 
  filter(sample_type == "Primary Tumor") %>% 
  pull(barcode)

expr_matrix <- as.matrix(counts_tpm[valid_rows, all_barcodes])
rownames(expr_matrix) <- gene_info$gene_name[valid_rows]

# 2. 定义 Signature (整合 Anergy 与 Exhaustion)
# 修正点：将变量名统一为 sig_list 以便后文引用
sig_list <- list(
  ASC_P_6 = c("ENG", "THY1", "NT5E", "PDGFRA", "VCAM1", "MCAM"),
  CD4T_Infiltration = c("CD4", "CD3D", "CD3E"),
  CD4T_Anergy = c("RNF128", "CBLB", "EGR2", "DGKA", "ITCH", "IZUMO1R"), 
  CD4T_Exhaustion = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4")
)

message(">>> 正在执行多维度 ssGSEA 打分 (约 30-60 秒)...")
# 修正点：确保使用 sig_list
if(packageVersion("GSVA") >= "1.50.0") {
  gsva_results <- gsva(ssgseaParam(expr_matrix, sig_list))
} else {
  gsva_results <- gsva(expr_matrix, sig_list, method="ssgsea")
}

# 3. 临床数据大融合
gsva_df <- as.data.frame(t(gsva_results))
gsva_df$barcode <- rownames(gsva_df)

brca_subtype <- TCGAquery_subtype(tumor = "brca") %>% select(patient, BRCA_Subtype_PAM50)

final_interaction_df <- as.data.frame(colData(brca_data)) %>%
  mutate(patient_id = substr(barcode, 1, 12)) %>%
  filter(barcode %in% all_barcodes) %>%
  left_join(brca_subtype, by = c("patient_id" = "patient")) %>%
  left_join(gsva_df, by = "barcode") %>%
  mutate(
    OS.time = ifelse(vital_status == "Dead", as.numeric(days_to_death), as.numeric(days_to_last_follow_up)),
    OS.status = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(OS.time) & OS.time > 0 & !is.na(BRCA_Subtype_PAM50))

# ------------------------------------------------------------------------------
# 4. 机制验证：在 Luminal B 中看 ASC_P 是否与 CD4T 失能相关
# ------------------------------------------------------------------------------
lumb_df <- final_interaction_df %>% filter(BRCA_Subtype_PAM50 == "LumB")

# 修正点：使用正确的对象名 lumb_df，并聚焦 CD4T_Anergy
cor_test <- cor.test(lumb_df$ASC_P_6, lumb_df$CD4T_Anergy, method = "spearman")
message(paste(">>> [机制验证] Luminal B 中 ASC_P 与 CD4T Anergy 相关性 Rho:", 
              round(cor_test$estimate, 3), " P-value:", round(cor_test$p.value, 6)))

# ------------------------------------------------------------------------------
# 5. 交互生存分析：定义“免疫逃逸指数” (ASC_P 高 且 CD4T 浸润低)
# ------------------------------------------------------------------------------
lumb_df <- lumb_df %>%
  mutate(
    ASC_Group = ifelse(ASC_P_6 > median(ASC_P_6), "ASC_High", "ASC_Low"),
    CD4_Group = ifelse(CD4T_Infiltration > median(CD4T_Infiltration), "CD4_High", "CD4_Low"),
    Combined_Group = paste0(ASC_Group, " / ", CD4_Group)
  )

fit_combined <- survfit(Surv(OS.time, OS.status) ~ Combined_Group, data = lumb_df)

# 修正点：调整 legend.labs 排序以匹配曲线
p_combined <- ggsurvplot(fit_combined, data = lumb_df, pval = TRUE,
                         title = "Luminal B: ASC_P & CD4T Infiltration Interaction",
                         palette = "jco", 
                         xlab = "Time (Days)", risk.table = TRUE)

print(p_combined)

# 建议设置较大的尺寸，以确保风险表（Risk Table）清晰可见
ggsave("Figure_Interaction_ASC_P_CD4T_LumB.pdf", 
       plot = p_combined$plot,      # 注意：survplot 对象需要提取 $plot
       width = 14, 
       height = 14, 
       device = "pdf")

# 如果您想把带风险表的完整图存下来，也可以用这个 survminer 自带的函数
pdf("Figure_2_Interaction_Full_with_Table.pdf", width = 10, height = 8)
print(p_combined)
dev.off()

message("🎉 保存成功！文件名：Figure_Interaction_ASC_P_CD4T_LumB.pdf")


# ------------------------------------------------------------------------------
# 6. 方式 1：针对四组分类的 Cox 回归（对应 KM 图）
#    这段代码会告诉您，“免疫逃逸组”比“免疫监视组”风险高多少倍。
# ------------------------------------------------------------------------------

# 加载生存分析核心包
library(survival)
library(survminer)

# 将 Combined_Group 转为 Factor，并设置 ASC_Low / CD4_High 为对照组 (Reference)
lumb_df$Combined_Group <- factor(lumb_df$Combined_Group)
lumb_df$Combined_Group <- relevel(lumb_df$Combined_Group, ref = "ASC_Low / CD4_High")

# 运行 Cox 模型
cox_categorical <- coxph(Surv(OS.time, OS.status) ~ Combined_Group, data = lumb_df)

# 查看详细结果 (关注 HR 和 P 值)
summary(cox_categorical)

# 绘制森林图（文章 Figure 2/3 常用）
p_forest_ASC_Low_High=ggforest(cox_categorical, data = lumb_df, main = "Hazard Ratio for ASC_P & CD4T Interaction")

ggsave("Figure_cox_categorical_Forest.pdf", 
       plot = p_forest_ASC_Low_High, 
       width = 10, 
       height =7, 
       device = "pdf")


# ------------------------------------------------------------------------------
# 6. 方式 2：基于原始评分的多元 Cox 回归
#   验证两个指标是否互不干扰，共同影响预后。
# ------------------------------------------------------------------------------
# 使用原始连续变量进行多元 Cox
# ASC_P_6 和 CD4T_Infiltration 的原始分值
cox_continuous <- coxph(Surv(OS.time, OS.status) ~ ASC_P_6 + CD4T_Infiltration, data = lumb_df)

# 查看结果
summary(cox_continuous)

# ==============================================================================
# 模块 10：绘制多元 Cox 回归森林图 (连续变量版)
# ==============================================================================
message(">>> 正在绘制方式 2 的多元 Cox 森林图...")

# 1. 运行多元 Cox 模型 (确保包含您的核心基质分和免疫浸润分)
cox_continuous <- coxph(Surv(OS.time, OS.status) ~ ASC_P_6 + CD4T_Infiltration, data = lumb_df)

# ------------------------------------------------------------------------------
# 方案 A：使用 survminer 自带的 ggforest (简单快速)
# ------------------------------------------------------------------------------
# 注意：对于连续变量，HR 代表分值每增加 1 个单位时死亡风险的变化
p_forest_simple <- ggforest(cox_continuous, 
                            data = lumb_df, 
                            main = "Multivariate Cox: ASC_P & CD4T (Continuous)",
                            fontsize = 1.0)
print(p_forest_simple)

# ------------------------------------------------------------------------------
# 方案 A：保存结果为 PDF
# ------------------------------------------------------------------------------
message(">>> 正在保存多元 Cox 森林图 (连续变量版) 为 PDF...")

# 建议使用较宽的比例（如 10x5 或 12x6），以确保左侧变量名和右侧数据不重叠
ggsave("Figure_Cox_Continuous_Standard_Forest.pdf", 
       plot = p_forest_simple, 
       width = 15, 
       height = 5, 
       device = "pdf")

# (可选) 同时保存一份 PNG 预览图，方便放入 PPT 汇报
ggsave("Figure_Cox_Continuous_Standard_Forest.png", 
       plot = p_forest_simple, 
       width = 15, 
       height = 5, 
       dpi = 300)

message("🎉 森林图已成功保存！")


# ==============================================================================
# 模块 11：Figure 2 核心机制验证 —— 相关性散点图 (Anergy & Infiltration)
# ==============================================================================
message(">>> 模块 11：正在生成 Figure 2A 和 2B 机制验证散点图...")

library(ggpubr)
library(cowplot)

# ------------------------------------------------------------------------------
# Figure 2A：ASC_P 与 CD4T Anergy (失能) 的正相关 (功能废掉)
# ------------------------------------------------------------------------------
p_2A <- ggscatter(lumb_df, x = "ASC_P_6", y = "CD4T_Anergy",
                  add = "reg.line",                         # 添加回归线
                  conf.int = TRUE,                          # 添加置信区间
                  cor.coef = TRUE,                          # 显示相关系数
                  cor.method = "spearman",                  # 使用 Spearman
                  xlab = "ASC_P Signature Score",
                  ylab = "CD4T Anergy Score",
                  title = "Fig 2A: ASC_P vs. CD4T Anergy",
                  color = "#E41A1C", shape = 21, fill = "lightgray", # 使用红色系强调风险
                  size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_2A

# ------------------------------------------------------------------------------
# Figure 2B：ASC_P 与 CD4T Infiltration (浸润) 的相关性 (物理阻碍/排斥)
# ------------------------------------------------------------------------------
p_2B <- ggscatter(lumb_df, x = "ASC_P_6", y = "CD4T_Infiltration",
                  add = "reg.line", 
                  conf.int = TRUE, 
                  cor.coef = TRUE, 
                  cor.method = "spearman",
                  xlab = "ASC_P Signature Score",
                  ylab = "CD4T Infiltration Score",
                  title = "Fig 2B: ASC_P vs. CD4T Infiltration",
                  color = "#377EB8", shape = 21, fill = "lightgray", # 使用蓝色系代表免疫浸润
                  size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_2B
# ------------------------------------------------------------------------------
# 保存与合成拼图
# ------------------------------------------------------------------------------
# 合并 2A 和 2B
figure_2_mechanisms <- plot_grid(p_2A, p_2B, 
                                 labels = c("A", "B"), 
                                 ncol = 2)

print(figure_2_mechanisms)

# 保存高分辨率 PDF
ggsave("Figure_2AB_Mechanisms_Correlation.pdf", 
       plot = figure_2_mechanisms, 
       width = 12, 
       height = 5, 
       device = "pdf")

message("🎉 Figure 2A (失能) 与 2B (浸润) 相关性验证图已保存！")

# 💡 导师小贴士：
# 1. 关注 Rho 值：2A 理论上应为显著正相关 (Rho > 0)，证明 ASC_P 诱导失能。
# 2. 关注 Rho 值：2B 如果为正相关，则证明 ASC_P 具有“诱导免疫耐受”效应，使T 细胞免疫功能失调进入。


# ==============================================================================
# 模块 12：免疫背景深度解析 —— CD4T 浸润与失能的相关性 (Luminal B)
# ==============================================================================
message(">>> 正在分析 CD4T 浸润与失能的内在联系...")

p_immune_internal <- ggscatter(lumb_df, x = "CD4T_Infiltration", y = "CD4T_Anergy",
                               add = "reg.line", 
                               conf.int = TRUE, 
                               cor.coef = TRUE, 
                               cor.method = "spearman",
                               xlab = "CD4T Infiltration Score (Quantity)",
                               ylab = "CD4T Anergy Score (Quality)",
                               title = "LumB: Immune Quantity vs. Dysfunction",
                               color = "#999999", shape = 21, fill = "white",
                               size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_immune_internal)

# 保存 PDF
ggsave("Figure_2C_Infiltration_vs_Anergy_Correlation.pdf", 
       plot = p_immune_internal, width = 5, height = 5)


# ==============================================================================
# 模块 13 (修正版)：GSEA 通路富集分析 —— 确保逻辑统一 (High vs Low)
# ==============================================================================
message(">>> 正在纠正比较逻辑：以 ASC_Low 为基准，分析 ASC_High 的通路活化...")

# 1. 关键修正：设置 Low 为对照组 (Reference Level)
# 这样 coef = 2 才是真正的 ASC_High 减去 ASC_Low
group_info <- factor(lumb_df$ASC_Group, levels = c("ASC_Low", "ASC_High"))

# 2. 重新运行 limma
design <- model.matrix(~ group_info)
fit <- lmFit(lumb_expr, design)
fit <- eBayes(fit)

# 这里的 logFC 正值代表在 ASC_High 中高表达
diff_res <- topTable(fit, coef = 2, number = Inf)
diff_res$gene <- rownames(diff_res)

# 3. 重新构建 Rank 列表
diff_res$stat <- sign(diff_res$logFC) * -log10(diff_res$P.Value)
gene_list <- diff_res$stat
names(gene_list) <- diff_res$gene
gene_list <- sort(gene_list, decreasing = TRUE)

# 4. 重新执行 GSEA
gsea_res_correct <- GSEA(gene_list, TERM2GENE = h_df, pvalueCutoff = 0.05, eps = 0)

# 5. 绘图：此时 Activated 就代表在 ASC_High 中真正活化的通路
p_gsea_correct <- dotplot(gsea_res_correct, showCategory = 20, split = ".sign") + 
  facet_grid(.~.sign) +
  scale_fill_gradientn(colours = c("#377EB8", "white", "#E41A1C")) + # 蓝色抑制，红色激活
  theme_bw() +
  labs(title = "GSEA: Hallmark Pathways (Activated in ASC_P High Group)",
       subtitle = "Reference: ASC_P Low | Comparison: High vs Low")

print(p_gsea_correct)

# 保存修正后的结果
ggsave("Figure_GSEA_Correct_Logic_High_vs_Low.pdf", 
       plot = p_gsea_correct, width = 10, height = 10)

# ==============================================================================
# 模块 14：通路核心基因热图展示 (纠正逻辑版)
# ==============================================================================
message(">>> 正在提取核心通路基因并准备热图数据...")

library(pheatmap)
library(RColorBrewer)

# 1. 确保比较逻辑正确：以 Low 为对照
lumb_df$ASC_Group <- factor(lumb_df$ASC_Group, levels = c("ASC_Low", "ASC_High"))

# 2. 选择我们要展示的通路 (从 msigdbr 获取基因列表)
target_pathways <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
                     "HALLMARK_TGF_BETA_SIGNALING",
                     "HALLMARK_IL6_JAK_STAT3_SIGNALING")

pathway_genes <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::filter(gs_name %in% target_pathways)

# 3. 准备表达矩阵 (挑选在通路中且在我们的表达矩阵中的基因)
# 为了热图美观，我们挑选每条通路中差异最显著的前 20 个基因
get_top_genes <- function(pathway_name, n = 20) {
  genes <- pathway_genes$gene_symbol[pathway_genes$gs_name == pathway_name]
  res <- diff_res %>% # 使用之前 limma 算的差异结果
    dplyr::filter(gene %in% genes) %>%
    dplyr::arrange(P.Value) %>%
    head(n) %>%
    pull(gene)
  return(res)
}

plot_genes <- unique(c(get_top_genes("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
                       get_top_genes("HALLMARK_TGF_BETA_SIGNALING"),
                       get_top_genes("HALLMARK_IL6_JAK_STAT3_SIGNALING")))

# 4. 提取并标准化数据
plot_matrix <- lumb_expr[plot_genes, ]
plot_matrix <- t(scale(t(plot_matrix))) # Z-score 标准化

# 5. 准备列注释 (样本信息)
annotation_col <- lumb_df %>%
  dplyr::select(barcode, ASC_Group, CD4T_Infiltration, OS.status) %>%
  dplyr::arrange(ASC_Group, desc(CD4T_Infiltration)) %>% # 按组别和浸润排序，使热图趋势明显
  as.data.frame()

rownames(annotation_col) <- annotation_col$barcode
annotation_col$barcode <- NULL
annotation_col$OS.status <- factor(annotation_col$OS.status, labels = c("Alive", "Dead"))

# 调整矩阵顺序以匹配注释
plot_matrix <- plot_matrix[, rownames(annotation_col)]

# 6. 绘图
message(">>> 正在绘制热图...")
p_heatmap <- pheatmap(plot_matrix,
                      annotation_col = annotation_col,
                      show_colnames = FALSE,
                      cluster_cols = FALSE, # 我们已经手动排序了，所以不聚类列，保持趋势
                      cluster_rows = TRUE,
                      color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
                      main = "Key Pathway Genes: ASC_P High vs Low (Luminal B)",
                      fontsize_row = 8,
                      gaps_col = table(annotation_col$ASC_Group)[1]) # 在两组之间加个缝隙

# 7. 保存结果
pdf("Figure_3_Pathway_Heatmap_EMT_TGFB.pdf", width = 10, height = 12)
print(p_heatmap)
dev.off()

message("🎉 热图已保存：Figure_3_Pathway_Heatmap_EMT_TGFB.pdf")

# ==============================================================================
# 模块 15：基于通路评分 (ssGSEA) 的热图展示 —— 简化逻辑，直击核心
# ==============================================================================
message(">>> 正在计算 Hallmark 通路评分并准备热图...")

library(GSVA)
library(pheatmap)
library(msigdbr)
library(dplyr)

# 1. 获取 Hallmark 基因集
h_list <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# 2. 执行通路水平的 ssGSEA 打分 (针对 Luminal B 表达矩阵)
# 这一步将“基因维度”降为“通路维度”
if(packageVersion("GSVA") >= "1.50.0") {
  pathway_ssgsea <- gsva(ssgseaParam(lumb_expr, h_list))
} else {
  pathway_ssgsea <- gsva(lumb_expr, h_list, method="ssgsea")
}

# 3. 挑选在 GSEA 中显著且与您假说最相关的核心通路
# 我们聚焦在免疫招募、基质重塑和免疫抑制这三个维度
selected_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", # 基质重塑 (您的 6-基因基础)
  "HALLMARK_TGF_BETA_SIGNALING",                # 诱导失能的核心 (Fig 2A 机制)
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",           # 免疫抑制环境
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",           # 炎症背景
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",         # 免疫招募信号 (解释 Fig 2B)
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_ANGIOGENESIS",                      # 血管生成 (与基质相关)
  "HALLMARK_HYPOXIA",                           # 缺氧环境
  "HALLMARK_COAGULATION"
)

# 4. 提取数据并进行 Z-score 标准化
plot_matrix_pathway <- pathway_ssgsea[selected_pathways, ]
plot_matrix_pathway <- t(scale(t(plot_matrix_pathway)))

# 5. 准备列注释 (样本信息) - 保持与之前一致的排序逻辑
annotation_col <- lumb_df %>%
  dplyr::select(barcode, ASC_Group, CD4T_Infiltration, OS.status) %>%
  dplyr::arrange(ASC_Group, desc(CD4T_Infiltration)) %>% 
  as.data.frame()

rownames(annotation_col) <- annotation_col$barcode
annotation_col$barcode <- NULL
annotation_col$OS.status <- factor(annotation_col$OS.status, labels = c("Alive", "Dead"))

# 匹配矩阵列顺序
plot_matrix_pathway <- plot_matrix_pathway[, rownames(annotation_col)]

# 6. 绘图：去掉前缀，让通路名称更清晰
rownames(plot_matrix_pathway) <- gsub("HALLMARK_", "", rownames(plot_matrix_pathway))

message(">>> 正在绘制通路水平热图...")
p_heatmap_pathway <- pheatmap(plot_matrix_pathway,
                              annotation_col = annotation_col,
                              show_colnames = FALSE,
                              cluster_cols = FALSE, 
                              cluster_rows = TRUE,
                              color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
                              main = "Pathway Activity: ASC_P High vs Low (Luminal B)",
                              fontsize_row = 10,
                              gaps_col = table(annotation_col$ASC_Group)[1])

# 7. 保存结果
pdf("Figure_3_Pathway_Score_Heatmap.pdf", width = 10, height = 7)
print(p_heatmap_pathway)
dev.off()

message("🎉 通路评分热图已保存：Figure_3_Pathway_Score_Heatmap.pdf")

# ==============================================================================
# 模块 16：升级版通路评分热图 —— 加入 CD4T_Anergy (失能) 分布
# ==============================================================================
message(">>> 正在生成包含失能评分 (Anergy) 的升级版热图...")

library(pheatmap)
library(dplyr)

# 1. 准备列注释 (加入 CD4T_Anergy)
# 保持之前的排序逻辑：优先按 ASC_Group 排序，组内按 CD4T_Infiltration 排序
annotation_col_v2 <- lumb_df %>%
  dplyr::select(barcode, ASC_Group, CD4T_Infiltration, CD4T_Anergy, OS.status) %>%
  dplyr::arrange(ASC_Group, desc(CD4T_Infiltration)) %>% 
  as.data.frame()

rownames(annotation_col_v2) <- annotation_col_v2$barcode
annotation_col_v2$barcode <- NULL
annotation_col_v2$OS.status <- factor(annotation_col_v2$OS.status, labels = c("Alive", "Dead"))

# 2. 匹配矩阵列顺序并进行重标准化
# 确保使用之前选定的 10 条核心通路 (selected_pathways)
plot_matrix_v2 <- pathway_ssgsea[selected_pathways, rownames(annotation_col_v2)]
plot_matrix_v2 <- t(scale(t(plot_matrix_v2))) 
rownames(plot_matrix_v2) <- gsub("HALLMARK_", "", rownames(plot_matrix_v2))

# 3. 绘图：此时顶部会有 ASC分组、浸润分、失能分 三重标注
message(">>> 正在绘制升级版热图...")
p_heatmap_v2 <- pheatmap(plot_matrix_v2,
                         annotation_col = annotation_col_v2,
                         show_colnames = FALSE,
                         cluster_cols = FALSE, # 关键：不聚类列，保持我们的排序逻辑
                         cluster_rows = TRUE,
                         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
                         main = "Pathway Activity & CD4T Status: ASC_P High vs Low",
                         fontsize_row = 10,
                         gaps_col = table(annotation_col_v2$ASC_Group)[1])

# 4. 保存结果
pdf("Figure_3_Pathway_Score_with_Anergy_Heatmap.pdf", width = 11, height = 8)
print(p_heatmap_v2)
dev.off()

message("🎉 升级版热图已保存：Figure_3_Pathway_Score_with_Anergy_Heatmap.pdf")


# ==============================================================================
# 模块 17：最终优化版热图 —— 按 CD4T_Anergy (失能评分) 排序
# ==============================================================================
message(">>> 正在生成基于 Anergy 排序的优化版热图...")

library(pheatmap)
library(dplyr)

# 1. 核心逻辑调整：按照 ASC_Group 分组，组内按照 CD4T_Anergy 从高到低排序
# 这样可以直观观察：通路活化是否与“失能程度”呈现完美的线性梯队
annotation_col_v3 <- lumb_df %>%
  dplyr::select(barcode, ASC_Group, CD4T_Anergy, CD4T_Infiltration, OS.status) %>%
  dplyr::arrange(ASC_Group, desc(CD4T_Anergy)) %>% # 改为按 Anergy 降序
  as.data.frame()

rownames(annotation_col_v3) <- annotation_col_v3$barcode
annotation_col_v3$barcode <- NULL
annotation_col_v3$OS.status <- factor(annotation_col_v3$OS.status, labels = c("Alive", "Dead"))

# 2. 匹配矩阵列顺序 (基于之前计算好的 pathway_ssgsea)
plot_matrix_v3 <- pathway_ssgsea[selected_pathways, rownames(annotation_col_v3)]
plot_matrix_v3 <- t(scale(t(plot_matrix_v3))) 
rownames(plot_matrix_v3) <- gsub("HALLMARK_", "", rownames(plot_matrix_v3))

# 3. 绘图
message(">>> 正在绘制 Anergy 排序热图...")
p_heatmap_v3 <- pheatmap(plot_matrix_v3,
                         annotation_col = annotation_col_v3,
                         show_colnames = FALSE,
                         cluster_cols = FALSE, # 关键：严格遵守我们的 Anergy 排序
                         cluster_rows = TRUE,
                         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
                         main = "Pathway Activity Sorted by CD4T Anergy Score (Luminal B)",
                         fontsize_row = 10,
                         gaps_col = table(annotation_col_v3$ASC_Group)[1])

# 4. 保存结果
pdf("Figure_3_Pathway_Sorted_by_Anergy.pdf", width = 11, height = 8)
print(p_heatmap_v3)
dev.off()

message("🎉 优化版热图已保存：Figure_3_Pathway_Sorted_by_Anergy.pdf")


# ==============================================================================
# 模块 18：补充证据 —— ASC_P 与 Treg 浸润的相关性分析 (针对 LumB)
# ==============================================================================
message(">>> 正在分析 ASC_P 与 Treg 浸润的相关性...")

# 1. 定义 Treg 核心 Signature (基于经典标志物 FOXP3, IL2RA, IKZF2, CCR4)
# 这是国际公认的识别 Treg 的最小基因集
sig_treg <- list(Treg_Infiltration = c("FOXP3", "IL2RA", "IKZF2", "CCR4"))

# 2. 执行 ssGSEA 打分
if(packageVersion("GSVA") >= "1.50.0") {
  treg_gsva <- gsva(ssgseaParam(expr_matrix[, lumb_df$barcode], sig_treg))
} else {
  treg_gsva <- gsva(expr_matrix[, lumb_df$barcode], sig_treg, method="ssgsea")
}

# 3. 将打分合并到 Luminal B 临床主表中
lumb_df$Treg_Infiltration <- as.numeric(treg_gsva[1, ])

# 4. 绘制相关性散点图 (Figure 2D)
library(ggpubr)
p_treg_cor <- ggscatter(lumb_df, x = "ASC_P_6", y = "Treg_Infiltration",
                        add = "reg.line", 
                        conf.int = TRUE, 
                        cor.coef = TRUE, 
                        cor.method = "spearman",
                        xlab = "ASC_P Signature Score",
                        ylab = "Treg Infiltration Score (FOXP3+)",
                        title = "Fig 2D: ASC_P vs. Treg Infiltration",
                        color = "#984EA3", shape = 21, fill = "lightgray", # 使用紫色代表抑制性 Treg
                        size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_treg_cor)

# 5. 保存结果
ggsave("Figure_2D_Treg_Correlation.pdf", plot = p_treg_cor, width = 5, height = 5)

message("🎉 Treg 相关性分析已保存：Figure_2D_Treg_Correlation.pdf")



# ==============================================================================
# 模块 18：补充证据 —— ASC_P 与 Treg 浸润的相关性分析 (针对 LumB)
# ==============================================================================
message(">>> 正在分析 ASC_P 与 Treg 浸润的相关性...")

# 1. 定义 Treg 核心 Signature (基于经典标志物 FOXP3, IL2RA, IKZF2, CCR4)
# 这是国际公认的识别 Treg 的最小基因集
sig_treg <- list(Treg_Infiltration = c("FOXP3", "IL2RA", "IKZF2", "CCR4"))

# 2. 执行 ssGSEA 打分
if(packageVersion("GSVA") >= "1.50.0") {
  treg_gsva <- gsva(ssgseaParam(expr_matrix[, lumb_df$barcode], sig_treg))
} else {
  treg_gsva <- gsva(expr_matrix[, lumb_df$barcode], sig_treg, method="ssgsea")
}

# 3. 将打分合并到 Luminal B 临床主表中
lumb_df$Treg_Infiltration <- as.numeric(treg_gsva[1, ])

# 4. 绘制相关性散点图 (Figure 2D)
library(ggpubr)
p_treg_cor <- ggscatter(lumb_df, x = "ASC_P_6", y = "Treg_Infiltration",
                        add = "reg.line", 
                        conf.int = TRUE, 
                        cor.coef = TRUE, 
                        cor.method = "spearman",
                        xlab = "ASC_P Signature Score",
                        ylab = "Treg Infiltration Score (FOXP3+)",
                        title = "Fig 2D: ASC_P vs. Treg Infiltration",
                        color = "#984EA3", shape = 21, fill = "lightgray", # 使用紫色代表抑制性 Treg
                        size = 2, alpha = 0.6) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_treg_cor)

# 5. 保存结果
ggsave("Figure_2D_Treg_Correlation.pdf", plot = p_treg_cor, width = 5, height = 5)

message("🎉 Treg 相关性分析已保存：Figure_2D_Treg_Correlation.pdf")


# ==============================================================================
# 模块 19：终极版“免疫抑制全景热图” —— 整合 Treg、Anergy 与核心通路
# ==============================================================================
message(">>> 正在生成终极版机制热图：按 Anergy 排序，整合 Treg 浸润...")

library(pheatmap)
library(dplyr)

# 1. 准备列注释：包含 ASC分组、失能评分(Anergy)、Treg浸润、总CD4T浸润和生存状态
# 核心排序逻辑：优先按 ASC_Group 区分，组内按 CD4T_Anergy 从高到低排列
annotation_col_final <- lumb_df %>%
  dplyr::select(barcode, ASC_Group, CD4T_Anergy, Treg_Infiltration, CD4T_Infiltration, OS.status) %>%
  dplyr::arrange(ASC_Group, desc(CD4T_Anergy)) %>% 
  as.data.frame()

rownames(annotation_col_final) <- annotation_col_final$barcode
annotation_col_final$barcode <- NULL
annotation_col_final$OS.status <- factor(annotation_col_final$OS.status, labels = c("Alive", "Dead"))

# 2. 准备矩阵数据 (使用之前 ssGSEA 计算的通路评分)
# 匹配最终确定的 10 条核心通路
plot_matrix_final <- pathway_ssgsea[selected_pathways, rownames(annotation_col_final)]

# Z-score 标准化：让不同通路的波动具有可比性
plot_matrix_final <- t(scale(t(plot_matrix_final))) 
rownames(plot_matrix_final) <- gsub("HALLMARK_", "", rownames(plot_matrix_final))

# 3. 设置注释条的颜色 (可选：让颜色更符合学术审美)
ann_colors = list(
  ASC_Group = c(ASC_Low = "#377EB8", ASC_High = "#E41A1C"),
  OS.status = c(Alive = "#BEBEBE", Dead = "black"),
  CD4T_Anergy = colorRampPalette(c("white", "#E41A1C"))(100),
  Treg_Infiltration = colorRampPalette(c("white", "#984EA3"))(100) # Treg 用紫色
)

# 4. 绘图：展示“基质-Treg-失能-通路”四位一体的逻辑
message(">>> 正在生成 PDF 结果...")
p_final <- pheatmap(plot_matrix_final,
                    annotation_col = annotation_col_final,
                    annotation_colors = ann_colors,
                    show_colnames = FALSE,
                    cluster_cols = FALSE, # 锁定我们的 Anergy 排序逻辑
                    cluster_rows = TRUE,
                    color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
                    main = "The Integrated Landscape of ASC_P Induced Immune Dysfunction",
                    fontsize_row = 10,
                    gaps_col = table(annotation_col_final$ASC_Group)[1])

# 5. 保存结果
pdf("Figure_3_Ultimate_Mechanism_Heatmap.pdf", width = 12, height = 8)
print(p_final)
dev.off()

message("🎉 终极版热图已生成：Figure_3_Ultimate_Mechanism_Heatmap.pdf")

save.image()


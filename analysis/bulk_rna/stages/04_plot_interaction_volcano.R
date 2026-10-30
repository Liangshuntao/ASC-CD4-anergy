# Figure 8B highlighted volcano source; visual edits not reproduced
# Executed by run_bulk_analysis.R; not a standalone entry point.

# SOURCE BEGIN 491:558
library(ggrepel)
library(ggplot2)
library(dplyr)

# 1. 准备 Top 15 基因 (保持原逻辑)
top_lab <- res_did %>%
  arrange(FDR) %>%
  slice(1:15) %>%
  pull(Gene)

# 2. 【关键修改】定义您想强制标记的基因 (手动添加 Nt5e 等)
my_interest_genes <- c("Nt5e", "Jag1", "Eng", "Has2", "Alcam","Emcn","Pecam1")

# 3. 合并标签列表 (Top 15 + 您感兴趣的基因，并去重)
final_labels <- unique(c(top_lab, my_interest_genes))

# 4. 准备绘图数据 (保持原逻辑)
dfv <- res_did %>%
  mutate(sig = case_when(
    FDR < 0.05 & logFC >=  0.5 ~ "Up",
    FDR < 0.05 & logFC <= -0.5 ~ "Down",
    TRUE ~ "ns"
  ))

# 5. 绘图
p_vol <- ggplot(dfv, aes(x = logFC, y = -log10(FDR))) +
  # (1) 基础点图层
  geom_point(aes(color = sig), alpha = 0.6, size = 1.6) +
  scale_color_manual(values = c("Up"="#D95F02", "Down"="#1B9E77", "ns"="grey85")) +
  
  # (2) 【新增】高亮重点基因的点 (加黑边、填黄色，醒目)
  geom_point(data = dfv %>% filter(Gene %in% my_interest_genes),
             color = "black",   # 边框黑色
             fill = "yellow",   # 填充黄色
             shape = 21,        # 空心圆模式(支持fill)
             size = 3.5,        # 稍微大一点
             stroke = 1) +      # 边框加粗
  
  # (3) 辅助线
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  
  # (4) 【修改】标签层 (使用合并后的 final_labels)
  ggrepel::geom_text_repel(
    data = dfv %>% filter(Gene %in% final_labels), # 这里改成了 final_labels
    aes(label = Gene),
    size = 4,                # 字稍微大一点
    fontface = "bold",       # 加粗
    box.padding = 0.5,       # 文本框边距
    point.padding = 0.5,     # 点边距
    max.overlaps = Inf,      # 允许重叠无限大(强制显示)
    min.segment.length = 0   # 总是显示连线
  ) +
  
  # (5) 主题修饰
  theme_bw() +
  ggtitle("ASC DID (P-specific coculture response): differential genes") +
  xlab("log2 fold change (interaction term)") +
  ylab(expression(-log[10](FDR))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        panel.grid.major = element_blank(), # 去除网格线更清爽
        panel.grid.minor = element_blank())

# 6. 显示并保存
print(p_vol)

ggsave("ASC_DID_volcano_highlight.3.pdf", p_vol, width = 10, height = 8)

# SOURCE END 491:558

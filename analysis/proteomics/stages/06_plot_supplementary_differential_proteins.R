# SOURCE BEGIN 208:264
# 4.4 火山图 + 热图（ASC）
colnames(res_ASC_anno)

p.asc.volcano=ggplot(res_ASC_anno,
       aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.7, size = 0.5) +
  scale_color_manual(values = c("ns" = "grey70", "sig" = "red")) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  ggtitle("Volcano - ASC_P vs ASC_N")

p.asc.volcano

ggsave(filename = "ASC protein Volcano.plot.pdf",
       plot = p.asc.volcano,
       width = 7,
       height =6 )



# 取 top 50 显著差异蛋白画热图
top50_asc <- res_ASC_anno %>%
  filter(sig == "sig") %>%
  arrange(FDR) %>%
  slice(1:min(50, n())) %>%
  pull(Accession)

mat_top50_asc <- expr_imp_asc[top50_asc, ]

library(pheatmap)

pdf("ASC_P vs ASC_N top50 DEPs.pdf", width = 4, height = 5)
pheatmap(mat_top50_asc,
         scale = "row",
         annotation_col = data.frame(condition = group_asc,
                                     row.names = colnames(mat_top50_asc)),
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = FALSE,
         main = "ASC_P vs ASC_N - top50 DEPs")


dev.off()

pdf("CD105+CD90+ ASC vs CD105-CD90+ ASC DEPs.pdf", width = 4, height = 5)
p.all.asc=pheatmap(expr_imp_asc,
         scale = "row",
         annotation_col = data.frame(condition = group_asc,
                                     row.names = colnames(mat_top50_asc)),
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = FALSE,
         main = "CD105+CD90+ ASC vs CD105-CD90+ ASC",
         border_color="white")
dev.off()


# SOURCE END 208:264

# SOURCE BEGIN 977:1030
# 5.4 火山图 + 热图（CD4）
p.cd4.volcano=ggplot(res_CD4_anno,
       aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.7, size = 0.5) +
  scale_color_manual(values = c("ns" = "grey70", "sig" = "red")) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed") +
  geom_hline(yintercept = -1*log10(0.05), linetype = "dashed") +
  theme_bw() +
  ggtitle("Volcano - CD4_P vs CD4_N")


p.cd4.volcano

ggsave(filename = "CD4 protein Volcano.plot.pdf",
       plot = p.cd4.volcano,
       width = 7,
       height =6 )





top50_cd4 <- res_CD4_anno %>%
  filter(sig == "sig") %>%
  arrange(FDR) %>%
  slice(1:min(50, n())) %>%
  pull(Accession)

mat_top50_cd4 <- expr_imp_cd4[top50_cd4, ]

pdf("CD4_P vs CD4_N top50 DEPs.pdf", width = 4, height = 5)
pheatmap(mat_top50_cd4,
         scale = "row",
         annotation_col = data.frame(condition = group_cd4,
                                     row.names = colnames(mat_top50_cd4)),
         treeheight_row = 0,
         treeheight_col = 0,
        
         show_rownames = FALSE,
         main = "CD4_P vs CD4_N - top50 DEPs")
dev.off()


pdf("CD4_P vs CD4_N DEPs.pdf", width = 4, height = 5)
p.all.asc=pheatmap(expr_imp_cd4,
                   scale = "row",
                   annotation_col = data.frame(condition = group_cd4,
                                               row.names = colnames(mat_top50_cd4)),
                   treeheight_row = 0,
                   treeheight_col = 0,
                   show_rownames = FALSE,
                   main = "CD4_P vs CD4_N",
                   border_color="white")
dev.off()

# SOURCE END 977:1030

# SOURCE BEGIN 147:169
# Principal component analysis by compartment
pca_asc <- prcomp(t(expr_imp_asc), scale. = FALSE)
pca_asc_df <- data.frame(
  PC1 = pca_asc$x[, 1],
  PC2 = pca_asc$x[, 2],
  sample_info_asc
)

p.asc.pca=ggplot(pca_asc_df,
       aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  theme_bw() +
  ggtitle("PCA - ASC only")

p.asc.pca

ggsave(filename = "ASC protein PCA.plot.pdf",
       plot =p.asc.pca,
       width = 8,
       height = 6,
       dpi = 300)


# SOURCE END 147:169

# SOURCE BEGIN 919:939
pca_cd4 <- prcomp(t(expr_imp_cd4), scale. = FALSE)
pca_cd4_df <- data.frame(
  PC1 = pca_cd4$x[, 1],
  PC2 = pca_cd4$x[, 2],
  sample_info_cd4
)

p.cd4.pca=ggplot(pca_cd4_df,
       aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  theme_bw() +
  ggtitle("PCA - CD4 only")

ggsave(filename = "CD4 protein PCA.plot.pdf",
       plot =p.cd4.pca,
       width = 8,
       height = 6,
       dpi = 300)



# SOURCE END 919:939

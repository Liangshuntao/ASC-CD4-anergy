#!/usr/bin/env Rscript
# Figure 4C-D replay from historical baseline inputs and saved differential results.
# Sources: 20250605 CD105/DEG.R (expression/heatmap blocks) and CD105.R (marker block).
# Does not refit DESeq2 or substitute the later 12-sample limma dataset.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || '--help' %in% args) {
  cat('Usage: Rscript plot_baseline_stromal_identity.R --counts BASELINE.csv --results HISTORICAL_RESULTS.csv --markers MARKERS.csv --out NEW_OUTPUT\n')
  quit(status=if ('--help' %in% args) 0L else 2L)
}
value <- function(k) {
  i <- which(args==k)
  if(length(i)!=1L || i==length(args)) stop('Required: ',k)
  args[i+1L]
}
counts <- normalizePath(value('--counts'),mustWork=TRUE)
results <- normalizePath(value('--results'),mustWork=TRUE)
markers <- normalizePath(value('--markers'),mustWork=TRUE)
out <- normalizePath(value('--out'),mustWork=FALSE)
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop('Use new or empty output')
if(tolower(out) %in% tolower(dirname(c(counts,results,markers)))) stop('Keep output separate from inputs')
suppressPackageStartupMessages({library(pheatmap);library(ggplot2)})
rpkm_matrix <- read.csv(counts)
rownames(rpkm_matrix) <- rpkm_matrix[,1]
rpkm_matrix <- rpkm_matrix[,-1]
stopifnot(ncol(rpkm_matrix)==6L)
colnames(rpkm_matrix) <- paste0(rep(c('CD105-MSC','CD105+MSC'),each=3),1:3)
res_df <- read.csv(results,check.names=FALSE)
stopifnot(all(c('gene','pvalue','log2FoldChange') %in% names(res_df)))
# Preserve actual historical pvalue filter (not the contradictory padj comment).
sig_res <- res_df[res_df$pvalue < 0.05 & abs(res_df$log2FoldChange) > 1,]
sig_rpkm_matrix <- rpkm_matrix[rownames(rpkm_matrix) %in% sig_res$gene,]
norm_sig_rpkm_matrix <- t(scale(t(sig_rpkm_matrix)))
apc_MSC <- read.csv(markers)
rownames(apc_MSC) <- apc_MSC[,1]
apc_MSC <- apc_MSC[,-1]
log_data <- log10(apc_MSC+1)
new_col_order <- c('CD105','CD54','CD166','CD73','CD90','CD29','CD44','CD106')
stopifnot(all(new_col_order %in% rownames(log_data)))
mat.1 <- log_data[match(new_col_order,rownames(log_data)),]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
p1 <- pheatmap(norm_sig_rpkm_matrix,scale='none',show_rownames=FALSE,
               show_colnames=TRUE,cluster_rows=TRUE,cluster_cols=TRUE,
               main='Heatmap of Differential Expressed Genes')
ggsave(file.path(out,'Figure4C_baseline_expression_heatmap.pdf'),plot=p1,width=5,height=8)
plot_MSC <- pheatmap(mat.1,cellwidth=20,cellheight=20,cluster_row=FALSE,cluster_cols=FALSE)
ggsave(file.path(out,'Figure4D_baseline_marker_heatmap.pdf'),plot=plot_MSC,width=6,height=8)
saveRDS(list(expression=norm_sig_rpkm_matrix,markers=mat.1),file.path(out,'plotted_matrices.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))

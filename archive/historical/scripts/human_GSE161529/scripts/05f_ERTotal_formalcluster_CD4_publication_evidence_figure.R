#!/usr/bin/env Rscript
# Publication evidence figure linking formal annotation and patient-level CD4 Tconv states to ASC_P-like CAF abundance.
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(ggplot2); library(patchwork); library(scales); library(digest) })
root <- "D:/Users/lstsw/Documents/Adipose_derived stromal_metabolic niches/adipose_stromal_metabolic_niche_lit/public_data_analysis/GSE161529_ASC_P_validation"
input_rds <- file.path(root,"outputs","ERTotal_cluster1_TNK_Harmony_optimized_formal_annotation","ERTotal_cluster1_Harmony_optimized_formal_annotation.rds")
patient_file <- file.path(root,"outputs","ERTotal_CD4_Tconv_Treg_formalcluster_LumB_anchor","ERTotal_formalcluster_patient_ASC_P_CAF_CD4_state_LumB_anchor.tsv")
assoc_file <- file.path(root,"outputs","ERTotal_CD4_Tconv_Treg_formalcluster_LumB_anchor","ERTotal_formalcluster_patient_level_ASC_P_CD4_LumB_associations.tsv")
cell_score_file <- file.path(root,"outputs","ERTotal_CD4_Tconv_Treg_formalcluster_LumB_anchor","ERTotal_formalcluster_CD4_Tconv_cell_state_scores.tsv")
annotation_file <- file.path(root,"tables","ERTotal_cluster1_optimized_formal_annotation.tsv")
output_dir <- file.path(root,"outputs","ERTotal_CD4_Tconv_Treg_formalcluster_LumB_anchor")
script_file <- file.path(root,"scripts","05f_ERTotal_formalcluster_CD4_publication_evidence_figure.R")
stopifnot(all(file.exists(c(input_rds,patient_file,assoc_file,cell_score_file,annotation_file))))
write_tsv <- function(x,n) write.table(x,file.path(output_dir,n),sep="\t",quote=FALSE,row.names=FALSE,na="")
obj <- readRDS(input_rds); tab <- read.delim(patient_file,check.names=FALSE,stringsAsFactors=FALSE); assoc <- read.delim(assoc_file,check.names=FALSE,stringsAsFactors=FALSE); cell_scores <- read.delim(cell_score_file,check.names=FALSE,stringsAsFactors=FALSE); ann <- read.delim(annotation_file,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(tab)==13L, identical(sort(unique(as.character(tab$nearest_centroid_reported))),c("Basal","LumA","LumB")), sum(tab$tconv_n)==1640L)
stopifnot(nrow(cell_scores)==1640L,all(c("cell_id","patient_id","ANERGY_ssGSEA","ACTIVATION_ssGSEA","EXHAUSTION_ssGSEA")%in%colnames(cell_scores)))
centroid_levels <- c("LumA","LumB","Basal")
centroid_cols <- c(LumA="#4C78A8",LumB="#E69F00",Basal="#7B3294")
centroid_shapes <- c(LumA=16,LumB=17,Basal=15)
program_cols <- c(ANERGY="#D55E00",ACTIVATION="#009E73",EXHAUSTION="#CC79A7")
tab$centroid <- factor(tab$nearest_centroid_reported,levels=centroid_levels)
tab$patient_short <- sub("^ERTotal__","",tab$patient_id)
tab$scoreable <- tab$tconv_n>0L
tab$reliable20 <- tab$tconv_n>=20L
tab$zero_ASC_candidate_zero_Tconv <- tab$asc_candidate_n==0L & tab$tconv_n==0L
tab <- tab[order(-as.integer(tab$zero_ASC_candidate_zero_Tconv),tab$asc_candidate_fraction,tab$patient_short),,drop=FALSE]
tab$patient_short <- factor(tab$patient_short,levels=tab$patient_short)
write_tsv(tab,"ERTotal_patient_CD4_Tconv_Fig11_ASC_order.tsv")

group_counts <- do.call(rbind,lapply(centroid_levels,function(g){x<-tab[tab$centroid==g,,drop=FALSE];data.frame(nearest_centroid=g,n_total=nrow(x),n_scoreable=sum(x$scoreable),n_ge20_Tconv=sum(x$reliable20),tconv_cells=sum(x$tconv_n),stringsAsFactors=FALSE)}))
write_tsv(group_counts,"ERTotal_nearest_centroid_group_counts.tsv")
metric_cols <- c(ANERGY="ANERGY_median",ACTIVATION="ACTIVATION_median",EXHAUSTION="EXHAUSTION_median")
long <- do.call(rbind,lapply(names(metric_cols),function(p)data.frame(patient_id=tab$patient_id,patient_short=tab$patient_short,centroid=tab$centroid,program=p,score=tab[[metric_cols[[p]]]],tconv_n=tab$tconv_n,reliable20=tab$reliable20,stringsAsFactors=FALSE)))
long$program <- factor(long$program,levels=rev(names(metric_cols))); long$centroid <- factor(long$centroid,levels=centroid_levels); long$score_z <- ave(long$score,long$program,FUN=function(x){z<-as.numeric(scale(x));z[!is.finite(z)]<-NA_real_;z})
write_tsv(long,"ERTotal_patient_CD4_Tconv_ssGSEA_by_nearest_centroid.tsv")
long$program_distribution <- factor(as.character(long$program),levels=names(metric_cols))
summary_rows <- list()
for(g in centroid_levels) for(p in names(metric_cols)) {
  ii <- long$centroid==g & long$program==p
  x <- long$score[ii]; ok <- is.finite(x)
  summary_rows[[length(summary_rows)+1L]] <- data.frame(
    nearest_centroid=g, program=p, n_scoreable=sum(ok), n_ge20=sum(long$reliable20[ii] & ok),
    median=if(any(ok)) median(x[ok]) else NA_real_,
    Q1=if(any(ok)) unname(quantile(x[ok],.25)) else NA_real_,
    Q3=if(any(ok)) unname(quantile(x[ok],.75)) else NA_real_, stringsAsFactors=FALSE)
}
summary_tab <- do.call(rbind,summary_rows)
write_tsv(summary_tab,"ERTotal_nearest_centroid_CD4_state_descriptive_summary.tsv")

# Cell-resolved heatmap data. Cells are shown descriptively; patient remains the inferential unit.
cell_scores$patient_short <- sub("^ERTotal__","",cell_scores$patient_id)
cell_scores <- cell_scores[order(match(cell_scores$patient_id,tab$patient_id),cell_scores$cell_id),,drop=FALSE]
cell_scores$cell_index <- ave(seq_len(nrow(cell_scores)),cell_scores$patient_id,FUN=seq_along)
cell_metric_cols <- c(ANERGY="ANERGY_ssGSEA",ACTIVATION="ACTIVATION_ssGSEA",EXHAUSTION="EXHAUSTION_ssGSEA")
cell_long <- do.call(rbind,lapply(names(cell_metric_cols),function(p)data.frame(patient_id=cell_scores$patient_id,patient_short=cell_scores$patient_short,cell_id=cell_scores$cell_id,cell_index=cell_scores$cell_index,program=p,score=cell_scores[[cell_metric_cols[[p]]]],stringsAsFactors=FALSE)))
missing_patients <- tab$patient_id[tab$tconv_n==0L]
if(length(missing_patients)) cell_long <- rbind(cell_long,do.call(rbind,lapply(missing_patients,function(pid)data.frame(patient_id=pid,patient_short=sub("^ERTotal__","",pid),cell_id=NA_character_,cell_index=1L,program=names(cell_metric_cols),score=NA_real_,stringsAsFactors=FALSE))))
patient_labels <- setNames(paste0(as.character(tab$patient_short),"|n=",tab$tconv_n),tab$patient_id)
cell_long$patient_label <- factor(patient_labels[cell_long$patient_id],levels=patient_labels[tab$patient_id])
cell_long$program <- factor(cell_long$program,levels=rev(names(cell_metric_cols)))
stopifnot(sum(!is.na(cell_long$cell_id))/3L==1640L,identical(as.integer(table(factor(cell_scores$patient_id,levels=tab$patient_id))),as.integer(tab$tconv_n)))
write_tsv(cell_long,"ERTotal_CD4_Tconv_per_cell_ssGSEA_heatmap_long.tsv")

marker_groups <- list(Cytotoxic_NK_gdT=c("NKG7","GNLY","KLRD1","PRF1","CD8A","CD8B","TRDC"),CD4_lineage=c("CD3D","CD3E","CD3G","CD4"),CD4_Tconv=c("CD40LG","IL7R","CCR7","TCF7","LEF1","KLRB1"),Treg=c("FOXP3","IL2RA","CTLA4","TIGIT","ENTPD1","CXCL13"),Epithelial=c("MUC1","TFF3","ESR1","CA12","AGR3"),Cycling=c("MKI67","TOP2A","TYMS"),Mast=c("TPSAB1","TPSB2","CPA3","MS4A2"))
marker_groups <- lapply(marker_groups,intersect,rownames(obj)); write_tsv(do.call(rbind,lapply(names(marker_groups),function(g)data.frame(marker_group=g,gene=marker_groups[[g]]))),"ERTotal_formalcluster_marker_dotplot_gene_panel.tsv")
cluster_cols <- c("0"="#F8766D","1"="#C49A00","2"="#53B400","3"="#00C094","4"="#00B6EB","5"="#A58AFF","6"="#FB61D7")
cluster_legend_labels <- c("0"="0  Cytotoxic/NK/gdT/CD8 T","1"="1  CD4 Tconv","2"="2  Treg","3"="3  Epithelial-like","4"="4  Cycling","5"="5  Mast","6"="6  Activated cytotoxic T/NK")
pA <- DimPlot(obj,reduction="umap_optimized",group.by="optimized_cluster",label=TRUE,repel=TRUE,label.size=4,pt.size=.22)+coord_fixed()+scale_color_manual(values=cluster_cols,breaks=names(cluster_legend_labels),labels=cluster_legend_labels)+labs(title="A  Formal 7-cluster map",color="Formal annotation")+theme_classic(base_size=10)+theme(plot.title=element_text(face="bold"),legend.position="bottom",legend.text=element_text(size=7))+guides(color=guide_legend(ncol=2,override.aes=list(size=3)))
pB <- DotPlot(obj,features=marker_groups,group.by="optimized_cluster",cols=c("grey92","#B2182B"),dot.scale=5)+RotatedAxis()+labs(title="B  Lineage-marker evidence",x=NULL,y="Optimized cluster")+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),axis.text.x=element_text(angle=55,hjust=1,size=7),legend.position="right")

pC_asc <- ggplot(tab,aes(patient_short,asc_candidate_fraction,fill=centroid))+geom_col(width=.78)+scale_fill_manual(values=centroid_cols,drop=FALSE,name="Nearest centroid")+scale_y_continuous(limits=c(0,max(tab$asc_candidate_fraction)*1.08),breaks=c(0,.1,.3,.5),expand=expansion(mult=c(0,.02)))+labs(title="C  Per-cell CD4 Tconv state landscape",subtitle="1,640 individual cells grouped by patient; patient remains the inferential unit",x=NULL,y="ASC_P-like\nCAF fraction")+theme_classic(base_size=8)+theme(plot.title=element_text(face="bold",size=11),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="right",plot.margin=margin(5,5,8,5))
pC_count <- ggplot(tab,aes(patient_short,tconv_n))+geom_col(fill="grey45",width=.78)+geom_hline(yintercept=20,linetype=2,color="#B2182B",linewidth=.45)+scale_y_continuous(trans=log1p_trans(),breaks=c(0,20,100,500))+labs(x=NULL,y="CD4 Tconv n\n(log1p)")+theme_classic(base_size=8)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=margin(8,5,0,5))
pC_heat <- ggplot(cell_long,aes(cell_index,program,fill=score))+geom_tile(width=1,height=.96)+facet_grid(~patient_label,scales="free_x",labeller=labeller(patient_label=function(x)sub("\\|","\n",x)))+scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,na.value="grey82",name="Per-cell\nssGSEA score")+scale_x_continuous(expand=c(0,0))+labs(x=NULL,y=NULL)+theme_minimal(base_size=8)+theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.text=element_text(size=6.2,face="bold"),strip.background=element_rect(fill="grey94",color="white"),panel.spacing.x=grid::unit(.7,"pt"),plot.margin=margin(0,5,5,5))
pC <- (pC_asc / plot_spacer() / pC_count / pC_heat)+plot_layout(heights=c(.62,.10,.48,1.15),guides="collect")

eligible <- tab[tab$eligible_joint,,drop=FALSE]
lab <- function(x,y){z<-assoc[assoc$x==x&assoc$y==y,,drop=FALSE];if(!nrow(z)||!is.finite(z$rho))"Descriptive only"else sprintf("Spearman rho = %.2f; P = %.3g; n = %d",z$rho,z$p_value,z$n)}
pD <- ggplot(eligible,aes(asc_candidate_fraction,ANERGY_median,color=centroid,shape=centroid))+geom_smooth(data=eligible,aes(x=asc_candidate_fraction,y=ANERGY_median,group=1),inherit.aes=FALSE,method="lm",formula=y~x,se=TRUE,level=.95,color=program_cols[["ANERGY"]],fill=program_cols[["ANERGY"]],alpha=.18,linewidth=.8)+geom_point(size=3.2)+scale_color_manual(values=centroid_cols)+scale_shape_manual(values=centroid_shapes)+labs(title="D  ASC_P-like CAFs and anergy",subtitle=lab("asc_candidate_fraction","ANERGY_median"),x="ASC_P-like candidates / CAF",y="Patient median ANERGY ssGSEA",color="Nearest centroid",shape="Nearest centroid")+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),legend.position="bottom")
pE <- ggplot(eligible,aes(asc_candidate_fraction,ACTIVATION_median,color=centroid,shape=centroid))+geom_smooth(data=eligible,aes(x=asc_candidate_fraction,y=ACTIVATION_median,group=1),inherit.aes=FALSE,method="lm",formula=y~x,se=TRUE,level=.95,color=program_cols[["ACTIVATION"]],fill=program_cols[["ACTIVATION"]],alpha=.18,linewidth=.8)+geom_point(size=3.2)+scale_color_manual(values=centroid_cols)+scale_shape_manual(values=centroid_shapes)+labs(title="E  ASC_P-like CAFs and activation",subtitle=lab("asc_candidate_fraction","ACTIVATION_median"),x="ASC_P-like candidates / CAF",y="Patient median ACTIVATION ssGSEA",color="Nearest centroid",shape="Nearest centroid")+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),legend.position="bottom")
pF <- ggplot(eligible,aes(asc_candidate_fraction,EXHAUSTION_median,color=centroid,shape=centroid))+geom_smooth(data=eligible,aes(x=asc_candidate_fraction,y=EXHAUSTION_median,group=1),inherit.aes=FALSE,method="lm",formula=y~x,se=TRUE,level=.95,color=program_cols[["EXHAUSTION"]],fill=program_cols[["EXHAUSTION"]],alpha=.18,linewidth=.8)+geom_point(size=3.2)+scale_color_manual(values=centroid_cols)+scale_shape_manual(values=centroid_shapes)+labs(title="F  ASC_P-like CAFs and exhaustion",subtitle=lab("asc_candidate_fraction","EXHAUSTION_median"),x="ASC_P-like candidates / CAF",y="Patient median EXHAUSTION ssGSEA",color="Nearest centroid",shape="Nearest centroid")+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),legend.position="bottom")
pG <- ggplot(eligible,aes(asc_candidate_fraction,treg_fraction_of_formal_CD4,color=centroid,shape=centroid))+geom_smooth(data=eligible,aes(x=asc_candidate_fraction,y=treg_fraction_of_formal_CD4,group=1),inherit.aes=FALSE,method="lm",formula=y~x,se=TRUE,level=.95,color="#7B3294",fill="#7B3294",alpha=.18,linewidth=.8)+geom_point(size=3.2)+scale_color_manual(values=centroid_cols)+scale_shape_manual(values=centroid_shapes)+labs(title="G  ASC_P-like CAFs and Treg fraction",subtitle=lab("asc_candidate_fraction","treg_fraction_of_formal_CD4"),x="ASC_P-like candidates / CAF",y="Treg / (Tconv + Treg)",color="Nearest centroid",shape="Nearest centroid")+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),legend.position="bottom")
pH <- ggplot(tab,aes(asc_candidate_fraction,tconv_n,color=centroid,shape=centroid))+geom_hline(yintercept=20,linetype=2,color="#B2182B",linewidth=.5)+geom_point(size=3.2)+scale_color_manual(values=centroid_cols)+scale_shape_manual(values=centroid_shapes)+scale_y_continuous(trans=log1p_trans(),breaks=c(0,1,5,20,100,500))+labs(title="H  Patient support across the ASC_P-like CAF axis",subtitle="All 13 patients; dashed line marks the prespecified 20-cell Tconv threshold",x="ASC_P-like candidates / CAF",y="Formal cluster 1 cells per patient (log1p)",color="Nearest centroid",shape="Nearest centroid")+theme_classic(base_size=9)+theme(plot.title=element_text(face="bold"),legend.position="bottom")

row1 <- (pA|pB) + plot_layout(widths=c(.82,1.45))
row3 <- ((pD|pE|pF|pG)+plot_layout(guides="collect")) & theme(legend.position="bottom")
fig <- (row1 / pC / row3 / pH) + plot_layout(heights=c(1.18,1.35,.84,.70))
png_primary <- file.path(output_dir,"Fig11_CD4_VERIFIED.png")
ggsave(png_primary,fig,width=18,height=19,dpi=320,bg="white",limitsize=FALSE)
pdf_primary <- file.path(output_dir,"Reviewer_Figure_ERTotal_CD4_Fig11_evidence_VERIFIED.pdf")
pdf_fallback <- file.path(output_dir,"Reviewer_Figure_ERTotal_CD4_Fig11_evidence_VERIFIED_v2.pdf")
pdf_written <- tryCatch({ggsave(pdf_primary,fig,width=18,height=19,bg="white",limitsize=FALSE);pdf_primary},error=function(e){ggsave(pdf_fallback,fig,width=18,height=19,bg="white",limitsize=FALSE);pdf_fallback})
writeLines(c("# Figure-generation run report",paste0("Formal RDS SHA256: ",digest(file=input_rds,algo="sha256")),paste0("Patient table SHA256: ",digest(file=patient_file,algo="sha256")),paste0("Cell-score table SHA256: ",digest(file=cell_score_file,algo="sha256")),paste0("Script SHA256: ",digest(file=script_file,algo="sha256")),paste0("PNG written: ",basename(png_primary)),paste0("PDF written: ",basename(pdf_written)),"Statistical unit for association panels: patient; cell-level scores are descriptive only.","Panel C displays all 1,640 formal cluster 1 cells using their raw per-cell ssGSEA scores, grouped by patient; patient labels and a count track report CD4 Tconv cell support.","Patients with zero ASC_P-like candidates and zero CD4 Tconv cells (0040 and 0042) are shown first; remaining patients are ordered by ascending ASC_P-like CAF fraction.","An explicit spacer separates the ASC_P-like CAF fraction and CD4 Tconv count tracks in panel C.","Association panels use the seven patients meeting the prespecified joint threshold of at least 20 formal CD4 Tconv cells and at least 20 CAFs.","Panels D-G show an OLS linear fit with 95% confidence band for visual reference; displayed correlation statistics remain patient-level Spearman results.","The primary figure includes ASC_P-like CAF associations with CD4 Tconv anergy, activation, exhaustion, and formal Treg fraction.","No PAM50 nearest-centroid supplementary distribution figure is generated and no between-centroid hypothesis test is reported."),file.path(output_dir,"ERTotal_CD4_Fig11_evidence_figure_RUN_REPORT.md"))
writeLines(capture.output(sessionInfo()),file.path(output_dir,"ERTotal_CD4_Fig11_evidence_figure_sessionInfo.txt"))

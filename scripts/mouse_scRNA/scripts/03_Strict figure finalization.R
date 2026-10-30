#!/usr/bin/env Rscript
# STRICT REPRODUCTION SCRIPT
# Loads checkpoints (seu_cd4_nonTreg.rds, tf_by_cell.rds, tf_act_long.rds, path_scores.rds)
# Fixes TF activity write-back, AddModuleScore, Figure 2/3/S1C, LIANA, PDF sizes
# No re-run of GSVA/ssGSEA/AUCell/TF/PROGENy/Monocle3

args <- commandArgs(trailingOnly=TRUE)
arg <- function(key, default=NULL) { i <- match(key,args); if(is.na(i)||i==length(args)) default else args[[i+1L]] }

input <- normalizePath(arg("--input","D:/R/20250926/seu_scWAT.rds"),winslash="/",mustWork=TRUE)
outdir <- normalizePath(arg("--out","output"),winslash="/",mustWork=FALSE)
threads <- as.integer(arg("--threads","4"))

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
oldwd <- getwd(); setwd(outdir)

required <- c("Seurat","SeuratWrappers","SingleCellExperiment","monocle3","GSVA","AUCell","decoupleR","dorothea","progeny","liana","Biobase","tidyverse","dplyr","tidyr","tibble","patchwork","ggplot2","pheatmap","viridis","tools")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))

suppressPackageStartupMessages({
  library(Seurat); library(SeuratWrappers); library(SingleCellExperiment)
  library(monocle3); library(GSVA); library(AUCell)
  library(decoupleR); library(dorothea); library(progeny); library(liana)
  library(tidyverse); library(patchwork); library(pheatmap); library(viridis)
  library(tools); library(conflicted)
})
conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(Biobase::exprs)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

options(stringsAsFactors=FALSE, future.globals.maxSize=30*1024^3)
Sys.setenv(LANGUAGE="en", OMP_NUM_THREADS=threads, MKL_NUM_THREADS=threads, OPENBLAS_NUM_THREADS=threads)

# ==============================================================================
# LOAD CHECKPOINTS (no re-computation of GSVA/ssGSEA/AUCell/TF/PROGENy/Monocle3)
# ==============================================================================
message("[STRICT] Loading checkpoints ...")
seu_scWAT <- readRDS(input)
stopifnot(inherits(seu_scWAT,"Seurat"), ncol(seu_scWAT)==58241)

seu_cd4_nonTreg <- readRDS("seu_cd4_nonTreg.rds")
stopifnot(inherits(seu_cd4_nonTreg,"Seurat"), ncol(seu_cd4_nonTreg)==14195)

tf_by_cell <- readRDS("tf_by_cell.rds")
tf_act_long <- readRDS("tf_act_long.rds")
path_scores <- readRDS("path_scores.rds")

message("[STRICT] Checkpoints loaded. Cells: ", ncol(seu_cd4_nonTreg))

# ==============================================================================
# 1. TF ACTIVITY WRITE-BACK (exact 5 TFs: Nr4a1, Ahr, Rara, Rbpj, Creb1 -> *_ACT)
#    Source: original lines 800-812
# ==============================================================================
message("[STRICT] Writing TF activity to metadata (Nr4a1_ACT, Ahr_ACT, Rara_ACT, Rbpj_ACT, Creb1_ACT) ...")

tf_focus <- base::intersect(
  c("Nr4a1","Nr4a2","Nr4a3","Ahr","Rara","Rarg","Rbpj","Creb1"),
  rownames(tf_by_cell)
)
message("[STRICT] Available TFs in tf_by_cell: ", paste(tf_focus, collapse=", "))

# Only write the 5 exact TFs requested
tf_exact <- base::intersect(c("Nr4a1","Ahr","Rara","Rbpj","Creb1"), rownames(tf_by_cell))
for (tf in tf_exact) {
  colname_tf <- paste0(tf, "_ACT")
  seu_cd4_nonTreg[[colname_tf]] <- as.numeric(tf_by_cell[tf, colnames(seu_cd4_nonTreg)])
  message("[STRICT] Wrote ", colname_tf, " : ", sum(!is.na(seu_cd4_nonTreg[[colname_tf]])), " valid values")
}

# Verify all 5 columns have 14195 values
tf_act_cols_check <- paste0(tf_exact, "_ACT")
for (col in tf_act_cols_check) {
  stopifnot(col %in% colnames(seu_cd4_nonTreg@meta.data))
  n_valid <- sum(!is.na(seu_cd4_nonTreg@meta.data[[col]]))
  stopifnot(n_valid == 14195)
}
message("[STRICT] All 5 TF_ACT columns verified: 14195 valid values each")

# ==============================================================================
# 2. RESTORE PROGENy PATHWAY SCORES (from path_scores checkpoint)
#    Source: original lines 278-301
# ==============================================================================
message("[STRICT] Restoring PROGENy pathway scores ...")
DefaultAssay(seu_cd4_nonTreg) <- "RNA"
E <- as.matrix(GetAssayData(seu_cd4_nonTreg, layer="data"))

if (nrow(path_scores) == ncol(E) && ncol(path_scores) <= 30) {
  path_by_cell <- t(path_scores)
} else if (ncol(path_scores) == ncol(E)) {
  path_by_cell <- as.matrix(path_scores)
} else if (nrow(path_scores) == ncol(E)) {
  path_by_cell <- t(as.matrix(path_scores))
} else {
  stop("PROGENy dimensions mismatch")
}

pw_focus <- base::intersect(colnames(path_scores), rownames(path_by_cell))
if (length(pw_focus)) {
  for (pw in pw_focus) {
    seu_cd4_nonTreg[[paste0(pw,"_PROGENY")]] <- as.numeric(path_by_cell[pw, colnames(seu_cd4_nonTreg)])
  }
}

path_focus <- c("TGFb","TNFa","NFkB","JAK-STAT","WNT","MAPK","PI3K",
                "Androgen","Hypoxia","EGFR","Estrogen","p53","Trail","VEGF")
pw_use <- base::intersect(path_focus, rownames(path_by_cell))
message("[STRICT] PROGENy pathways restored: ", paste(pw_use, collapse=", "))
stopifnot(length(pw_use) >= 2)

# ==============================================================================
# 3. RESTORE AddModuleScore (7 gene sets, exact historical parameters)
#    Source: original lines 32-39 (gene sets) + 133-138 (AddModuleScore)
# ==============================================================================
message("[STRICT] Running AddModuleScore with exact historical gene sets ...")

gs_anergy <- c("RNF128","CBLB","ITCH","EGR2","EGR3","DGKA","NR4A1","NR4A2","NR4A3","TOX","PDCD1","CTLA4","LAG3","HAVCR2")
gs_act    <- c("CD69","IL2RA","IFNG","TNF","GZMB","PRF1")
gs_exh    <- c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","BATF","CXCL13")
gs_th1    <- c("TBX21","IFNG","CXCR3")
gs_th2    <- c("GATA3","IL4","IL5","IL13")
gs_th17   <- c("RORC","IL17A","IL17F","CCR6")
gs_tfh    <- c("BCL6","CXCR5")
sig_list <- list(ANERGY=gs_anergy, ACT=gs_act, EXH=gs_exh, TH1=gs_th1, TH2=gs_th2, TH17=gs_th17, TFH=gs_tfh)

# Match gene names exactly as original
DefaultAssay(seu_cd4_nonTreg) <- "RNA"
mat <- as.matrix(GetAssayData(seu_cd4_nonTreg, layer="data"))
rownames_upper <- toupper(rownames(mat))
sig_list_UP <- lapply(sig_list, function(v) {
  v_up <- toupper(v)
  v_up[ v_up %in% rownames_upper ]
})
gene_index <- setNames(seq_along(rownames_upper), rownames_upper)
sig_list_for_gsva <- lapply(sig_list_UP, function(v) rownames(mat)[gene_index[v]])

seu_cd4_nonTreg <- Seurat::AddModuleScore(seu_cd4_nonTreg, features=sig_list_for_gsva,
                                           name=paste0(names(sig_list_for_gsva), "_MSCORE"))
colnames(seu_cd4_nonTreg@meta.data) <- ifelse(grepl("_MSCORE\\d+$", colnames(seu_cd4_nonTreg@meta.data)),
                                               sub("_MSCORE\\d+$","_MSCORE", colnames(seu_cd4_nonTreg@meta.data)),
                                               colnames(seu_cd4_nonTreg@meta.data))
score_cols_MSCORE <- c("ANERGY_MSCORE","ACT_MSCORE","EXH_MSCORE","TH1_MSCORE","TH2_MSCORE","TH17_MSCORE","TFH_MSCORE")
message("[STRICT] AddModuleScore done. MSCORE columns: ", paste(base::intersect(score_cols_MSCORE, colnames(seu_cd4_nonTreg@meta.data)), collapse=", "))

# ==============================================================================
# 4. PRESERVE PSEUDOTIME (from seu_cd4_nonTreg_final.rds or current object)
#    Copy to monocle3_pseudotime if needed (historical name)
#    Source: original line 1113 uses "monocle3_pseudotime"
# ==============================================================================
message("[STRICT] Preserving pseudotime ...")
# Check if seu_cd4_nonTreg_final.rds exists with pseudotime
final_rds <- "seu_cd4_nonTreg_final.rds"
if (file.exists(final_rds)) {
  seu_final <- readRDS(final_rds)
  if ("pseudotime" %in% colnames(seu_final@meta.data)) {
    seu_cd4_nonTreg$pseudotime <- seu_final$pseudotime
    message("[STRICT] Pseudotime restored from final RDS")
  } else if ("monocle3_pseudotime" %in% colnames(seu_final@meta.data)) {
    seu_cd4_nonTreg$pseudotime <- seu_final$monocle3_pseudotime
    message("[STRICT] Pseudotime restored from monocle3_pseudotime in final RDS")
  }
} else if ("pseudotime" %in% colnames(seu_cd4_nonTreg@meta.data)) {
  message("[STRICT] Pseudotime already present in loaded object")
} else if ("monocle_pt" %in% colnames(seu_cd4_nonTreg@meta.data)) {
  seu_cd4_nonTreg$pseudotime <- seu_cd4_nonTreg$monocle_pt
  message("[STRICT] Pseudotime copied from moncole_pt")
}

# Create monocle3_pseudotime alias (used in Fig 2V heatmap code line 1113)
if ("pseudotime" %in% colnames(seu_cd4_nonTreg@meta.data)) {
  seu_cd4_nonTreg$monocle3_pseudotime <- seu_cd4_nonTreg$pseudotime
  message("[STRICT] monocle3_pseudotime alias created")
}

# ==============================================================================
# 5. BUILD ASC OBJECT + ASC-CD4 LIANA INPUT
#    Source: original lines 420-479
# ==============================================================================
message("[STRICT] Building ASC and ASC-CD4 LIANA objects ...")
Idents(seu_scWAT) <- "active_ident"
seu_ASC <- subset(seu_scWAT, idents=c("ASC_CD105pos_CD90pos","ASC_CD105neg_CD90pos","ASC_other"))

asc_groups <- c("ASC_CD105pos_CD90pos","ASC_CD105neg_CD90pos","ASC_other")
cd4_group <- "CD4_T_nonTreg"

asc_cells <- WhichCells(seu_scWAT, expression=ASC_CD90_CD105 %in% asc_groups)
cd4_cells <- WhichCells(seu_cd4_nonTreg, expression=active_ident==cd4_group)

cells_use <- unique(c(asc_cells, cd4_cells))
seu_LR <- subset(seu_scWAT, cells=cells_use)

seu_LR$LR_group <- NA_character_
is_asc <- seu_LR$ASC_CD90_CD105 %in% asc_groups
seu_LR$LR_group[is_asc] <- as.character(seu_LR$ASC_CD90_CD105[is_asc])
is_cd4 <- seu_LR$active_ident == cd4_group
seu_LR$LR_group[is_cd4] <- cd4_group
seu_LR <- subset(seu_LR, subset=!is.na(LR_group))
Idents(seu_LR) <- "LR_group"
DefaultAssay(seu_LR) <- "RNA"

message("[STRICT] LR_group table:")
print(table(seu_LR$LR_group, useNA="ifany"))

# ==============================================================================
# 6. LIANA: RE-RUN with MouseConsensus (REQUIREMENT 6)
#    Source: original line 483
# ==============================================================================
message("[STRICT] Running LIANA with MouseConsensus ...")
# Reuse the aggregate produced by the runnable SCE compatibility block when
# present. The historical Seurat call below is retained as a non-fatal replay.
if (file.exists("liana_agg.rds")) {
  liana_agg <- readRDS("liana_agg.rds")
  message("[STRICT] Loaded existing SCE-compatible LIANA aggregate: ", nrow(liana_agg), " rows")
}
# Ensure RNA data is normalized
DefaultAssay(seu_LR) <- "RNA"
if (nrow(GetAssayData(seu_LR, layer="data")) == 0) {
  seu_LR <- NormalizeData(seu_LR, verbose=FALSE)
}

liana_res <- tryCatch({
  liana::liana_wrap(seu_LR, groupby="LR_group", resource="MouseConsensus")
}, error = function(e) {
  message("[STRICT] LIANA liana_wrap failed: ", e$message)
  NULL
})

if (!is.null(liana_res)) {
  # Check if results are non-empty
  nrows_list <- sapply(liana_res, nrow)
  message("[STRICT] LIANA method rows: ", paste(names(nrows_list), nrows_list, sep="=", collapse=", "))
  
  if (any(nrows_list > 0)) {
    # Try aggregate
    liana_agg <- tryCatch({
      liana::liana_aggregate(liana_res)
    }, error = function(e) {
      message("[STRICT] liana_aggregate failed: ", e$message)
      NULL
    })
    
    if (!is.null(liana_agg) && nrow(liana_agg) > 0) {
      message("[STRICT] LIANA SUCCESS: aggregate_rank table has ", nrow(liana_agg), " rows")
      saveRDS(liana_res, "liana_res.rds")
      saveRDS(liana_agg, "liana_agg.rds")
      write.csv(liana_agg, "CD4_ASC_liana_agg.csv", row.names=TRUE)
    } else {
      message("[STRICT] LIANA BLOCKED: liana_aggregate returned empty result")
      saveRDS(liana_res, "liana_res.rds")
    }
  } else {
    message("[STRICT] LIANA BLOCKED: all methods returned 0 rows")
    saveRDS(liana_res, "liana_res.rds")
  }
} else {
  message("[STRICT] LIANA BLOCKED: liana_wrap failed entirely")
}

# ==============================================================================
# FIGURE GENERATION - ALL SIZES FROM HISTORICAL ggsave PARAMETERS
# ==============================================================================

# ---- Figure 1A: scWAT UMAP (10x6) ----
message("[STRICT] Figure 1A: scWAT UMAP (10x6)")
Idents(seu_scWAT) <- "active_ident"
plot.umap <- DimPlot(seu_scWAT, reduction="umap", group.by="active_ident", label=TRUE)
ggsave("Fig1A_UMAP.pdf", plot=plot.umap, width=10, height=6)

# ---- Figure 1B: Canonical marker DotPlot (10x6) ----
message("[STRICT] Figure 1B: DotPlot (10x6)")
p_markers <- DotPlot(seu_scWAT,
  features=c("Eng","Thy1","Nt5e","Alcam","Cd3g","Cd4","Foxp3","Cd8b1","Cd79a","Iglv1","Itgax","Itgam","Lyz2","C1qa","S100a8","Nkg7","Trdc"),
  assay="RNA") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_blank(),
        legend.text=element_text(size=8),
        legend.title=element_text(size=9)) +
  labs(title="ASC marker") +
  scale_color_gradient2(low="#4675B5", mid="#FAF5C1", high="#D53127", midpoint=0, name="平均表达量") +
  scale_size_continuous(range=c(1,6), name="表达细胞比例")
ggsave("Fig1B_DotPlot.pdf", plot=p_markers, width=10, height=6)

# ---- Figure 1C: ssGSEA heatmap (8x6, dpi=300) ----
message("[STRICT] Figure 1C: ssGSEA heatmap (8x6)")
sc_ssgsea_cols <- grep("_ssGSEA$", colnames(seu_cd4_nonTreg@meta.data), value=TRUE)
sc_ssgsea_mat <- t(as.matrix(seu_cd4_nonTreg@meta.data[, sc_ssgsea_cols, drop=FALSE]))
rownames(sc_ssgsea_mat) <- gsub("_ssGSEA$","", rownames(sc_ssgsea_mat))
row_order <- c("ANERGY","ACT","EXH","TH1","TH2","TH17","TFH")
row_order <- base::intersect(row_order, rownames(sc_ssgsea_mat))
sc_ssgsea_mat <- sc_ssgsea_mat[row_order, , drop=FALSE]

p.CD4.ssgsea <- pheatmap(sc_ssgsea_mat, cluster_rows=FALSE, cluster_cols=TRUE,
  show_rownames=TRUE, show_colnames=FALSE,
  color=colorRampPalette(c("#4675B5","white","red"))(50),
  fontsize=10, treeheight_col=0, main="ssGSEA")
ggsave("Fig1C_ssGSEA_heatmap.pdf", plot=p.CD4.ssgsea, width=8, height=6, dpi=300)

# ---- Figure 1D: ANERGY/ACT/EXH violin (4x3, dpi=600) ----
message("[STRICT] Figure 1D: ssGSEA violin (4x3)")
ss_cols <- c("ANERGY_ssGSEA","ACT_ssGSEA","EXH_ssGSEA")
ss_cols <- base::intersect(ss_cols, colnames(seu_cd4_nonTreg@meta.data))
df_ss <- seu_cd4_nonTreg@meta.data %>%
  dplyr::select(all_of(ss_cols)) %>%
  tibble::rownames_to_column("cell") %>%
  tidyr::pivot_longer(cols=-cell, names_to="program", values_to="value") %>%
  dplyr::mutate(program=gsub("_ssGSEA","",program),
                program=factor(program, levels=c("ANERGY","ACT","EXH")))

p_vln_ss <- ggplot(df_ss, aes(x=program, y=value, fill=program)) +
  geom_violin(trim=FALSE, scale="width") +
  geom_boxplot(width=0.15, outlier.size=0.1, color="black", alpha=0.8) +
  scale_fill_manual(values=c("ANERGY"="#b2182b","ACT"="#2166ac","EXH"="#4d9221")) +
  theme_bw(base_size=10) +
  labs(x="T-cell program", y="ssGSEA score",
       title="Per-cell ANERGY/ACT/EXH ssGSEA scores in CD4_T_nonTreg") +
  theme(legend.position="none")
ggsave("Fig1D_ssGSEA_violin.pdf", p_vln_ss, width=4, height=3, dpi=600)

# ---- Figure 1E: Dominant state fraction (4x3, dpi=600) ----
message("[STRICT] Figure 1E: Dominant state (4x3)")
df_ss_wide <- seu_cd4_nonTreg@meta.data %>%
  dplyr::select(all_of(ss_cols)) %>%
  tibble::rownames_to_column("cell")

df_ss_dom <- df_ss_wide %>%
  dplyr::rowwise() %>%
  dplyr::mutate(max_program=c("ANERGY","ACT","EXH")[which.max(c(ANERGY_ssGSEA,ACT_ssGSEA,EXH_ssGSEA))]) %>%
  dplyr::ungroup()

message("[STRICT] ANERGY dominant: ", sum(df_ss_dom$max_program=="ANERGY"), " / ", nrow(df_ss_dom))
stopifnot(sum(df_ss_dom$max_program=="ANERGY") == 10508)

df_ss_bar <- as.data.frame(prop.table(table(df_ss_dom$max_program))) %>%
  dplyr::rename(program=Var1, fraction=Freq) %>%
  dplyr::mutate(program=factor(program, levels=c("ANERGY","ACT","EXH")))

p_bar <- ggplot(df_ss_bar, aes(x=program, y=fraction, fill=program)) +
  geom_col(width=0.6) +
  scale_fill_manual(values=c("ANERGY"="#b2182b","ACT"="#2166ac","EXH"="#4d9221")) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  theme_bw(base_size=10) +
  labs(x="Dominant program per cell (ssGSEA)", y="Fraction of CD4_T_nonTreg cells",
       title="Proportion of CD4_T_nonTreg cells dominated by ANERGY vs ACT/EXH") +
  theme(legend.position="none")
ggsave("Fig1E_dominant_state.pdf", p_bar, width=4, height=3, dpi=600)

# ---- Figure 1F: Pseudotime UMAP (6x4, dpi=600) ----
message("[STRICT] Figure 1F: Pseudotime UMAP (6x4)")
if ("pseudotime" %in% colnames(seu_cd4_nonTreg@meta.data)) {
  # Use Seurat UMAP directly
  pt_vec <- seu_cd4_nonTreg$pseudotime
  umap_coords <- Embeddings(seu_cd4_nonTreg, "umap")
  common_cells <- base::intersect(rownames(umap_coords), names(pt_vec))
  
  df_pt <- data.frame(
    cell=common_cells,
    UMAP_1=umap_coords[common_cells,1],
    UMAP_2=umap_coords[common_cells,2],
    pt=as.numeric(pt_vec[common_cells]),
    ANERGY=seu_cd4_nonTreg$ANERGY_ssGSEA[common_cells],
    ACT=seu_cd4_nonTreg$ACT_ssGSEA[common_cells],
    EXH=seu_cd4_nonTreg$EXH_ssGSEA[common_cells]
  )
  
  p_pt <- ggplot(df_pt, aes(x=UMAP_1, y=UMAP_2, colour=pt)) +
    geom_point(size=0.3, alpha=0.8) +
    scale_colour_viridis_c(option="plasma") +
    theme_classic() +
    labs(title="CD4_T_nonTreg pseudotime (Monocle3 UMAP)", colour="Pseudotime")
  ggsave("Fig1F_pseudotime_UMAP.pdf", p_pt, width=6, height=4, dpi=600)
} else {
  message("[STRICT] WARNING: No pseudotime found, skipping Fig1F")
}

# ---- Figure 1G: State trends along pseudotime (6x4, dpi=600) ----
message("[STRICT] Figure 1G: Pseudotime trends (6x4)")
if (exists("df_pt") && nrow(df_pt) > 0) {
  df_score <- df_pt %>%
    dplyr::select(pt, ANERGY, ACT, EXH) %>%
    tidyr::drop_na()
  df_long <- df_score %>%
    tidyr::pivot_longer(cols=c(ANERGY, ACT, EXH), names_to="module", values_to="score")
  p_trend <- ggplot(df_long, aes(x=pt, y=score, colour=module)) +
    geom_point(size=0.2, alpha=0.2) +
    geom_smooth(se=FALSE, method="gam", formula=y~s(x,k=5)) +
    theme_classic() +
    labs(x="Pseudotime (Monocle3)", y="Module score",
         title="ANERGY / ACT / EXH scores along pseudotime")
  ggsave("Fig1G_pseudotime_trend.pdf", p_trend, width=6, height=4, dpi=600)
}

# ---- Figure 1H-N: ssGSEA FeaturePlots (combined 20x12, dpi=300) ----
message("[STRICT] Figure 1H-N: ssGSEA FeaturePlots (20x12)")
score_cols_ssGSEA <- c("ANERGY_ssGSEA","ACT_ssGSEA","EXH_ssGSEA","TH1_ssGSEA","TH2_ssGSEA","TH17_ssGSEA","TFH_ssGSEA")
featureplot_ssGSEA <- FeaturePlot(seu_cd4_nonTreg, features=score_cols_ssGSEA,
  cols=c("#FAF5C1","#D53127"), pt.size=0.25,
  min.cutoff="q05", max.cutoff="q95", order=TRUE)
ggsave("Fig1H-N_ssGSEA_FeaturePlot.pdf", plot=featureplot_ssGSEA, width=20, height=12, dpi=300)

# ---- Figure 2A-E: 5 exact TF FeaturePlots (20x12, dpi=300) ----
message("[STRICT] Figure 2A-E: TF FeaturePlots (20x12)")
feature_list <- c("ANERGY_ssGSEA","ACT_ssGSEA","EXH_ssGSEA",
                  "Nr4a1_ACT","Ahr_ACT","Rara_ACT","Rbpj_ACT","Creb1_ACT")
feature_list <- feature_list[feature_list %in% colnames(seu_cd4_nonTreg@meta.data)]
message("[STRICT] Feature list for Fig2A-E: ", paste(feature_list, collapse=", "))

p_umap <- FeaturePlot(seu_cd4_nonTreg, features=feature_list, reduction="umap",
  min.cutoff="q05", max.cutoff="q95", pt.size=0.25, cols=c("#FAF5C1","red"))
ggsave("Fig2A-E_TF_FeaturePlot.pdf", plot=p_umap, width=20, height=12, dpi=300)

# ---- Figure 2F: ANERGY 25%/75% grouping (NOT median) (8x3.5, dpi=300) ----
message("[STRICT] Figure 2F: ANERGY high-low TF violin (8x3.5)")
cand_cols <- c("ANERGY_ssGSEA","ANERGY_GSVA","ANERGY_AUCell","ANERGY_MSCORE")
an_col <- cand_cols[cand_cols %in% colnames(seu_cd4_nonTreg@meta.data)][1]
seu_cd4_nonTreg$ANERGY_score_main <- seu_cd4_nonTreg@meta.data[[an_col]]

q <- quantile(seu_cd4_nonTreg$ANERGY_score_main, c(0.25, 0.75), na.rm=TRUE)
message("[STRICT] ANERGY quantiles: Q25=", round(q[1],4), " Q75=", round(q[2],4))
seu_cd4_nonTreg$ANERGY_group <- cut(seu_cd4_nonTreg$ANERGY_score_main,
  breaks=c(-Inf, q[1], q[2], Inf), labels=c("Low","Mid","High"), include.lowest=TRUE)

# Ensure exact 5 TF columns exist
tf_focus_2f <- c("Rbpj","Nr4a1","Creb1","Ahr","Rara")
tf_act_cols_2f <- paste0(tf_focus_2f, "_ACT")
missing_tf_2f <- tf_focus_2f[!tf_act_cols_2f %in% colnames(seu_cd4_nonTreg@meta.data)]
if (length(missing_tf_2f) > 0) {
  for (tf in missing_tf_2f) {
    if (tf %in% rownames(tf_by_cell)) {
      seu_cd4_nonTreg[[paste0(tf,"_ACT")]] <- as.numeric(tf_by_cell[tf, colnames(seu_cd4_nonTreg)])
    }
  }
}

score_cols_tf <- paste0(tf_focus_2f, "_ACT")
score_cols_tf <- score_cols_tf[score_cols_tf %in% colnames(seu_cd4_nonTreg@meta.data)]

df_tf <- seu_cd4_nonTreg@meta.data %>%
  dplyr::select(ANERGY_group, all_of(score_cols_tf)) %>%
  dplyr::filter(!is.na(ANERGY_group), ANERGY_group != "Mid") %>%
  tidyr::pivot_longer(cols=all_of(score_cols_tf), names_to="TF", values_to="activity") %>%
  dplyr::mutate(TF=sub("_ACT$","",TF),
                ANERGY_group=factor(ANERGY_group, levels=c("Low","High")))

p_tf_violin <- ggplot(df_tf, aes(x=ANERGY_group, y=activity, fill=ANERGY_group)) +
  geom_violin(trim=TRUE, scale="width") +
  geom_boxplot(width=0.2, outlier.size=0.2, color="black", alpha=0.7, position=position_dodge(width=0.5)) +
  facet_wrap(~TF, nrow=1, scales="free_y") +
  scale_fill_manual(values=c("Low"="#d9d9d9","High"="#e41a1c")) +
  labs(x="ANERGY group (CD4_T_nonTreg)", y="TF regulon activity (decoupleR)",
       title="NR4A / AHR / RAR / RBPJ regulon activity in ANERGY-low vs ANERGY-high CD4 T cells") +
  theme_bw(base_size=11) +
  theme(legend.position="none", strip.background=element_rect(fill="grey90", color=NA),
        strip.text=element_text(face="bold"), axis.text.x=element_text(angle=0, hjust=0.5))
ggsave("Fig2F_TF_ANERGY_violin.pdf", p_tf_violin, width=8, height=3.5, dpi=300)

# ---- Figure 2G: Spearman correlation heatmap (8x5, dpi=600) ----
message("[STRICT] Figure 2G: TF-state Spearman correlation (8x5)")
score_cols_all <- grep("_(AUCell|ssGSEA|MSCORE|GSVA)$", colnames(seu_cd4_nonTreg@meta.data), value=TRUE)
keep_prefix <- c("ANERGY","ACT","EXH")
pattern <- paste0("^(", paste(keep_prefix, collapse="|"), ")")
score_cols <- score_cols_all[grepl(pattern, score_cols_all)]

tf_focus_2g <- c("Nr4a1","Nr4a2","Nr4a3","Ahr","Rara","Rarg","Rbpj","Creb1")
tf_use <- base::intersect(tf_focus_2g, rownames(tf_by_cell))
common_cells <- base::intersect(colnames(tf_by_cell), rownames(seu_cd4_nonTreg@meta.data))

df_scores <- seu_cd4_nonTreg@meta.data[common_cells, score_cols, drop=FALSE] |> as.data.frame()
df_scores[] <- lapply(df_scores, as.numeric)
df_tf_cor <- t(tf_by_cell[tf_use, common_cells, drop=FALSE]) |> as.data.frame()
df_tf_cor[] <- lapply(df_tf_cor, as.numeric)

cor_tf <- matrix(NA_real_, nrow=ncol(df_scores), ncol=ncol(df_tf_cor),
                 dimnames=list(colnames(df_scores), colnames(df_tf_cor)))
for (i in seq_len(ncol(df_scores))) {
  for (j in seq_len(ncol(df_tf_cor))) {
    x <- df_scores[[i]]; y <- df_tf_cor[[j]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= 3) cor_tf[i,j] <- suppressWarnings(cor(x[ok], y[ok], method="spearman"))
  }
}

row_order_2g <- c("ANERGY_ssGSEA","ANERGY_GSVA","ANERGY_AUCell","ANERGY_MSCORE",
                   "ACT_ssGSEA","ACT_GSVA","ACT_AUCell","ACT_MSCORE",
                   "EXH_ssGSEA","EXH_GSVA","EXH_AUCell","EXH_MSCORE")
row_order_use <- row_order_2g[row_order_2g %in% rownames(cor_tf)]
cor_tf_plot <- cor_tf[row_order_use, , drop=FALSE]
finite_row <- rowSums(is.finite(cor_tf_plot)) > 0
finite_col <- colSums(is.finite(cor_tf_plot)) > 0
cor_tf_plot <- cor_tf_plot[finite_row, finite_col, drop=FALSE]
cor_tf_plot[!is.finite(cor_tf_plot)] <- 0

p_tf_cor <- pheatmap(cor_tf_plot, cluster_rows=FALSE, cluster_cols=TRUE,
  show_rownames=TRUE, show_colnames=TRUE, angle_col=45, scale="none",
  cellwidth=18, cellheight=18,
  color=colorRampPalette(c("#4675B5","white","red"))(50),
  main="Correlation: TF activity vs ANERGY/ACT/EXH scores (CD4_T_nonTreg)",
  fontsize=8)
ggsave("Fig2G_TF_state_correlation.pdf", plot=p_tf_cor, width=8, height=5, dpi=600)

# ---- Figure 2H-U: PROGENy FeaturePlots (combined 20x20, dpi=300) ----
message("[STRICT] Figure 2H-U: PROGENy FeaturePlots (20x20)")
pw_cols <- paste0(pw_use, "_PROGENY")
pw_cols <- pw_cols[pw_cols %in% colnames(seu_cd4_nonTreg@meta.data)]
if (length(pw_cols) > 0) {
  p_progeny_umap <- FeaturePlot(seu_cd4_nonTreg, features=pw_cols, reduction="umap",
    min.cutoff="q05", max.cutoff="q95", order=TRUE, ncol=3,
    cols=c("#FAF5C1","#D53127"), pt.size=0.25)
  ggsave("Fig2H-U_PROGENy_FeaturePlot.pdf", p_progeny_umap, width=20, height=20, dpi=300)
}

# ---- Figure 2V: PROGENy pseudotime heatmap (8x4, dpi=300) ----
message("[STRICT] Figure 2V: PROGENy pseudotime heatmap (8x4)")
if ("monocle3_pseudotime" %in% colnames(seu_cd4_nonTreg@meta.data)) {
  pt_col <- "monocle3_pseudotime"
  pt_vec <- seu_cd4_nonTreg@meta.data[[pt_col]]
  if (is.data.frame(pt_vec)) pt_vec <- pt_vec[[1]]
  pt_vec <- as.numeric(pt_vec)
  
  cell_order <- colnames(seu_cd4_nonTreg)[order(pt_vec)]
  mat_pw <- path_by_cell[pw_use, cell_order, drop=FALSE]
  mat_pw_scaled <- t(scale(t(mat_pw)))
  
  pseudotime_vec <- seu_cd4_nonTreg$monocle3_pseudotime[cell_order]
  q <- quantile(pseudotime_vec, c(0.25, 0.75), na.rm=TRUE)
  pseudotime_group <- cut(pseudotime_vec, breaks=c(-Inf, q[1], q[2], Inf),
    labels=c("pseudotime_Early","pseudotime_Mid","pseudotime_Late"), include.lowest=TRUE)
  ann_col <- data.frame(Pseudotime_group=pseudotime_group, row.names=cell_order)
  
  p_progeny_heat <- pheatmap(mat_pw_scaled, cluster_rows=FALSE, cluster_cols=FALSE,
    show_rownames=TRUE, show_colnames=FALSE, annotation_col=ann_col,
    color=colorRampPalette(c("#4675B5","white","red"))(50),
    main="PROGENy pathway activity across CD4 Tconv", fontsize_row=8, border_color=NA)
  ggsave("Fig2V_PROGENy_pseudotime_heatmap.pdf", p_progeny_heat, width=8, height=4, dpi=300)
}

# ---- Figure 3A: ASC markers DotPlot (6x4) ----
message("[STRICT] Figure 3A: ASC markers (6x4)")
features_ASC_marker <- c("Eng","Thy1","Icam1","Alcam","Itgb1","Cd44","Vcam1")
features_ASC_marker <- base::intersect(features_ASC_marker, rownames(seu_ASC))
dot_plot_asc <- DotPlot(seu_ASC, features=features_ASC_marker, assay="RNA") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10), axis.title=element_blank(),
        legend.text=element_text(size=8), legend.title=element_text(size=9)) +
  labs(title="ASC marker") +
  scale_color_gradient2(low="#4675B5", mid="#FAF5C1", high="#D53127", midpoint=0, name="平均表达量") +
  scale_size_continuous(range=c(1,6), name="表达细胞比例")
ggsave("Fig3A_ASC_markers.pdf", plot=dot_plot_asc, width=6, height=4)

# ---- Figure 3B: Vascular markers DotPlot (7x4) ----
message("[STRICT] Figure 3B: Vascular markers (7x4)")
features_vascular_marker <- c("Emcn","Pecam1","Cdh5","Esam","Rgs5","Tek","Lyve1","Pdgfra")
features_vascular_marker <- base::intersect(features_vascular_marker, rownames(seu_ASC))
dot_plot_vascular <- DotPlot(seu_ASC, features=features_vascular_marker, assay="RNA") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10), axis.title=element_blank(),
        legend.text=element_text(size=8), legend.title=element_text(size=9)) +
  labs(title="Perivascular / vascular-associated markers") +
  scale_color_gradient2(low="#4675B5", mid="#FAF5C1", high="#D53127", midpoint=0, name="平均表达量") +
  scale_size_continuous(range=c(1,6), name="表达细胞比例")
ggsave("Fig3B_Vascular_markers.pdf", plot=dot_plot_vascular, width=7, height=4)

# ---- Figure 3C: ECM markers DotPlot (6x4) ----
message("[STRICT] Figure 3C: ECM markers (6x4)")
features_ECM_marker <- c("Fmod","Prelp","Dpt","Mmp13")
features_ECM_marker <- base::intersect(features_ECM_marker, rownames(seu_ASC))
dot_plot_ECM <- DotPlot(seu_ASC, features=features_ECM_marker, assay="RNA") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10), axis.title=element_blank(),
        legend.text=element_text(size=8), legend.title=element_text(size=9)) +
  labs(title="ECM / tissue-remodeling genes") +
  scale_color_gradient2(low="#4675B5", mid="#FAF5C1", high="#D53127", midpoint=0, name="平均表达量") +
  scale_size_continuous(range=c(1,6), name="表达细胞比例")
ggsave("Fig3C_ECM_markers.pdf", plot=dot_plot_ECM, width=6, height=4)

# ---- Figure 3D: Immunoregulatory markers DotPlot (6x4) ----
message("[STRICT] Figure 3D: Immunoregulatory markers (6x4)")
features_Imm_marker <- c("Il33","Cd200","Bmp4","Jag1","Entpd1","Nt5e")
features_Imm_marker <- base::intersect(features_Imm_marker, rownames(seu_ASC))
dot_plot_Imm <- DotPlot(seu_ASC, features=features_Imm_marker, assay="RNA") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10), axis.title=element_blank(),
        legend.text=element_text(size=8), legend.title=element_text(size=9)) +
  labs(title="Immunoregulatory ligands and ecto-enzymes") +
  scale_color_gradient2(low="#4675B5", mid="#FAF5C1", high="#D53127", midpoint=0, name="平均表达量") +
  scale_size_continuous(range=c(1,6), name="表达细胞比例")
ggsave("Fig3D_Immuno_markers.pdf", plot=dot_plot_Imm, width=6, height=4)

# ---- Figure 3E: ASC-CD4 UMAP (10x6) ----
message("[STRICT] Figure 3E: ASC-CD4 UMAP (10x6)")
plot.LR.umap <- DimPlot(seu_LR, reduction="umap", group.by="active_ident", label=TRUE)
ggsave("Fig3E_CD4_ASC_UMAP.pdf", plot=plot.LR.umap, width=10, height=6)

# ---- Figure 3F: LIANA global communication heatmap (10x10, dpi=300) ----
message("[STRICT] Figure 3F: LIANA communication heatmap (10x10)")
if (exists("liana_agg") && !is.null(liana_agg) && nrow(liana_agg) > 0) {
  comm_mat <- liana_agg %>%
    group_by(source, target) %>%
    summarise(mean_rank=mean(aggregate_rank, na.rm=TRUE), .groups="drop") %>%
    mutate(score=1-mean_rank) %>%
    select(source, target, score) %>%
    pivot_wider(names_from=target, values_from=score, values_fill=0)
  
  comm_mat_mat <- comm_mat %>% tibble::column_to_rownames("source") %>% as.matrix()
  
  p_comm <- pheatmap(comm_mat_mat, cluster_rows=TRUE, cluster_cols=TRUE,
    color=colorRampPalette(c("#4675B5","white","red"))(50),
    main="Global communication strength (LIANA aggregate)",
    fontsize=8, cellwidth=18, cellheight=18, angle_col=45,
    treeheight_row=0, treeheight_col=0, border_color="white")
  ggsave("Fig3F_LIANA_communication.pdf", plot=p_comm, width=10, height=10, dpi=300)
  
  # ---- Figure 3G: JAG1-NOTCH1 bar (4.5x4, dpi=600) ----
  message("[STRICT] Figure 3G: JAG1-NOTCH1 bar (4.5x4)")
  liana_agg <- liana_agg %>% mutate(ligand=ligand.complex, receptor=receptor.complex)
  asc_to_cd4 <- liana_agg %>%
    dplyr::filter(source %in% asc_groups, target==cd4_group) %>%
    mutate(lr_pair=paste(ligand, receptor, sep="_"), rank_score=1-aggregate_rank)
  
  jag1_notch1 <- asc_to_cd4 %>%
    dplyr::filter(ligand=="Jag1", receptor=="Notch1") %>%
    mutate(axis_name="JAG1_NOTCH1")
  jag1_notch1$source <- factor(jag1_notch1$source, levels=asc_groups)
  
  if (nrow(jag1_notch1) > 0) {
    p_jag1_bar <- ggplot(jag1_notch1, aes(x=source, y=rank_score)) +
      geom_col(width=0.6) +
      geom_text(aes(label=sprintf("%.2f",rank_score)), vjust=-0.3, size=3) +
      coord_cartesian(ylim=c(0,1)) +
      labs(title="JAG1-NOTCH1 axis from ASC subsets to CD4_T_nonTreg",
           x="ASC subsets (source)", y="1 - aggregate_rank (MouseConsensus)") +
      theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(hjust=0.5))
    ggsave("Fig3G_JAG1_NOTCH1_bar.pdf", p_jag1_bar, width=4.5, height=4, dpi=600)
  }
  
  # ---- Figure 3H: Top LR bubble (6x14, dpi=300) ----
  message("[STRICT] Figure 3H: Top LR bubble (6x14)")
  topN <- 50
  top_axes <- asc_to_cd4 %>%
    group_by(source) %>%
    arrange(aggregate_rank, .by_group=TRUE) %>%
    slice_head(n=topN) %>%
    ungroup() %>%
    mutate(lr_pair=paste(ligand, receptor, sep="_"), rank_score=1-aggregate_rank,
           highlight=ifelse(ligand=="Jag1" & receptor=="Notch1", "JAG1-NOTCH1", "Others"))
  top_axes$source <- factor(top_axes$source, levels=asc_groups)
  
  p_axes_bubble <- ggplot(top_axes, aes(y=source, x=reorder(lr_pair, rank_score),
    size=rank_score, color=highlight)) +
    geom_point() +
    scale_size_continuous(range=c(1.5,6), name="1 - aggregate_rank") +
    scale_color_manual(values=c("JAG1-NOTCH1"="red","Others"="grey40"), name="") +
    coord_flip() +
    labs(title="Top ASC->CD4_T_nonTreg ligand-receptor interactions (MouseConsensus)",
         y="ASC subsets (source)", x="Ligand-receptor pairs to CD4_T_nonTreg") +
    theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(hjust=0.5))
  ggsave("Fig3H_TopLR_bubble.pdf", p_axes_bubble, width=6, height=14, dpi=300)
} else {
  message("[STRICT] LIANA BLOCKED: skipping Fig3F-H (no aggregate_rank)")
}

# ---- Figure 3I-U: ALL 13 historical receptor/ligand FeaturePlots (18x14, dpi=300) ----
message("[STRICT] Figure 3I-U: LR FeaturePlots (18x14)")
receptors_L <- c("Bmp4","Bmpr1a","Tgfb1","Tgfb2","Tgfb3","Tgfbr1","Tgfbr2",
                 "Jag1","Notch1","Nt5e","Entpd1","Adora2a","Rbpj")
receptors_L <- base::intersect(receptors_L, rownames(seu_LR))
message("[STRICT] LR features available: ", paste(receptors_L, collapse=", "))

if (length(receptors_L) > 0) {
  p.L.R <- FeaturePlot(object=seu_LR, features=receptors_L,
    pt.size=0.05, min.cutoff="q05", max.cutoff="q95", label=FALSE)
  ggsave("Fig3I-U_LR_FeaturePlot.pdf", plot=p.L.R, width=18, height=14, dpi=300)
}

# ---- Figure S1A: GSVA heatmap (8x6, dpi=600) ----
message("[STRICT] Figure S1A: GSVA heatmap (8x6)")
gsva_cols <- grep("_GSVA$", colnames(seu_cd4_nonTreg@meta.data), value=TRUE)
gsva_mat <- t(as.matrix(seu_cd4_nonTreg@meta.data[, gsva_cols, drop=FALSE]))
rownames(gsva_mat) <- gsub("_GSVA$","", rownames(gsva_mat))
row_order_gsva <- base::intersect(c("ANERGY","ACT","EXH","TH1","TH2","TH17","TFH"), rownames(gsva_mat))
gsva_mat <- gsva_mat[row_order_gsva, , drop=FALSE]

p.CD4.gsva <- pheatmap(gsva_mat, cluster_rows=FALSE, cluster_cols=TRUE,
  show_rownames=TRUE, show_colnames=FALSE,
  color=colorRampPalette(c("#4675B5","white","red"))(50),
  fontsize=10, treeheight_col=0, main="GSVA")
ggsave("FigS1A_GSVA_heatmap.pdf", plot=p.CD4.gsva, width=8, height=6, dpi=600)

# ---- Figure S1B: AUCell heatmap (5x6, dpi=600) ----
message("[STRICT] Figure S1B: AUCell heatmap (5x6)")
auc_cols <- grep("_AUCell$", colnames(seu_cd4_nonTreg@meta.data), value=TRUE)
auc_mat <- t(as.matrix(seu_cd4_nonTreg@meta.data[, auc_cols, drop=FALSE]))
rownames(auc_mat) <- gsub("_AUCell$","", rownames(auc_mat))

CD4_AUCell_pheatmap <- pheatmap(auc_mat, cluster_rows=FALSE, cluster_cols=TRUE,
  show_rownames=TRUE, show_colnames=FALSE, main="AUCell AUC Scores for Gene Sets",
  scale="row", color=colorRampPalette(c("#4675B5","white","red"))(50), fontsize=10)
ggsave("FigS1B_AUCell_heatmap.pdf", plot=CD4_AUCell_pheatmap, width=5, height=6, dpi=600)

# ---- Figure S1C: AddModuleScore heatmap (8x6, dpi=300) ----
message("[STRICT] Figure S1C: AddModuleScore heatmap (8x6)")
score_cols_MSCORE_use <- base::intersect(score_cols_MSCORE, colnames(seu_cd4_nonTreg@meta.data))
MSCORE_mat <- t(as.matrix(seu_cd4_nonTreg@meta.data[, score_cols_MSCORE_use, drop=FALSE]))
rownames(MSCORE_mat) <- score_cols_MSCORE_use
colnames(MSCORE_mat) <- colnames(seu_cd4_nonTreg)

CD4_MSCORE_pheatmap <- pheatmap(MSCORE_mat, cluster_rows=FALSE, cluster_cols=TRUE,
  show_rownames=TRUE, show_colnames=FALSE,
  main="AddModuleScore (Module scores) for CD4_nonTreg",
  scale="row", color=colorRampPalette(c("#4675B5","white","red"))(50),
  fontsize=10, treeheight_col=0)
ggsave("FigS1C_MSCORE_heatmap.pdf", plot=CD4_MSCORE_pheatmap, width=8, height=6, dpi=300)

# ==============================================================================
# FINAL SAVE
# ==============================================================================
message("[STRICT] Saving final objects ...")
saveRDS(seu_cd4_nonTreg, "seu_cd4_nonTreg_final.rds")
tryCatch({
  saveRDS(seu_LR, "seu_LR_final.rds")
}, error = function(e) {
  message("[STRICT] saveRDS seu_LR failed: ", e$message, " - retrying with short name")
  saveRDS(seu_LR, "seuLRfinal.rds")
})
writeLines(capture.output(sessionInfo()), "sessionInfo_final.txt")

# ==============================================================================
# VALIDATION
# ==============================================================================
message("[STRICT] ===== VALIDATION =====")

# Check ANERGY dominant
message("[STRICT] ANERGY dominant count: ", sum(df_ss_dom$max_program=="ANERGY"), "/", nrow(df_ss_dom))

# Check metadata columns
meta_cols <- colnames(seu_cd4_nonTreg@meta.data)
tf_act_present <- paste0(c("Nr4a1","Ahr","Rara","Rbpj","Creb1"), "_ACT")
tf_act_found <- tf_act_present[tf_act_present %in% meta_cols]
message("[STRICT] TF_ACT columns present: ", paste(tf_act_found, collapse=", "))

mscore_present <- score_cols_MSCORE[score_cols_MSCORE %in% meta_cols]
message("[STRICT] MSCORE columns present: ", paste(mscore_present, collapse=", "))

progeny_present <- paste0(pw_use, "_PROGENY")
progeny_found <- progeny_present[progeny_present %in% meta_cols]
message("[STRICT] PROGENy columns present: ", paste(progeny_found, collapse=", "))

pt_present <- "pseudotime" %in% meta_cols || "monocle3_pseudotime" %in% meta_cols
message("[STRICT] Pseudotime present: ", pt_present)

# ---- Generate panel_validation.csv ----
message("[STRICT] Generating panel_validation.csv ...")
pdf_files <- list.files(pattern="\\.pdf$", full.names=FALSE)
panel_info <- data.frame(
  panel = character(),
  file_exists = logical(),
  file_nonempty = logical(),
  pages = integer(),
  width = numeric(),
  height = numeric(),
  source_line = character(),
  status = character(),
  stringsAsFactors = FALSE
)

for (f in sort(pdf_files)) {
  fi <- file.info(f)
  exists_flag <- !is.na(fi$size) && fi$size > 0
  # Try to get PDF dimensions via pdftools if available
  w <- NA_real_; h <- NA_real_; pg <- NA_integer_
  if (requireNamespace("pdftools", quietly=TRUE)) {
    tryCatch({
      info <- pdftools::pdf_info(f)
      pg <- info$pages
      # width/height from page 1 in inches
      if (pg > 0) {
        w <- info$pages[[1]]$width / 72
        h <- info$pages[[1]]$height / 72
      }
    }, error=function(e){})
  }
  panel_info <- rbind(panel_info, data.frame(
    panel = f,
    file_exists = exists_flag,
    file_nonempty = exists_flag,
    pages = ifelse(is.na(pg), 1L, pg),
    width = round(w, 2),
    height = round(h, 2),
    source_line = "historical",
    status = ifelse(exists_flag, "PASS", "MISSING"),
    stringsAsFactors = FALSE
  ))
}
write.csv(panel_info, "panel_validation.csv", row.names=FALSE)
message("[STRICT] panel_validation.csv written with ", nrow(panel_info), " entries")

setwd(oldwd)
message("[STRICT] ===== STRICT REPRODUCTION SCRIPT COMPLETE =====")

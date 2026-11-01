# Explicit execution interface; statistical calculations below are unchanged.
if ("--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat("Usage: Rscript SCRIPT --input ANNOTATED_OBJECT.rds --out OUTPUT_DIRECTORY [--threads 4] [--preflight] [--resume]\n")
  quit(status = 0L)
}
if (!all(c("--input", "--out") %in% commandArgs(trailingOnly = TRUE)))
  stop("Both --input and --out must be supplied explicitly")
#!/usr/bin/env Rscript
# RUNNABLE SUBMISSION VERSION
# Historical analysis and plotting statements are retained from the historical mouse source (see archive provenance).
# Only dependency order and documented execution-interface fixes differ.
args <- commandArgs(trailingOnly=TRUE)
arg <- function(key, default=NULL) { i <- match(key,args); if(is.na(i)||i==length(args)) default else args[[i+1L]] }
input <- normalizePath(arg("--input",NULL),winslash="/",mustWork=TRUE)
outdir <- normalizePath(arg("--out","output"),winslash="/",mustWork=FALSE)
threads <- as.integer(arg("--threads","4"))
preflight <- "--preflight" %in% args
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
required <- c("Seurat","SeuratWrappers","SingleCellExperiment","monocle3","GSVA","AUCell","decoupleR","dorothea","progeny","liana","Biobase","conflicted","tidyverse","dplyr","tidyr","tibble","patchwork","ggplot2","pheatmap","viridis")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing packages (install outside this script): ",paste(missing,collapse=", "))
suppressPackageStartupMessages({library(Seurat);library(SeuratWrappers);library(SingleCellExperiment);library(monocle3);library(GSVA);library(AUCell);library(decoupleR);library(dorothea);library(progeny);library(liana);library(tidyverse);library(patchwork);library(pheatmap);library(viridis);library(conflicted)})
conflicted::conflicts_prefer(base::intersect)
options(stringsAsFactors=FALSE,future.globals.maxSize=30*1024^3)
Sys.setenv(LANGUAGE="en",OMP_NUM_THREADS=threads,MKL_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads)
seu_scWAT <- readRDS(input)
stopifnot(inherits(seu_scWAT,"Seurat"),ncol(seu_scWAT)==58241,length(unique(seu_scWAT$sample_name))==13,"umap" %in% names(seu_scWAT@reductions))
stopifnot(all(c("active_ident","ASC_CD90_CD105","sample_name","intervention_group") %in% colnames(seu_scWAT@meta.data)))
input_counts <- table(seu_scWAT$active_ident)
stopifnot(input_counts[["CD4_T_nonTreg"]]==14195,input_counts[["ASC_CD105pos_CD90pos"]]==277,input_counts[["ASC_CD105neg_CD90pos"]]==1752,input_counts[["ASC_other"]]==4641)
oldwd <- getwd(); setwd(outdir)
write.csv(as.data.frame(input_counts),"input_cell_counts.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"sessionInfo.txt")
if(preflight) { message("PREFLIGHT PASS"); setwd(oldwd); quit(save="no",status=0) }

# ==============================================================================
# FIGURE 1 UPSTREAM | state scores core | original lines 996-1056 | gene sets, CD4 subset, GSVA and ssGSEA
# ==============================================================================
# Define state signatures
gs_anergy <- c("RNF128","CBLB","ITCH","EGR2","EGR3","DGKA","NR4A1","NR4A2","NR4A3","TOX","PDCD1","CTLA4","LAG3","HAVCR2")
gs_act    <- c("CD69","IL2RA","IFNG","TNF","GZMB","PRF1")
gs_exh    <- c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","BATF","CXCL13")
gs_th1    <- c("TBX21","IFNG","CXCR3")
gs_th2    <- c("GATA3","IL4","IL5","IL13")
gs_th17   <- c("RORC","IL17A","IL17F","CCR6")
gs_tfh    <- c("BCL6","CXCR5")
sig_list <- list(ANERGY=gs_anergy, ACT=gs_act, EXH=gs_exh, TH1=gs_th1, TH2=gs_th2, TH17=gs_th17, TFH=gs_tfh)

# Use active_ident as the grouping column
stopifnot("active_ident" %in% colnames(seu_scWAT@meta.data))
seu_scWAT$active_ident <- factor(trimws(seu_scWAT$active_ident))
Idents(seu_scWAT) <- "active_ident"

# Extract log-normalized expression from the data layer
DefaultAssay(seu_scWAT) <- "RNA"

# Input expression should be normalized before signature scoring

# Select CD4_T_nonTreg cells using active_ident
seu_cd4_nonTreg <- subset(seu_scWAT, subset = active_ident == "CD4_T_nonTreg")

# Extract the CD4_T_nonTreg expression matrix

mat<- as.matrix(GetAssayData(seu_cd4_nonTreg, layer = "data"))  # genes x cells



# Match gene symbols in uppercase
rownames_upper <- toupper(rownames(mat))
sig_list_UP <- lapply(sig_list, function(v) {
  v_up <- toupper(v)
  v_up[ v_up %in% rownames_upper ]
})
gene_index <- setNames(seq_along(rownames_upper), rownames_upper)



# Calculate state scores
# 4.1 GSVA / ssGSEA
message("[GSVA] running gsva + ssgsea ...")
sig_list_for_gsva <- lapply(sig_list_UP, function(v) rownames(mat)[gene_index[v]])

sc_gsva   <- GSVA::gsva(mat, sig_list_for_gsva, method="gsva", kcdf="Gaussian", mx.diff=TRUE, parallel.sz=10)

sc_ssgsea <- GSVA::gsva(mat, sig_list_for_gsva, method="ssgsea", kcdf="Gaussian", parallel.sz=10)

for (nm in rownames(sc_gsva)) {
  seu_cd4_nonTreg[[paste0(nm,"_GSVA")]]   <- as.numeric(sc_gsva[nm, colnames(seu_cd4_nonTreg)])
  seu_cd4_nonTreg[[paste0(nm,"_ssGSEA")]] <- as.numeric(sc_ssgsea[nm, colnames(seu_cd4_nonTreg)])
}

saveRDS(seu_cd4_nonTreg, file="seu_cd4_nonTreg.rds") # Save the processed object

seu_cd4_nonTreg=readRDS("seu_cd4_nonTreg.rds" ) 




# Visualize signature scores

# ==============================================================================
# FIGURE 1 UPSTREAM | AUCell | original lines 1165-1196 | AUCell ranking and score columns; reused by Figure S1B
# ==============================================================================
# 4.2 AUCell
message("[AUCell] building rankings + AUC ...")
# Convert signature genes to initial-capital mouse symbols
library(tools)

sig_list_UP_capitalized <- lapply(sig_list_UP, function(v) {
  # Capitalize the first character only
  v_cap <- toTitleCase(tolower(v))
  return(v_cap)
})

# Inspect converted symbols
print(sig_list_UP_capitalized)

# Check the output
# Examples of converted symbols: Rnf128, Cblb and Itch

rankings <- AUCell::AUCell_buildRankings(mat, nCores = max(1, parallel::detectCores()-2), plotStats = FALSE)

aucMaxRank <- ceiling(0.05 * nrow(rankings))  # top 5%

cells_AUC  <- AUCell::AUCell_calcAUC(sig_list_UP_capitalized, rankings, aucMaxRank = aucMaxRank)

auc_mat    <- as.matrix(AUCell::getAUC(cells_AUC))

for (nm in rownames(auc_mat)) {
  seu_cd4_nonTreg[[paste0(nm,"_AUCell")]] <- as.numeric(auc_mat[nm, colnames(seu_cd4_nonTreg)])
}

saveRDS(seu_cd4_nonTreg, file="seu_cd4_nonTreg.rds") # Save the processed object

seu_cd4_nonTreg=readRDS("seu_cd4_nonTreg.rds" ) 

# ==============================================================================
# FIGURE 1 UPSTREAM | AddModuleScore | original lines 1249-1256 | module score columns; reused by Figure S1C
# ==============================================================================
message("[Seurat] AddModuleScore ...")
seu_cd4_nonTreg <- Seurat::AddModuleScore(seu_cd4_nonTreg, features = sig_list_for_gsva,
                                    name = paste0(names(sig_list_for_gsva), "_MSCORE"))
colnames(seu_cd4_nonTreg@meta.data) <- ifelse(grepl("_MSCORE\\d+$", colnames(seu_cd4_nonTreg@meta.data)),
                                        sub("_MSCORE\\d+$","_MSCORE", colnames(seu_cd4_nonTreg@meta.data)),
                                        colnames(seu_cd4_nonTreg@meta.data))

score_cols_MSCORE <- c("ANERGY_MSCORE", "ACT_MSCORE", "EXH_MSCORE", "TH1_MSCORE", "TH2_MSCORE", "TH17_MSCORE", "TFH_MSCORE")

# ==============================================================================
# COMPUTE | TF and PROGENy | original lines 1693-1850 | DoRothEA/decoupleR and PROGENy
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(tibble)
  library(decoupleR)    # TF regulon activity
  library(progeny)      # Pathway activity
  library(pheatmap)
  library(ggplot2); library(patchwork)
})

# Check prerequisites
stopifnot(inherits(seu_cd4_nonTreg, "Seurat"))
DefaultAssay(seu_cd4_nonTreg) <- "RNA"

E <- as.matrix(GetAssayData(seu_cd4_nonTreg, layer = "data"))  # genes x cells
stopifnot(nrow(E) > 100, ncol(E) > 50)

# Construct mouse DoRothEA regulons
# Use the reference supplied by the dorothea package
regulons <- NULL
if (requireNamespace("dorothea", quietly = TRUE)) {
  data("dorothea_mm", package = "dorothea")
  regulons <- dorothea::dorothea_mm %>%
    dplyr::filter(confidence %in% c("A","B","C")) %>%           # Filter by regulon confidence
    dplyr::transmute(
      tf        = as.character(tf),
      target    = as.character(target),
      mor       = ifelse(mor %in% c(1,-1), mor, 1),             # Regulatory sign; replace values outside 1/-1 with 1
      likelihood= if ("likelihood" %in% names(.)) likelihood else 1
    )
} else {
  stop("找不到包 `dorothea`。请先安装：BiocManager::install('dorothea')")
}

# Retain targets present in the expression matrix
regulons <- regulons %>% dplyr::filter(target %in% rownames(E))
stopifnot(nrow(regulons) > 1)


# Calculate TF activity using decoupleR::run_wmean

tf_act_long <- decoupleR::run_wmean(
  mat        = E,
  network    = regulons,
  .source    = "tf",
  .target    = "target",
  .mor       = "mor",
  .likelihood= "likelihood",
  minsize    = 5         # Require at least five matched targets per TF
)

saveRDS(tf_act_long, file="tf_act_long.rds") 

tf_act_long=readRDS("tf_act_long.rds" ) 

# Handle column-name differences across decoupleR versions

score_col     <- c("score","statistic")[c("score","statistic") %in% names(tf_act_long)][1]
condition_col <- c("condition","sample")[c("condition","sample") %in% names(tf_act_long)][1]
stopifnot(length(score_col)==1, length(condition_col)==1)

# Reshape activity values to a TF-by-cell matrix
# Calculate TF activity
tf_act_long <- decoupleR::run_wmean(
  mat        = E,
  network    = regulons,
  .source    = "tf",
  .target    = "target",
  .mor       = "mor",
  .likelihood= "likelihood",
  minsize    = 5
)

# Select one statistic per TF-cell pair
# Inspect available statistics
table(tf_act_long$statistic)
tf_act_long <- tf_act_long %>%
  dplyr::filter(statistic == "corr_wmean")  # Use the selected corr_wmean statistic

# Retain one value per TF-cell pair
score_col     <- c("score","statistic")[c("score","statistic") %in% names(tf_act_long)][1]
condition_col <- c("condition","sample")[c("condition","sample") %in% names(tf_act_long)][1]
stopifnot(length(score_col)==1, length(condition_col)==1)

tf_by_cell <- tf_act_long %>%
  dplyr::select(source, !!rlang::sym(condition_col), !!rlang::sym(score_col)) %>%
  dplyr::rename(
    TF   = source,
    cell = !!rlang::sym(condition_col),
    score= !!rlang::sym(score_col)
  ) %>%
  dplyr::mutate(score = as.numeric(score)) %>%  # Coerce scores to numeric
  tidyr::pivot_wider(
    names_from  = cell,
    values_from = score
  ) %>%
  tibble::column_to_rownames("TF") %>%
  as.matrix()

# Check matrix dimensions
dim(tf_by_cell)
head(rownames(tf_by_cell))
head(colnames(tf_by_cell))

table(tf_act_long$statistic)



saveRDS(tf_by_cell, file="tf_by_cell.rds") 

tf_by_cell=readRDS("tf_by_cell.rds" ) 



# Calculate mouse PROGENy pathway activity
# Returned matrix orientation can vary across versions
path_scores <- progeny(
  E,
  organism  = "Mouse",
  top       = 500,     # Use the top 500 footprint genes per pathway
  perm      = 1,       # Permutation parameter
  z_scores  = TRUE,    # Return z-scores
  scale     = TRUE     # Standardize input expression by gene
)

dim(path_scores)
head(rownames(path_scores))
head(colnames(path_scores))
colnames(path_scores)
saveRDS(path_scores, file="path_scores.rds") 

path_scores=readRDS("path_scores.rds" ) 
 

# Orient the matrix as pathways by cells
if (nrow(path_scores) == ncol(E) && ncol(path_scores) <= 30) {
  # Transpose when the inferred orientation is reversed
  path_by_cell <- t(path_scores)
} else if (ncol(path_scores) == ncol(E)) {
  # Rows represent pathways and columns represent cells
  path_by_cell <- as.matrix(path_scores)
} else if (nrow(path_scores) == ncol(E)) {
  path_by_cell <- t(as.matrix(path_scores))
} else {
  stop("Progeny 返回维度与细胞数不匹配，请 str(path_scores) 检查。")
}


# Store pathway scores in cell metadata
pw_focus <- intersect(colnames(path_scores), rownames(path_by_cell))

 colnames(path_scores)

if (length(pw_focus)) {
  for (pw in pw_focus) {
    seu_cd4_nonTreg[[paste0(pw,"_PROGENY")]] <- as.numeric(path_by_cell[pw, colnames(seu_cd4_nonTreg)])
  }
}


# ==============================================================================
# COMPUTE | PROGENy order | original lines 1958-1966 | historical pathway order used by Fig. 2V
# ==============================================================================
# Select PROGENy pathways for display
path_focus <- c(
  "TGFb","TNFa","NFkB","JAK-STAT","WNT","MAPK","PI3K",
  "Androgen","Hypoxia","EGFR","Estrogen","p53","Trail","VEGF"
)

pw_use <- intersect(path_focus, rownames(path_by_cell))
pw_use
stopifnot(length(pw_use) >= 2)

# ==============================================================================
# COMPUTE | Monocle3 | original lines 3903-3998 | trajectory, historical root rule and pseudotime
# ==============================================================================
# Load dependencies
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)  # as.cell_data_set
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})


# Construct the CD4 lineage object from seu_scWAT

# Set grouping identities from active_ident
stopifnot("active_ident" %in% colnames(seu_scWAT@meta.data))
Idents(seu_scWAT) <- "active_ident"

# Specify CD4 lineage cluster labels

cd4_lineage_levels <- c(
  "CD4_T_nonTreg"
)
# Retain cluster labels present in the object
cd4_lineage_levels <- intersect(cd4_lineage_levels, levels(seu_scWAT$active_ident))
stopifnot("CD4_T_nonTreg" %in% cd4_lineage_levels)

# Subset CD4 lineage cells
cd4_lineage_cells <- colnames(seu_scWAT)[seu_scWAT$active_ident %in% cd4_lineage_levels]
length(cd4_lineage_cells)

seu_cd4_lineage <- subset(seu_scWAT, cells = cd4_lineage_cells)
seu_cd4_lineage$active_ident <- factor(seu_cd4_lineage$active_ident)
Idents(seu_cd4_lineage) <- "active_ident"

# Ensure PCA and UMAP embeddings are available
DefaultAssay(seu_cd4_lineage) <- "RNA"

if (nrow(SeuratObject::GetAssayData(seu_cd4_lineage, assay="RNA", layer="scale.data")) == 0) {
  # Normalize expression if the required layer is absent
  seu_cd4_lineage <- NormalizeData(seu_cd4_lineage, verbose = FALSE)
  seu_cd4_lineage <- FindVariableFeatures(seu_cd4_lineage, verbose = FALSE)
  seu_cd4_lineage <- ScaleData(seu_cd4_lineage, verbose = FALSE)
}

if (!"pca" %in% names(seu_cd4_lineage@reductions)) {
  seu_cd4_lineage <- RunPCA(seu_cd4_lineage, npcs = 30, verbose = FALSE)
}
if (!"umap" %in% names(seu_cd4_lineage@reductions)) {
  seu_cd4_lineage <- RunUMAP(seu_cd4_lineage, dims = 1:20, verbose = FALSE)
}

# Convert Seurat to a monocle3 cell_data_set and learn the trajectory

cds_cd4 <- as.cell_data_set(seu_cd4_lineage)

# Carry active_ident into the trajectory object
colData(cds_cd4)$active_ident <- seu_cd4_lineage$active_ident

# Use the Seurat UMAP embedding in monocle3
cds_cd4@int_colData$reducedDims$UMAP <- 
  seu_cd4_lineage@reductions$umap@cell.embeddings

# Cluster cells and learn the principal graph
cds_cd4 <- cluster_cells(cds_cd4, reduction_method = "UMAP")
cds_cd4 <- learn_graph(cds_cd4)

# Select root cells and calculate pseudotime

cd4_meta_cds <- as.data.frame(colData(cds_cd4))

# Prefer explicit activated or naive CD4 cluster labels for root selection
root_cells <- character(0)

if ("CD4_T_act" %in% cd4_meta_cds$active_ident) {
  root_cells <- rownames(cd4_meta_cds)[cd4_meta_cds$active_ident == "CD4_T_act"]
} else if ("CD4_T_naive" %in% cd4_meta_cds$active_ident) {
  root_cells <- rownames(cd4_meta_cds)[cd4_meta_cds$active_ident == "CD4_T_naive"]
}

# Otherwise select high-ACT, low-ANERGY cells as roots
if (length(root_cells) == 0) {
  # This fallback requires ACT_GSVA and ANERGY_GSVA scores
  if (all(c("ACT_GSVA","ANERGY_GSVA") %in% colnames(seu_cd4_nonTreg@meta.data))) {
    tmp <- seu_cd4_nonTreg@meta.data
    tmp$cell <- rownames(tmp)
    tmp$ACT_minus_ANERGY <- tmp$ACT_GSVA - tmp$ANERGY_GSVA
    thr <- quantile(tmp$ACT_minus_ANERGY, 0.95, na.rm = TRUE)
    root_from_nt <- tmp$cell[tmp$ACT_minus_ANERGY >= thr]
    root_cells <- intersect(root_from_nt, colnames(cds_cd4))
  }
}

stopifnot("没有找到合适的 root 细胞，请检查 active_ident 或 ANERGY/ACT 分数。" = length(root_cells) > 0)

# Calculate pseudotime
cds_cd4 <- order_cells(cds_cd4, root_cells = root_cells)
pt_all  <- pseudotime(cds_cd4)   # Named numeric vector：names = cell barcodes

# ==============================================================================
# COMPUTE | ASC object | original lines 633-635 | three historical ASC labels
# ==============================================================================
seu_ASC <- subset(seu_scWAT, idents = c("ASC_CD105pos_CD90pos",
                                        "ASC_CD105neg_CD90pos",
                                        "ASC_other"))

# ==============================================================================
# COMPUTE | ASC-CD4 LIANA input | original lines 2740-2792 | historical ASC-CD4 subset and LR_group
# ==============================================================================
# Inspect metadata columns
colnames(seu_scWAT@meta.data)
table(seu_scWAT$ASC_CD90_CD105)
table(seu_scWAT$active_ident)

# Define the retained cell populations
# Include the selected ASC groups and CD4_T_nonTreg cells
asc_groups <- c("ASC_CD105pos_CD90pos", "ASC_CD105neg_CD90pos","ASC_other")  
cd4_group  <- "CD4_T_nonTreg"

# Identify the selected cells in the original object
# Select ASC cells using ASC_CD90_CD105
asc_cells <- WhichCells(
  seu_scWAT,
  expression = ASC_CD90_CD105 %in% asc_groups
)

# Select CD4_T_nonTreg cells using active_ident
cd4_cells <- WhichCells(
  seu_cd4_nonTreg,
  expression = active_ident == cd4_group
)

length(asc_cells); length(cd4_cells)

# Combine cell names and construct the subset
cells_use <- unique(c(asc_cells, cd4_cells))

seu_LR <- subset(seu_scWAT, cells = cells_use)
seu_LR

# Define ligand-receptor analysis groups
seu_LR$LR_group <- NA_character_

# Use ASC_CD90_CD105 labels for ASC groups
is_asc <- seu_LR$ASC_CD90_CD105 %in% asc_groups
seu_LR$LR_group[is_asc] <- as.character(seu_LR$ASC_CD90_CD105[is_asc])

# Use active_ident for the CD4_T_nonTreg group
is_cd4 <- seu_LR$active_ident == cd4_group
seu_LR$LR_group[is_cd4] <- cd4_group

table(seu_LR$LR_group, useNA = "ifany")

# Remove cells without an LR_group label
seu_LR <- subset(seu_LR, subset = !is.na(LR_group))

# Set active identities
Idents(seu_LR) <- "LR_group"

# Ensure the RNA assay has normalized expression
DefaultAssay(seu_LR) <- "RNA"
# Use seu_LR for the selected ASC and CD4 ligand-receptor analyses

# [INTERFACE FIX LIANA-02] SeuratObject 5 compatibility: pass the unchanged
# expression matrix and LR_group labels through SingleCellExperiment.
conflicted::conflicts_prefer(Biobase::exprs)
sce_LR <- as.SingleCellExperiment(seu_LR, assay = "RNA")
SingleCellExperiment::colLabels(sce_LR) <- factor(
  SingleCellExperiment::colData(sce_LR)$LR_group,
  levels = c("ASC_CD105pos_CD90pos", "ASC_CD105neg_CD90pos", "ASC_other", "CD4_T_nonTreg")
)
liana_methods <- c("sca", "cellphonedb", "natmi", "connectome", "logfc")
liana_res <- setNames(lapply(liana_methods, function(method_name) {
  liana::liana_wrap(
    sce_LR,
    method = method_name,
    resource = "MouseConsensus",
    idents_col = "LR_group",
    min_cells = 10,
    .simplify = TRUE
  )
}), liana_methods)
liana_agg <- liana::liana_aggregate(liana_res)
saveRDS(liana_res,"liana_res.Rds")
saveRDS(liana_agg,"liana_agg.rds")
write.csv(liana_agg,"CD4_ASC_liana_agg.csv",row.names=TRUE)

# ==============================================================================
# FIGURE 1 | A | original lines 95-97 | whole scWAT UMAP
# ==============================================================================
plot.umap = DimPlot(seu_scWAT, reduction = "umap",group.by = "active_ident",label = TRUE)
plot.umap
ggsave("UMAP.pdf", plot = plot.umap, width = 10, height = 6)

# ==============================================================================
# FIGURE 1 | B | original lines 198-228 | canonical marker DotPlot
# ==============================================================================
p_markers=DotPlot(seu_scWAT, 
                  features = c("Eng","Thy1","Nt5e","Alcam","Cd3g","Cd4","Foxp3","Cd8b1","Cd79a","Iglv1","Itgax","Itgam","Lyz2","C1qa","S100a8","Nkg7","Trdc"),
                  assay = "RNA"
             ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) +
  labs(
    title = "ASC marker"
  ) +
  scale_color_gradient2(
    low = "#4675B5", 
    mid = "#FAF5C1", 
    high = "#D53127",
    midpoint = 0,
    name = "平均表达量"
  ) +
  scale_size_continuous(
    range = c(1, 6),
    name = "表达细胞比例"
  )



p_markers

ggsave("all_cluster_marker.dot_plot.pdf", plot = p_markers, width = 10, height = 6)

# ==============================================================================
# FIGURE 1 | C | original lines 1058-1086 | ssGSEA heatmap
# ==============================================================================
library(pheatmap)
p.CD4.ssgsea=pheatmap(sc_ssgsea,cluster_rows = T, cluster_cols = T, show_rownames = T, show_colnames = F,main = "ssGSEA")
p.CD4.ssgsea
p.CD4.gsva=pheatmap(sc_gsva, cluster_rows = T, cluster_cols = T, show_rownames = T, show_colnames = F,main = "GSVA")

# Set display order
# Specify row order
row_order <- c("ANERGY",  "ACT","EXH","TH1", "TH2", "TH17", "TFH")

# Align rows with the requested order
sc_ssgsea <- sc_ssgsea[row_order, , drop = FALSE]

# Draw the ordered heatmap
p.CD4.ssgsea <- pheatmap(
  sc_ssgsea, 
  cluster_rows = F,     # Cluster rows
  cluster_cols = TRUE,     # Cluster columns
  show_rownames = TRUE,    # Show row labels
  show_colnames = FALSE,   # Hide column labels
  color = colorRampPalette(c("#4675B5",  "white", "red"))(50),  # Color scale
  fontsize = 10,
  treeheight_col = 0,
  main = "ssGSEA"          # Heatmap title
)
# Display the heatmap
p.CD4.ssgsea

# Save the heatmap as PDF
ggsave("CD4_ssGSEA_heatmap.pdf", plot = p.CD4.ssgsea, width = 8, height = 6, dpi = 300)

# ==============================================================================
# FIGURE 1 | D | original lines 1425-1456 | ANERGY ACT EXH violin
# ==============================================================================
# Select ssGSEA score columns
ss_cols <- c("ANERGY_ssGSEA","ACT_ssGSEA","EXH_ssGSEA")
ss_cols <- intersect(ss_cols, colnames(seu_cd4_nonTreg@meta.data))

df_ss <- seu_cd4_nonTreg@meta.data %>%
  dplyr::select(all_of(ss_cols)) %>%
  tibble::rownames_to_column("cell") %>%
  tidyr::pivot_longer(cols = -cell,
                      names_to = "program",
                      values_to = "value") %>%
  dplyr::mutate(
    program = gsub("_ssGSEA","", program),
    program = factor(program, levels = c("ANERGY","ACT","EXH"))
  )

p_vln_ss <- ggplot(df_ss, aes(x = program, y = value, fill = program)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.1, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("ANERGY" = "#b2182b",
                               "ACT"    = "#2166ac",
                               "EXH"    = "#4d9221")) +
  theme_bw(base_size = 10) +
  labs(
    x = "T-cell program",
    y = "ssGSEA score",
    title = "Per-cell ANERGY/ACT/EXH ssGSEA scores in CD4_T_nonTreg"
  ) +
  theme(legend.position = "none")

p_vln_ss
ggsave("CD4_ANERGY_ACT_EXH_ssGSEA_violin.pdf",
       p_vln_ss, width = 4, height = 3, dpi = 600)

# ==============================================================================
# FIGURE 1 | E | original lines 1458-1502 | dominant state fraction
# ==============================================================================
# Calculate cell-state proportions

# Use ssGSEA scores for dominant-state assignment
ss_cols <- c("ANERGY_ssGSEA","ACT_ssGSEA","EXH_ssGSEA")
ss_cols <- intersect(ss_cols, colnames(seu_cd4_nonTreg@meta.data))

df_ss_wide <- seu_cd4_nonTreg@meta.data %>%
  dplyr::select(all_of(ss_cols)) %>%
  tibble::rownames_to_column("cell")

# Assign the highest-scoring program to each cell
df_ss_dom <- df_ss_wide %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    max_program = c("ANERGY","ACT","EXH")[
      which.max(c(ANERGY_ssGSEA, ACT_ssGSEA, EXH_ssGSEA))
    ]
  ) %>%
  dplyr::ungroup()

table(df_ss_dom$max_program)
prop.table(table(df_ss_dom$max_program))

# Plot state proportions as bars
df_ss_bar <- as.data.frame(prop.table(table(df_ss_dom$max_program))) %>%
  dplyr::rename(program = Var1, fraction = Freq) %>%
  dplyr::mutate(program = factor(program, levels = c("ANERGY","ACT","EXH")))

p_bar <- ggplot(df_ss_bar, aes(x = program, y = fraction, fill = program)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("ANERGY" = "#b2182b",
                               "ACT"    = "#2166ac",
                               "EXH"    = "#4d9221")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 10) +
  labs(
    x = "Dominant program per cell (ssGSEA)",
    y = "Fraction of CD4_T_nonTreg cells",
    title = "Proportion of CD4_T_nonTreg cells dominated by ANERGY vs ACT/EXH"
  ) +
  theme(legend.position = "none")

p_bar
ggsave("CD4_ANERGY_ACT_EXH_dominant_fraction_ssGSEA.pdf",
       p_bar, width = 4, height = 3, dpi = 600)

# ==============================================================================
# FIGURE 1 | F | original lines 4344-4392 | Monocle3 pseudotime UMAP
# ==============================================================================
# Extract UMAP and pseudotime and align cells with Seurat

library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)

# Retrieve UMAP coordinates and pseudotime from monocle3
umap_cds <- reducedDims(cds_cd4)$UMAP  # matrix: cells x 2
pt_vec   <- pseudotime(cds_cd4)       # Named numeric vector

# Align by cell name
common_cells <- intersect(rownames(umap_cds), colnames(seu_cd4_nonTreg))
length(common_cells)   # Check the CD4_T_nonTreg cell count

umap_use <- umap_cds[common_cells, , drop = FALSE]
pt_use   <- as.numeric(pt_vec[common_cells])

# Store pseudotime in the Seurat object
seu_cd4_nonTreg$monocle_pt <- NA_real_
seu_cd4_nonTreg$monocle_pt[common_cells] <- pt_use

# Plot UMAP colored by pseudotime

# Construct the ggplot data frame
df_pt <- data.frame(
  cell    = common_cells,
  UMAP_1  = umap_use[common_cells, 1],
  UMAP_2  = umap_use[common_cells, 2],
  pt      = pt_use,
  ident   = seu_cd4_nonTreg$active_ident[common_cells],
  ANERGY  = seu_cd4_nonTreg$ANERGY_ssGSEA[common_cells],
  ACT     = seu_cd4_nonTreg$ACT_ssGSEA[common_cells],
  EXH     = seu_cd4_nonTreg$EXH_ssGSEA[common_cells]
)

# Color UMAP points by pseudotime
p_pt <- ggplot(df_pt, aes(x = UMAP_1, y = UMAP_2, colour = pt)) +
  geom_point(size = 0.3, alpha = 0.8) +
  scale_colour_viridis_c(option = "plasma") +
  theme_classic() +
  labs(
    title = "CD4_T_nonTreg pseudotime (Monocle3 UMAP)",
    colour = "Pseudotime"
  )

p_pt
ggsave("CD4_T_nonTreg_UMAP_pseudotime_manual.pdf",
       p_pt, width = 6, height = 4, dpi = 600)

# ==============================================================================
# FIGURE 1 | G | original lines 4394-4422 | state trends along pseudotime
# ==============================================================================
# Examine ANERGY, ACT and EXH scores along pseudotime

# Retain cells with complete scores
df_score <- df_pt %>%
  select(pt, ANERGY, ACT, EXH) %>%
  tidyr::drop_na()

# Reshape scores to long format for faceted plots
df_long <- df_score %>%
  tidyr::pivot_longer(
    cols      = c(ANERGY, ACT, EXH),
    names_to  = "module",
    values_to = "score"
  )

# Fit smooth score trends along pseudotime
p_trend <- ggplot(df_long, aes(x = pt, y = score, colour = module)) +
  geom_point(size = 0.2, alpha = 0.2) +
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, k = 5)) +
  theme_classic() +
  labs(
    x = "Pseudotime (Monocle3)",
    y = "Module score",
    title = "ANERGY / ACT / EXH scores along pseudotime"
  )

p_trend
ggsave("CD4_T_nonTreg_module_scores_vs_pseudotime.pdf",
       p_trend, width = 6, height = 4, dpi = 600)

# ==============================================================================
# FIGURE 1 | H-N | original lines 1102-1145 | historical combined ssGSEA FeaturePlot
# ==============================================================================
# Plot GSVA and ssGSEA scores as violin plots
library(Seurat)
library(ggplot2)

score_cols <- c("ANERGY_GSVA", "ACT_GSVA", "EXH_GSVA", "TH1_GSVA", "TH2_GSVA", "TH17_GSVA", "TFH_GSVA",
                "ANERGY_ssGSEA", "ACT_ssGSEA", "EXH_ssGSEA", "TH1_ssGSEA", "TH2_ssGSEA", "TH17_ssGSEA", "TFH_ssGSEA")

score_cols_GSVA <- c("ANERGY_GSVA", "ACT_GSVA", "EXH_GSVA", "TH1_GSVA", "TH2_GSVA", "TH17_GSVA", "TFH_GSVA")

score_cols_ssGSEA <- c("ANERGY_ssGSEA", "ACT_ssGSEA", "EXH_ssGSEA", "TH1_ssGSEA", "TH2_ssGSEA", "TH17_ssGSEA", "TFH_ssGSEA")

# Compute the correlation matrix
score_data <- seu_cd4_nonTreg@meta.data[, score_cols]
corr_matrix <- cor(score_data, method = "spearman", use = "pairwise.complete.obs")

# Plot the correlation heatmap
library(pheatmap)
p_corr_GSVA_ssGSEA=pheatmap(corr_matrix, 
                            color = colorRampPalette(c("#4675B5",  "white", "red"))(50),  # Color scale
                            fontsize = 10,
                            treeheight_col = 0,
                            treeheight_row  = 0,
                            main = "Correlation of Functional Scores in CD4_T_nonTreg", 
                            
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean")

p_corr_GSVA_ssGSEA

ggsave("CD4_corr_GSVA_ssGSEA_heatmap.pdf", plot = p_corr_GSVA_ssGSEA, width = 6, height = 6, dpi = 600)

# Plot the existing CD4_T_nonTreg UMAP
featureplot_ssGSEA=FeaturePlot(seu_cd4_nonTreg, 
                               features = score_cols_ssGSEA, 
                               cols = c("#FAF5C1", "#D53127"), 
                               pt.size = 0.25, 
                               
                               min.cutoff = "q05", 
                               max.cutoff = "q95", 
                               order = TRUE)
featureplot_ssGSEA
ggsave("CD4_ssGSEA_FeaturePlot.pdf", 
       plot = featureplot_ssGSEA,
       width = 20, height = 12, dpi = 300)

# ==============================================================================
# FIGURE 2 | A-E | original lines 2255-2313 | historical combined TF FeaturePlot
# ==============================================================================
# Display ANERGY, selected TF activities and NOTCH_PROGENY on UMAP

# UMAP overlays for state and activity scores

library(Seurat)
library(ggplot2)
library(patchwork)
library(tibble)
library(dplyr)

DefaultAssay(seu_cd4_nonTreg) <- "RNA"

# Store the selected TF activities in metadata

tf_focus <- intersect(
  c("Nr4a1","Nr4a2","Nr4a3","Ahr","Rara","Rarg","Rbpj","Creb1"),
  rownames(tf_by_cell)
)

if (length(tf_focus) > 0) {
  for (tf in tf_focus) {
    colname_tf <- paste0(tf, "_ACT")
    if (!colname_tf %in% colnames(seu_cd4_nonTreg@meta.data)) {
      seu_cd4_nonTreg[[colname_tf]] <- as.numeric(tf_by_cell[tf, colnames(seu_cd4_nonTreg)])
    }
  }
}


# Specify features to display

feature_list <- c(
  # State-score columns
  "ANERGY_ssGSEA", "ACT_ssGSEA", "EXH_ssGSEA",
  # Selected TF activity columns
  "Nr4a1_ACT",
  "Ahr_ACT","Rara_ACT","Rbpj_ACT","Creb1_ACT"
)

feature_list <- feature_list[feature_list %in% colnames(seu_cd4_nonTreg@meta.data)]
feature_list

# Plot UMAP feature maps

p_umap <- FeaturePlot(
  seu_cd4_nonTreg,
  features   = feature_list,
  reduction  = "umap",
  min.cutoff = "q05",
  max.cutoff = "q95",
  pt.size = 0.25, 
  cols       = c("#FAF5C1","red")
) 
p_umap

ggsave(
  filename = "CD4_T_nonTreg_UMAP_ANERGY_TF_NOTCH.pdf",
  plot     = p_umap,
  width    = 20, height = 12, dpi = 300)

# ==============================================================================
# FIGURE 2 | F | original lines 2315-2438 | ANERGY high-low TF activity
# ==============================================================================
# Compare TF activities by anergy-score group
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# Select the ANERGY score used for stratification
cand_cols <- c("ANERGY_ssGSEA", "ANERGY_GSVA", "ANERGY_AUCell", "ANERGY_MSCORE")
an_col    <- cand_cols[cand_cols %in% colnames(seu_cd4_nonTreg@meta.data)][1]
stopifnot(length(an_col) == 1)

seu_cd4_nonTreg$ANERGY_score_main <- seu_cd4_nonTreg@meta.data[[an_col]]

# Define Low, Mid and High groups using score quartiles
q <- quantile(seu_cd4_nonTreg$ANERGY_score_main, c(0.25, 0.75), na.rm = TRUE)
seu_cd4_nonTreg$ANERGY_group <- cut(
  seu_cd4_nonTreg$ANERGY_score_main,
  breaks = c(-Inf, q[1], q[2], Inf),
  labels = c("Low", "Mid", "High"),
  include.lowest = TRUE
)


# Selected transcription factors
tf_focus <- c("Rbpj", "Nr4a1", "Creb1", "Ahr", "Rara")

# Identify existing activity columns and those to import from tf_by_cell
tf_act_cols <- paste0(tf_focus, "_ACT")
missing_tf  <- tf_focus[!tf_act_cols %in% colnames(seu_cd4_nonTreg@meta.data)]

if (length(missing_tf) > 0) {
  stopifnot(exists("tf_by_cell"))  # Use the activity matrix generated by DoRothEA/decoupleR
  
  common_cells <- intersect(colnames(tf_by_cell), colnames(seu_cd4_nonTreg))
  
  for (tf in missing_tf) {
    if (!tf %in% rownames(tf_by_cell)) next
    
    v <- as.numeric(tf_by_cell[tf, common_cells])
    
    # Align activity values to Seurat cell order
    tmp <- rep(NA_real_, ncol(seu_cd4_nonTreg))
    names(tmp) <- colnames(seu_cd4_nonTreg)
    tmp[common_cells] <- v
    
    seu_cd4_nonTreg[[paste0(tf, "_ACT")]] <- tmp
  }
}


# Retain available TF columns
score_cols_tf <- paste0(tf_focus, "_ACT")
score_cols_tf <- score_cols_tf[score_cols_tf %in% colnames(seu_cd4_nonTreg@meta.data)]

# Compare ANERGY Low and High groups, excluding Mid
df_tf <- seu_cd4_nonTreg@meta.data %>%
  dplyr::select(ANERGY_group, all_of(score_cols_tf)) %>%
  dplyr::filter(!is.na(ANERGY_group), ANERGY_group != "Mid") %>%
  tidyr::pivot_longer(
    cols      = all_of(score_cols_tf),
    names_to  = "TF",
    values_to = "activity"
  ) %>%
  dplyr::mutate(
    TF           = sub("_ACT$", "", TF),
    ANERGY_group = factor(ANERGY_group, levels = c("Low", "High"))
  )

# Plot each TF in a separate comparison panel
p_tf_violin <- ggplot(df_tf, aes(x = ANERGY_group, y = activity, fill = ANERGY_group)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.2, outlier.size = 0.2, color = "black",
               alpha = 0.7, position = position_dodge(width = 0.5)) +
  facet_wrap(~ TF, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Low" = "#d9d9d9", "High" = "#e41a1c")) +
  labs(
    x     = "ANERGY group (CD4_T_nonTreg)",
    y     = "TF regulon activity (decoupleR)",
    title = "NR4A / AHR / RAR / RBPJ regulon activity in ANERGY-low vs ANERGY-high CD4 T cells"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_tf_violin

ggsave("Fig2A_TF_activity_ANERGY_high_vs_low.pdf",
       p_tf_violin, width = 8, height = 3.5, dpi = 300)


# Calculate TF activity comparison P values
library(dplyr)

# TF activity columns
tf_features <- c("Rbpj_ACT","Nr4a1_ACT","Creb1_ACT","Ahr_ACT","Rara_ACT")

meta <- seu_cd4_nonTreg@meta.data

# Retain cells in the Low or High group
df_tf <- meta %>%
  filter(!is.na(ANERGY_group)) %>%
  select(ANERGY_group, all_of(tf_features))

# Apply Wilcoxon rank-sum tests to High versus Low groups
tf_pvals <- sapply(tf_features, function(tf) {
  g1 <- df_tf %>% filter(ANERGY_group == "Low")  %>% pull(tf)
  g2 <- df_tf %>% filter(ANERGY_group == "High") %>% pull(tf)
  wilcox.test(g2, g1)$p.value   # High versus Low comparison
})

tf_padj <- p.adjust(tf_pvals, method = "BH")

tf_stats <- data.frame(
  TF    = tf_features,
  p_val = tf_pvals,
  p_adj = tf_padj,
  stringsAsFactors = FALSE
)
tf_stats
write.csv(x =tf_stats,file = "TF_stats.csv",row.names = T)

# ==============================================================================
# FIGURE 2 | G | original lines 2142-2253 | TF-state correlation
# ==============================================================================
# Spearman correlations between state scores and TF activities
# Identify state-score columns
# Collect AUCell, GSVA, ssGSEA and MSCORE columns
score_cols_all <- grep(
  "_(AUCell|ssGSEA|MSCORE|GSVA)$",
  colnames(seu_cd4_nonTreg@meta.data),
  value = TRUE
)

# Retain ANERGY, ACT and EXH prefixes
keep_prefix <- c("ANERGY", "ACT", "EXH")
pattern <- paste0("^(", paste(keep_prefix, collapse = "|"), ")")
score_cols <- score_cols_all[ grepl(pattern, score_cols_all) ]

length(score_cols)
score_cols

# Select transcription factors for correlation analysis
# Candidates include RBPJ, AHR, RARA/RARG, CREB and NR4A factors
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

tf_focus <- c("Nr4a1","Nr4a2","Nr4a3",
              "Ahr","Rara","Rarg",
              "Rbpj","Creb1")

tf_use <- intersect(tf_focus, rownames(tf_by_cell))
common_cells <- intersect(colnames(tf_by_cell),
                          rownames(seu_cd4_nonTreg@meta.data))

# State-score matrix: cells by scores
df_scores <- seu_cd4_nonTreg@meta.data[common_cells, score_cols, drop = FALSE] |>
  as.data.frame()
df_scores[] <- lapply(df_scores, as.numeric)

# TF activity matrix: cells by TFs
df_tf <- t(tf_by_cell[tf_use, common_cells, drop = FALSE]) |>
  as.data.frame()
df_tf[] <- lapply(df_tf, as.numeric)

# Compute state-score-by-TF Spearman correlations
cor_tf <- matrix(
  NA_real_,
  nrow = ncol(df_scores),
  ncol = ncol(df_tf),
  dimnames = list(colnames(df_scores), colnames(df_tf))
)
for (i in seq_len(ncol(df_scores))) {
  for (j in seq_len(ncol(df_tf))) {
    x <- df_scores[[i]]
    y <- df_tf[[j]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= 3) {
      cor_tf[i, j] <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
    }
  }
}

cor_tf

# Order rows by ANERGY, ACT and EXH

row_order <- c(
  "ANERGY_ssGSEA", "ANERGY_GSVA", "ANERGY_AUCell", "ANERGY_MSCORE",
  "ACT_ssGSEA",    "ACT_GSVA",    "ACT_AUCell",    "ACT_MSCORE",
  "EXH_ssGSEA",    "EXH_GSVA",    "EXH_AUCell",    "EXH_MSCORE"
)
row_order_use <- row_order[row_order %in% rownames(cor_tf)]
cor_tf_plot   <- cor_tf[row_order_use, , drop = FALSE]

# Remove nonfinite-only rows/columns and replace remaining nonfinite values with zero

finite_row <- rowSums(is.finite(cor_tf_plot)) > 0
finite_col <- colSums(is.finite(cor_tf_plot)) > 0
cor_tf_plot <- cor_tf_plot[finite_row, finite_col, drop = FALSE]
cor_tf_plot[!is.finite(cor_tf_plot)] <- 0

cor_tf_plot

write.csv(cor_tf_plot, "CD4_TF_vs_state_ANERGY_ACT_EXH_scores_correlation.csv", row.names = FALSE)

saveRDS(cor_tf_plot, file="cor_tf_plot.rds") # Save the processed object


# Plot the correlation heatmap
library(pheatmap)
p_tf_cor <- pheatmap::pheatmap(
  cor_tf_plot,
  cluster_rows  = FALSE,
  cluster_cols  = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col     = 45,
  scale         = "none",
  cellwidth = 18,     # Heatmap cell height
  cellheight = 18,      # Heatmap cell width
  color         = colorRampPalette(c("#4675B5","white","red"))(50),
  main          = "Correlation: TF activity vs ANERGY/ACT/EXH scores (CD4_T_nonTreg)",
  fontsize      = 8
)

p_tf_cor

ggsave(
  "CD4_TF_vs_state_scores_correlation_mouseTF.pdf",
  plot   = p_tf_cor,
  width  = 8,
  height = 5,
  dpi    = 600
)

# ==============================================================================
# FIGURE 2 | H-U | original lines 1871-1885 | historical combined PROGENy FeaturePlot
# ==============================================================================
# Visualization and correlation analysis
p_progeny_umap <- FeaturePlot(
  seu_cd4_nonTreg,
  features    = pw_cols,
  reduction   = "umap",      # Specify the UMAP reduction name
  min.cutoff  = "q05",
  max.cutoff  = "q95",
  order       = TRUE,
  ncol        = 3,
  cols = c("#FAF5C1", "#D53127"), 
  pt.size = 0.25
)
p_progeny_umap
ggsave("Fig2H_CD4_PROGENY_pathways_UMAP.pdf",
       p_progeny_umap, width = 20, height = 20, dpi = 300)

# ==============================================================================
# FIGURE 2 | V | original lines 1888-1948 | PROGENy heatmap ordered by pseudotime
# ==============================================================================
library(pheatmap)
pt_col <- "monocle3_pseudotime"   # Specify the pseudotime metadata column

# Extract pseudotime from metadata
pt_vec <- seu_cd4_nonTreg@meta.data[[pt_col]]

# Use the first column if returned as a data frame
if (is.data.frame(pt_vec)) {
  pt_vec <- pt_vec[[1]]
}

# Coerce pseudotime to numeric
pt_vec <- as.numeric(pt_vec)

# Check missing values
summary(pt_vec)

# Order cells by pseudotime
cell_order <- colnames(seu_cd4_nonTreg)[order(pt_vec)]
length(cell_order)
head(cell_order)


mat_pw <- path_by_cell[pw_use, cell_order, drop = FALSE]

# Apply gene-wise z-score scaling
mat_pw_scaled <- t(scale(t(mat_pw)))  # Row-wise z-score

# Construct cell annotations
# Use ANERGY_ssGSEA quartile groups
if ("monocle3_pseudotime" %in% colnames(seu_cd4_nonTreg@meta.data)) {
  pseudotime_vec <- seu_cd4_nonTreg$monocle3_pseudotime[cell_order]
  q <- quantile(pseudotime_vec, c(0.25, 0.75), na.rm = TRUE)
  pseudotime_group <- cut(
    pseudotime_vec,
    breaks = c(-Inf, q[1], q[2], Inf),
    labels = c("pseudotime_Eraly","pseudotime_Mid","pseudotime_Late"),
    include.lowest = TRUE
  )
  ann_col <- data.frame(Pseudotime_group = pseudotime_group)
  rownames(ann_col) <- cell_order
} else {
  ann_col <- NULL
}

p_progeny_heat <- pheatmap(
  mat_pw_scaled,
  cluster_rows   = F,
  cluster_cols   = F,
  show_rownames  = TRUE,
  show_colnames  = FALSE,
  annotation_col = ann_col,
  color          = colorRampPalette(c("#4675B5", "white", "red"))(50),
  main           = "PROGENy pathway activity across CD4 Tconv",
  fontsize_row   = 8,
  border_color   = NA
)

ggsave("Fig2H_CD4_PROGENY_pathways_heatmap.pdf",
       p_progeny_heat, width = 8, height = 4, dpi = 300)


# ==============================================================================
# FIGURE 3 | A | original lines 430-470 | ASC markers
# ==============================================================================

features_ASC_marker <- c(
  ## 1) Core ASC / stromal & adhesion markers
  "Eng",    # CD105
  "Thy1",   # CD90
  "Icam1",  # CD54
  "Alcam",  # CD166
  "Itgb1",  # CD29
  "Cd44",
  "Vcam1")  # CD106

dot_plot_asc <- DotPlot(
  seu_ASC,
  features = features_ASC_marker,
  assay = "RNA"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) +
  labs(
    title = "ASC marker"
  ) +
  scale_color_gradient2(
    low = "#4675B5", 
    mid = "#FAF5C1", 
    high = "#D53127",
    midpoint = 0,
    name = "平均表达量"
  ) +
  scale_size_continuous(
    range = c(1, 6),
    name = "表达细胞比例"
  )

print(dot_plot_asc)

ggsave("ASC_marker.dot_plot.20251209.1.pdf", plot = dot_plot_asc, width = 6, height = 4)

# ==============================================================================
# FIGURE 3 | B | original lines 472-512 | vascular markers
# ==============================================================================
features_vascular_marker <- c(  
  ## 2) Perivascular / vascular-associated markers
  "Emcn",
  "Pecam1", # CD31
  "Cdh5",   # VE-cadherin
  "Esam",
  "Rgs5",
  "Tek",    # Tie2
  "Lyve1",  # lymphatic / vascular-associated
  "Pdgfra") # ASC/stromal markers displayed after vascular markers

dot_plot_vascular <- DotPlot(
  seu_ASC,
  features = features_vascular_marker,
  assay = "RNA"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) +
  labs(
    title = "Perivascular / vascular-associated markers"
  ) +
  scale_color_gradient2(
    low = "#4675B5", 
    mid = "#FAF5C1", 
    high = "#D53127",
    midpoint = 0,
    name = "平均表达量"
  ) +
  scale_size_continuous(
    range = c(1, 6),
    name = "表达细胞比例"
  )

print(dot_plot_vascular)

ggsave("Vascular_marker.dot_plot.20251209.1.pdf", plot = dot_plot_vascular, width = 7, height = 4)

# ==============================================================================
# FIGURE 3 | C | original lines 514-553 | ECM markers
# ==============================================================================



features_ECM_marker <- c(    
  ## 3) ECM / tissue-remodeling genes
  "Fmod",
  "Prelp",
  "Dpt",
  "Mmp13")

dot_plot_ECM <- DotPlot(
  seu_ASC,
  features = features_ECM_marker,
  assay = "RNA"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) +
  labs(
    title = "ECM / tissue-remodeling genes"
  ) +
  scale_color_gradient2(
    low = "#4675B5", 
    mid = "#FAF5C1", 
    high = "#D53127",
    midpoint = 0,
    name = "平均表达量"
  ) +
  scale_size_continuous(
    range = c(1, 6),
    name = "表达细胞比例"
  )

print(dot_plot_ECM)

ggsave("ECM_marker.dot_plot.20251209.1.pdf", plot = dot_plot_ECM, width = 6, height = 4)

# ==============================================================================
# FIGURE 3 | D | original lines 555-597 | immunoregulatory markers
# ==============================================================================



features_Imm_marker <- c(     
  ## 4) Immunoregulatory ligands and ecto-enzymes
  "Il33",
  "Cd200",
  "Bmp4",
  "Jag1",
  "Entpd1", # CD39
  "Nt5e"    # CD73
)

dot_plot_Imm <- DotPlot(
  seu_ASC,
  features = features_Imm_marker,
  assay = "RNA"
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) +
  labs(
    title = "Immunoregulatory ligands and ecto-enzymes"
  ) +
  scale_color_gradient2(
    low = "#4675B5", 
    mid = "#FAF5C1", 
    high = "#D53127",
    midpoint = 0,
    name = "平均表达量"
  ) +
  scale_size_continuous(
    range = c(1, 6),
    name = "表达细胞比例"
  )

print(dot_plot_Imm)

ggsave("Imm_marker.dot_plot.20251209.1.pdf", plot = dot_plot_Imm, width = 6, height = 4)

# ==============================================================================
# FIGURE 3 | E | original lines 3463-3469 | ASC-CD4 UMAP objects
# ==============================================================================
plot.LR.umap = DimPlot(seu_LR, reduction = "umap",group.by = "active_ident",label = TRUE)
plot.LR.umap

plot.CD4.umap = DimPlot(seu_cd4_nonTreg, reduction = "umap",group.by = "active_ident",label = TRUE)

plot.CD4.umap
ggsave("CD4_UMAP.pdf", plot = plot.CD4.umap, width = 6.5141, height = 3.2754)
ggsave("CD4_ASC_UMAP.pdf", plot = plot.LR.umap, width = 10, height = 6)

# ==============================================================================
# FIGURE 3 | F | original lines 2928-2971 | global LIANA communication strength
# ==============================================================================
# Aggregate communication scores across all groups
library(dplyr)
library(tidyr)
library(pheatmap)

# Compute mean communication scores by source-target pair
comm_mat <- liana_agg %>%
  group_by(source, target) %>%
  summarise(
    mean_rank = mean(aggregate_rank, na.rm = TRUE),  # Smaller aggregate ranks indicate stronger support
    .groups   = "drop"
  ) %>%
  mutate(
    score = 1 - mean_rank  # Transform ranks so larger values indicate stronger support
  ) %>%
  select(source, target, score) %>%
  pivot_wider(names_from = target, values_from = score, values_fill = 0)

comm_mat_mat <- comm_mat %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# Plot the communication heatmap
p_comm <- pheatmap(
  comm_mat_mat,
  cluster_rows = T,
  cluster_cols = T,
  color        = colorRampPalette(c("#4675B5", "white", "red"))(50),
  main         = "Global communication strength (LIANA aggregate)",
  fontsize     = 8,
  cellwidth = 18,     # Heatmap cell height
  cellheight = 18,      # Heatmap cell width
  angle_col    = 45,
  treeheight_row=0,
  treeheight_col=0,
  border_color = "white"
)
p_comm

ggsave("Global communication strength (LIANA aggregate).pdf",
       plot =p_comm ,
       width = 10,
       height =10 ,
       dpi =300 )

# ==============================================================================
# FIGURE 3 | G-H | original lines 3494-3610 | JAG1-NOTCH1 and top LR plots using precomputed liana_agg
# ==============================================================================
asc_groups <- c("ASC_CD105pos_CD90pos",
                "ASC_CD105neg_CD90pos",
                "ASC_other")
cd4_group  <- "CD4_T_nonTreg"

# Aggregate LIANA results and add ligand/receptor columns

liana_agg <- liana_aggregate(liana_res)

# MouseConsensus output provides ligand.complex and receptor.complex
# Copy these fields to ligand and receptor for downstream operations
liana_agg <- liana_agg %>%
  mutate(
    ligand   = ligand.complex,
    receptor = receptor.complex
  )

# Select ASC-to-CD4_T_nonTreg ligand-receptor pairs

asc_to_cd4 <- liana_agg %>%
  dplyr::filter(
    source %in% asc_groups,
    target == cd4_group
  ) %>%
  mutate(
    lr_pair    = paste(ligand, receptor, sep = "_"),
    rank_score = 1 - aggregate_rank   # Transform aggregate ranks to a larger-is-stronger score
  )

# Compare predicted JAG1-NOTCH1 communication across ASC groups

jag1_notch1 <- asc_to_cd4 %>%
  dplyr::filter(ligand == "Jag1", receptor == "Notch1") %>%
  mutate(axis_name = "JAG1_NOTCH1")

# Inspect represented groups
jag1_notch1 %>% dplyr::select(source, target, rank_score)

# Specify the ASC group order on the x-axis
jag1_notch1$source <- factor(jag1_notch1$source, levels = asc_groups)

p_jag1_bar <- ggplot(jag1_notch1,
                     aes(x = source, y = rank_score)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rank_score)),
            vjust = -0.3, size = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "JAG1–NOTCH1 axis from ASC subsets to CD4_T_nonTreg",
    x     = "ASC subsets (source)",
    y     = "1 - aggregate_rank (MouseConsensus)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )

p_jag1_bar
ggsave("Fig3D_JAG1_NOTCH1_ASC_to_CD4_bar.pdf",
       p_jag1_bar, width = 4.5, height = 4, dpi = 600)

# Manuscript-related communication summary
## “JAG1–NOTCH1 as a top candidate axis emanating from ASC_CD105pos_CD90pos
##  towards CD4_T_nonTreg, with markedly higher consensus scores than the same
##  axis originating from ASC_CD105neg_CD90pos or ASC_other.”

# Plot top ligand-receptor pairs and highlight JAG1-NOTCH1

topN <- 50   # Select the top N ligand-receptor pairs per ASC group by aggregate_rank

top_axes <- asc_to_cd4 %>%
  group_by(source) %>%
  arrange(aggregate_rank, .by_group = TRUE) %>%
  slice_head(n = topN) %>%
  ungroup() %>%
  mutate(
    lr_pair    = paste(ligand, receptor, sep = "_"),
    rank_score = 1 - aggregate_rank,
    highlight  = ifelse(ligand == "Jag1" & receptor == "Notch1",
                        "JAG1–NOTCH1", "Others")
  )

top_axes$source <- factor(top_axes$source, levels = asc_groups)

p_axes_bubble <- ggplot(
  top_axes,
  aes(y     = source,
      x     = reorder(lr_pair, rank_score),
      size  = rank_score,
      color = highlight)
) +
  geom_point() +
  scale_size_continuous(range = c(1.5, 6),
                        name  = "1 - aggregate_rank") +
  scale_color_manual(values = c("JAG1–NOTCH1" = "red",
                                "Others"       = "grey40"),
                     name = "") +
  coord_flip() +
  labs(
    title = "Top ASC→CD4_T_nonTreg ligand–receptor interactions (MouseConsensus)",
    y     = "ASC subsets (source)",
    x     = "Ligand–receptor pairs to CD4_T_nonTreg"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )

p_axes_bubble

ggsave("Fig3D_TOPLR_ASC_to_CD4_bubble.2.pdf",
       p_axes_bubble, 
       width = 6, 
       height = 14, 
       dpi = 300)

# ==============================================================================
# FIGURE 3 | I-U | original lines 3424-3459 | historical combined LR FeaturePlot
# ==============================================================================

t_receptors <- c("Bmpr1a","Bmpr1b","Bmpr2","Tgfbr1","Tgfbr2","Pdcd1","Havcr2","Notch1","Notch2","Met",
                 "Ahr","Ptger2","Ptger4","Adora2a")

receptors_L <- c("Bmp4","Bmpr1a",
                 "Tgfb1","Tgfb2","Tgfb3","Tgfbr1","Tgfbr2",
                 "Jag1","Notch1",
                 "Nt5e","Entpd1","Adora2a",
                 "Rbpj")

p.l.R=FeaturePlot(object = seu_scWAT,
            features = receptors_L,
           # cols = c("#4675B5","white","red"), 
            pt.size = 0.25, 
            
            min.cutoff = "q05", 
            max.cutoff = "q95", 
            label = F)
p.l.R

p.L.R=FeaturePlot(object = seu_LR,
                  features = receptors_L,
                  # cols = c("#4675B5","white","red"), 
                  pt.size = 0.05, 
                  
                  min.cutoff = "q05", 
                  max.cutoff = "q95", 
                  label = F)
p.L.R

# Save as PDF
ggsave("CD4_ASC_LRaxes_FeaturePlot.pdf",
       plot   = p.L.R,
       width  = 18,
       height = 14,
       dpi    = 300)

# ==============================================================================
# FIGURE S1 | A | original lines 1088-1099 | GSVA heatmap
# ==============================================================================
p.CD4.gsva=pheatmap(sc_gsva,
                    cluster_rows = F, 
                    cluster_cols = T, 
                    show_rownames = T, 
                    show_colnames = F,
                    color = colorRampPalette(c("#4675B5",  "white", "red"))(50),  # Color scale
                    fontsize = 10,
                    treeheight_col = 0,
                    main = "GSVA")
p.CD4.gsva

ggsave("CD4_GSVA_heatmap.pdf", plot = p.CD4.gsva, width = 8, height = 6, dpi = 600)

# ==============================================================================
# FIGURE S1 | B | original lines 1223-1246 | AUCell heatmap
# ==============================================================================
# Plot gene-set AUC heatmap
library(pheatmap)

# Extract AUC values
auc_mat <- as.matrix(AUCell::getAUC(cells_AUC))

# Draw the heatmap
CD4_AUCell_pheatmap=pheatmap(auc_mat,
         cluster_rows = F,   # Cluster rows
         cluster_cols = TRUE,   # Cluster columns
         show_rownames = TRUE,  # Show gene-set labels
         show_colnames = FALSE, # Hide cell labels
         main = "AUCell AUC Scores for Gene Sets",
         scale = "row",         # Standardize rows
         color = colorRampPalette(c("#4675B5",  "white", "red"))(50),  # Color scale
         fontsize = 10)         # Font size


ggsave("CD4_AUCell_pheatmap.pdf", 
       plot = CD4_AUCell_pheatmap,
       width = 5, height = 6, dpi = 600)




# ==============================================================================
# FIGURE S1 | C | original lines 1272-1305 | AddModuleScore heatmap
# ==============================================================================
# Plot AddModuleScore results
library(pheatmap)

# Check that score columns exist in metadata
score_cols_MSCORE <- intersect(score_cols_MSCORE, colnames(seu_cd4_nonTreg@meta.data))

# Extract the module-by-cell score matrix
MSCORE_mat <- t(as.matrix(seu_cd4_nonTreg@meta.data[, score_cols_MSCORE, drop = FALSE]))
rownames(MSCORE_mat) <- score_cols_MSCORE          # Rows are modules
colnames(MSCORE_mat) <- colnames(seu_cd4_nonTreg)  # Columns are cell barcodes

# Draw the heatmap
CD4_MSCORE_pheatmap <- pheatmap(
  MSCORE_mat,
  cluster_rows   = FALSE,  # Preserve the selected module order
  cluster_cols   = TRUE,   # Cluster columns
  show_rownames  = TRUE,   # Show module labels
  show_colnames  = FALSE,  # Hide cell labels
  main           = "AddModuleScore (Module scores) for CD4_nonTreg",
  scale          = "row",  # Standardize each module by row
  color          = colorRampPalette(c("#4675B5", "white", "red"))(50),
  fontsize       = 10,
  treeheight_col = 0 # Set dendrogram heights with treeheight_row and treeheight_col
)



# Save as PDF
ggsave("CD4_MSCORE_pheatmap.pdf",
       plot   = CD4_MSCORE_pheatmap,
       width  =8,
       height = 6,
       dpi    = 300)


saveRDS(seu_cd4_nonTreg,"seu_cd4_nonTreg_final.rds")
saveRDS(seu_LR,"seu_LR_final.rds")
writeLines(capture.output(sessionInfo()),"sessionInfo_final.txt")
setwd(oldwd)
message("RUNNABLE MANUSCRIPT SCRIPT COMPLETE")

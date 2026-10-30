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
## ========= 1) 定义签名 =========
gs_anergy <- c("RNF128","CBLB","ITCH","EGR2","EGR3","DGKA","NR4A1","NR4A2","NR4A3","TOX","PDCD1","CTLA4","LAG3","HAVCR2")
gs_act    <- c("CD69","IL2RA","IFNG","TNF","GZMB","PRF1")
gs_exh    <- c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","BATF","CXCL13")
gs_th1    <- c("TBX21","IFNG","CXCR3")
gs_th2    <- c("GATA3","IL4","IL5","IL13")
gs_th17   <- c("RORC","IL17A","IL17F","CCR6")
gs_tfh    <- c("BCL6","CXCR5")
sig_list <- list(ANERGY=gs_anergy, ACT=gs_act, EXH=gs_exh, TH1=gs_th1, TH2=gs_th2, TH17=gs_th17, TFH=gs_tfh)

## ========= 2) 以 active_ident 作为分群列 =========
stopifnot("active_ident" %in% colnames(seu_scWAT@meta.data))
seu_scWAT$active_ident <- factor(trimws(seu_scWAT$active_ident))
Idents(seu_scWAT) <- "active_ident"

## ========= 3) 取表达矩阵（log-normalized “data” 槽）=========
DefaultAssay(seu_scWAT) <- "RNA"

# 若未 NormalizeData()，先运行：seu_scWAT <- NormalizeData(seu_scWAT)

# 根据 active_ident 提取 CD4_T_nonTreg 细胞
seu_cd4_nonTreg <- subset(seu_scWAT, subset = active_ident == "CD4_T_nonTreg")

# 获取 CD4_T_nonTreg 的表达矩阵

mat<- as.matrix(GetAssayData(seu_cd4_nonTreg, layer = "data"))  # genes x cells

###  mat <- as.matrix(GetAssayData(seu_scWAT, layer = "data"))  # genes x cells 提取总的

# 大写匹配，提升鲁棒性
rownames_upper <- toupper(rownames(mat))
sig_list_UP <- lapply(sig_list, function(v) {
  v_up <- toupper(v)
  v_up[ v_up %in% rownames_upper ]
})
gene_index <- setNames(seq_along(rownames_upper), rownames_upper)



## ========= 4) 评分 =========
# 4.1 GSVA / ssGSEA
message("[GSVA] running gsva + ssgsea ...")
sig_list_for_gsva <- lapply(sig_list_UP, function(v) rownames(mat)[gene_index[v]])

sc_gsva   <- GSVA::gsva(mat, sig_list_for_gsva, method="gsva", kcdf="Gaussian", mx.diff=TRUE, parallel.sz=10)

sc_ssgsea <- GSVA::gsva(mat, sig_list_for_gsva, method="ssgsea", kcdf="Gaussian", parallel.sz=10)

for (nm in rownames(sc_gsva)) {
  seu_cd4_nonTreg[[paste0(nm,"_GSVA")]]   <- as.numeric(sc_gsva[nm, colnames(seu_cd4_nonTreg)])
  seu_cd4_nonTreg[[paste0(nm,"_ssGSEA")]] <- as.numeric(sc_ssgsea[nm, colnames(seu_cd4_nonTreg)])
}

saveRDS(seu_cd4_nonTreg, file="seu_cd4_nonTreg.rds") ####保存分析后的数据，方便下次直接读取。

seu_cd4_nonTreg=readRDS("seu_cd4_nonTreg.rds" ) 




### ========可视化=========== ###

# ==============================================================================
# FIGURE 1 UPSTREAM | AUCell | original lines 1165-1196 | AUCell ranking and score columns; reused by Figure S1B
# ==============================================================================
# 4.2 AUCell
message("[AUCell] building rankings + AUC ...")
# 将 sig_list_UP 中的基因符号转换为首字母大写
library(tools)

sig_list_UP_capitalized <- lapply(sig_list_UP, function(v) {
  # 转换为只有首字母大写
  v_cap <- toTitleCase(tolower(v))
  return(v_cap)
})

# 查看转换后的结果
print(sig_list_UP_capitalized)

# 检查输出
# 转换后的基因符号应该只有首字母大写，例如 "Rnf128", "Cblb", "Itch", ...

rankings <- AUCell::AUCell_buildRankings(mat, nCores = max(1, parallel::detectCores()-2), plotStats = FALSE)

aucMaxRank <- ceiling(0.05 * nrow(rankings))  # top 5%

cells_AUC  <- AUCell::AUCell_calcAUC(sig_list_UP_capitalized, rankings, aucMaxRank = aucMaxRank)

auc_mat    <- as.matrix(AUCell::getAUC(cells_AUC))

for (nm in rownames(auc_mat)) {
  seu_cd4_nonTreg[[paste0(nm,"_AUCell")]] <- as.numeric(auc_mat[nm, colnames(seu_cd4_nonTreg)])
}

saveRDS(seu_cd4_nonTreg, file="seu_cd4_nonTreg.rds") ####保存分析后的数据，方便下次直接读取。

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

## ===== 0) 前提检查 =====
stopifnot(inherits(seu_cd4_nonTreg, "Seurat"))
DefaultAssay(seu_cd4_nonTreg) <- "RNA"

E <- as.matrix(GetAssayData(seu_cd4_nonTreg, layer = "data"))  # genes x cells
stopifnot(nrow(E) > 100, ncol(E) > 50)

## ===== 1) 获取/构建 DoRothEA regulons（小鼠）=====
###### 优先使用 dorothea 包内置（更稳、离线可用）#####
regulons <- NULL
if (requireNamespace("dorothea", quietly = TRUE)) {
  data("dorothea_mm", package = "dorothea")
  regulons <- dorothea::dorothea_mm %>%
    dplyr::filter(confidence %in% c("A","B","C")) %>%           # 置信度筛选
    dplyr::transmute(
      tf        = as.character(tf),
      target    = as.character(target),
      mor       = ifelse(mor %in% c(1,-1), mor, 1),             # 激活/抑制；缺失置 1
      likelihood= if ("likelihood" %in% names(.)) likelihood else 1
    )
} else {
  stop("找不到包 `dorothea`。请先安装：BiocManager::install('dorothea')")
}

# 仅保留命中的 target
regulons <- regulons %>% dplyr::filter(target %in% rownames(E))
stopifnot(nrow(regulons) > 1)


## ===== 2) 计算 TF 活性（decoupleR::run_wmean）=====
# 提醒：wmean 快、稳；viper 更鲁棒但慢（可换 run_viper）
tf_act_long <- decoupleR::run_wmean(
  mat        = E,
  network    = regulons,
  .source    = "tf",
  .target    = "target",
  .mor       = "mor",
  .likelihood= "likelihood",
  minsize    = 5         # 每个 TF 至少命中 5 个靶
)

saveRDS(tf_act_long, file="tf_act_long.rds") 

tf_act_long=readRDS("tf_act_long.rds" ) 

# 不同 decoupleR 版本列名略有差异：score/statistic；condition/sample

score_col     <- c("score","statistic")[c("score","statistic") %in% names(tf_act_long)][1]
condition_col <- c("condition","sample")[c("condition","sample") %in% names(tf_act_long)][1]
stopifnot(length(score_col)==1, length(condition_col)==1)

# 转宽：TF x Cell 矩阵
#计算 TF 活性（decoupleR::run_wmean）=====
tf_act_long <- decoupleR::run_wmean(
  mat        = E,
  network    = regulons,
  .source    = "tf",
  .target    = "target",
  .mor       = "mor",
  .likelihood= "likelihood",
  minsize    = 5
)

# 选定一个 statistic，避免 TF×cell 多行
# 建议先看一眼有哪些：
table(tf_act_long$statistic)
tf_act_long <- tf_act_long %>%
  dplyr::filter(statistic == "corr_wmean")  # 或 "norm_wmean"

# 后面和你原来的逻辑类似，只不过现在每个 TF×cell 只有 1 行了
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
  dplyr::mutate(score = as.numeric(score)) %>%  # 保险起见压成 numeric
  tidyr::pivot_wider(
    names_from  = cell,
    values_from = score
  ) %>%
  tibble::column_to_rownames("TF") %>%
  as.matrix()

# 检查行/列
dim(tf_by_cell)
head(rownames(tf_by_cell))
head(colnames(tf_by_cell))

table(tf_act_long$statistic)



saveRDS(tf_by_cell, file="tf_by_cell.rds") 

tf_by_cell=readRDS("tf_by_cell.rds" ) 



## ===== 3) PROGENy 通路活性（小鼠）=====
# progeny 返回：样本 x 通路 或 通路 x 样本（不同版本会变）
path_scores <- progeny(
  E,
  organism  = "Mouse",
  top       = 500,     # 每条通路用前 500 个 footprint 基因
  perm      = 1,       # 不做置换，加速
  z_scores  = TRUE,    # 输出 z-score
  scale     = TRUE     # 对输入按基因标准化
)

dim(path_scores)
head(rownames(path_scores))
head(colnames(path_scores))
colnames(path_scores)
saveRDS(path_scores, file="path_scores.rds") 

path_scores=readRDS("path_scores.rds" ) 
 

# 统一为 Pathway x Cell
if (nrow(path_scores) == ncol(E) && ncol(path_scores) <= 30) {
  # 这种分布不合理，通常说明方向相反；下面强转
  path_by_cell <- t(path_scores)
} else if (ncol(path_scores) == ncol(E)) {
  # 行通常是 pathway，列是 cell
  path_by_cell <- as.matrix(path_scores)
} else if (nrow(path_scores) == ncol(E)) {
  path_by_cell <- t(as.matrix(path_scores))
} else {
  stop("Progeny 返回维度与细胞数不匹配，请 str(path_scores) 检查。")
}


# 把通路写入 meta.data
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
## 想看的 PROGENy 通路（根据你前面分析用的那批）
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
## ========== 0. 加载依赖 ==========
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)  # as.cell_data_set
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})


## ========== 1. 从 seu_scWAT 构建 CD4 谱系对象 ==========

# 1.1 用 active_ident 作为分群列（和你之前保持一致）
stopifnot("active_ident" %in% colnames(seu_scWAT@meta.data))
Idents(seu_scWAT) <- "active_ident"

# 1.2 定义 CD4 谱系里要包含的簇名字
#    这里给一个模板，你可以根据自己的数据在这个向量里增减
cd4_lineage_levels <- c(
  "CD4_T_nonTreg"
)
# 只保留确实存在的那些
cd4_lineage_levels <- intersect(cd4_lineage_levels, levels(seu_scWAT$active_ident))
stopifnot("CD4_T_nonTreg" %in% cd4_lineage_levels)

# 1.3 抽取 CD4 谱系细胞
cd4_lineage_cells <- colnames(seu_scWAT)[seu_scWAT$active_ident %in% cd4_lineage_levels]
length(cd4_lineage_cells)

seu_cd4_lineage <- subset(seu_scWAT, cells = cd4_lineage_cells)
seu_cd4_lineage$active_ident <- factor(seu_cd4_lineage$active_ident)
Idents(seu_cd4_lineage) <- "active_ident"

## ========== 2. 确保有 PCA + UMAP ==========
DefaultAssay(seu_cd4_lineage) <- "RNA"

if (nrow(SeuratObject::GetAssayData(seu_cd4_lineage, assay="RNA", layer="scale.data")) == 0) {
  # 如果之前没标准化过，这里兜一下底
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

## ========== 3. Seurat → monocle3 的 cell_data_set，并学习轨迹 ==========

cds_cd4 <- as.cell_data_set(seu_cd4_lineage)

# 把 active_ident 带进去（方便之后 root / 可视化）
colData(cds_cd4)$active_ident <- seu_cd4_lineage$active_ident

# 用 Seurat 的 UMAP 作为 monocle3 的 UMAP
cds_cd4@int_colData$reducedDims$UMAP <- 
  seu_cd4_lineage@reductions$umap@cell.embeddings

# 用 monocle3 自己的聚类 & learn_graph
cds_cd4 <- cluster_cells(cds_cd4, reduction_method = "UMAP")
cds_cd4 <- learn_graph(cds_cd4)

## ========== 4. 选择轨迹起点（root）并计算伪时间 ==========

cd4_meta_cds <- as.data.frame(colData(cds_cd4))

# 4.1 首选：如果有显式的 "CD4_T_act" / "CD4_T_naive" 簇，用它们当 root
root_cells <- character(0)

if ("CD4_T_act" %in% cd4_meta_cds$active_ident) {
  root_cells <- rownames(cd4_meta_cds)[cd4_meta_cds$active_ident == "CD4_T_act"]
} else if ("CD4_T_naive" %in% cd4_meta_cds$active_ident) {
  root_cells <- rownames(cd4_meta_cds)[cd4_meta_cds$active_ident == "CD4_T_naive"]
}

# 4.2 如果没有明确 ACT/naive 簇，则退而求其次：用 ANERGY/ACT 打分自动找一批高 ACT 低 ANERGY 的 root
if (length(root_cells) == 0) {
  # 需要保证你在 seu_cd4_nonTreg 里确实有 ACT_GSVA / ANERGY_GSVA
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

# 4.3 计算伪时间
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
##### 0) 前提：看一下 meta 列
colnames(seu_scWAT@meta.data)
table(seu_scWAT$ASC_CD90_CD105)
table(seu_scWAT$active_ident)

## 1) 定义要保留的群体名称 --------------------------
# 按你现在的说法：两个 ASC 组 + 1 个 CD4_T_nonTreg
asc_groups <- c("ASC_CD105pos_CD90pos", "ASC_CD105neg_CD90pos","ASC_other")  # 如有第三个 ASC 组，直接在这里加
cd4_group  <- "CD4_T_nonTreg"

## 2) 在原始对象中找这些细胞 ------------------------
# ASC 细胞（基于 ASC_CD90_CD105 列）
asc_cells <- WhichCells(
  seu_scWAT,
  expression = ASC_CD90_CD105 %in% asc_groups
)

# CD4_T_nonTreg 细胞（基于 active_ident 列）
cd4_cells <- WhichCells(
  seu_cd4_nonTreg,
  expression = active_ident == cd4_group
)

length(asc_cells); length(cd4_cells)

## 3) 合并细胞名 & 构建子对象 ------------------------
cells_use <- unique(c(asc_cells, cd4_cells))

seu_LR <- subset(seu_scWAT, cells = cells_use)
seu_LR

## 4) 定义 LR 分组列：ASC 子群 + CD4_T_nonTreg ----------
seu_LR$LR_group <- NA_character_

# ASC 用 ASC_CD90_CD105 的标签
is_asc <- seu_LR$ASC_CD90_CD105 %in% asc_groups
seu_LR$LR_group[is_asc] <- as.character(seu_LR$ASC_CD90_CD105[is_asc])

# CD4_T_nonTreg 用 active_ident 标记
is_cd4 <- seu_LR$active_ident == cd4_group
seu_LR$LR_group[is_cd4] <- cd4_group

table(seu_LR$LR_group, useNA = "ifany")

# 去掉任何 LR_group 仍为 NA 的细胞（理论上应该没有）
seu_LR <- subset(seu_LR, subset = !is.na(LR_group))

# 设置 Idents
Idents(seu_LR) <- "LR_group"

# 确认 RNA assay 已标准化
DefaultAssay(seu_LR) <- "RNA"
# 只包含 ASC 各亚群 + CD4_T_nonTreg 的新对象 seu_LR，后面所有 LR 分析都在它上面做。

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

######调整展示顺序
# 设置自定义的行名顺序
row_order <- c("ANERGY",  "ACT","EXH","TH1", "TH2", "TH17", "TFH")

# 确保行名顺序符合期望的顺序
sc_ssgsea <- sc_ssgsea[row_order, , drop = FALSE]

# 绘制热图，行名按照设置的顺序排列
p.CD4.ssgsea <- pheatmap(
  sc_ssgsea, 
  cluster_rows = F,     # 对行进行聚类
  cluster_cols = TRUE,     # 对列进行聚类
  show_rownames = TRUE,    # 显示行名
  show_colnames = FALSE,   # 不显示列名
  color = colorRampPalette(c("#4675B5",  "white", "red"))(50),  # 色标设置
  fontsize = 10,
  treeheight_col = 0,
  main = "ssGSEA"          # 热图标题
)
# 输出热图
p.CD4.ssgsea

# 保存热图为 PDF 文件
ggsave("CD4_ssGSEA_heatmap.pdf", plot = p.CD4.ssgsea, width = 8, height = 6, dpi = 300)

# ==============================================================================
# FIGURE 1 | D | original lines 1425-1456 | ANERGY ACT EXH violin
# ==============================================================================
# 只选 ssGSEA 列
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
###计算细胞比例

# 用 ssGSEA 作为例子
ss_cols <- c("ANERGY_ssGSEA","ACT_ssGSEA","EXH_ssGSEA")
ss_cols <- intersect(ss_cols, colnames(seu_cd4_nonTreg@meta.data))

df_ss_wide <- seu_cd4_nonTreg@meta.data %>%
  dplyr::select(all_of(ss_cols)) %>%
  tibble::rownames_to_column("cell")

# 计算每个细胞的“max program”
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

# 画成条形图
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
####1. 从 cds_cd4 提取 UMAP & pseudotime，并对齐到 Seurat

library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)

## 1) 从 monocle3 取 UMAP 坐标 & 伪时间
umap_cds <- reducedDims(cds_cd4)$UMAP  # matrix: cells x 2
pt_vec   <- pseudotime(cds_cd4)       # Named numeric vector

# 确保顺序按 cell name 对齐
common_cells <- intersect(rownames(umap_cds), colnames(seu_cd4_nonTreg))
length(common_cells)   # 应该是 CD4_T_nonTreg 的细胞数

umap_use <- umap_cds[common_cells, , drop = FALSE]
pt_use   <- as.numeric(pt_vec[common_cells])

## 2) 把 pseudotime 写回 Seurat 对象
seu_cd4_nonTreg$monocle_pt <- NA_real_
seu_cd4_nonTreg$monocle_pt[common_cells] <- pt_use

## 2. 手工画 “UMAP + pseudotime” 散点图

## 3) 组装一个 data.frame 用于 ggplot
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

## 4) 伪时间着色 UMAP
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
## 3. 沿伪时间看 ANERGY / ACT / EXH 的趋势

## 5) 只保留有完整分数的细胞
df_score <- df_pt %>%
  select(pt, ANERGY, ACT, EXH) %>%
  tidyr::drop_na()

## 6) 长格式：方便用 facet 画多条曲线
df_long <- df_score %>%
  tidyr::pivot_longer(
    cols      = c(ANERGY, ACT, EXH),
    names_to  = "module",
    values_to = "score"
  )

## 7) 伪时间上的平滑曲线（GAM 或 loess）
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
# 小提琴图展示 GSVA 和 ssGSEA 分数
library(Seurat)
library(ggplot2)

score_cols <- c("ANERGY_GSVA", "ACT_GSVA", "EXH_GSVA", "TH1_GSVA", "TH2_GSVA", "TH17_GSVA", "TFH_GSVA",
                "ANERGY_ssGSEA", "ACT_ssGSEA", "EXH_ssGSEA", "TH1_ssGSEA", "TH2_ssGSEA", "TH17_ssGSEA", "TFH_ssGSEA")

score_cols_GSVA <- c("ANERGY_GSVA", "ACT_GSVA", "EXH_GSVA", "TH1_GSVA", "TH2_GSVA", "TH17_GSVA", "TFH_GSVA")

score_cols_ssGSEA <- c("ANERGY_ssGSEA", "ACT_ssGSEA", "EXH_ssGSEA", "TH1_ssGSEA", "TH2_ssGSEA", "TH17_ssGSEA", "TFH_ssGSEA")

# 计算相关性矩阵
score_data <- seu_cd4_nonTreg@meta.data[, score_cols]
corr_matrix <- cor(score_data, method = "spearman", use = "pairwise.complete.obs")

# 使用 pheatmap 绘制相关性热图
library(pheatmap)
p_corr_GSVA_ssGSEA=pheatmap(corr_matrix, 
                            color = colorRampPalette(c("#4675B5",  "white", "red"))(50),  # 色标设置
                            fontsize = 10,
                            treeheight_col = 0,
                            treeheight_row  = 0,
                            main = "Correlation of Functional Scores in CD4_T_nonTreg", 
                            
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean")

p_corr_GSVA_ssGSEA

ggsave("CD4_corr_GSVA_ssGSEA_heatmap.pdf", plot = p_corr_GSVA_ssGSEA, width = 6, height = 6, dpi = 600)

# 如果seu_cd4_nonTreg中已有UMAP降维结果，直接可视化
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
#######在 CD4_T_nonTreg UMAP 上叠加 ANERGY + 5 个 TF 活性 + NOTCH_PROGENY#############

## ========= 2) UMAP 叠加 ANERGY + 5 个 TF + NOTCH_PROGENY =========

library(Seurat)
library(ggplot2)
library(patchwork)
library(tibble)
library(dplyr)

DefaultAssay(seu_cd4_nonTreg) <- "RNA"

## 2.1 把 5 个 TF 活性写回 meta.data（如果之前没写过可以执行这段）

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


## 2.4 需要展示的 features 列表

feature_list <- c(
  # 状态分数（你也可以改/追加 GSVA/MSCORE 等）
  "ANERGY_ssGSEA", "ACT_ssGSEA", "EXH_ssGSEA",
  # 5个 TF 活性
  "Nr4a1_ACT",
  "Ahr_ACT","Rara_ACT","Rbpj_ACT","Creb1_ACT"
)

feature_list <- feature_list[feature_list %in% colnames(seu_cd4_nonTreg@meta.data)]
feature_list

## 2.5 UMAP 特征图

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
##TF和ANERGY 分数分析
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1) 从 meta.data 里挑一个你最想用来分层的 ANERGY 分数
cand_cols <- c("ANERGY_ssGSEA", "ANERGY_GSVA", "ANERGY_AUCell", "ANERGY_MSCORE")
an_col    <- cand_cols[cand_cols %in% colnames(seu_cd4_nonTreg@meta.data)][1]
stopifnot(length(an_col) == 1)

seu_cd4_nonTreg$ANERGY_score_main <- seu_cd4_nonTreg@meta.data[[an_col]]

# 2) 用四分位数把细胞分成 Low / Mid / High 三档
q <- quantile(seu_cd4_nonTreg$ANERGY_score_main, c(0.25, 0.75), na.rm = TRUE)
seu_cd4_nonTreg$ANERGY_group <- cut(
  seu_cd4_nonTreg$ANERGY_score_main,
  breaks = c(-Inf, q[1], q[2], Inf),
  labels = c("Low", "Mid", "High"),
  include.lowest = TRUE
)


# 你关注的 TF
tf_focus <- c("Rbpj", "Nr4a1", "Creb1", "Ahr", "Rara")

# 看看哪些 *_ACT 列已经存在，哪些需要从 tf_by_cell 里补
tf_act_cols <- paste0(tf_focus, "_ACT")
missing_tf  <- tf_focus[!tf_act_cols %in% colnames(seu_cd4_nonTreg@meta.data)]

if (length(missing_tf) > 0) {
  stopifnot(exists("tf_by_cell"))  # 确保你运行 DoRothEA/decoupleR 时生成了这个矩阵
  
  common_cells <- intersect(colnames(tf_by_cell), colnames(seu_cd4_nonTreg))
  
  for (tf in missing_tf) {
    if (!tf %in% rownames(tf_by_cell)) next
    
    v <- as.numeric(tf_by_cell[tf, common_cells])
    
    # 对齐到 Seurat 对象的细胞顺序
    tmp <- rep(NA_real_, ncol(seu_cd4_nonTreg))
    names(tmp) <- colnames(seu_cd4_nonTreg)
    tmp[common_cells] <- v
    
    seu_cd4_nonTreg[[paste0(tf, "_ACT")]] <- tmp
  }
}


# 真实存在的 TF 列
score_cols_tf <- paste0(tf_focus, "_ACT")
score_cols_tf <- score_cols_tf[score_cols_tf %in% colnames(seu_cd4_nonTreg@meta.data)]

# 只看 ANERGY 低 & 高两组（去掉中间档）
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

# 画图：每个 TF 一格，比较 Low vs High
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


###计算每个TFTF 在 ANERGY Low vs High 之间的 p 值
library(dplyr)

# 你的 5 个 TF 活性列
tf_features <- c("Rbpj_ACT","Nr4a1_ACT","Creb1_ACT","Ahr_ACT","Rara_ACT")

meta <- seu_cd4_nonTreg@meta.data

# 只保留有分组的细胞（Low/High）
df_tf <- meta %>%
  filter(!is.na(ANERGY_group)) %>%
  select(ANERGY_group, all_of(tf_features))

# 对每个 TF 做 Wilcoxon 秩和检验（High vs Low）
tf_pvals <- sapply(tf_features, function(tf) {
  g1 <- df_tf %>% filter(ANERGY_group == "Low")  %>% pull(tf)
  g2 <- df_tf %>% filter(ANERGY_group == "High") %>% pull(tf)
  wilcox.test(g2, g1)$p.value   # 高 vs 低
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
## ===== 5) 状态分数 vs TF 活性相关性（Spearman）=====
## ===== A) 定位状态分数列（和 PROGENy 那段保持一致）=====
# 1) 所有状态评分列（AUCell / GSVA / ssGSEA / MSCORE）
score_cols_all <- grep(
  "_(AUCell|ssGSEA|MSCORE|GSVA)$",
  colnames(seu_cd4_nonTreg@meta.data),
  value = TRUE
)

# 2) 只保留 ANERGY / ACT / EXH 这些前缀（和你现在设定一致）
keep_prefix <- c("ANERGY", "ACT", "EXH")
pattern <- paste0("^(", paste(keep_prefix, collapse = "|"), ")")
score_cols <- score_cols_all[ grepl(pattern, score_cols_all) ]

length(score_cols)
score_cols

## ===== B) 选 TF：可以先聚焦一批感兴趣的 TF =====
## 你之前关注：RBPJ/NOTCH, AHR, RARA/RARG, CREB/NR4A1-3 等
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

## 状态分数：cells × state_scores
df_scores <- seu_cd4_nonTreg@meta.data[common_cells, score_cols, drop = FALSE] |>
  as.data.frame()
df_scores[] <- lapply(df_scores, as.numeric)

## TF 活性：cells × TF
df_tf <- t(tf_by_cell[tf_use, common_cells, drop = FALSE]) |>
  as.data.frame()
df_tf[] <- lapply(df_tf, as.numeric)

## Spearman 相关：模块 × TF
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

## 4.6 按你想要的 ANERGY/ACT/EXH 顺序排行

row_order <- c(
  "ANERGY_ssGSEA", "ANERGY_GSVA", "ANERGY_AUCell", "ANERGY_MSCORE",
  "ACT_ssGSEA",    "ACT_GSVA",    "ACT_AUCell",    "ACT_MSCORE",
  "EXH_ssGSEA",    "EXH_GSVA",    "EXH_AUCell",    "EXH_MSCORE"
)
row_order_use <- row_order[row_order %in% rownames(cor_tf)]
cor_tf_plot   <- cor_tf[row_order_use, , drop = FALSE]

## 4.7 去掉全 NA 行/列，并把残余 NaN/Inf 设成 0，避免 hclust 报错

finite_row <- rowSums(is.finite(cor_tf_plot)) > 0
finite_col <- colSums(is.finite(cor_tf_plot)) > 0
cor_tf_plot <- cor_tf_plot[finite_row, finite_col, drop = FALSE]
cor_tf_plot[!is.finite(cor_tf_plot)] <- 0

cor_tf_plot

write.csv(cor_tf_plot, "CD4_TF_vs_state_ANERGY_ACT_EXH_scores_correlation.csv", row.names = FALSE)

saveRDS(cor_tf_plot, file="cor_tf_plot.rds") ####保存分析后的数据，方便下次直接读取。


## 4.8 画热图
library(pheatmap)
p_tf_cor <- pheatmap::pheatmap(
  cor_tf_plot,
  cluster_rows  = FALSE,
  cluster_cols  = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col     = 45,
  scale         = "none",
  cellwidth = 18,     # 设定cell 的大小
  cellheight = 18,      # 设定cell 的大小
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
######################可视化及相关性分析########################################
p_progeny_umap <- FeaturePlot(
  seu_cd4_nonTreg,
  features    = pw_cols,
  reduction   = "umap",      # 你 UMAP 的名字如果不是 "umap"，改这里
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
pt_col <- "monocle3_pseudotime"   # 或者你实际用的列名，比如 "pseudotime"

# 1) 从 meta.data 里取出伪时间这一列
pt_vec <- seu_cd4_nonTreg@meta.data[[pt_col]]

# 2) 如果不小心是 data.frame/tibble，就取第一列
if (is.data.frame(pt_vec)) {
  pt_vec <- pt_vec[[1]]
}

# 3) 强制转成数值型
pt_vec <- as.numeric(pt_vec)

# 4) 检查一下有没有 NA
summary(pt_vec)

# 5) 按伪时间排序细胞
cell_order <- colnames(seu_cd4_nonTreg)[order(pt_vec)]
length(cell_order)
head(cell_order)


mat_pw <- path_by_cell[pw_use, cell_order, drop = FALSE]

# 可选：z-score 按行标准化，让色阶更好看
mat_pw_scaled <- t(scale(t(mat_pw)))  # 行 z-score

# 细胞注释：比如按 ANERGY 高/低，或伪时间分组
# 这里举例按 ANERGY_ssGSEA 四分位分组
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
  "Pdgfra") # ASC/stromal marker，放在血管组末尾作为过渡

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
##### 1. 全局通讯强度矩阵（所有 group 间的整体强度）#####
library(dplyr)
library(tidyr)
library(pheatmap)

## 1.1 计算 source × target 的平均通讯强度
comm_mat <- liana_agg %>%
  group_by(source, target) %>%
  summarise(
    mean_rank = mean(aggregate_rank, na.rm = TRUE),  # 小=强
    .groups   = "drop"
  ) %>%
  mutate(
    score = 1 - mean_rank  # 转成“越大越强”，直观一点
  ) %>%
  select(source, target, score) %>%
  pivot_wider(names_from = target, values_from = score, values_fill = 0)

comm_mat_mat <- comm_mat %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

## 1.2 画 heatmap
p_comm <- pheatmap(
  comm_mat_mat,
  cluster_rows = T,
  cluster_cols = T,
  color        = colorRampPalette(c("#4675B5", "white", "red"))(50),
  main         = "Global communication strength (LIANA aggregate)",
  fontsize     = 8,
  cellwidth = 18,     # 设定cell 的大小
  cellheight = 18,      # 设定cell 的大小
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

## ===== 1. 聚合 LIANA 结果，并补 ligand / receptor 列 =====

liana_agg <- liana_aggregate(liana_res)

# MouseConsensus 聚合表里只有 ligand.complex / receptor.complex
# 为了方便后面 dplyr 操作，这里显式复制成 ligand / receptor
liana_agg <- liana_agg %>%
  mutate(
    ligand   = ligand.complex,
    receptor = receptor.complex
  )

## ===== 2. 取 ASC → CD4_T_nonTreg 的所有 LR  =====

asc_to_cd4 <- liana_agg %>%
  dplyr::filter(
    source %in% asc_groups,
    target == cd4_group
  ) %>%
  mutate(
    lr_pair    = paste(ligand, receptor, sep = "_"),
    rank_score = 1 - aggregate_rank   # 把 “越小越好” 的 rank 变成 “越大越强”
  )

## ===== 3. 抽取 JAG1–NOTCH1 轴，并画一个柱状图比较三个 ASC 群 =====

jag1_notch1 <- asc_to_cd4 %>%
  dplyr::filter(ligand == "Jag1", receptor == "Notch1") %>%
  mutate(axis_name = "JAG1_NOTCH1")

# 看一下三组有没有都出现
jag1_notch1 %>% dplyr::select(source, target, rank_score)

# 为了在 x 轴上按你习惯的顺序排 ASC 三群
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

## 这张图对应你文中的句子：
## “JAG1–NOTCH1 as a top candidate axis emanating from ASC_CD105pos_CD90pos
##  towards CD4_T_nonTreg, with markedly higher consensus scores than the same
##  axis originating from ASC_CD105neg_CD90pos or ASC_other.”

## ===== 4. 辅助图：各 ASC 群对 CD4 的 top LR + 高亮 JAG1–NOTCH1  =====

topN <- 50   # 每个 ASC 群取前 N 条 LR（按 aggregate_rank）

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

## 保存 PDF
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
                    color = colorRampPalette(c("#4675B5",  "white", "red"))(50),  # 色标设置
                    fontsize = 10,
                    treeheight_col = 0,
                    main = "GSVA")
p.CD4.gsva

ggsave("CD4_GSVA_heatmap.pdf", plot = p.CD4.gsva, width = 8, height = 6, dpi = 600)

# ==============================================================================
# FIGURE S1 | B | original lines 1223-1246 | AUCell heatmap
# ==============================================================================
# 绘制基因集 AUC 值的热图
library(pheatmap)

# 提取 AUC 值矩阵
auc_mat <- as.matrix(AUCell::getAUC(cells_AUC))

# 绘制热图
CD4_AUCell_pheatmap=pheatmap(auc_mat,
         cluster_rows = F,   # 聚类行
         cluster_cols = TRUE,   # 聚类列
         show_rownames = TRUE,  # 显示行名（基因集）
         show_colnames = FALSE, # 不显示列名（细胞）
         main = "AUCell AUC Scores for Gene Sets",
         scale = "row",         # 对行进行标准化
         color = colorRampPalette(c("#4675B5",  "white", "red"))(50),  # 色标设置
         fontsize = 10)         # 字体大小


ggsave("CD4_AUCell_pheatmap.pdf", 
       plot = CD4_AUCell_pheatmap,
       width = 5, height = 6, dpi = 600)




# ==============================================================================
# FIGURE S1 | C | original lines 1272-1305 | AddModuleScore heatmap
# ==============================================================================
# 绘制 AddModuleScore 的热图
library(pheatmap)

# 确保打分列都存在于 meta.data
score_cols_MSCORE <- intersect(score_cols_MSCORE, colnames(seu_cd4_nonTreg@meta.data))

# 提取 MSCORE 矩阵：行 = 基因集（模块），列 = 细胞
MSCORE_mat <- t(as.matrix(seu_cd4_nonTreg@meta.data[, score_cols_MSCORE, drop = FALSE]))
rownames(MSCORE_mat) <- score_cols_MSCORE          # 行名：各模块
colnames(MSCORE_mat) <- colnames(seu_cd4_nonTreg)  # 列名：细胞 barcodes

# 绘制热图
CD4_MSCORE_pheatmap <- pheatmap(
  MSCORE_mat,
  cluster_rows   = FALSE,  # 不对基因集（模块）做聚类，保持设定顺序
  cluster_cols   = TRUE,   # 对列聚类
  show_rownames  = TRUE,   # 显示行名（模块）
  show_colnames  = FALSE,  # 不显示列名（细胞）
  main           = "AddModuleScore (Module scores) for CD4_nonTreg",
  scale          = "row",  # 对每个模块按行标准化
  color          = colorRampPalette(c("#4675B5", "white", "red"))(50),
  fontsize       = 10,
  treeheight_col = 0 #聚类树高度（占的空间）：treeheight_row, treeheight_col
)



# 保存为 PDF
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

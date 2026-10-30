# Software and execution prerequisites

Install packages separately, using compatible original environments/reference versions.

| Module | Package groups |
|---|---|
| Bulk primary | tidyverse, limma, tibble, dplyr; full plots add pheatmap, msigdbr, ggplot2, ggrepel |
| Bulk supplement | limma, ggplot2, openxlsx, msigdbr, ggrepel |
| Mouse | Seurat, SeuratObject, SeuratWrappers, SingleCellExperiment, monocle3, GSVA, AUCell, decoupleR, dorothea, progeny, liana, Biobase, conflicted, tidyverse, patchwork, pheatmap, viridis; see stage checks |
| Human adipose preprocessing | Matrix, data.table, dplyr, ggplot2, patchwork; QC adds SingleCellExperiment, SummarizedExperiment, scDblFinder, tibble; integration adds Seurat, harmony, clustree, ggraph |
| Human adipose figure | Seurat, SeuratObject, Matrix, ggplot2, patchwork, dplyr, scales, digest, yaml |
| Human breast | Seurat, SeuratObject, Matrix, GSVA, ggplot2, patchwork, scales, digest, jsonlite, ggtext, png, grid |
| Proteomic replay | limma, dplyr, tibble; plots add ggplot2, pheatmap |
| Clinical | GSVA, survival, survminer, dplyr, ggpubr, ggplot2, patchwork; TCGA adds SummarizedExperiment, cowplot, pheatmap, clusterProfiler, msigdbr, limma, tidyr, enrichplot |

[Historical versions](archive/historical/environment/software_versions.md)
distinguish original records from later reproduction runs. Unknown historical
versions are not inferred from the current workstation. The retained TCGA pathway
stage requires GSVA::ssgseaParam; the older mouse workflow uses the legacy API.
Separate environments may be necessary.

Reference gene sets are inputs: supply original versions and existing package
resources/caches. Recommended scripts initiate no explicit downloads; this does
not guarantee that third-party package internals never access the network.

Local clinical preflight found missing survminer and ggpubr; neither was installed.
Windows locale/Unicode warnings were observed. Numerical table comparisons passed,
but pixel-identical plot typography is not claimed.

# Software requirements

Install the packages required by the selected module before execution.

| Module | Package groups |
|---|---|
| Bulk primary | tidyverse, limma, tibble, dplyr; full plots add pheatmap, msigdbr, ggplot2, ggrepel |
| Bulk supplement | limma, ggplot2, openxlsx, msigdbr, ggrepel |
| Baseline stromal identity (Figure 4C–D) | pheatmap, ggplot2 |
| Mouse | Seurat, SeuratObject, SeuratWrappers, SingleCellExperiment, monocle3, GSVA, AUCell, decoupleR, dorothea, progeny, liana, Biobase, conflicted, tidyverse, patchwork, pheatmap, viridis; see stage checks |
| Human adipose preprocessing | Matrix, data.table, dplyr, ggplot2, patchwork; QC adds SingleCellExperiment, SummarizedExperiment, scDblFinder, tibble; integration adds Seurat, harmony, clustree, ggraph |
| Human adipose figure | Seurat, SeuratObject, Matrix, ggplot2, patchwork, dplyr, scales, digest, yaml |
| Human breast | Seurat, SeuratObject, Matrix, GSVA, ggplot2, patchwork, scales, digest, jsonlite, ggtext, png, grid, ragg |
| Proteomic replay | limma, dplyr, tibble; plots add ggplot2, pheatmap |
| Clinical | GSVA, survival, survminer, dplyr, ggpubr, ggplot2, patchwork; TCGA adds SummarizedExperiment, cowplot, pheatmap, clusterProfiler, msigdbr, limma, tidyr, enrichplot |

## Recorded versions

| Scope | Software / reference versions |
|---|---|
| Single-cell reproduction environment | R 4.3.3; Seurat 5.4.0; SeuratObject 5.3.0; monocle3 1.4.26; GSVA 1.50.5; AUCell 1.24.0; decoupleR 2.9.7; dorothea 1.14.1; progeny 1.24.0; liana 0.1.14 |
| Original mouse communication analysis | Seurat 5.3.0 |
| Human ER-positive integration | harmony 1.2.4 |
| PRJNA932038 count generation | Cell Ranger 10.1.0; GRCh38-2024-A reference; intronic reads included |
| Bulk supplementary enrichment | msigdbr 26.1.0; MSigDB 2026.1 |

Versions not listed above were not established from the available records.
The listed versions describe module-specific environments, not a single
environment verified for all workflows.

## Compatibility notes

- The TCGA pathway stage requires the GSVA API providing `GSVA::ssgseaParam`;
  the mouse workflow uses the legacy GSVA API. Separate R package environments
  may be necessary; GSVA 1.50.5 is not a universal requirement for both workflows.
- Bulk Hallmark analyses use human Hallmark sets mapped to mouse gene symbols,
  not mouse-native Hallmark sets. Offline execution requires the corresponding
  MSigDB and ortholog-mapping resources in the `msigdbr` cache; see the
  [Bulk RNA instructions](analysis/bulk_rna/README.md).
- Baseline identity plotting reads a saved differential-expression table and
  does not require DESeq2 to fit a new model.
- Cell Ranger is required for upstream count generation, not for the R workflows
  starting from existing count matrices or annotated objects.

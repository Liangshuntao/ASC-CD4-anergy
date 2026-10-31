# Human breast adipose candidate evidence

Figure S7 is generated from the annotated PRJNA932038 Seurat object and its
manually reviewed cluster map.

```
Rscript analysis/human_adipose/plot_asc_like_candidate_evidence.R --input /path/annotated.rds --annotation-map /path/mapping.csv --out /path/new_output --preflight
Rscript analysis/human_adipose/plot_asc_like_candidate_evidence.R --input /path/annotated.rds --annotation-map /path/mapping.csv --out /path/new_output
```

The input object must contain the study annotations and embeddings used for
plotting. The cluster map is provided in metadata/reviewed_cluster_map.csv.

See [DATA_ACCESS.md](../../DATA_ACCESS.md) for input availability.
The staged preprocessing interface is described in PREPROCESSING.md, followed by
explicit import of the reviewed cluster map. Dependencies: Seurat, SeuratObject, Matrix, ggplot2, patchwork,
dplyr, scales, digest and yaml. No packages are installed by this script.
Preflight checks input paths and package dependencies.
Historical output names are retained to preserve the existing figure/data mapping.

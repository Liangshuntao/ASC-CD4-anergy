# PRJNA932038 preprocessing and reviewed annotation

The source chain for the final annotated object is the original staged workflow
(step_01/02/03 and reviewed Step04 mapping), not the separate exploratory SIMPLE
scDblFinder script. The latter remains historical only. Original config is archived;
seed, sample mapping, QC formulas, Harmony dimensions and locked resolution 0.4
are unchanged. Human review remains explicit, and no approval file is generated.

Run one stage at a time, using a new output directory each time:

```
Rscript analysis/human_adipose/run_preprocessing.R --stage read --counts-dir /path/six_cellranger_sample_folders --out /path/read_output --preflight
Rscript analysis/human_adipose/run_preprocessing.R --stage read --counts-dir /path/six_cellranger_sample_folders --out /path/read_output
Rscript analysis/human_adipose/run_preprocessing.R --stage qc --previous /path/read_output/objects/PRJNA932038_preQC.rds --reviewed-input --out /path/qc_output
Rscript analysis/human_adipose/run_preprocessing.R --stage integrate --previous /path/qc_output/objects/PRJNA932038_QC_singlets.rds --reviewed-input --out /path/integration_output
```

The input folder contains s1_SRR23526669 through s6_SRR23526664, each with a
filtered_feature_bc_matrix subfolder. The six sample/patient/tissue labels and
QC settings are in config/preprocessing_parameters.R. Cell Ranger inputs were
generated with v10.1.0, GRCh38-2024-A and intronic reads included. The R runner
does not download reads, run Cell Ranger or install any dependencies.

QC/doublet filtering is sample-wise. Integration retains the original candidate
resolution grid for audit and the previously selected 0.4 cluster assignment.
--reviewed-input acknowledges review of the supplied checkpoint; the runner does
not bypass or modify the original project's approval chain. It is a publication
replay interface, not a controller for the historical gated project.

Final cell identities are reviewed external annotation inputs, not automatically
inferred by the runner. The figure-replay script independently applies and checks
the reviewed cluster map. Use the final annotated checkpoint for exact figure
replay. Recreating manual annotation from the clustered object uses the supplied
mapping only; the annotation's biological validity is a PI decision. Original
marker/SingleR exploratory aids remain in the historical archive and are not
silently promoted to final labels.

To import that exact mapping into the clustered checkpoint:

```
Rscript analysis/human_adipose/apply_reviewed_annotation.R --input /path/integration_output/objects/PRJNA932038_global_harmony_clustered.rds --map analysis/human_adipose/metadata/reviewed_cluster_map.csv --reviewed-map --out /path/new_annotated.rds
```

The utility imports only the two reviewed label columns and requires an exact
cluster-count match; it does not alter assays, embeddings or cluster assignments.

Dependencies: read stage: Matrix, dplyr, ggplot2, patchwork, data.table;
QC: additionally SingleCellExperiment, SummarizedExperiment, scDblFinder, tibble;
integration: Seurat, harmony, clustree, ggraph and the plotting dependencies.
The full 63,450-nucleus processing was not rerun during code curation.

# PRJNA932038 preprocessing and reviewed annotation

The workflow reads Cell Ranger count matrices, performs sample-wise quality
control and doublet filtering, and integrates samples with Harmony. Final cell
identities are assigned using the supplied manually reviewed cluster-to-cell-type
mapping, rather than inferred automatically.

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

## Parameters and input checks

| Argument | Meaning |
|---|---|
| `--stage read` | Read the six Cell Ranger count-matrix folders. |
| `--stage qc` | Apply sample-wise QC and doublet filtering to the preceding checkpoint. |
| `--stage integrate` | Integrate QC-filtered cells and generate clustering outputs. |
| `--counts-dir` | Parent directory containing the six sample folders; required for `read`. |
| `--previous` | RDS checkpoint from the preceding stage; required for `qc` and `integrate`. |
| `--reviewed-input` | Confirm that the preceding checkpoint and its QC summaries have been inspected; required for `qc` and `integrate`. This flag does not perform QC review. |
| `--out` | New or empty output directory, separate from the input location. |
| `--preflight` | Check required paths and packages without reading count matrices or running the analysis. |

The configuration sets the random seed to 1234, uses dimensions 1–30 for
integration, and groups samples by `sample_id` for Harmony correction. The
selected global and stromal clustering resolutions are both 0.4. Sample labels
and QC thresholds are specified in `config/preprocessing_parameters.R`.

## Manual annotation

Inspect cluster markers before applying cell-type labels. The supplied mapping
records the manual annotations used for the study; it should only be applied to
the corresponding clustering output. Import it with:

```
Rscript analysis/human_adipose/apply_reviewed_annotation.R --input /path/integration_output/objects/PRJNA932038_global_harmony_clustered.rds --map analysis/human_adipose/metadata/reviewed_cluster_map.csv --reviewed-map --out /path/new_annotated.rds
```

For this command, `--input` is the clustered Seurat RDS, `--map` is the
cluster-to-cell-type CSV, and `--out` is a new annotated RDS file. The required
`--reviewed-map` flag confirms that the supplied manual mapping has been checked;
it does not infer or validate cell identities.

The utility imports `celltype_v2` and `major_celltype_v2` and requires an exact
cluster-count match. It does not alter assays, embeddings or cluster assignments,
and refuses to overwrite an existing output file. Use the final annotated object
and the same mapping for figure generation.

Dependencies: read stage: Matrix, dplyr, ggplot2, patchwork, data.table;
QC: additionally SingleCellExperiment, SummarizedExperiment, scDblFinder, tibble;
integration: Seurat, harmony, clustree, ggraph and the plotting dependencies.

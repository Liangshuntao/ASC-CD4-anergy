# Human breast cohort evidence workflow

This workflow replays the final GSE161529 cohort figures, patient-level PAM50
resemblance and ER-positive CD4-state associations. It is not a new cell-identity
discovery pipeline. Approved annotation tables, PAM50 centroids and the explicit
sample/patient map used in the source project are supplied under metadata/.

```
Rscript analysis/human_breast/run_breast_analysis.R --data-dir /path/author_seurat_objects --normal-object /path/SeuratObject_NormTotal_annotated.rds --formal-object /path/reviewed_ertotal_tnk.rds --out /path/new_output --preflight
Rscript analysis/human_breast/run_breast_analysis.R --data-dir /path/author_seurat_objects --normal-object /path/SeuratObject_NormTotal_annotated.rds --formal-object /path/reviewed_ertotal_tnk.rds --out /path/new_output
```

Input data-dir contains SeuratObject_NormTotal.rds, SeuratObject_ERTotal.rds,
SeuratObject_HER2.rds and SeuratObject_TNBC.rds from the original author release.
--normal-object explicitly supplies the originally used annotated normal object;
--formal-object supplies the reviewed ER-positive seven-cluster T/NK checkpoint.
These are immutable inputs. Raw author matrices and processed checkpoints are not
redistributed here; source links are in the original availability documentation.

The runner copies only small reviewed metadata/reference tables into a new output
workspace, never modifying input objects. It then runs: four cohort candidate
figures, final panel refinements, patient pseudobulk PAM50 resemblance, formal-cluster
CD4 Tconv/Treg scores/associations, the final association figure and cohort assembly.
No exploratory cluster-map search or alternative CD4 definition is substituted.

Patient/sample units, malignant epithelial labels, minimum-cell thresholds,
normalization, centroid correlations, state signatures, multiple-testing corrections
and plotting calculations remain exactly those in the selected final source.
PAM50 resemblance is exploratory, not a clinical subtype diagnosis. Cells are not
treated as biological replicates for patient association tests.

Input/output discovery and internal helper references have been made portable.
There are no automatic package installations or reference downloads. Dependencies
are checked explicitly. Large single-cell computations were not rerun solely for
code cleanup; inspect validation records for the actual extent of testing.

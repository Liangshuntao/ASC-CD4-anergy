# GSE183288 manuscript single-cell analysis code

This directory is the curated submission package for the single-cell analyses shown in manuscript Figures 1–3 and Figure S1.

## Scope and provenance

- Canonical input: the manually annotated Seurat object `seu_scWAT.rds` (58,241 cells; 13 samples).
- The input object is not redistributed because of its size. Supply its local path with `--input`.
- `scripts/01_EXACT historical extraction.R` is the read-only provenance archive extracted from the historical analysis script.
- `scripts/02_Runnable analysis.R` retains the historical calculations and panel parameters, with only documented execution-interface and dependency-order repairs.
- `scripts/03_Strict figure finalization.R` is the CodeBuddy-validated checkpoint replay used to generate the final non-LIANA panels.
- `scripts/04_LIANA Figure3F-H.R` is the CodeBuddy-validated SeuratObject 5-compatible LIANA implementation using `SingleCellExperiment` and `MouseConsensus`.

No FASTQ processing, Cell Ranger processing, de novo annotation, or replacement biological analysis is performed here.

## Recommended execution

Run from a clean R session with the package versions listed in `validation/sessionInfo*.txt`:

```powershell
Rscript --vanilla run_submission.R `
  --input D:/R/20250926/seu_scWAT.rds `
  --out output `
  --threads 4
```

For a quick input and dependency check:

```powershell
Rscript --vanilla run_submission.R --input D:/R/20250926/seu_scWAT.rds --out output --preflight
```

The runner never writes to the input directory. Existing output directories are rejected unless `--resume` is supplied.

## Important execution note

The historical core calculation is computationally expensive because it converts a large expression layer to a dense matrix and runs cell-level GSVA/ssGSEA. Preserve the generated checkpoint files when resuming. The final LIANA step uses `SingleCellExperiment`; this is an execution compatibility repair for SeuratObject 5 and does not change the specified `MouseConsensus` resource, groups, or five LIANA methods.

## Scientific limitations retained from the manuscript analysis

- The historical main analysis pools SC, TC, SH, and TH samples.
- Cell-level Wilcoxon tests and correlations are descriptive/exploratory and do not model biological samples as independent replicates.
- LIANA is pooled across cells and does not model sample-to-sample variation.
- The Monocle3 root uses the historical ACT-minus-ANERGY rule and is therefore not independent of the interpreted state scores.
- CD4 and ASC annotations are inherited manual labels; their original label-generation code was not recovered.
- The atlas UMAP is inherited, and the recovered code does not establish a fully independent secondary CD4 reclustering workflow.

See `provenance/INTERFACE_FIXES.md`, `provenance/METHOD_AUDIT.md`, and `provenance/figure_panel_traceability.csv` for the complete audit trail.


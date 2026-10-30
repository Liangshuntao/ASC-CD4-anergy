# Bulk RNA analysis

## Supported source and execution order

The original workflow audit identifies the post-line-155 limma interaction workflow
as the final model; the first exploratory t-test workflow is excluded. Source
ranges and permitted interface edits are recorded in provenance/source_blocks.csv.
All computation expressions in extracted blocks are retained.

The runner requires an explicit expression CSV and a new output directory:

```
Rscript analysis/bulk_rna/run_bulk_analysis.R --input /path/expression.csv --out /path/new_output --preflight
Rscript analysis/bulk_rna/run_bulk_analysis.R --input /path/expression.csv --out /path/new_output --stage model
Rscript analysis/bulk_rna/run_bulk_analysis.R --input /path/expression.csv --out /path/new_output --stage all
```

`model` runs only the primary preparation, limma fit and original table exports.
`all` runs stages 01 through 07 in a shared environment, retaining the sensitivity
model as a separate stage. A new output directory is required; input and output
must be separate. No package is installed automatically.

Input columns: Gene ID, Gene Symbol, and the 12 supplied Normalized read count
columns for ASC_N_1..3, ASC_P_1..3, after_ASC_N_1..3 and after_ASC_P_1..3.
Use the supplied fractional normalized counts, not substituted raw integer counts.
Statistical unit: sample column, three per subtype/timepoint. The original model
is unpaired and has no donor or batch covariate; suffixes do not establish pairing.
The primary design, contrasts, log2(x+1) transformation, median aggregation and
BH thresholds are unchanged. The existing cell-cycle score enters only the
separate historical sensitivity model. Never silently substitute it as primary.

## Dependencies and outputs

Model: tidyverse, limma, tibble, dplyr. Full plots/sensitivity additionally require
pheatmap, msigdbr, ggplot2 and ggrepel. Install compatible versions separately;
environment evidence remains in archive/historical/environment/software_versions.md.
msigdbr resources must already be available according to that package's requirements;
the runner is not permission to fetch new reference datasets.

Original output basenames are retained to keep figure/source correspondence.
Primary tables include ASC_RNA_Pre_P_vs_N_all.csv, ASC_RNA_Post_P_vs_N_all.csv,
and ASC_RNA_DID_P_specific_change_all.csv. Other stages produce the original
heatmap, highlighted volcano, CAMERA and sensitivity plot files.

## Current validation boundary

The current source/input reproduces Figure 8A variance (89.26% / 4.42%); an earlier
audit's 90.64% / 4.43% values are not current evidence. This
runner does not claim pixel-identical Figure 8 reproduction, and image text/legend edits
are not treated as computation changes. Reviewer-response scripts and remaining
modules are still being curated. For reviewer-response analyses, use the separately
curated run_bulk_supplement.R with --input, --did-table and a new --out directory.
The --did-table input is ASC_RNA_DID_P_specific_change_all.csv from the primary
model. Stages run within-subset limma contrasts, overlap/ORA summary, then
volcano/CAMERA figures. Their calculations are unchanged; runtime comparison of
the full supplement remains pending. No automatic installation was added.

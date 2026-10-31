# Bulk RNA analysis

## Model and execution order

The primary analysis uses a limma linear model with ASC subtype, coculture
timepoint and their interaction (`~ 0 + subtype * timepoint`). Contrasts compare
ASC_P with ASC_N before and after coculture and estimate the difference in
coculture-associated changes between the two subtypes (difference-in-differences).

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
BH thresholds are unchanged. Cell-cycle adjustment is performed as a separate
sensitivity analysis.

## Dependencies and outputs

Model: tidyverse, limma, tibble, dplyr. Full plots/sensitivity additionally require
pheatmap, msigdbr, ggplot2 and ggrepel. Install compatible versions separately;
see [DEPENDENCIES.md](../../DEPENDENCIES.md) for recorded software versions.

Hallmark enrichment uses the MSigDB human Hallmark collection (`H`), mapped to
mouse gene symbols through `msigdbr(species = "Mus musculus", category = "H")`;
these are not mouse-native Hallmark sets. The supplementary analysis records
`msigdbr` 26.1.0 and MSigDB 2026.1. Use the corresponding reference release when
reproducing those results, since gene-set membership can change between releases.

For offline execution, the required MSigDB and ortholog-mapping resources must
already be present in the cache used by the installed `msigdbr` version and
accessible to the R session. Installing the package alone does not ensure that
these resources are cached. Check that the Hallmark call above succeeds before
running the enrichment stages. With an empty cache, `msigdbr` may attempt to
retrieve resources online; the workflow does not bundle the reference cache.

Original output basenames are retained to keep figure/source correspondence.
Primary tables include ASC_RNA_Pre_P_vs_N_all.csv, ASC_RNA_Post_P_vs_N_all.csv,
and ASC_RNA_DID_P_specific_change_all.csv. Other stages produce the original
heatmap, highlighted volcano, CAMERA and sensitivity plot files.

## Supplementary response analysis

For the within-subset responses shown in Figure S5, use
run_bulk_supplement.R with --input, --did-table and a new
--out directory.
The --did-table input is ASC_RNA_DID_P_specific_change_all.csv from the primary
model. Stages run within-subset limma contrasts, overlap/ORA summary, then
volcano/CAMERA figures.

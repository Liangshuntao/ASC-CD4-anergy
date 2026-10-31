# Figure 4C–D source data

| File | Contents |
|---|---|
| baseline_expression.csv | Six-sample baseline stromal expression table for Figure 4C. |
| baseline_differential_results.csv | Saved differential-expression results used to select genes for Figure 4C. |
| baseline_marker_expression.csv | Marker-expression table for Figure 4D. |

From the repository root, run:

```sh
Rscript --vanilla analysis/bulk_rna/plot_baseline_stromal_identity.R --counts data/figure4_baseline/baseline_expression.csv --results data/figure4_baseline/baseline_differential_results.csv --markers data/figure4_baseline/baseline_marker_expression.csv --out /path/to/new_output
```

The original values and column headers are retained. This baseline dataset is
separate from the coculture RNA-seq dataset. See the
[analysis notes](../../analysis/bulk_rna/BASELINE_IDENTITY.md) for transformations
and gene-selection criteria.

Creator: Shuntao Liang. These data are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

If you use these data in your research, please cite the associated publication.

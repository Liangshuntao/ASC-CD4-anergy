# Baseline stromal identity: Figure 4C–D

Run `plot_baseline_stromal_identity.R --help` for explicit inputs. This is a
separate six-column baseline dataset, not the 12-sample coculture matrix used
for Figure 8. Required originals are `20250604 cd105 cd90 count.csv`, saved
`differential_expression_results CD105+ VS CD105-.csv`, and `MSC marker.csv`.
These tables are included under [data/figure4_baseline](../../data/figure4_baseline/README.md).

The historical DEG.R rounds the matrix and fits DESeq2; this plotting entry
uses its saved result table without refitting the model. The
historical heatmap selects pvalue < 0.05 and absolute log2 fold change > 1
(unadjusted P), and row-standardizes the
original expression values. The marker panel uses log10(value + 1), without
row scaling. Marker ordering is CD105/CD54/CD166/CD73/CD90/CD29/CD44/CD106.
Do not replace these choices with the later limma workflow.

Final figure assembly removes dendrogram/axis decorations and changes marker
labels to gene symbols. Such layout edits are separate from the calculation.

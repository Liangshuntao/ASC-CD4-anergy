# ASC–CD4 anergy: analysis code

Analysis and plotting code accompanying the manuscript:

**Perivascular-like CD105⁺CD90⁺ adipose stromal cells promote anergy-like CD4⁺ T-cell hyporesponsiveness**

The repository includes mouse and human single-cell analyses, bulk RNA-seq and
proteomic analyses, and breast cancer cohort analyses. The scripts support the
computational analyses in Figures 1–4 and 8–14 and Supplementary Figures S1–S3
and S5–S9, as detailed below. Coverage is panel-specific.

| Workflow | Entry | Figures |
|---|---|---|
| Mouse single-cell | analysis/mouse_scrna/run_mouse_analysis.R | 1–3, S1–S3 |
| Baseline stromal identity | analysis/bulk_rna/plot_baseline_stromal_identity.R | 4C–D |
| Bulk RNA-seq | analysis/bulk_rna/run_bulk_analysis.R | 8 |
| Bulk response context | analysis/bulk_rna/run_bulk_supplement.R | S5 |
| Proteomics | analysis/proteomics/run_proteomics_replay.R | 9, S6 |
| Human breast | analysis/human_breast/run_breast_analysis.R | 10, 13, S8, S9 |
| Human adipose | analysis/human_adipose/plot_asc_like_candidate_evidence.R | S7 |
| Clinical cohorts | analysis/clinical_cohorts/run_clinical_analysis.R | 11, 12, 14 |

Read each module README for execution order. Entries support --help. Example:

```sh
Rscript --vanilla analysis/bulk_rna/run_bulk_analysis.R --input /path/to/expression.csv --out /path/to/new_output --stage all
```

Inputs and outputs must be explicit. Packages are checked, not automatically
installed. See [DEPENDENCIES.md](DEPENDENCIES.md) and
[DATA_ACCESS.md](DATA_ACCESS.md). Human adipose preprocessing is described in
analysis/human_adipose/PREPROCESSING.md.

## Figure correspondence

[Figure_script_mapping.csv](Figure_script_mapping.csv) records checked panel
relationships and remaining verification limits. Historical output basenames
are not authoritative final panel IDs.
Composite arrangement, selected subplots and typography may be manual.
Figure 4C–D uses a separate baseline dataset and the entry
`analysis/bulk_rna/plot_baseline_stromal_identity.R`; see its
[input and historical-method notes](analysis/bulk_rna/BASELINE_IDENTITY.md).
Other Figure 4 panels, Figures 5–7, S4, S10, the graphical abstract and Figure 9D
are not claimed as R-generated experimental or schematic outputs.

## Testing scope

Scripts were checked for R syntax and entry-point help. Selected bulk RNA,
proteomic and figure-source tables were compared with saved results. The full
single-cell and clinical workflows were not rerun end to end. Required inputs
and their availability are listed in [DATA_ACCESS.md](DATA_ACCESS.md).

## Citation

If you use this code in your research, please cite the associated publication.

## Contact

Shuntao Liang — [lstswjs@163.com](mailto:lstswjs@163.com).

## Licensing

The study-specific code is released under the [MIT License](LICENSE).
The [Figure 4C–D source data](data/figure4_baseline/README.md) are licensed under
CC BY 4.0.
Third-party tools and reference resources retain their own licensing terms.

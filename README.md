# ASC CD4 T cell analysis code

Curated analysis interfaces for the ASC–CD4 T-cell study. Start under `analysis/`;
do not execute the historical archive as a single pipeline.

| Analysis | Recommended entry | Instructions |
|---|---|---|
| Mouse single-cell | `analysis/mouse_scrna/run_mouse_analysis.R` | [Mouse workflow](analysis/mouse_scrna/README.md) |
| Bulk primary and sensitivity | `analysis/bulk_rna/run_bulk_analysis.R` | [Bulk workflow](analysis/bulk_rna/README.md) |
| Bulk within-subset responses | `analysis/bulk_rna/run_bulk_supplement.R` | [Bulk workflow](analysis/bulk_rna/README.md) |
| Human adipose preprocessing | `analysis/human_adipose/run_preprocessing.R` | [Preprocessing](analysis/human_adipose/PREPROCESSING.md) |
| Human adipose evidence figure | `analysis/human_adipose/plot_asc_like_candidate_evidence.R` | [Figure replay](analysis/human_adipose/README.md) |
| Human breast cohorts | `analysis/human_breast/run_breast_analysis.R` | [Breast workflow](analysis/human_breast/README.md) |
| Proteomic checkpoint replay | `analysis/proteomics/run_proteomics_replay.R` | [Proteomic workflow](analysis/proteomics/README.md) |
| TCGA and METABRIC | `analysis/clinical_cohorts/run_clinical_analysis.R` | [Clinical workflow](analysis/clinical_cohorts/README.md) |

Each entry supports `--help`. Supply input/output paths explicitly and use a new
output directory. Packages are checked, not installed. Manual annotations and
checkpoint review are explicit inputs, not inferred anew. Read
[DEPENDENCIES.md](DEPENDENCIES.md) and [DATA_ACCESS.md](DATA_ACCESS.md).

## Source and figures

[Figure_script_mapping.csv](Figure_script_mapping.csv) maps all 14 main and 10
supplementary figures to workflows or identifies experimental/manual panels.
Manuscript panels may be selected and assembled from combined original plots.

[archive/historical](archive/historical) preserves original source bytes,
including exploratory and superseded code. Those files are NOT recommended
execution entries. The initial bulk t-test is excluded from the curated primary
limma workflow; the final METABRIC DSS branch replaces the earlier OS branch as
the recommended source. `provenance/` records source blocks and interface edits.

## Validation and limits

Seven bulk/proteomic result tables matched stored historical numeric results
exactly. Bulk PCA variance also matches the current manuscript. Source-transform,
syntax, entry-point, path and annotation-import checks are in `validation/`.

These checks do not establish full raw-data-to-paper reproducibility. Large
single-cell and clinical analyses were not completely rerun. Reviewed annotations,
processed checkpoints and original reference versions remain required. Proteomic
replay avoids a new random MinProb imputation. See
[CURATION_STATUS.md](CURATION_STATUS.md) for exact validation scope.

The score historically called TIDE is a custom expression-derived composite,
not official TIDE. Its median labels are not validated treatment-response
predictions. Original variable/output names remain for traceability; manuscript
terminology must be reconciled separately.

## Licensing

No open-source licence is assigned pending author/institution confirmation.
Third-party tools/reference resources retain their original licences/attribution.
Raw sequencing data and large expression/checkpoint objects are not redistributed.

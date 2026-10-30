# ASC CD4 T cell anergy analysis code

Code curation in progress, 30 October 2026. This is not the final submission release.

## Start here

Recommended entry points are under `analysis/`, with explicit input/output arguments.
Historical scripts are under `archive/historical/` and are NOT recommended entry points.
Do not run the historical bulk script from beginning to end: its first t-test
section is exploratory and is not the final limma/DID workflow.

- [Bulk RNA instructions](analysis/bulk_rna/README.md): extracted primary model,
  plots and sensitivity stages. Three primary result tables were regenerated and
  matched stored historical numerical results exactly. Current Figure 8A variance
  was also reproduced (89.26% / 4.42%). Full visual equivalence is not claimed.
- [Mouse single-cell instructions](analysis/mouse_scrna/README.md): functionally
  named stages, explicit input/output arguments and updated runner references.
  Interface transformations and syntax were checked; computations were not rerun.
- [Proteomic replay](analysis/proteomics/README.md): final protein/pathway models
  replay from the historical processed checkpoint; four result tables reproduce
  exactly. Original stochastic preprocessing is not silently rerun.
- [Clinical cohorts](analysis/clinical_cohorts/README.md): final TCGA blocks and
  the METABRIC DSS branch, with explicit input files and no automatic download.
  Full numerical clinical replay has not been validated in this environment.
- Human breast-cohort replay and earlier human-adipose preprocessing remain under
  active source reconciliation; those historical files are not recommended entries.
- Bulk reviewer-supplement scripts now have an explicit three-stage runner;
  computation bodies are unchanged, with runtime-result comparison still pending.
- [Human adipose figure replay](analysis/human_adipose/README.md): the final
  annotated-object figure script now accepts explicit object, mapping and output
  paths. This does not replace the still-unfinished preprocessing curation.

See [CURATION_STATUS.md](CURATION_STATUS.md) for outstanding work and validation scope.

## Status and scope

This package contains recovered custom R scripts for mouse single-cell analyses,
human public-cohort validation, stromal bulk RNA-seq, proteomics and clinical-cohort
analyses. It is not a new software package or a verified end-to-end reanalysis.
See Figure_script_mapping.csv for coverage and unresolved attribution.
See archive/historical/ARCHIVE_STATUS.md for checks on the initial archive and
provenance/ for curation-specific source mappings.

The original scripts are preserved without changing analysis logic. Some retain
historical workstation paths, interactive operations, automatic installation or
download commands, intermediate-object dependencies and approval checks. Inspect
them and configure inputs before use; do not source all files sequentially.
No analyses, installations or downloads were executed during packaging.

## Modules and execution

- analysis/mouse_scrna: curated runner; starts from a manually
  annotated Seurat object rather than raw reads. Original annotation-generation
  code was not recovered. Historical panel mapping is provided in provenance/.
- archive/historical/scripts/human_PRJNA932038: preprocessing and revised reviewer-figure code.
  Intermediate annotation object and original project configuration are required;
  this is not yet a standalone pipeline. No controlled-access data are included.
- archive/historical/scripts/human_GSE161529: annotation, ASC_P-like candidate plots, patient-level
  PAM50 resemblance and ER-positive CD4-state analyses. Uses the authors' public
  Seurat resources and separately generated intermediate annotations.
- analysis/bulk_rna: extracted final primary limma and sensitivity code. Supply the
  original normalized-expression table and required differential-analysis tables.
- analysis/proteomics: provenance-linked models and plots replayed from processed
  historical state, rather than a new random imputation.
- analysis/clinical_cohorts: provenance-linked TCGA analyses and final METABRIC
  DSS branch. Original median-OS/exploratory branches remain historical only.

## Data and environment

Mouse source accession: GSE183288. Human validation: PRJNA932038 and GSE161529.
GSE161529 author objects: https://figshare.com/articles/dataset/17058077
TCGA input: TCGA-BRCA through the GDC; retain the data-access conditions.
Raw data, expression matrices and large intermediate objects are not redistributed.
New bulk/proteomics data accessions and remaining processed-input access routes
must be reconciled with the manuscript before claiming complete availability.
Software evidence and workflow boundaries are in archive/historical/environment/software_versions.md.

## Custom score caveat

The score labelled TIDE in the historical TCGA script is a custom combination of
expression-derived components, not the official TIDE algorithm. Median-based
labels are not validated treatment-response predictions. Manuscript terminology
and the repository description must be reconciled before final release.

## Licensing and reuse

No open-source licence has been assigned pending author/institution confirmation.
Third-party tools and reference resources retain their original licences and
should be obtained from their original sources. Third-party author scripts are
not redistributed in this package.

# ASC CD4 T cell anergy analysis code

Initial code archive assembled on 30 October 2026. This is not the final submission release.

## Status and scope

This package contains recovered custom R scripts for mouse single-cell analyses,
human public-cohort validation, stromal bulk RNA-seq, proteomics and clinical-cohort
analyses. It is not a new software package or a verified end-to-end reanalysis.
See Figure_script_mapping.csv for coverage and unresolved attribution.
See ARCHIVE_STATUS.md for the checks performed and the limits of this initial archive.

The original scripts are preserved without changing analysis logic. Some retain
historical workstation paths, interactive operations, automatic installation or
download commands, intermediate-object dependencies and approval checks. Inspect
them and configure inputs before use; do not source all files sequentially.
No analyses, installations or downloads were executed during packaging.

## Modules and execution

- scripts/mouse_scRNA: existing curated runner and README; starts from a manually
  annotated Seurat object rather than raw reads. Original annotation-generation
  code was not recovered. Historical panel mapping is provided in provenance/.
- scripts/human_PRJNA932038: preprocessing and revised reviewer-figure code.
  Intermediate annotation object and original project configuration are required;
  this is not yet a standalone pipeline. No controlled-access data are included.
- scripts/human_GSE161529: annotation, ASC_P-like candidate plots, patient-level
  PAM50 resemblance and ER-positive CD4-state analyses. Uses the authors' public
  Seurat resources and separately generated intermediate annotations.
- scripts/bulk_RNA: historical analysis and reviewer-response scripts. Supply the
  original normalized-expression table and required differential-analysis tables.
- scripts/proteomics: historical interactive source only; not validated for batch
  execution. Includes installation calls and sequential interactive corrections.
- scripts/clinical_cohorts: historical TCGA scripts. These may download public
  data. METABRIC-specific source attribution remains unresolved.

## Data and environment

Mouse source accession: GSE183288. Human validation: PRJNA932038 and GSE161529.
GSE161529 author objects: https://figshare.com/articles/dataset/17058077
TCGA input: TCGA-BRCA through the GDC; retain the data-access conditions.
Raw data, expression matrices and large intermediate objects are not redistributed.
New bulk/proteomics data accessions and remaining processed-input access routes
must be reconciled with the manuscript before claiming complete availability.
Software evidence and workflow boundaries are in environment/software_versions.md.

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

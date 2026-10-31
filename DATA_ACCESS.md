# Data inputs and availability

This repository contains analysis code, annotation/reference tables and the
[Figure 4C–D source data](data/figure4_baseline/README.md). Bulk RNA-seq and
proteomic data are supplied separately in the manuscript supplementary-data
package and are not duplicated here. Large annotated single-cell objects and
prepared clinical objects are not included; the required inputs are listed below.

## Input inventory

| Workflow | Source | Required inputs |
|---|---|---|
| Mouse | GEO GSE183288 | Original annotated seu_scWAT.rds, inherited labels/embeddings |
| Human adipose | PRJNA932038 | Six Cell Ranger folders; reviewed annotated object/map for figure replay |
| Human breast | GSE161529; [author resource](https://figshare.com/articles/dataset/17058077) | Four author Seurat objects, annotated normal object and reviewed formal ER-positive T/NK object; seven small mapping/reference tables are included |
| Bulk RNA | Newly generated study data | Original fractional normalized-count CSV; primary DID table for supplementary context |
| Baseline stromal identity (4C–D) | Separate early six-sample study dataset | Three tables included under [data/figure4_baseline](data/figure4_baseline/README.md); see [baseline analysis instructions](analysis/bulk_rna/BASELINE_IDENTITY.md) |
| Proteomics | Newly generated study data | Saved processed checkpoint exported using export_proteomics_checkpoint.R, including the imputed expression and pathway scores |
| TCGA | TCGA-BRCA/GDC | Original GDC-prepared SummarizedExperiment (tpm_unstrand) and exact subtype table |
| METABRIC | Original cBioPortal brca_metabric export | Clinical patient/sample files and Illumina microarray expression |

HRA004141 is not used in these workflows.

For TCGA, the required exports are the original `brca_data` object as RDS and
`brca_subtype` as TSV. A newly downloaded cohort may differ in sample membership
or processing and is not interchangeable with this prepared input.

The supplied human-breast centroid table is the genefu 2.34.0 pam50.scale reference;
source attribution remains in the PAM50 stage. Reviewed annotations are fixed
inputs, not newly inferred clinical identities. Third-party resources retain their licences.

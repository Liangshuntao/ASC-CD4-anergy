# Data inputs and access boundaries

This repository provides code. An accession alone does not replace the reviewed
annotations and processed checkpoints required for exact replay.

| Workflow | Source | Required inputs |
|---|---|---|
| Mouse | GEO GSE183288 | Original annotated seu_scWAT.rds, inherited labels/embeddings |
| Human adipose | PRJNA932038 | Six Cell Ranger folders; reviewed annotated object/map for figure replay |
| Human breast | GSE161529; [author resource](https://figshare.com/articles/dataset/17058077) | Four author Seurat objects, annotated normal object and reviewed formal ER-positive T/NK object; seven small mapping/reference tables are included |
| Bulk RNA | Newly generated study data | Original fractional normalized-count CSV; primary DID table for supplementary context |
| Proteomics | Newly generated study data | Historical processed checkpoint exported by the whitelisted utility, not a fresh random imputation |
| TCGA | TCGA-BRCA/GDC | Original GDC-prepared SummarizedExperiment (tpm_unstrand) and exact subtype table |
| METABRIC | Original cBioPortal brca_metabric export | Clinical patient/sample files and Illumina microarray expression |

New-experiment deposition identifiers and processed-object access routes must be
reported separately in the manuscript. No accession or availability of an unuploaded
object is invented here. HRA004141 was requested but is not an input to these workflows.

For TCGA, export the unchanged brca_data object as RDS and brca_subtype as TSV from
the original environment. Do not silently substitute a newly downloaded cohort.
The large historical TCGA workspace was not loaded during curation.

The supplied human-breast centroid table is the genefu 2.34.0 pam50.scale reference;
source attribution remains in the PAM50 stage. Reviewed annotations are fixed
inputs, not newly inferred clinical identities. Third-party resources retain their licences.

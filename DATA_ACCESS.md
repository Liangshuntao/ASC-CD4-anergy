# Input access and reconstruction boundaries

Public source accessions do not replace the annotation and intermediate files
needed by the recovered scripts. An accession alone does not establish exact
reproducibility.

| Module | Source | Additional requirement |
|---|---|---|
| Mouse single-cell | GEO GSE183288 | Canonical annotated seu_scWAT.rds and annotation reconstruction/access route |
| Human adipose | PRJNA932038 | Six sample count matrices, project config and step04 annotation object |
| Human breast cohorts | GSE161529; Figshare 17058077 | Author Seurat objects, reviewed mappings, PAM50 reference and generated intermediates |
| Bulk RNA-seq | Newly generated study data | Original expression table, metadata and model-result inputs; repository accession to confirm |
| Proteomics | Newly generated study data | ASC CD4 protein.csv and final preprocessing/model provenance; accession to confirm |
| Clinical cohorts | TCGA-BRCA; METABRIC | Exact clinical/subtype data, cohort filters, METABRIC script/input provenance |

HRA004141 is not packaged and must not be described as used or available through
this repository merely because access was requested.

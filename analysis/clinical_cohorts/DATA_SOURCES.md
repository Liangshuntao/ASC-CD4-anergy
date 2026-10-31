# Clinical cohort data sources

## METABRIC

Study ID: `brca_metabric` (Breast Cancer, METABRIC, Nature 2012 and Nature
Communications 2016).

Obtain the clinical patient, clinical sample and microarray expression files from
the [official study](https://www.cbioportal.org/study/summary?id=brca_metabric)
or its [DataHub directory](https://github.com/cBioPortal/datahub/tree/master/public/brca_metabric).
Retain the filenames shown below and pass their parent directory using
`--metabric-dir`.

| File | SHA-256 of the study input |
|---|---|
| data_clinical_patient.txt | 2a6bbadc6e4b2b07c0505559d438da959fc5c9aba55c7194abc394e2d6454995 |
| data_clinical_sample.txt | 9a2291e065e2a0a5337b4eaeb09836e70b1bc3e1ed17714a013b53515697e6dc |
| data_mrna_illumina_microarray.txt | b8aa84ea4d205483e4b2c2794a5eac6e73d4644e0f13e74c940f8e3d81d93dd6 |

The original download revision was not recorded. These hashes identify the exact
input files; a newer download may differ.

The study-specific [licence](https://github.com/cBioPortal/datahub/blob/master/public/brca_metabric/LICENSE)
is ODbL 1.0. Attribute cBioPortal and the original studies: Curtis et al., Nature
2012; Pereira et al., Nature Communications 2016; Rueda et al., Nature 2019.
The code's MIT licence does not replace these data terms.

## TCGA-BRCA

The source is the [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/),
project TCGA-BRCA. This workflow requires the prepared expression object and
PAM50 table used for the study. Re-downloading current GDC data may not recreate
the historical sample set. See [DATA_ACCESS.md](../../DATA_ACCESS.md) for the
processed-input availability status.

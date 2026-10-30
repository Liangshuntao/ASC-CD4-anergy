# Software versions for manuscript preparation

Prepared: 30 October 2026. Versions below are taken from archived workflow records, not from a newly queried installation. This is a version-reporting resource; it does not establish that every manuscript figure was regenerated in the same environment.

## Versions supported by archived records

| Software | Version | Applicable workflow and scope |
|---|---|---|
| R | 4.3.3 | Archived mouse communication, single-cell reproduction, human single-cell/nucleus and bulk-supplement records; not yet a verified global version for all historical TCGA/proteomic computations |
| Seurat | 5.3.0 | Historical mouse ASC–CD4 communication analysis |
| Seurat | 5.4.0 | Later mouse single-cell reproduction and human single-cell/nucleus workflows |
| SeuratObject | 5.3.0 | Later reproduction and human workflows |
| monocle3 | 1.4.26 | Mouse single-cell reproducibility run; earliest historical trajectory-generation version remains unconfirmed |
| GSVA | 1.50.5 | Mouse reproducibility and ER-positive CD4 functional-state scoring records |
| AUCell | 1.24.0 | Mouse single-cell reproducibility record |
| decoupleR | 2.9.7 | Mouse single-cell reproducibility record |
| dorothea | 1.14.1 | Mouse single-cell reproducibility record |
| progeny | 1.24.0 | Mouse single-cell reproducibility record |
| liana | 0.1.14 | Historical mouse communication and later reproducibility records |
| harmony | 1.2.4 | GSE161529 ER-positive T/NK integration records |
| scDblFinder | 1.16.0 | Human-validation workflow session record at ASC_P_human_validation_simple_scDblFinder/logs/sessionInfo.txt |
| SingleR | 2.4.1 | Human annotation package-version record at step04/attempt_008/tables/step04_package_versions.csv; checked alongside locked configuration, not inferred from the current workstation |
| celldex | 1.12.0 | Same human annotation package-version record |
| limma | 3.58.1 | Available in archived mouse and human workflow environments; exact version for historical bulk/proteomics models still requires workflow-specific evidence |
| Cell Ranger | 10.1.0 | All six PRJNA932038 sample-level web_summary.html reports, with count commands |

PRJNA932038 transcriptome reference: **GRCh38-2024-A**, with **intronic reads included**. A reference resource is not a software package version.

The R package `cellranger` 1.1.0 appearing in an R sessionInfo is unrelated to the 10x Genomics Cell Ranger 10.1.0 pipeline; do not confuse the two.

## English text supported by the records

The following can be used in a reproducibility/software subsection, retaining the distinction between historical and later runs:

> The archived single-cell reproducibility workflow was run in R v4.3.3 using Seurat v5.4.0, SeuratObject v5.3.0, monocle3 v1.4.26, GSVA v1.50.5, AUCell v1.24.0, decoupleR v2.9.7, dorothea v1.14.1, progeny v1.24.0, and liana v0.1.14. Earlier mouse cell–cell communication analyses used Seurat v5.3.0 and liana v0.1.14. Human ER-positive T-cell integration and state-scoring workflows used Seurat v5.4.0, harmony v1.2.4, and GSVA v1.50.5.

For the human single-nucleus preprocessing subsection:

> Publicly available sequencing reads from PRJNA932038 were reprocessed using Cell Ranger v10.1.0 and the GRCh38-2024-A reference, with intronic reads included. Downstream human annotation records documented scDblFinder v1.16.0, SingleR v2.4.1, and celldex v1.12.0 in their respective workflows.

Do not replace “reproducibility workflow” with “all analyses” without confirming which outputs are used in the final submission.

## Items not yet suitable for an unqualified version statement

| Item | Current status |
|---|---|
| R 4.3.2 | No supporting execution record found in the inspected workflows |
| genefu 2.34.0 | Listed in the manuscript/QA; exact PAM50 execution-version evidence not located |
| imputeLCMD | Exact original proteomic run version not recovered |
| clusterProfiler / org.Mm.eg.db | Locked configuration lists clusterProfiler 4.10.1, but corresponding proteomic execution versions are not closed |
| TCGAbiolinks / survival / survminer | Scripts located; original TCGA/METABRIC run versions not established. A package loaded by a different workflow is not evidence for these models |
| FlowJo | Requires the actual cytometry project/application metadata |
| GSE161529 author reference v3.0.0 | Reference-resource claim, not a Cell Ranger pipeline version; not newly reverified here |

## Custom R score previously labelled TIDE

The PI has confirmed that the final manuscript uses the custom score implemented in `D:/R/20260319 TCGA/TCGA_BRCA_ASC.R`. This is **custom R code, without a verified author-assigned script version**; it must not be assigned an official TIDE software version. The source file and its archived date/hash can identify this code. The interpretation and naming issue is separate from this software inventory and remains unresolved by a version table.

Source evidence index: `version_evidence.json` and `E1_software_version_findings.md` in this directory, plus the workflow-specific records cited above. No original manuscript or analysis script was modified.

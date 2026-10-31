# Clinical cohort analyses

Use run_clinical_analysis.R with --stage tcga or --stage metabric.
The METABRIC analysis uses disease-specific survival (DSS) and the joint
top/bottom quartile subset of ASC_P and CD4 T-cell scores.

```
Rscript analysis/clinical_cohorts/run_clinical_analysis.R --stage tcga --tcga-object /path/prepared_brca.rds --subtypes /path/brca_subtypes.tsv --out /path/new_output --preflight
Rscript analysis/clinical_cohorts/run_clinical_analysis.R --stage metabric --metabric-dir /path/brca_metabric --out /path/new_output --preflight
```

Remove --preflight to execute. TCGA requires the original GDC-prepared
SummarizedExperiment (tpm_unstrand assay) and the exact subtype table used at the
time of analysis. Do not substitute a newly downloaded cohort and assume identical
membership. The runner performs no GDC download or subtype API request.
METABRIC requires data_clinical_patient.txt, data_clinical_sample.txt and
data_mrna_illumina_microarray.txt from the original cBioPortal study export.
See [DATA_SOURCES.md](DATA_SOURCES.md) for official access routes, input hashes
and data licensing.

Original signatures, thresholds, cohort filters, time/event definitions, score
normalization and survival formulas are retained. The METABRIC branch preserves
the historical OS_MONTHS/VITAL_STATUS definition used for DSS. No redefinition
of endpoint, covariate adjustment or cutoff optimization is introduced.

The score historically labelled TIDE is a custom expression composite, not the
official TIDE implementation. The median-based True/False labels are retained
for source traceability, not validated clinical response predictions. Plot titles
use neutral terms.

Dependencies are checked before execution. GSVA API selection checks for
ssgseaParam availability. See [DEPENDENCIES.md](../../DEPENDENCIES.md) for
software versions and compatibility requirements.

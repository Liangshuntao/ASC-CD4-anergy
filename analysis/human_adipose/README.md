# Human breast adipose candidate evidence

Recommended figure replay begins with the audited annotated PRJNA932038 Seurat
object and its reviewed cluster map. It is not a FASTQ-to-annotation pipeline.

```
Rscript analysis/human_adipose/plot_asc_like_candidate_evidence.R --input /path/annotated.rds --annotation-map /path/mapping.csv --out /path/new_output --preflight
Rscript analysis/human_adipose/plot_asc_like_candidate_evidence.R --input /path/annotated.rds --annotation-map /path/mapping.csv --out /path/new_output
```

All final candidate, denominator, marker-count and nested ENG/THY1 definitions,
metadata assertions, embeddings and plot calculations are unchanged. Changes
are limited to explicit file paths, dependency checks and input manifest paths.
The object must contain the original audited fields; the script does not invent
annotations or silently regenerate a different embedding. Original source and
preprocessing are kept under archive/historical/scripts/human_PRJNA932038.
The reviewed mapping CSV is available there under steps/step04/attempt_011/tables/.

This covers the previously generated reviewer evidence used for Figure S7; source
access/reconstruction for the annotated object and preprocessing curation still
need closure. Dependencies: Seurat, SeuratObject, Matrix, ggplot2, patchwork,
dplyr, scales, digest and yaml. No packages are installed by this script.
Preflight verifies paths/dependencies, not scientific validity or all object metadata.
Historical output names are retained to preserve the existing figure/data mapping.

# Submission scope

The submission package contains only manuscript-relevant single-cell code and its audit trail. Exploratory scripts, failed retry scripts, batch wrappers, large RDS checkpoints, generated PDF figures, and `Rplots.pdf` are intentionally excluded.

The exact historical extraction is retained separately from runnable code so that interface repairs cannot be mistaken for original analysis statements. All accepted repairs are limited to input/output handling, dependency order, supported object accessors, namespace conflicts, and the Seurat-to-SingleCellExperiment compatibility path required by LIANA under SeuratObject 5.


# Mouse single cell workflow

Recommended entry point: run_mouse_analysis.R. Historical source remains in
archive/historical/scripts/mouse_scRNA and is not the recommended entry point.

```
Rscript analysis/mouse_scrna/run_mouse_analysis.R --input /path/seu_scWAT.rds --out /path/new_output --threads 4 --preflight
Rscript analysis/mouse_scrna/run_mouse_analysis.R --input /path/seu_scWAT.rds --out /path/new_output --threads 4
```

Stages run in this order: compute signatures/activities/trajectory, finalize figure
panels from checkpoints, then LIANA ASC-to-CD4 communication. --resume retains the
original checkpoint logic. The runner rejects nonempty output unless --resume is
specified. No package installation is performed. Dependencies and original version
records remain in archive/historical/scripts/mouse_scRNA/README.md and validation/.

The required input is the historically annotated Seurat object, not a raw matrix.
Its annotation metadata, sample labels and inherited embeddings are required.
The stage scripts preserve historical state gene sets, algorithm parameters,
Monocle3 root rule, cell counts, tests and pooling. Some analyses are descriptive
cell-level comparisons, not independent biological-sample inference. Do not treat
the interface cleanup as a new scientific validation or independent reanalysis.

The current curation changes explicit input requirements, names and internal
stage references only. See provenance/mouse_interface_changes.json. No expensive
single-cell computations were rerun during curation. The current figure mapping
records the selected LR/BMP/TGF-beta/adenosine panels; typography and composite
layout edits are not claimed to be reproduced pixel-identically.

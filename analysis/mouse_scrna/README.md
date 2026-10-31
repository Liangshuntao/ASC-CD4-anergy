# Mouse single cell workflow

Recommended entry point: run_mouse_analysis.R. Historical exploratory scripts
are excluded from this figure-code package.

```
Rscript analysis/mouse_scrna/run_mouse_analysis.R --input /path/seu_scWAT.rds --out /path/new_output --threads 4 --preflight
Rscript analysis/mouse_scrna/run_mouse_analysis.R --input /path/seu_scWAT.rds --out /path/new_output --threads 4
```

Stages run in this order: compute signatures/activities/trajectory, finalize figure
panels from checkpoints, then LIANA ASC-to-CD4 communication. --resume retains the
original checkpoint logic. The runner rejects nonempty output unless --resume is
specified. No package installation is performed. See ../../DEPENDENCIES.md for
package requirements and historical-versus-reproduction version boundaries.

The required input is the historically annotated Seurat object, not a raw matrix.
Its annotation metadata, sample labels and inherited embeddings are required.
The stage scripts preserve historical state gene sets, algorithm parameters,
Monocle3 root rule, cell counts, tests and pooling. Some analyses are descriptive
cell-level comparisons rather than biological-sample-level tests.

Figure 2H uses the correlation matrix between state scores and pathway activities.
The figure mapping identifies the selected LR/BMP/TGF-beta/adenosine panels.

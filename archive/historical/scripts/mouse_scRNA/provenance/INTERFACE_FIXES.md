# Runnable submission: documented changes

Authoritative historical archive: `GSE183288_manuscript_single_cell_EXACT_EXTRACT.R`  
Runnable submission version: `GSE183288_manuscript_single_cell_RUNNABLE.R`

No gene set, biological label, score definition, cutoff, statistical test, trajectory-root rule, LIANA resource, plot color, plot dimension, point size, heatmap setting, or manuscript conclusion was intentionally changed.

| ID | Type | Change | Reason |
|---|---|---|---|
| RUN-01 | Input interface | Added `--input`, defaulting to `D:/R/20250926/seu_scWAT.rds`. | Remove dependence on the current working directory. |
| RUN-02 | RDS interface | Removed historical `load("seu_scWAT.rds")`; retained `readRDS()` semantics. | The canonical object was written with `saveRDS()` and cannot be opened with `load()`. |
| RUN-03 | Output interface | Added `--out`; all historical relative output names are written below this directory. | Preserve `D:/R/20250926` as read-only. |
| RUN-04 | Environment | Removed `pak::cache_clean()`, `pak::pkg_install()` and runtime `install.packages()`/`BiocManager::install()`. Added a fail-fast package check. | Submission code must not mutate the package library or require network access. |
| RUN-05 | Dependency order | Moved GSVA/ssGSEA, AUCell and AddModuleScore into the Figure 1 upstream scoring phase; their S1 heatmaps remain in Figure S1. | These scores are part of the Figure 1 analysis and are dependencies of later Figure 2 calculations. |
| RUN-06 | Dependency order | Moved TF/PROGENy, Monocle3, ASC-CD4 object construction and LIANA aggregation before their first consuming Figure panel. | The historical interactive workspace supplied these objects implicitly. |
| RUN-07 | Seurat interface | Replaced the invalid `names(RNA@data)` test with a read-only `GetAssayData(..., layer="scale.data")` row-count test. | Check the same intended condition using an explicit supported accessor. |
| RUN-08 | LIANA interface | Treat `liana_wrap()` as the method-result list and call `liana_aggregate()` before tidy filtering. | Current LIANA 0.1.14 does not support filtering the unaggregated list as a tibble. |
| RUN-09 | Figure order | Figure 3E saves only its historical `CD4_ASC_UMAP.pdf`; the historical composite requiring the later `p.L.R` object is not executed there. | Keep Figure display order without introducing a forward workspace dependency. |
| RUN-10 | Validation | Added object class, UMAP, metadata, 58,241-cell, 13-sample, CD4 and ASC count checks plus `--preflight`. | Fail before expensive computation when the canonical input is wrong. |
| RUN-11 | Seurat v5 interface | Replaced defunct `GetAssayData(..., slot="data")` calls with `layer="data"`. | SeuratObject 5.0+ rejects the old `slot` argument; the requested normalized-expression layer is unchanged. |
| RUN-12 | LIANA input object | Convert the unchanged ASC/CD4 Seurat subset to `SingleCellExperiment`, set `colLabels` from `LR_group`, and run the same `MouseConsensus` analysis with `sca`, `cellphonedb`, `natmi`, `connectome`, and `logfc`. | LIANA's Seurat accessor path calls the removed `GetAssayData(slot=)` interface under SeuratObject 5; the SCE path was independently executed and validated. |
| RUN-13 | LIANA command interface | Replaced the CodeBuddy test machine's absolute input/output paths with `--input` and `--out`. | Permit clean-session execution outside the original workstation directory without changing analysis content. |

The historical cell-level Wilcoxon/correlation tests, pooled LIANA design, mixed SC/TC/SH/TH scope, and ACT-minus-ANERGY Monocle root remain unchanged. They are methodological caveats, not interface fixes.

## Verification status (2026-07-03)

- Static R parsing: PASS.
- Clean-session `--preflight`: PASS; canonical Seurat class, UMAP, metadata, 58,241 cells, 13 samples, 14,195 CD4 Tconv and ASC counts 277/1,752/4,641 matched.
- Full run: INCOMPLETE, not failed analytically. The external runner stopped after its 3,600-second limit while historical ssGSEA was at iteration 24. GSVA had completed. The log reported 10 workers and a 1.8-GiB sparse-to-dense allocation.
- No Figure output from this timed run may be labelled reproduced. A longer uninterrupted run is required.

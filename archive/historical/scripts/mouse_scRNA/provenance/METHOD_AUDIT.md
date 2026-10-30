# Method audit — do not silently modify strict mode

- The canonical object pools SC/TC/SH/TH while the manuscript uses “steady-state” language.
- `active_ident` and `ASC_CD90_CD105` are recovered manual annotations; their generation code is unavailable.
- Legacy CD4 trajectory uses the inherited atlas UMAP rather than a clearly documented secondary reclustering.
- The Monocle3 root is the top 5% of `ACT_GSVA - ANERGY_GSVA`, creating outcome-informed directionality.
- ANERGY high/low Wilcoxon tests and TF/pathway correlations treat cells as observations; inference is pseudoreplicated.
- Pooled LIANA does not model mouse-to-mouse variability.
- The original exploratory script installed packages during execution, used incorrect `load(.rds)`, depended on workspace objects and did not consistently lock seeds.

Strict mode reproduces these choices so they can be audited. Any corrected analysis must be a separately named version and must not overwrite strict outputs.


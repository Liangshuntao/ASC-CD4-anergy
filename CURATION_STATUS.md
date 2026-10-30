# Submission code curation status

The full requested curation is IN PROGRESS, not complete.

## Completed and checked

- Historical source was moved to archive/historical without changing source bytes.
  All 53 source-file SHA256 hashes match the initial archive manifest.
- Bulk: seven functional stage filenames and an explicit input/output runner;
  the exploratory t-test is not in the recommended workflow. Duplicate result
  exports and nonselected plotting variants remain only in historical source.
- Ten bulk source blocks match their original parsed R expressions after the
  documented I/O substitutions. No statistical formula or parameter was changed.
- The primary bulk model was run on the existing 12-sample matrix in an isolated
  validation output. All seven numerical fields for all genes in each of three
  historical result tables matched exactly; significance labels also matched.
- Mouse: the runner and three stages were functionally renamed. Internal stage
  references were updated and private default input paths removed. Exact allowed
  text transformations were checked. The full single-cell workflow was not rerun.
- Seventeen recommended R files parse. Windows locale/Unicode warnings persist;
  parsing is not runtime validation of every stage.
- The final PRJNA932038 candidate-evidence figure has an explicit input object,
  annotation map and separate output interface. Paths/dependencies preflight passed;
  preprocessing curation and complete figure replay validation remain pending.

## Work still required

- Current Figure 8A variance has now been checked: the current input and unchanged
  PCA code yield 89.26% / 4.42%, matching the figure. Earlier audit values were stale.
- Reviewer-response bulk scripts and execution order are curated; runtime-result
  comparison and complete dependency/input checks remain pending.
- Curate human single-cell workflows with reviewed annotation/config dependencies.
- Recover final proteomic blocks rather than promote the interactive historical script.
- Reconcile clinical and METABRIC source blocks and the custom score naming issue.
- Complete exact final-panel lineage, input routes, package requirements and
  interface tests for all modules before a final submission release.

The bulk numerical comparison verifies the primary model only. It does not prove
visual equivalence of the current manuscript panels, full pathway reproducibility,
or validity of unrelated analyses. Original calculations and study design remain
unchanged, including their existing limitations.

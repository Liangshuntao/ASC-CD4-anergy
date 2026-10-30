# Code curation and validation scope

Updated 31 October 2026. Code curation is complete locally; remote publication
verification is recorded separately in the local project handoff.

## Changes

- Historical source is separated under archive/historical, with original source
  bytes and Git history preserved.
- Six modules have functional English filenames, explicit input/output interfaces,
  dependency checks and stage order. Internal calls and the 24-figure map use new paths.
- Recommended source contains no package installation or explicit download calls.
- Source records document path, encoding, matrix-access and diagnostic-label repairs.
  Scientific signatures, units, contrasts, thresholds and calculations were not redesigned.

## Checked evidence

- Three bulk primary tables: all seven numeric fields and significance labels
  matched historical results exactly. PCA: 89.26% / 4.42%, matching current Figure 8A.
- Four proteomic protein/pathway tables: all seven numeric fields matched exactly.
  Main and supplementary plot stages executed from the saved processed checkpoint.
- Bulk/protein/clinical blocks passed expression comparisons after listed edits.
  Mouse/human source bodies passed exact interface-transformation checks.
- All recommended R files parse. Nine main/import interfaces passed --help tests.
  Mouse preflight, bulk supplementary first-stage preflight, human-adipose preflight,
  human-breast preflight and helper/config smoke checks passed. Annotation import
  preserves cell order and rejects incomplete/duplicate/count-mismatched maps.
- Small validation summaries are under validation/; no expression matrix is included.

## Limits

The complete large single-cell/clinical calculations were not rerun. Clinical
preflight found missing local survminer and ggpubr; no installation was performed.
Exact original input/reference versions and reviewed annotations remain necessary.
Windows locale/Unicode warnings and separately edited manuscript layouts mean
pixel-identical figure reproduction is not certified.

The raw proteomic script had an assignment error and unrecorded imputation RNG;
the recommended replay explicitly uses the original processed checkpoint.
The custom clinical score is not official TIDE; median labels are not validated
therapy-response predictions. Curation does not erase these evidence boundaries.

No final scientific submission tag or Release was created by this task.

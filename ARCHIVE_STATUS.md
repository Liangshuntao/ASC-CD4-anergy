# Initial archive status

This initial public archive is provided for code inspection. It is not tagged
v1.0-submission and does not establish full end-to-end reproducibility.

## Checks performed

- 53 source files were collected, including 41 R scripts; source hashes are in code_manifest.csv.
- All 41 R scripts passed R parse() syntax checks. Windows locale and character
  conversion warnings occurred; syntax success does not establish runtime success.
- A credential-pattern scan found no matches. This is not a comprehensive security audit.
- An initial mapping covers 14 main and 10 supplementary figures. Some exact
  panel-to-script relationships remain unresolved and are labelled accordingly.
- No analysis code, package installation or data download was executed during packaging.
- Raw sequencing data, expression matrices, large intermediate objects, private
  reviewer correspondence and local packaging logs are not included.

## Known limitations

- Historical paths and intermediate-object dependencies remain in the source scripts.
- Mouse single-cell analysis starts from inherited manual annotations; the original
  annotation-generation workflow was not recovered.
- Human workflows require original project inputs, mappings and intermediate objects.
- Proteomic source is historical interactive code, not a validated batch pipeline.
- Exact METABRIC source attribution and several input access routes remain unresolved.
- The score called TIDE in the historical clinical script is a custom score, not
  official TIDE; manuscript terminology needs reconciliation before a final release.
- No open-source licence is assigned pending author/institution confirmation.

The initial archive is not a promise that all manuscript results can be regenerated
from the public raw datasets with a single command. Subsequent corrections should
be documented in new commits; original scientific calculations are not changed here.

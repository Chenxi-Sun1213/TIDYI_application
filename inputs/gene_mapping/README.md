# Frozen replay3 mouse gene mapping

Except for Objective 2, formal mouse symbol-to-Ensembl conversion reads
`mouse_symbol_to_ensembl.csv` through the shared `gene_mapping.py` module.
The file is an exact frozen copy of the archived MyGene cache at
`/scratch/chenxi.sun/.tidyi_gene_cache/mouse_symbol_to_ensembl.csv`.
It uses the original columns `symbol` and `ensembl_gene_clean`.

The project owner selected the replay3 policy on 2026-10-06:

- Use the archived cache alone, without AER/BioMart supplementation or online queries.
- Translate symbols with a valid mouse Ensembl ID; retain unresolved symbols.
- Preserve existing valid input Ensembl IDs and retain all expression features.
- Make duplicate mapped names unique; do not merge their expression counts.

The file contains 29,857 unique symbols: 25,428 valid mouse mappings,
4,428 empty entries, and one non-Ensembl value (`Scaper -> 244891`).
The loader ignores that invalid value, retaining `Scaper` as its original symbol,
matching replay3. The source annotation release was not recorded.
`mapping_provenance.json` records the source and file checksum.

Objective 1, the AER preparation in Objective 5, and the Objective 6 support
pipeline share this file. AER retains existing valid original Ensembl IDs and
uses the cache only where an ID is absent. Objectives 3 and 6 consume saved rate
tables in their normal master runs; this mapping update does not regenerate those
inputs. Arigoni human genes retain their original IDs.

Objective 2 was restored on 2026-10-12 to
`inputs/classifier_data/mart_export_gene_names.txt` and the historical rule
that retains only genes present in that table before rate inference.

`build_reference.py` copies the archived cache and rebuilds its provenance. It is
an optional maintenance utility, never called automatically by a master.
Historical `reference_outputs/` retain their original gene-mapping results.

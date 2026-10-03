# Frozen mouse gene mapping

All handover analyses that translate mouse gene symbols read
`mouse_symbol_to_ensembl.csv` through the shared `gene_mapping.py` module.
Analysis runs do not query a service, update a cache, or download annotations.

The reference contains 41,870 symbols: 39,948 have a resolved mouse Ensembl ID,
1,884 have no valid candidate, and 38 remain ambiguous. The columns are
`symbol`, `ensembl_gene`, `status`, `sources`, and `candidate_ensembl_ids`.
An empty `ensembl_gene` means the symbol must not be translated automatically.

The project owner confirmed this resolution policy on 2026-10-03:

1. Keep a unique valid result from the archived TIDYI mouse mapping cache.
2. For other symbols, supplement from the AER original 10x feature metadata and
   the retained BioMart human–mouse orthologue export. Resolve only when the
   local sources provide exactly one distinct valid mouse ID.
3. Preserve existing input Ensembl IDs. Do not replace them using symbol lookup.
4. Leave unresolved symbols unchanged in Objectives 1 and 6. Objective 2 keeps
   its existing mapped-gene filtering step; AER keeps only valid mouse IDs and
   preserves its uniqueness and coverage checks.

There are 162 symbols with multiple candidate IDs across sources. The 124 with
a unique legacy-cache result retain that result; the remaining 38 are left
unresolved. Candidate IDs are retained in the CSV for inspection. The legacy
`Scaper` entry contained the Entrez ID `244891`; that invalid Ensembl value was
excluded and replaced by the unambiguous local Ensembl candidate.

The source annotation release was not recorded in the supplied files. This is
a frozen reference assembled from existing local sources, not a claim about
the latest Ensembl release. `mapping_provenance.json` records source paths and
SHA-256 checksums. Running an analysis requires only the frozen CSV, not the
original files used to build it. `build_reference.py` is an optional rebuild
script and is never called by a master.

`mapping_coverage.csv` records a metadata-only audit of all 37 ASTRID H5AD
inputs, before downstream filtering: Mukin BM 98.49%, Mukin skin 98.12%, and
TMS BM 93.74% of feature names resolve. These percentages describe gene-name
coverage, not cell retention or classifier marker coverage.

Objective 2 previously used only the BioMart orthologue subset. The unified
reference can retain additional genes or choose different IDs, so rerun rate
values may differ. Saved notebook outputs and `reference_outputs/` remain
historical results and have not been regenerated.

Offline regression checks, without fitting TIDYI or running the six analyses:

```bash
cd /scratch/chenxi.sun/RA_project_handover
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
```

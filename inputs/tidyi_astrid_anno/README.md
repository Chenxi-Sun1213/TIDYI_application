# ASTRID-annotation-based TIDYI rates

TIDYI rates inferred from curated data based on ASTRID annotation (pooled sample).

These five historical CSVs are the default rate inputs for Objective 6 model
fitting. Cells are pooled within each comparison group and tissue before
cell-type-level rate inference; the rows are not per-animal estimates.
The ASTRID comparison workflow includes its saved cell-type correction rules.

The five files cover Mukin young SPF versus old SPF, young SPF versus young GF,
old SPF versus old GF, and TMS male versus female in young and old animals.
Each CSV contains 13 columns, including conditions, tissue, cell type,
proliferation/death/growth rates, turnover, and cell percentage.

The files in [tidyi_refine](../tidyi_refine/README.md) use an author-annotation-based
lineage. Shared filenames and columns do not make the two rate sets
interchangeable.

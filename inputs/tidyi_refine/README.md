# Author-annotation-based TIDYI rates

TIDYI rates inferred from curated data based on author annotation.

These five historical CSVs are the default saved rate inputs for Objective 3
when `--use-precomputed-rates` is selected. A full Objective 3 run regenerates
rates under `03_death_proliferation/outputs/rate_generation/ty_df/` instead.

The historical pipeline starts from `tissue_wise_knn` annotations in the curated
Mukin/TMS H5ADs. Mukin additionally applies comparison-specific constrained kNN
relabeling and entropy filtering; TMS does not receive that additional step.
The author-annotation-based description therefore identifies the lineage,
rather than implying that the original labels were used without processing.

Cells are pooled within each comparison group and tissue before cell-type-level
rate inference. The five files cover Mukin young SPF versus old SPF, young SPF
versus young GF, old SPF versus old GF, and TMS male versus female in young and
old animals. Each CSV contains 13 columns, including conditions, tissue, cell
type, proliferation/death/growth rates, turnover, and cell percentage.

The files in [tidyi_astrid_anno](../tidyi_astrid_anno/README.md) use a different
annotation lineage and are not interchangeable with these tables.

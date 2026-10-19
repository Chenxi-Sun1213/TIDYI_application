# Objective 3 rate-generation inputs

The Objective 3 master uses these versioned sample-statistics tables together
with the authoritative external single-cell matrices:

- `mukin_sample_stas.csv`: copied from
  `thesis_project/TIDYI/TIDYI_pattern/mukin_sample_stas.csv`.
- `tms_sample_stas.csv`: copied from
  `thesis_project/TIDYI/TIDYI_pattern/tms_sample_stas.csv`.
- Mukin H5AD (external, not moved):
  `thesis_project/mukin_qcdata/mukin_data_tissuewise_knn.h5ad`.
- TMS H5AD (external, not moved):
  `thesis_project/tms_qcdata/tms_10x_tissuewise_knn.h5ad`.

The master invokes `support/tidyi_comparison_analysis.py` for the five formal
comparisons. It starts from `tissue_wise_knn`, uses the `raw` layer, and
applies an 80-cell minimum separately within each comparison group, tissue, and
fine `analysis_celltype`. The later `celltye_cross-tissue` aggregation is not
part of the fixed Objective 3 rate lineage. Gene conversion uses the handover's
frozen offline mouse symbol-to-Ensembl mapping.

For Mukin only, four state-like query labels undergo 10-neighbor transfer in 30
PCA dimensions. A transferred cell is retained only when normalized entropy is
strictly below 0.3 and its target is in the query label's whitelist. TMS does
not receive this additional entropy filter.

The generated rate CSVs are written below
`outputs/rate_generation/ty_df/`; this directory is generated output, not an
authoritative input.

## Reference CV generation

The master also invokes `support/reference_cv_generation.py`. It uses the
same saved QC sample tables to select Mukin 2m SPF and TMS 1m male reference
samples. Reference rates come from `live-dead_annotation` in the original
external matrices `thesis_project/mukin_qcdata/mukin_data.h5ad` and
`thesis_project/tms_qcdata/tms_10x.h5ad`. Historical label merges and a 100-cell
minimum apply per sample and cell type. The saved rate-group sets establish
this cutoff; the historical notebook's current default of 10 differs.

The tissuewise H5ADs listed above separately supply `tissue_wise_knn` cell
counts for the historical matching rule. Their count labels are not substituted
for the reference rate labels. The two CV tables are regenerated under
`outputs/reference_cv/` using the TIDYI selected by the master, together with
sample rates, exclusions, and a manifest. Frozen CV tables under
`inputs/tidyi_statistics/` remain for explicit precomputed compatibility runs.

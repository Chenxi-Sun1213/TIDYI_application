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

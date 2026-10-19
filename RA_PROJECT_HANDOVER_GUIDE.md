# RA project handover guide

> Prepared: 2026-10-02  
> Scope: six formal RA analyses

## 1. Handover rule

Each formal objective has exactly one `master_*.py` entry point. The master
declares its inputs, selects the audited notebook cells or standalone analysis,
creates outputs, verifies that every declared output was refreshed, and writes
a JSON run manifest.

The original notebooks are retained as analysis sources. Exploratory cancer,
developmental, and zero-death branches outside the six objectives are not part
of this handover.

## 2. Formal objectives

| ID | Objective | Master entry point |
|---|---|---|
| 1 | TIDYI rate reproducibility across replicates | `01_rate_reproducibility/master_rate_reproducibility.py` |
| 2 | Mouse-human skin rate comparison | `02_mouse_human_skin/master_mouse_human_skin.py` |
| 3 | Death-proliferation dominance | `03_death_proliferation/master_death_proliferation_dominance.py` |
| 4 | Arigoni SciData 2024 cancer validation | `04_arigoni_cancer/master_arigoni_scidata_2024_validation.py` |
| 5 | AER developmental validation | `05_aer_developmental/master_aer_developmental_validation.py` |
| 6 | Systems biology model fitting | `06_systems_biology/master_systems_biology_model_fitting.py` |

## 3. Data locations and lineage

### 3.1 Shared large data

| Dataset | Authoritative location | Used by | Status |
|---|---|---|---|
| Arigoni BE1run12 raw 10x matrices | `/srv/mfs/hausserlab/data/ARigoni_Scidata_2024/BE1run12` | 4 | already present |
| AER reconstructed raw counts | `/srv/mfs/hausserlab/data/AER2025_Juliane/aer_six_samples_raw_counts.h5ad` | 5 | copied and SHA-256 verified |
| AER post-entropy checkpoint | `/srv/mfs/hausserlab/data/AER2025_Juliane/aer_post_entropy.h5ad` | 5 | copied and SHA-256 verified |

The AER post-entropy checkpoint uses `cell_type` as the downstream annotation
column. It was produced after QC, manual ASTRID cluster review, per-sample kNN
transfer, and entropy filtering at `< 0.3` for transferred query cells.

### 3.2 Mukin and TMS data retained in place

At the project owner's request, the following data are not moved in this phase:

```text
/scratch/chenxi.sun/thesis_project/tms_qcdata/AllRef_BM_0204
/scratch/chenxi.sun/thesis_project/mukin_qcdata/AllRef_BM_0204
/scratch/chenxi.sun/thesis_project/mukin_qcdata/AllRef_skin_0204
/scratch/chenxi.sun/thesis_project/mukin_qcdata/mukin_data_tissuewise_knn.h5ad
/scratch/chenxi.sun/thesis_project/tms_qcdata/tms_10x_tissuewise_knn.h5ad
```

Objectives 1 and 2 read one `*_outASTRID.h5ad` per sample. Objective 1 uses
`celltype_hb` for TMS/Mukin BM and `SingleR_CellType` for Mukin skin. Objective
2 starts from `SingleR_CellType` and constructs `kc_human_validation` using the
documented Krt10/Krt15 rule. Objective 3 now regenerates its five pooled-condition
rate tables from the two tissuewise H5ADs, starting from `tissue_wise_knn`;
Mukin then applies comparison-specific constrained kNN relabeling. Rates remain
grouped by fine `analysis_celltype`; the later cross-tissue aggregation is not
used for Objective 3. Its reference CV producer separately infers sample-level
rates from the original `mukin_data.h5ad` and `tms_10x.h5ad` annotations and uses
the tissuewise H5ADs for cell-count matching. Both producers are called directly
by the master; historical notebooks are not executed. Objective 6
uses the same tissuewise H5ADs for tissuewise abundance and the ASTRID BM
objects for mouse-level abundance.

### 3.3 Lightweight inputs copied into the handover

`inputs/` contains two distinct five-table TIDYI rate sets:
[author-annotation-based rates](inputs/tidyi_refine/README.md) for Objective 3
compatibility runs and [ASTRID-annotation-based pooled rates](inputs/tidyi_astrid_anno/README.md)
for Objective 6. It also contains the shared offline gene mapping and the
Objective 2 classifier gene-name map. Reference CV tables are generated under
`03_death_proliferation/outputs/reference_cv/`, rather than bundled in `inputs/`.
Objective 3 also contains its historical rate producer, adapted to the shared
frozen offline gene mapping, plus exact copies of the two small sample-statistics
tables; the large Mukin/TMS H5ADs remain at their authoritative paths.

### 3.4 Shared offline mouse gene mapping

Except for Objective 2, mouse symbol-to-Ensembl conversion uses
`inputs/gene_mapping/mouse_symbol_to_ensembl.csv`, loaded by `gene_mapping.py`.
The owner selected the replay3 archived-cache-only policy on 2026-10-06.
The file is an exact copy of the old cache with 25,428 valid mouse mappings;
there is no AER/BioMart supplementation. Unresolved symbols are retained in the
expression matrices. Source checksums and the resolution policy are
documented in `inputs/gene_mapping/README.md` and `mapping_provenance.json`.

Objective 1, Objective 3 rate generation, AER preparation in
Objective 5, and the Objective 6 support script share this reference. All mouse
symbol mapping in these paths is offline. Existing AER mouse IDs and Arigoni
human IDs retain their input identity. Objective 6 continues to use saved rate
tables. Missing reference files fail locally instead of triggering an online
lookup.

On 2026-10-12, Objective 2 was restored to its historical mapping and filtering:
`inputs/classifier_data/mart_export_gene_names.txt` maps mouse symbols, and only
genes present in that table are retained before TIDYI inference. The local table
is byte-identical to the historical installed classifier-data table. This
restores the historical normalization denominator. Other objectives retain the
replay3 policy.

## 4. Analysis units and denominators

| Objective | Analysis unit / denominator |
|---|---|
| 1 | TIDYI is run per animal-cell type after dataset-specific filtering; replicate comparisons use matched rate rows. Minimum group size is 100 cells. |
| 2 | Mouse keratinocyte TIDYI groups are mouse × `kc_human_validation`; minimum group size follows the retained notebook logic. |
| 3 | Rate generation pools cells by comparison group × tissue and calculates cell-type percentages before the per-group 80-cell filter. Mukin query labels require normalized entropy `< 0.3` and whitelist approval; TMS has no analogous entropy filter. Each downstream row is a pooled tissue-cell-type contrast. Method 3 uses `Δgrowth = Δproliferation - Δdeath`; the regulation threshold is `2 ln(2)/(24×3)`. |
| 4 | Each TIDYI group is an Arigoni cancer cell line; minimum group size is 100 cells. |
| 5 | TIDYI groups are sample × `cell_type`; groups must contain more than 50 cells. AER Method 3 comparisons are stage/genotype contrasts; its weight uses Mukin typical CVs and unweighted mean rates across the available AER sample-cell-type rows. |
| 6 | Abundance is dataset × mouse × tissue × cell type; model fitting uses eight pooled conditions per selected cell type. |

## 5. Outputs and references

Each objective has two separate directories:

- `outputs/`: empty before the first handover run; the master writes here.
- `reference_outputs/`: pre-handover results copied for comparison.

Never point a master at `reference_outputs/`. A successful run must generate
all outputs declared near the top of the master and write a run manifest in
`outputs/`.

The pre-handover references do not contain `tidyi_runtime_provenance.json`,
because the old notebook runs did not create it. The Objective 4 consolidated
main figure is generated by its master.

## 6. Manual execution

Use the existing shared-server environment and run one master at a time:

```bash
python "/scratch/chenxi.sun/RA_project_handover/01_rate_reproducibility/master_rate_reproducibility.py"
python "/scratch/chenxi.sun/RA_project_handover/02_mouse_human_skin/master_mouse_human_skin.py"
python "/scratch/chenxi.sun/RA_project_handover/03_death_proliferation/master_death_proliferation_dominance.py"
python "/scratch/chenxi.sun/RA_project_handover/04_arigoni_cancer/master_arigoni_scidata_2024_validation.py"
python "/scratch/chenxi.sun/RA_project_handover/05_aer_developmental/master_aer_developmental_validation.py"
python "/scratch/chenxi.sun/RA_project_handover/06_systems_biology/master_systems_biology_model_fitting.py"
```

Objectives 1, 2, 3, 4, and 5 import TIDYI from the active environment unless
`--tidyi-source-root` is provided. Objective 3 generates its five rate tables
under `03_death_proliferation/outputs/rate_generation/ty_df/` and recomputes its
reference CV tables under `03_death_proliferation/outputs/reference_cv/` before
running the dominance analysis. Its explicit `--use-precomputed-rates`
compatibility mode skips both producers and reads `--tidyi-input-root` plus
`--reference-cv-root` (the latter defaults to its generated `outputs/reference_cv/`
directory). Objective 5 reads the typical-CV table from the same Objective 3
output directory, so run Objective 3 first or transfer those generated CV
outputs with the project. Objective 6 continues to use copied rate tables with
the same 13-column contract.

## 7. Known boundaries

- Mukin/TMS large inputs still use their existing absolute paths by decision.
- AER starts from the retained post-entropy checkpoint and does not rerun Cell
  Ranger, ASTRID, manual cluster review, or kNN label transfer.
- The deleted CCLE and GSE157220/Mix-seq raw H5AD files belonged to historical,
  out-of-scope cancer analyses and are not inputs to Objective 4.
- A real replacement TIDYI should be connected through the existing command
  arguments. Do not modify downstream scientific logic unless its observed API
  or schema actually differs.

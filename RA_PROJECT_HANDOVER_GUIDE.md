# RA project handover guide

> Prepared: 2026-10-02  
> Scope: six formal RA analyses  
> Execution status: assembled and statically checked; not run from this directory

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
documented Krt10/Krt15 rule. Objective 6 uses `tissue_wise_knn` in the two
tissuewise H5ADs for tissuewise abundance and the ASTRID BM objects for
mouse-level abundance.

### 3.3 Lightweight inputs copied into the handover

`inputs/` contains the five legacy TIDYI rate tables for objectives 3 and 6,
the two Method 3 CV tables, and the classifier gene-name map. These copies keep
the formal masters independent of the old RA-project directory while
preserving the original input content.

## 4. Analysis units and denominators

| Objective | Analysis unit / denominator |
|---|---|
| 1 | TIDYI is run per animal-cell type after dataset-specific filtering; replicate comparisons use matched rate rows. Minimum group size is 100 cells. |
| 2 | Mouse keratinocyte TIDYI groups are mouse × `kc_human_validation`; minimum group size follows the retained notebook logic. |
| 3 | Each row is a condition-level pooled tissue-cell-type contrast. Method 3 fractions use `Δgrowth = Δproliferation - Δdeath`; the regulation threshold is `2 ln(2)/(24×3)`. |
| 4 | Each TIDYI group is an Arigoni cancer cell line; minimum group size is 100 cells. |
| 5 | TIDYI groups are sample × `cell_type`; groups must contain more than 50 cells. AER Method 3 comparisons are stage/genotype contrasts. |
| 6 | Abundance is dataset × mouse × tissue × cell type; model fitting uses eight pooled conditions per selected cell type. |

## 5. Outputs and references

Each objective has two separate directories:

- `outputs/`: empty before the first handover run; the master writes here.
- `reference_outputs/`: pre-handover results copied for comparison.

Never point a master at `reference_outputs/`. A successful run must generate
all outputs declared near the top of the master and write a run manifest in
`outputs/`.

The pre-handover references do not contain `tidyi_runtime_provenance.json`,
because the old notebook runs did not create it. Objective 4 also has no
reference copy of the newly declared consolidated main figure. These files must
be created by the first master run, but their prior absence is not a failed
reference-copy check.

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

Objectives 1, 2, 4, and 5 import TIDYI from the selected Jupyter kernel unless
`--tidyi-source-root` is provided. Objectives 3 and 6 default to the copied
legacy rate tables and accept `--tidyi-input-root` for a replacement set with
the same 13-column contract.

## 7. Acceptance checklist

For each objective:

1. Run only its master entry point from a fresh shell or kernel environment.
2. Confirm that every declared input exists before computation starts.
3. Confirm every declared output exists, is non-empty, and has a new timestamp.
4. Inspect the JSON manifest for input paths, TIDYI version labels, annotation
   columns, thresholds, and output checksums.
5. Compare output tables against `reference_outputs/` using their biological
   keys and numeric values, not row order alone.
6. Compare every formal figure visually with its reference counterpart.
7. Record any expected difference caused by an explicitly changed TIDYI model.
8. Do not delete old source data or the historical project until all six
   objectives have passed.

## 8. Known boundaries

- The new masters have not yet been executed; current validation is static.
- Mukin/TMS large inputs still use their existing absolute paths by decision.
- AER starts from the retained post-entropy checkpoint and does not rerun Cell
  Ranger, ASTRID, manual cluster review, or kNN label transfer.
- The deleted CCLE and GSE157220/Mix-seq raw H5AD files belonged to historical,
  out-of-scope cancer analyses and are not inputs to Objective 4.
- A real replacement TIDYI should be connected through the existing command
  arguments. Do not modify downstream scientific logic unless its observed API
  or schema actually differs.

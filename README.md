# TIDYI application handover

Six RA analyses, each with one `master_*.py` entry point. Masters run the
required notebook cells or Python workflows, check declared outputs, and
write run manifests.

| Project | Entry point | Analysis |
|---|---|---|
| 01 | [Rate reproducibility](01_rate_reproducibility/master_rate_reproducibility.py) | Mukin/TMS replicate rate comparisons |
| 02 | [Mouse–human skin](02_mouse_human_skin/master_mouse_human_skin.py) | Keratinocyte rate validation |
| 03 | [Death–proliferation dominance](03_death_proliferation/master_death_proliferation_dominance.py) | H5AD → rates and reference CVs → Method 3 ([details](03_death_proliferation/README.md)) |
| 04 | [Arigoni cancer](04_arigoni_cancer/master_arigoni_scidata_2024_validation.py) | Arigoni SciData 2024 cell-line validation |
| 05 | [AER development](05_aer_developmental/master_aer_developmental_validation.py) | Post-entropy checkpoint → BioLUCID, TIDYI, Method 3 and figures ([details](05_aer_developmental/README.md)) |
| 06 | [Systems biology](06_systems_biology/master_systems_biology_model_fitting.py) | Abundance and model fitting |

## Run

Run one master at a time in the existing `scx_py310` environment. For example,
from this directory:

```bash
conda activate scx_py310
python 01_rate_reproducibility/master_rate_reproducibility.py --kernel-name scx_py310
```

Notebook-based masters accept `--kernel-name`; all masters provide `--help`.
Objective 3 regenerates rates and reference CVs by default;
`--use-precomputed-rates` selects its saved-input compatibility mode.
Objective 5 reads the CV table generated in `03_death_proliferation/outputs/reference_cv/`;
run Objective 3 first or transfer its generated CV outputs with the project.
Use `--tidyi-source-root` to select another TIDYI implementation for Objectives
01–05, with version/model/formula labels recorded in the manifest.

## Files and data

- Keep notebooks and helpers: they contain calculations used by the masters,
  including `01_rate_reproducibility/rates_reproduce_test.ipynb`.
- Shared `inputs/` contains frozen mappings and two distinct legacy rate sets:
  [author-annotation-based rates](inputs/tidyi_refine/README.md) for Objective 3
  compatibility runs, and [ASTRID-annotation-based pooled rates](inputs/tidyi_astrid_anno/README.md)
  for Objective 6. Mouse mapping uses the frozen replay3 cache; Objective 2
  retains its legacy mart mapping and gene filter.
- Large inputs remain in `thesis_project` and `/srv/mfs/hausserlab/data/`.
  See [the handover guide](RA_PROJECT_HANDOVER_GUIDE.md) for paths and provenance.
- `outputs/` holds generated results and run records; `reference_outputs/`
  holds historical comparison results and must remain unchanged. Both are
  Git-ignored and need separate transfer when handing over existing results.

After code changes, compare fresh results with references and review figures;
static checks alone do not establish scientific reproducibility.

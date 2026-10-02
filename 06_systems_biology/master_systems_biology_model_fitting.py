#!/usr/bin/env python3
"""Master entry point for Objective 6: systems biology model fitting.

Legacy run using the current pooled TIDYI rate directory:

    python master_systems_biology_model_fitting.py

Run with rate CSVs produced by a colleague's replacement TIDYI:

    python master_systems_biology_model_fitting.py \
        --tidyi-input-root /path/to/new/ty_df \
        --tidyi-version-label <version-or-commit> \
        --tidyi-model-label <model-name> \
        --tidyi-formula-label <formula-name>

The replacement directory must contain the same five filenames and the same
13-column TIDYI contract used by the legacy analysis.

Inputs:
- ``inputs/tidyi_astrid_anno`` (copied pooled-condition rate CSVs).
- thesis Mukin/TMS tissuewise H5AD files.
- thesis Mukin/TMS AllRef_BM_0204 ASTRID directories.
- local support/tidyi_astrid_comparison_analysis.py, which contains the TMS prefix and Mukin
  RBC-contamination correction rules used for the abundance table.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


RA_ROOT = Path(__file__).resolve().parents[1]
if str(RA_ROOT) not in sys.path:
    sys.path.insert(0, str(RA_ROOT))

from handover_utils import SourceReplacement, run_notebook_master  # noqa: E402


SOURCE_NOTEBOOK = Path(__file__).with_name("system_biology.ipynb")
OUTPUT_DIR = Path(__file__).parent / "outputs"
LEGACY_TIDYI_INPUT_ROOT = RA_ROOT / "inputs/tidyi_astrid_anno"
PIPELINE_HELPER = Path(__file__).parent / "support/tidyi_astrid_comparison_analysis.py"
MUKIN_TISSUEWISE_H5AD = Path(
    "/scratch/chenxi.sun/thesis_project/mukin_qcdata/mukin_data_tissuewise_knn.h5ad"
)
TMS_TISSUEWISE_H5AD = Path(
    "/scratch/chenxi.sun/thesis_project/tms_qcdata/tms_10x_tissuewise_knn.h5ad"
)
MUKIN_ASTRID_ROOT = Path(
    "/scratch/chenxi.sun/thesis_project/mukin_qcdata/AllRef_BM_0204"
)
TMS_ASTRID_ROOT = Path(
    "/scratch/chenxi.sun/thesis_project/tms_qcdata/AllRef_BM_0204"
)

RATE_FILENAMES = (
    "mukin_young_spf_vs_old_spf_tidyi.csv",
    "mukin_young_spf_vs_young_gf_tidyi.csv",
    "mukin_old_spf_vs_old_gf_tidyi.csv",
    "tms_young_male_vs_young_female_tidyi.csv",
    "tms_old_male_vs_old_female_tidyi.csv",
)

SELECTED_CELLS = (2, 4, 12, 13, 15, 17, 18, 21, 23, 26)

REFINED_CELLTYPE_SUFFIXES = (
    "b_pre",
    "dendritic",
    "monocyte_pro",
    "myeloid_progenitor",
    "myelomonocyte",
    "t_cd8",
)
EXPECTED_OUTPUTS = (
    OUTPUT_DIR / "astrid_mouse_cell_abundance.csv",
    OUTPUT_DIR / "abundance_variation_Marrow_Mukin.png",
    OUTPUT_DIR / "abundance_variation_Marrow_TMS.png",
    OUTPUT_DIR / "step1_scatter_Marrow_combined.png",
    OUTPUT_DIR / "M4_grid_Marrow_combined_results.csv",
    OUTPUT_DIR / "M23_grid_Marrow_combined_results.csv",
    OUTPUT_DIR / "M23_refined_exponential_results.csv",
    OUTPUT_DIR / "M23_refined_exponential_profiles.csv",
    *[
        OUTPUT_DIR / f"M23_refined_{figure_kind}_{celltype}.png"
        for celltype in REFINED_CELLTYPE_SUFFIXES
        for figure_kind in ("profiles", "fits", "parameter_grids", "valleys")
    ],
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kernel-name", default="python3")
    parser.add_argument(
        "--tidyi-input-root",
        type=Path,
        default=LEGACY_TIDYI_INPUT_ROOT,
        help="Directory containing the five required pooled-condition TIDYI CSVs.",
    )
    parser.add_argument("--tidyi-version-label", default="legacy-rate-files")
    parser.add_argument("--tidyi-model-label", default="legacy-model")
    parser.add_argument("--tidyi-formula-label", default="legacy-formula")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tidyi_input_root = args.tidyi_input_root.expanduser().resolve()
    rate_inputs = tuple(tidyi_input_root / name for name in RATE_FILENAMES)
    replacements = {
        2: (
            SourceReplacement(
                "/scratch/chenxi.sun/thesis_project/TIDYI/TIDYI_pattern/ty_pipeline_scripts",
                str(PIPELINE_HELPER.parent),
            ),
        ),
        4: (
            SourceReplacement(
                "/scratch/chenxi.sun/TIDYI_Astrid-anno/ty_df",
                str(tidyi_input_root),
            ),
        ),
        12: (
            SourceReplacement(
                "/scratch/chenxi.sun/RA project/system_biology/outputs/astrid_mouse_cell_abundance.csv",
                str(OUTPUT_DIR / "astrid_mouse_cell_abundance.csv"),
            ),
        ),
        13: (
            SourceReplacement(
                "/scratch/chenxi.sun/RA project/system_biology/outputs",
                str(OUTPUT_DIR),
                expected_count=2,
            ),
        ),
        **{
            cell: (
                SourceReplacement(
                    "/scratch/chenxi.sun/TIDYI_Astrid-anno/ty_df",
                    str(tidyi_input_root),
                ),
                SourceReplacement(
                    "/scratch/chenxi.sun/RA project/system_biology/outputs",
                    str(OUTPUT_DIR),
                ),
            )
            for cell in (21, 23, 26)
        },
    }
    run_notebook_master(
        master_name="Objective 6 - systems biology model fitting",
        source_notebook=SOURCE_NOTEBOOK,
        cell_numbers=SELECTED_CELLS,
        expected_outputs=EXPECTED_OUTPUTS,
        authoritative_inputs=(
            *rate_inputs,
            MUKIN_TISSUEWISE_H5AD,
            TMS_TISSUEWISE_H5AD,
            MUKIN_ASTRID_ROOT,
            MUKIN_ASTRID_ROOT / "BM_combined_ASTRID_results.csv",
            TMS_ASTRID_ROOT,
            PIPELINE_HELPER,
        ),
        manifest_path=OUTPUT_DIR / "master_systems_biology_model_fitting_run_manifest.json",
        kernel_name=args.kernel_name,
        source_replacements=replacements,
        run_metadata={
            "tidyi_input_root": str(tidyi_input_root),
            "tidyi_version_label": args.tidyi_version_label,
            "tidyi_model_label": args.tidyi_model_label,
            "tidyi_formula_label": args.tidyi_formula_label,
            "required_tidyi_filenames": list(RATE_FILENAMES),
            "rate_annotation_source": "ASTRID comparison workflow with celltype_hb corrections",
            "abundance_unit": "dataset x mouse x tissue x cell_type",
            "model_input_unit": "pooled condition x cell_type",
        },
    )


if __name__ == "__main__":
    main()

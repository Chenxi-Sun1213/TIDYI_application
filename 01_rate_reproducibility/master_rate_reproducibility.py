#!/usr/bin/env python3
"""Master entry point for Objective 1: TIDYI rate reproducibility.

Run with the legacy TIDYI installed in the selected Jupyter kernel:

    python master_rate_reproducibility.py

Run against a colleague's replacement TIDYI checkout:

    python master_rate_reproducibility.py \
        --tidyi-source-root /path/containing/the/tidyi/package \
        --tidyi-version-label <version-or-commit> \
        --tidyi-model-label <model-name> \
        --tidyi-formula-label <formula-name>

External inputs:
- thesis_project/tms_qcdata/AllRef_BM_0204
  Produced by: thesis TMS ASTRID workflow.
- thesis_project/mukin_qcdata/AllRef_BM_0204
  Produced by: thesis Mukin bone-marrow ASTRID workflow.
- thesis_project/mukin_qcdata/AllRef_skin_0204
  Produced by: thesis Mukin skin ASTRID workflow.
- BM_combined_ASTRID_results.csv in the Mukin bone-marrow directory
  Produced by: the same ASTRID workflow; used for RBC-contamination review.

Support analysis assembled here:
- rates_reproduce_test.ipynb, selected one-based cells listed in SELECTED_CELLS.

The source notebook historically wrote to locations inside the old RA project.
This handover master redirects every declared output to this objective's
``outputs`` directory without changing calculations.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


RA_ROOT = Path(__file__).resolve().parents[1]
if str(RA_ROOT) not in sys.path:
    sys.path.insert(0, str(RA_ROOT))

from handover_utils import (  # noqa: E402
    SourceReplacement,
    configure_tidyi_source,
    run_notebook_master,
    tidyi_runtime_capture_cell,
)


SOURCE_NOTEBOOK = Path(__file__).with_name("rates_reproduce_test.ipynb")
OUTPUT_DIR = Path(__file__).parent / "outputs"
TIDYI_RESULTS = OUTPUT_DIR / "tidyi_result_df.csv"
TIDYI_RUNTIME_PROVENANCE = OUTPUT_DIR / "tidyi_runtime_provenance.json"

# Explicitly resolves the notebook's former hidden execution order.
SELECTED_CELLS = (2, 4, 5, 8, 11, 12, 14, 17, 18, 19, 20, 21)

EXPECTED_OUTPUT_NAMES = (
    "mukin_bm_age_19m_sex_F_microbiota_SPF_tissue_BM_SPF_19m_2_vs_SPF_19m_1.png",
    "mukin_bm_age_2m_sex_F_microbiota_GF_tissue_BM_GF_2m_3_vs_GF_2m_2.png",
    "mukin_bm_age_2m_sex_F_microbiota_SPF_tissue_BM_SPF_2m_4_vs_SPF_2m_5.png",
    "mukin_sk_age_19m_sex_F_microbiota_GF_tissue_SK_GF_19m_1_vs_GF_19m_3.png",
    "mukin_sk_age_19m_sex_F_microbiota_SPF_tissue_SK_SPF_19m_2_vs_SPF_19m_3.png",
    "mukin_sk_age_2m_sex_F_microbiota_GF_tissue_SK_GF_2m_3_vs_GF_2m_1.png",
    "mukin_sk_age_2m_sex_F_microbiota_SPF_tissue_SK_SPF_2m_4_vs_SPF_2m_5.png",
    "tms_bm_age_18m_sex_M_microbiota_nan_tissue_BM_18-M-52_vs_18-M-53.png",
    "tms_bm_age_1m_sex_M_microbiota_nan_tissue_BM_1-M-62_vs_1-M-63.png",
    "tms_bm_age_21m_sex_F_microbiota_nan_tissue_BM_21-F-54_vs_21-F-55.png",
    "tms_bm_age_24m_sex_M_microbiota_nan_tissue_BM_24-M-58_vs_24-M-60.png",
    "tms_bm_age_30m_sex_M_microbiota_nan_tissue_BM_30-M-2_vs_30-M-5.png",
    "tms_bm_age_3m_sex_F_microbiota_nan_tissue_BM_3-F-57_vs_3-F-56.png",
    "pairwise_rate_comparisons.csv",
    "pairwise_rate_pearson_correlations.csv",
    "rate_reproducibility_connected_dot_plot.pdf",
    "rate_reproducibility_permutation_null_grid.pdf",
    "mukin_sk_keratinocyte_spf_gf_death_vs_prolif_scatter.pdf",
    "tidyi_runtime_provenance.json",
)
EXPECTED_OUTPUTS = (TIDYI_RESULTS, *[OUTPUT_DIR / name for name in EXPECTED_OUTPUT_NAMES])

TMS_ASTRID_DIR = Path("/scratch/chenxi.sun/thesis_project/tms_qcdata/AllRef_BM_0204")
MUKIN_BM_ASTRID_DIR = Path(
    "/scratch/chenxi.sun/thesis_project/mukin_qcdata/AllRef_BM_0204"
)
MUKIN_SKIN_ASTRID_DIR = Path(
    "/scratch/chenxi.sun/thesis_project/mukin_qcdata/AllRef_skin_0204"
)
MUKIN_BM_ASTRID_RESULTS = MUKIN_BM_ASTRID_DIR / "BM_combined_ASTRID_results.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kernel-name", default="python3")
    parser.add_argument(
        "--tidyi-source-root",
        type=Path,
        default=None,
        help="Optional directory containing the replacement tidyi package.",
    )
    parser.add_argument("--tidyi-version-label", default="legacy-environment")
    parser.add_argument("--tidyi-model-label", default="legacy-model")
    parser.add_argument("--tidyi-formula-label", default="legacy-formula")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_tidyi_source(args.tidyi_source_root)

    old_output_root = "/scratch/chenxi.sun/RA project/rates reproduce/rep_align_scatterplot copy"
    new_output_root = str(OUTPUT_DIR)
    old_tidyi_results = "/scratch/chenxi.sun/RA project/rates reproduce/tidyi_result_df.csv"
    replacements = {
        14: (SourceReplacement(old_tidyi_results, str(TIDYI_RESULTS)),),
        18: (SourceReplacement(old_tidyi_results, str(TIDYI_RESULTS)),),
        19: (SourceReplacement(old_output_root, new_output_root),),
        20: (SourceReplacement(old_output_root, new_output_root),),
        21: (
            SourceReplacement(old_output_root, new_output_root),
            SourceReplacement(old_tidyi_results, str(TIDYI_RESULTS)),
        ),
    }

    run_notebook_master(
        master_name="Objective 1 - TIDYI rate reproducibility",
        source_notebook=SOURCE_NOTEBOOK,
        cell_numbers=SELECTED_CELLS,
        expected_outputs=EXPECTED_OUTPUTS,
        authoritative_inputs=(
            TMS_ASTRID_DIR,
            MUKIN_BM_ASTRID_DIR,
            MUKIN_SKIN_ASTRID_DIR,
            MUKIN_BM_ASTRID_RESULTS,
        ),
        manifest_path=OUTPUT_DIR / "master_rate_reproducibility_run_manifest.json",
        kernel_name=args.kernel_name,
        source_replacements=replacements,
        append_cells=(tidyi_runtime_capture_cell(TIDYI_RUNTIME_PROVENANCE),),
        run_metadata={
            "tidyi_source_root": (
                str(args.tidyi_source_root.resolve()) if args.tidyi_source_root else None
            ),
            "tidyi_version_label": args.tidyi_version_label,
            "tidyi_model_label": args.tidyi_model_label,
            "tidyi_formula_label": args.tidyi_formula_label,
            "annotation_columns": {
                "tms_bm": "celltype_hb",
                "mukin_bm": "celltype_hb",
                "mukin_skin": "SingleR_CellType",
            },
            "minimum_cells_per_sample_celltype": 100,
        },
    )


if __name__ == "__main__":
    main()

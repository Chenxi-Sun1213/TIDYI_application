#!/usr/bin/env python3
"""Master entry point for Objective 4: Arigoni SciData 2024 validation.

Only the Arigoni SciData 2024 BE1run12 dataset is in scope. CCLE, Mix-seq, and
the historical zero-death diagnostics are intentionally excluded.

Legacy run:

    python master_arigoni_scidata_2024_validation.py

Replacement TIDYI run:

    python master_arigoni_scidata_2024_validation.py \
        --tidyi-source-root /path/containing/the/tidyi/package \
        --tidyi-version-label <version-or-commit> \
        --tidyi-model-label <model-name> \
        --tidyi-formula-label <formula-name>

External input:
- /srv/mfs/hausserlab/data/ARigoni_Scidata_2024/BE1run12
  Produced externally from the Arigoni Scientific Data 2024 release.

Support analysis assembled here:
- Arigoni_Scidata_2024.ipynb, selected one-based cells in SELECTED_CELLS.
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


SOURCE_NOTEBOOK = Path(__file__).with_name("Arigoni_Scidata_2024.ipynb")
DATA_DIR = Path("/srv/mfs/hausserlab/data/ARigoni_Scidata_2024/BE1run12")
OUTPUT_DIR = Path(__file__).parent / "outputs"
TIDYI_RUNTIME_PROVENANCE = OUTPUT_DIR / "tidyi_runtime_provenance.json"
MAIN_FIGURE = OUTPUT_DIR / "arigoni_tidyi_rates_by_cell_line.png"
SELECTED_CELLS = (2, 4, 7, 9, 11, 13, 15, 16, 17)

EXPECTED_OUTPUTS = (
    OUTPUT_DIR / "tidyi_cell_df.csv",
    OUTPUT_DIR / "tidyi_results_df.csv",
    MAIN_FIGURE,
    TIDYI_RUNTIME_PROVENANCE,
)


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
    save_main_figure = (
        "fig.savefig(OUTPUT_DIR / 'arigoni_tidyi_rates_by_cell_line.png', "
        "dpi=300, bbox_inches='tight')"
    )
    replacements = {
        4: (
            SourceReplacement(
                "/scratch/chenxi.sun/RA project/cancer pattern/output/Arigoni_Scidata_2024",
                str(OUTPUT_DIR),
            ),
        ),
    }
    run_notebook_master(
        master_name="Objective 4 - Arigoni SciData 2024 cancer cell-line validation",
        source_notebook=SOURCE_NOTEBOOK,
        cell_numbers=SELECTED_CELLS,
        expected_outputs=EXPECTED_OUTPUTS,
        authoritative_inputs=(DATA_DIR, DATA_DIR / "CRL5868" / "metrics_summary.csv"),
        manifest_path=OUTPUT_DIR / "master_arigoni_scidata_2024_run_manifest.json",
        kernel_name=args.kernel_name,
        source_replacements=replacements,
        append_cells=(
            save_main_figure,
            tidyi_runtime_capture_cell(TIDYI_RUNTIME_PROVENANCE),
        ),
        run_metadata={
            "tidyi_source_root": (
                str(args.tidyi_source_root.resolve()) if args.tidyi_source_root else None
            ),
            "tidyi_version_label": args.tidyi_version_label,
            "tidyi_model_label": args.tidyi_model_label,
            "tidyi_formula_label": args.tidyi_formula_label,
            "dataset_scope": "Arigoni SciData 2024 BE1run12 only",
            "grouping_column": "cell_line",
            "minimum_cells_per_tidyi_group": 100,
        },
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Master entry point for Objective 2: mouse-human skin rate comparison.

Legacy run:

    python master_mouse_human_skin.py

Replacement TIDYI run:

    python master_mouse_human_skin.py \
        --tidyi-source-root /path/containing/the/tidyi/package \
        --tidyi-version-label <version-or-commit> \
        --tidyi-model-label <model-name> \
        --tidyi-formula-label <formula-name>

External inputs:
- thesis_project/mukin_qcdata/AllRef_skin_0204
  Produced by: thesis Mukin skin ASTRID workflow.
- classifier_data/mart_export_gene_names.txt
  Produced externally with the installed classifier-data package.

The source annotation is ``SingleR_CellType``. The additional
``kc_human_validation`` label is produced here from Krt10/Krt15 expression and
must not be described as a manually reviewed ASTRID cluster mapping.
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


SOURCE_NOTEBOOK = Path(__file__).with_name("human_skin_validation_inmice.ipynb")
OUTPUT_DIR = Path(__file__).parent / "outputs"
TIDYI_RUNTIME_PROVENANCE = OUTPUT_DIR / "tidyi_runtime_provenance.json"
SELECTED_CELLS = (2, 4, 6, 7, 8, 9, 11, 13, 14, 15)

EXPECTED_OUTPUT_NAMES = (
    "mukin_skin_umap_annotation_leiden.png",
    "mukin_skin_umap_krt15_krt10.png",
    "tidyi_rates_kc_human_validation.csv",
    "tidyi_rates_kc_human_validation_scatter_GF.png",
    "tidyi_rates_kc_human_validation_scatter_SPF.png",
    "tidyi_runtime_provenance.json",
)
EXPECTED_OUTPUTS = tuple(OUTPUT_DIR / name for name in EXPECTED_OUTPUT_NAMES)

MUKIN_SKIN_ASTRID_DIR = Path(
    "/scratch/chenxi.sun/thesis_project/mukin_qcdata/AllRef_skin_0204"
)
CLASSIFIER_GENE_MAP = RA_ROOT / "inputs/classifier_data/mart_export_gene_names.txt"


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
    replacements = {
        8: (
            SourceReplacement(
                "/scratch/chenxi.sun/RA project/rates reproduce/skin_validation_output",
                str(OUTPUT_DIR),
            ),
        ),
        13: (
            SourceReplacement(
                "/scratch/common/miniconda3/envs/scx_py310/lib/python3.10/"
                "site-packages/classifier_data/mart_export_gene_names.txt",
                str(CLASSIFIER_GENE_MAP),
            ),
        ),
    }
    run_notebook_master(
        master_name="Objective 2 - mouse-human skin rate comparison",
        source_notebook=SOURCE_NOTEBOOK,
        cell_numbers=SELECTED_CELLS,
        expected_outputs=EXPECTED_OUTPUTS,
        authoritative_inputs=(MUKIN_SKIN_ASTRID_DIR, CLASSIFIER_GENE_MAP),
        manifest_path=OUTPUT_DIR / "master_mouse_human_skin_run_manifest.json",
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
            "source_annotation_column": "SingleR_CellType",
            "validation_annotation_column": "kc_human_validation",
            "keratinocyte_rule": "Krt10 > Krt15 - 1.2 => differentiated; otherwise progenitor",
        },
    )


if __name__ == "__main__":
    main()

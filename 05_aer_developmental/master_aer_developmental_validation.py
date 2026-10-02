#!/usr/bin/env python3
"""Master entry point for Objective 5: AER developmental validation.

This handover master assembles the active downstream AER analysis from the two
verified H5AD inputs stored under the shared ``AER2025_Juliane`` data directory.
The audited notebook and helpers are retained beside this master.

Legacy run:

    python master_aer_developmental_validation.py

Replacement TIDYI run:

    python master_aer_developmental_validation.py \
        --tidyi-source-root /path/containing/the/tidyi/package \
        --tidyi-version-label <version-or-commit> \
        --tidyi-model-label <model-name> \
        --tidyi-formula-label <formula-name>

External inputs and producers:
- /srv/mfs/hausserlab/data/AER2025_Juliane/aer_six_samples_raw_counts.h5ad
  Produced by: AER2025_Juliane/scanpy_analysis/run_scanpy_reconstruction.py.
- /srv/mfs/hausserlab/data/AER2025_Juliane/aer_post_entropy.h5ad
  Produced by: AER_DISTILL_analysis.ipynb after QC, ASTRID manual mapping,
  per-sample kNN transfer, and entropy filtering.
- TIDYI_statistics typical/reference CV tables
  Produced by: the TIDYI statistics workflow; used for Method 3.
- inputs/distill_celltype_percentages_by_sample.tsv
  Produced by: the AER final-export/composition stage; used by the publication
  figure script.

Formal outputs are limited to local output sections 03, 05, all eight retained
files in 07 (four primary plus four reference-sensitivity files), and
pub_ready_fig. Upstream QC and annotation checkpoints are not regenerated.
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


PROJECT_ROOT = Path(__file__).parent
AER_DATA_ROOT = Path("/srv/mfs/hausserlab/data/AER2025_Juliane")
SOURCE_NOTEBOOK = PROJECT_ROOT / "AER_DISTILL_analysis.ipynb"
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
INPUT_ROOT = PROJECT_ROOT / "inputs"
RAW_RECONSTRUCTION = AER_DATA_ROOT / "aer_six_samples_raw_counts.h5ad"
POST_ENTROPY_H5AD = AER_DATA_ROOT / "aer_post_entropy.h5ad"
COMPOSITION_TABLE = INPUT_ROOT / "distill_celltype_percentages_by_sample.tsv"
TYPICAL_CV = RA_ROOT / "inputs/tidyi_statistics/stats_notebook_aging_typical_cv_values.csv"
REFERENCE_CV = (
    RA_ROOT
    / "inputs/tidyi_statistics/stats_notebook_ref_cv_separate_by_dataset_celltype.csv"
)
SENSITIVITY_SCRIPT = Path(__file__).with_name(
    "aer_method3_reference_sensitivity.py"
)
TIDYI_RUNTIME_PROVENANCE = OUTPUT_ROOT / "05_tidyi/tidyi_runtime_provenance.json"

# The retained post-entropy checkpoint makes upstream QC/ASTRID/kNN cells
# unnecessary here. Cell 19 loads that checkpoint when RUN_KNN remains false.
SELECTED_CELLS = (2, 4, 19, 21, 22, 24, 26, 27, 33)

EXPECTED_OUTPUTS = (
    OUTPUT_ROOT / "03_biolucid/celltype_by_sample_input_counts.csv",
    OUTPUT_ROOT / "03_biolucid/global_results.csv",
    OUTPUT_ROOT / "03_biolucid/per_sample_results.csv",
    OUTPUT_ROOT / "03_biolucid/biolucid_scatter.png",
    OUTPUT_ROOT / "05_tidyi/gene_mapping_report.csv",
    OUTPUT_ROOT / "05_tidyi/tidyi_genes_used.txt",
    OUTPUT_ROOT / "05_tidyi/sample_celltype_group_counts.csv",
    OUTPUT_ROOT / "05_tidyi/sample_celltype_rates.csv",
    OUTPUT_ROOT / "05_tidyi/tidyi_rates.png",
    TIDYI_RUNTIME_PROVENANCE,
    OUTPUT_ROOT
    / "07_growth_rate_decomposition/aer_method3_aer_global_mean_regulation_magnitude_filter.csv",
    OUTPUT_ROOT
    / "07_growth_rate_decomposition/aer_method3_aer_global_mean_unfiltered_heatmap.png",
    OUTPUT_ROOT
    / "07_growth_rate_decomposition/aer_method3_aer_global_mean_magnitude_filtered_heatmap.png",
    OUTPUT_ROOT
    / "07_growth_rate_decomposition/aer_method3_aer_global_mean_magnitude_filter_histogram_5bins.png",
    OUTPUT_ROOT
    / "07_growth_rate_decomposition/aer_method3_regulation_magnitude_filter.csv",
    OUTPUT_ROOT / "07_growth_rate_decomposition/aer_method3_unfiltered_heatmap.png",
    OUTPUT_ROOT
    / "07_growth_rate_decomposition/aer_method3_magnitude_filtered_heatmap.png",
    OUTPUT_ROOT
    / "07_growth_rate_decomposition/aer_method3_magnitude_filter_histogram_5bins.png",
    OUTPUT_ROOT / "pub_ready_fig/aer_tidyi_rates_pub_ready.pdf",
    OUTPUT_ROOT / "pub_ready_fig/aer_tidyi_rates_pub_ready.svg",
    OUTPUT_ROOT / "pub_ready_fig/celltype_composition_pub_ready.pdf",
    OUTPUT_ROOT / "pub_ready_fig/celltype_composition_pub_ready.svg",
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
    parser.add_argument(
        "--postprocess-python",
        default=sys.executable,
        help="Python executable used for make_pub_ready_figures.py.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_tidyi_source(args.tidyi_source_root)
    replacements = {
        2: (
            SourceReplacement(
                'PROJECT_ROOT = Path("/scratch/chenxi.sun/AER_DISTILL")',
                f'PROJECT_ROOT = Path({str(PROJECT_ROOT)!r})',
            ),
        ),
        4: (
            SourceReplacement(
                'INPUT_H5AD = Path(\n    "/scratch/chenxi.sun/AER2025_Juliane/scanpy_analysis/outputs/"\n    "h5ad/aer_six_samples_raw_counts.h5ad"\n)',
                f'INPUT_H5AD = Path({str(RAW_RECONSTRUCTION)!r})',
            ),
            SourceReplacement(
                'ASTRID_MARKER_FORMULA_SOURCE = Path(\n    "/scratch/chenxi.sun/thesis_project/ASTRID_pkg/astrid/data/"\n    "ExpectedCellTypesMarkersMouse.csv"\n)',
                'ASTRID_MARKER_FORMULA_SOURCE = PROJECT_ROOT / "inputs" / "ExpectedCellTypesMarkersMouse.csv"',
            ),
            SourceReplacement(
                'for required_path in [\n    INPUT_H5AD,\n    ASTRID_SCRIPT,\n    ASTRID_REFERENCE,\n    ASTRID_MARKER_FORMULA,\n]:',
                'for required_path in [INPUT_H5AD, ASTRID_MARKER_FORMULA]:',
            ),
        ),
        19: (
            SourceReplacement(
                'POST_ENTROPY_H5AD = OUT["restore"] / "aer_post_entropy.h5ad"',
                f'POST_ENTROPY_H5AD = Path({str(POST_ENTROPY_H5AD)!r})',
            ),
            SourceReplacement(
                'transfer_summary = pd.read_csv(OUT["restore"] / "knn_transfer_summary.csv")',
                'transfer_summary = pd.read_csv(PROJECT_ROOT / "inputs" / "knn_transfer_summary.csv")',
            ),
            SourceReplacement(
                'entropy_qc_by_sample = pd.read_csv(OUT["restore"] / "entropy_qc_by_sample.csv")',
                'entropy_qc_by_sample = pd.read_csv(PROJECT_ROOT / "inputs" / "entropy_qc_by_sample.csv")',
            ),
            SourceReplacement(
                'entropy_qc_by_celltype = pd.read_csv(\n        OUT["restore"] / "entropy_qc_by_sample_celltype.csv"\n    )',
                'entropy_qc_by_celltype = pd.read_csv(\n        PROJECT_ROOT / "inputs" / "entropy_qc_by_sample_celltype.csv"\n    )',
            ),
        ),
        22: (SourceReplacement("RUN_BIOLUCID = False", "RUN_BIOLUCID = True"),),
        24: (
            SourceReplacement(
                "RUN_TIDYI_PREPARATION = False", "RUN_TIDYI_PREPARATION = True"
            ),
        ),
        26: (
            SourceReplacement("RUN_TIDYI_TRAINING = False", "RUN_TIDYI_TRAINING = True"),
        ),
        27: (
            SourceReplacement(
                "RUN_TIDYI_APPLICATION = False", "RUN_TIDYI_APPLICATION = True"
            ),
        ),
        33: (
            SourceReplacement(
                'METHOD3_TYPICAL_CV_PATH = Path(\n    "/scratch/chenxi.sun/TIDYI_statistics/results/aging_paired_ttest/"\n    "effect_size_typical_cv/stats_notebook_aging_typical_cv_values.csv"\n)',
                f'METHOD3_TYPICAL_CV_PATH = Path({str(TYPICAL_CV)!r})',
            ),
            SourceReplacement(
                'METHOD3_REFERENCE_CV_DETAIL_PATH = Path(\n    "/scratch/chenxi.sun/TIDYI_statistics/results/ref_cv/"\n    "stats_notebook_ref_cv_separate_by_dataset_celltype.csv"\n)',
                f'METHOD3_REFERENCE_CV_DETAIL_PATH = Path({str(REFERENCE_CV)!r})',
            ),
        ),
    }
    run_notebook_master(
        master_name="Objective 5 - AER developmental validation",
        source_notebook=SOURCE_NOTEBOOK,
        cell_numbers=SELECTED_CELLS,
        expected_outputs=EXPECTED_OUTPUTS,
        authoritative_inputs=(
            RAW_RECONSTRUCTION,
            POST_ENTROPY_H5AD,
            COMPOSITION_TABLE,
            TYPICAL_CV,
            REFERENCE_CV,
            PROJECT_ROOT / "aer_distill_helpers.py",
            PROJECT_ROOT / "make_pub_ready_figures.py",
            SENSITIVITY_SCRIPT,
        ),
        manifest_path=OUTPUT_ROOT / "master_aer_developmental_validation_run_manifest.json",
        kernel_name=args.kernel_name,
        source_replacements=replacements,
        append_cells=(tidyi_runtime_capture_cell(TIDYI_RUNTIME_PROVENANCE),),
        post_commands=(
            (
                args.postprocess_python,
                str(SENSITIVITY_SCRIPT),
                "--primary-table",
                str(
                    OUTPUT_ROOT
                    / "07_growth_rate_decomposition/aer_method3_aer_global_mean_regulation_magnitude_filter.csv"
                ),
                "--output-dir",
                str(OUTPUT_ROOT / "07_growth_rate_decomposition"),
            ),
            (args.postprocess_python, str(PROJECT_ROOT / "make_pub_ready_figures.py")),
        ),
        run_metadata={
            "tidyi_source_root": (
                str(args.tidyi_source_root.resolve()) if args.tidyi_source_root else None
            ),
            "tidyi_version_label": args.tidyi_version_label,
            "tidyi_model_label": args.tidyi_model_label,
            "tidyi_formula_label": args.tidyi_formula_label,
            "post_entropy_annotation_column": "cell_type",
            "manual_astrid_annotation_column": "astrid_reannotation",
            "entropy_threshold": 0.3,
            "tidyi_minimum_cells_exclusive": 50,
            "formal_output_sections": [
                "03_biolucid",
                "05_tidyi",
                "07_primary_and_reference_sensitivity",
                "pub_ready_fig",
            ],
        },
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Master entry point for Objective 3: death-proliferation dominance.

This master contains only the minimum computation needed to reproduce the six
confirmed regulation-magnitude-filter outputs. Scatter plots, effect-size plots,
growth-noise filtering, and the older decomposition figures are out of scope.

Legacy run:

    python master_death_proliferation_dominance.py

Run with rate CSVs produced by a colleague's replacement TIDYI:

    python master_death_proliferation_dominance.py \
        --tidyi-input-root /path/to/new/ty_df \
        --tidyi-version-label <version-or-commit> \
        --tidyi-model-label <model-name> \
        --tidyi-formula-label <formula-name>

Legacy rate producer and annotation lineage:
- Producer: thesis_project/TIDYI/TIDYI_pattern/ty_pipeline_scripts/
  tidyi_comparison_analysis.py.
- Run evidence: TIDYI_refine/logs/all_comparisons_refined.log.
- Base annotation: tissue_wise_knn in the Mukin/TMS tissuewise H5AD files.
- Mukin additionally used comparison-specific constrained kNN relabeling before
  TIDYI; these rates are not a manually reviewed ASTRID-label product.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import PercentFormatter  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402


RA_ROOT = Path(__file__).resolve().parents[1]
if str(RA_ROOT) not in sys.path:
    sys.path.insert(0, str(RA_ROOT))

from handover_utils import (  # noqa: E402
    TIDYI_RATE_COLUMNS,
    file_record,
    require_paths,
    snapshot_outputs,
    validate_refreshed_outputs,
)


LEGACY_TIDYI_INPUT_ROOT = RA_ROOT / "inputs/tidyi_refine"
TYPICAL_CV_PATH = (
    RA_ROOT / "inputs/tidyi_statistics/stats_notebook_aging_typical_cv_values.csv"
)
REFERENCE_CV_DETAIL_PATH = (
    RA_ROOT
    / "inputs/tidyi_statistics/stats_notebook_ref_cv_separate_by_dataset_celltype.csv"
)
OUTPUT_DIR = Path(__file__).parent / "outputs"
REGULATION_MAGNITUDE_THRESHOLD = 2 * np.log(2) / (24 * 3)

RATE_FILENAMES = (
    "mukin_young_spf_vs_old_spf_tidyi.csv",
    "mukin_young_spf_vs_young_gf_tidyi.csv",
    "mukin_old_spf_vs_old_gf_tidyi.csv",
    "tms_young_male_vs_young_female_tidyi.csv",
    "tms_old_male_vs_old_female_tidyi.csv",
)

ATLAS_TISSUES = {
    "Mukin": ["Marrow", "Epi-e", "Spleen"],
    "TMS": ["Marrow", "Spleen"],
}
ATLAS_DATASET = {"Mukin": "mukin", "TMS": "tms"}
CONTRASTS = (
    {
        "atlas": "Mukin",
        "contrast": "Young GF vs Young SPF",
        "condition_a": "young_gf",
        "condition_b": "young_spf",
        "source_a": "mukin_young_spf_vs_young_gf",
        "source_b": "mukin_young_spf_vs_young_gf",
    },
    {
        "atlas": "Mukin",
        "contrast": "Old SPF vs Young SPF",
        "condition_a": "old_spf",
        "condition_b": "young_spf",
        "source_a": "mukin_young_spf_vs_old_spf",
        "source_b": "mukin_young_spf_vs_old_spf",
    },
    {
        "atlas": "Mukin",
        "contrast": "Old GF vs Young GF",
        "condition_a": "old_gf",
        "condition_b": "young_gf",
        "source_a": "mukin_old_spf_vs_old_gf",
        "source_b": "mukin_young_spf_vs_young_gf",
    },
    {
        "atlas": "Mukin",
        "contrast": "Old GF vs Old SPF",
        "condition_a": "old_gf",
        "condition_b": "old_spf",
        "source_a": "mukin_old_spf_vs_old_gf",
        "source_b": "mukin_old_spf_vs_old_gf",
    },
    {
        "atlas": "TMS",
        "contrast": "Young female vs Young male",
        "condition_a": "young_female",
        "condition_b": "young_male",
        "source_a": "tms_young_male_vs_young_female",
        "source_b": "tms_young_male_vs_young_female",
    },
    {
        "atlas": "TMS",
        "contrast": "Old male vs Young male",
        "condition_a": "old_male",
        "condition_b": "young_male",
        "source_a": "tms_old_male_vs_old_female",
        "source_b": "tms_young_male_vs_young_female",
    },
    {
        "atlas": "TMS",
        "contrast": "Old female vs Young female",
        "condition_a": "old_female",
        "condition_b": "young_female",
        "source_a": "tms_old_male_vs_old_female",
        "source_b": "tms_young_male_vs_young_female",
    },
    {
        "atlas": "TMS",
        "contrast": "Old female vs Old male",
        "condition_a": "old_female",
        "condition_b": "old_male",
        "source_a": "tms_old_male_vs_old_female",
        "source_b": "tms_old_male_vs_old_female",
    },
)

EXPECTED_OUTPUTS = tuple(
    OUTPUT_DIR / name
    for name in (
        "mukin_regulation_magnitude_filter.csv",
        "mukin_method3_magnitude_filtered_heatmap.png",
        "mukin_method3_magnitude_filter_histograms.png",
        "tms_regulation_magnitude_filter.csv",
        "tms_method3_magnitude_filtered_heatmap.png",
        "tms_method3_magnitude_filter_histograms.png",
    )
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tidyi-input-root",
        type=Path,
        default=LEGACY_TIDYI_INPUT_ROOT,
        help="Directory containing the five required TIDYI rate CSVs.",
    )
    parser.add_argument("--tidyi-version-label", default="legacy-rate-files")
    parser.add_argument("--tidyi-model-label", default="legacy-model")
    parser.add_argument("--tidyi-formula-label", default="legacy-formula")
    return parser.parse_args()


def load_rate_tables(input_root: Path) -> tuple[dict[str, pd.DataFrame], tuple[Path, ...]]:
    input_paths = tuple(input_root / filename for filename in RATE_FILENAMES)
    require_paths((*input_paths, TYPICAL_CV_PATH, REFERENCE_CV_DETAIL_PATH))
    tables: dict[str, pd.DataFrame] = {}
    for path in input_paths:
        table = pd.read_csv(path)
        missing = set(TIDYI_RATE_COLUMNS).difference(table.columns)
        if missing:
            raise KeyError(f"{path} is missing TIDYI columns: {sorted(missing)}")
        tables[path.stem.removesuffix("_tidyi")] = table
    return tables, input_paths


def group_table(
    tables: dict[str, pd.DataFrame], source_name: str, atlas: str, group: str
) -> pd.DataFrame:
    table = tables[source_name]
    return table[
        table["analysis_dataset"].eq(ATLAS_DATASET[atlas])
        & table["comparison_group"].eq(group)
        & table["tissue"].isin(ATLAS_TISSUES[atlas])
    ].copy()


def build_paired_rates(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for spec in CONTRASTS:
        left = group_table(tables, spec["source_a"], spec["atlas"], spec["condition_a"])
        right = group_table(tables, spec["source_b"], spec["atlas"], spec["condition_b"])
        keep = ["tissue", "celltype", "prolif_rate", "death_rate"]
        paired = left[keep].merge(
            right[keep],
            on=["tissue", "celltype"],
            how="inner",
            suffixes=("_a", "_b"),
            validate="one_to_one",
        )
        paired.insert(0, "atlas", spec["atlas"])
        paired.insert(1, "contrast", spec["contrast"])
        paired.insert(2, "condition_a", spec["condition_a"])
        paired.insert(3, "condition_b", spec["condition_b"])
        rows.append(paired)
    paired_rates = pd.concat(rows, ignore_index=True)
    key = ["atlas", "tissue", "celltype", "contrast"]
    if paired_rates.duplicated(key).any():
        raise ValueError(f"Unexpected duplicate paired-rate keys: {key}")
    return paired_rates


def calculate_method3(paired_rates: pd.DataFrame) -> pd.DataFrame:
    typical_cv_table = pd.read_csv(TYPICAL_CV_PATH)
    reference_cv_detail = pd.read_csv(REFERENCE_CV_DETAIL_PATH)
    typical_cv = {
        (row["dataset"], row["rate"]): float(row["typical_cv"])
        for _, row in typical_cv_table.iterrows()
    }
    cv_source = {
        (row["dataset"], row["rate"]): (row["cv_source_dataset"], row["cv_rate"])
        for _, row in typical_cv_table.iterrows()
    }

    global_mean_rate: dict[tuple[str, str], float] = {}
    for atlas in ("Mukin", "TMS"):
        for rate in ("prolif_rate", "death_rate"):
            source_dataset, source_rate = cv_source[(atlas, rate)]
            source_rows = reference_cv_detail[
                reference_cv_detail["dataset"].eq(source_dataset)
                & reference_cv_detail["rate"].eq(source_rate)
                & reference_cv_detail["cv_status"].eq("computed")
            ]
            if source_rows.empty:
                raise ValueError(
                    f"No computed reference CV rows for {atlas} {rate}: "
                    f"{source_dataset} {source_rate}"
                )
            global_mean_rate[(atlas, rate)] = float(source_rows["mean_rate"].mean())

    result = paired_rates.copy()
    result["delta_prolif"] = result["prolif_rate_a"] - result["prolif_rate_b"]
    result["delta_death"] = result["death_rate_a"] - result["death_rate_b"]
    result["death_component"] = -result["delta_death"]
    result["delta_growth"] = result["delta_prolif"] + result["death_component"]
    result["row_label"] = result["tissue"] + " | " + result["celltype"]

    def method3_weight(atlas: str) -> float:
        numerator = typical_cv[(atlas, "death_rate")] * global_mean_rate[
            (atlas, "death_rate")
        ]
        denominator = typical_cv[(atlas, "prolif_rate")] * global_mean_rate[
            (atlas, "prolif_rate")
        ]
        if denominator <= 0:
            raise ValueError(f"Non-positive Method 3 denominator for {atlas}")
        return numerator / denominator

    result["k_m3"] = result["atlas"].map(
        {atlas: method3_weight(atlas) for atlas in ("Mukin", "TMS")}
    )
    denominator = result["death_component"].abs() + result["k_m3"] * result[
        "delta_prolif"
    ].abs()
    result["death_fraction_m3"] = np.where(
        denominator > 0,
        result["death_component"].abs() / denominator,
        np.nan,
    )
    if not np.allclose(
        result["delta_growth"],
        result["delta_prolif"] - result["delta_death"],
        equal_nan=True,
    ):
        raise AssertionError("delta growth invariant failed")
    available = result["death_fraction_m3"].dropna()
    if not available.between(0, 1).all():
        raise AssertionError("Method 3 death fractions must be between zero and one")
    return result


def magnitude_table(method3: pd.DataFrame, atlas: str) -> pd.DataFrame:
    result = method3[method3["atlas"].eq(atlas)].copy()
    result["abs_delta_proliferation"] = result["delta_prolif"].abs()
    result["abs_delta_death"] = result["delta_death"].abs()
    result["regulation_magnitude"] = (
        result["abs_delta_proliferation"] + result["abs_delta_death"]
    )
    result["regulation_threshold"] = REGULATION_MAGNITUDE_THRESHOLD
    result["passes_magnitude_filter"] = (
        result["regulation_magnitude"] >= REGULATION_MAGNITUDE_THRESHOLD
    )
    result["death_fraction_m3_magnitude_filtered"] = result[
        "death_fraction_m3"
    ].where(result["passes_magnitude_filter"])
    return result


def save_magnitude_outputs(table: pd.DataFrame, atlas: str) -> None:
    prefix = atlas.lower()
    output_columns = [
        "atlas",
        "tissue",
        "celltype",
        "contrast",
        "condition_a",
        "condition_b",
        "delta_prolif",
        "delta_death",
        "abs_delta_proliferation",
        "abs_delta_death",
        "regulation_magnitude",
        "regulation_threshold",
        "passes_magnitude_filter",
        "death_fraction_m3",
        "death_fraction_m3_magnitude_filtered",
    ]
    table[output_columns].to_csv(
        OUTPUT_DIR / f"{prefix}_regulation_magnitude_filter.csv", index=False
    )

    tissue_rank = {value: index for index, value in enumerate(ATLAS_TISSUES[atlas])}
    row_order = (
        table[["tissue", "celltype", "row_label"]]
        .drop_duplicates()
        .assign(tissue_rank=lambda frame: frame["tissue"].map(tissue_rank))
        .sort_values(["tissue_rank", "celltype"])["row_label"]
        .tolist()
    )
    column_order = [
        spec["contrast"] for spec in CONTRASTS if spec["atlas"] == atlas
    ]
    original = table.pivot_table(
        index="row_label",
        columns="contrast",
        values="death_fraction_m3",
        aggfunc="first",
    ).reindex(index=row_order, columns=column_order)
    filtered = table.pivot_table(
        index="row_label",
        columns="contrast",
        values="death_fraction_m3_magnitude_filtered",
        aggfunc="first",
    ).reindex(index=row_order, columns=column_order)
    availability = original.notna().astype(float)
    annotations = filtered.applymap(lambda value: "" if pd.isna(value) else f"{value:.2f}")

    figure_height = max(5.2, 0.22 * len(filtered) + 1.5)
    figure_width = max(6.2, 1.25 * len(column_order) + 3.0)
    figure, axis = plt.subplots(
        figsize=(figure_width, figure_height), constrained_layout=True
    )
    sns.heatmap(
        availability,
        ax=axis,
        cmap=ListedColormap(["#D9D9D9", "#FFFFFF"]),
        vmin=0,
        vmax=1,
        cbar=False,
        linewidths=0.35,
        linecolor="white",
    )
    sns.heatmap(
        filtered,
        mask=filtered.isna(),
        ax=axis,
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
        linewidths=0,
        annot=annotations,
        fmt="",
        annot_kws={"fontsize": 5.8},
        cbar_kws={"label": "Death contribution fraction"},
    )
    axis.set_title(
        "Method 3 death contribution after regulation-magnitude filtering\n"
        f"{atlas}; colored cells satisfy |Δβ| + |Δα| ≥ "
        f"{REGULATION_MAGNITUDE_THRESHOLD:.4f} h⁻¹",
        fontsize=10,
        pad=28,
    )
    axis.legend(
        handles=[
            Patch(
                facecolor="white",
                edgecolor="#BDBDBD",
                label=(
                    f"Filtered: |Δβ| + |Δα| < "
                    f"{REGULATION_MAGNITUDE_THRESHOLD:.4f} h⁻¹"
                ),
            ),
            Patch(
                facecolor="#D9D9D9",
                edgecolor="#BDBDBD",
                label="Fraction unavailable",
            ),
        ],
        loc="lower left",
        bbox_to_anchor=(0, 1.0),
        ncol=2,
        frameon=False,
        fontsize=7,
    )
    axis.set_xlabel("Contrast: condition A vs condition B", fontsize=8)
    axis.set_ylabel("Tissue | celltype", fontsize=8)
    axis.tick_params(axis="x", labelrotation=35, labelsize=8)
    axis.tick_params(axis="y", labelsize=6)
    figure.savefig(
        OUTPUT_DIR / f"{prefix}_method3_magnitude_filtered_heatmap.png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close(figure)

    all_fractions = table["death_fraction_m3"].dropna()
    retained_fractions = table["death_fraction_m3_magnitude_filtered"].dropna()
    histogram_rows = (
        ("Before filter", all_fractions, "#9ECAE1"),
        ("After filter", retained_fractions, "#4C78A8"),
    )
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(7.2, 5.4),
        sharex=True,
        sharey="col",
        constrained_layout=True,
    )
    for row_index, (label, values, color) in enumerate(histogram_rows):
        for column_index, number_of_bins in enumerate((5, 10)):
            axis = axes[row_index, column_index]
            bins = np.linspace(0, 1, number_of_bins + 1)
            axis.hist(
                values,
                bins=bins,
                weights=np.ones(len(values)) / len(values),
                color=color,
                edgecolor="white",
                linewidth=0.5,
            )
            axis.set_title(
                f"{label}; {number_of_bins} bins (n={len(values)})", fontsize=9
            )
            axis.set_xlim(0, 1)
            axis.set_xticks(np.linspace(0, 1, 6))
            axis.tick_params(axis="both", labelsize=7)
            axis.yaxis.set_major_formatter(PercentFormatter(1.0))
            axis.spines[["top", "right"]].set_visible(False)
            if row_index == 1:
                sns.rugplot(
                    x=values, ax=axis, height=0.05, color="black", linewidth=0.7
                )
    for axis in axes[:, 0]:
        axis.set_ylabel("Proportion of entries", fontsize=8)
    for axis in axes[1, :]:
        axis.set_xlabel("Death contribution fraction (Method 3)", fontsize=8)
    figure.suptitle(
        f"{atlas} Method 3 fraction distribution before and after\n"
        "regulation-magnitude filtering",
        fontsize=10,
    )
    figure.savefig(
        OUTPUT_DIR / f"{prefix}_method3_magnitude_filter_histograms.png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close(figure)


def main() -> None:
    args = parse_args()
    input_root = args.tidyi_input_root.expanduser().resolve()
    tables, rate_paths = load_rate_tables(input_root)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pre_run_output_records = [
        file_record(path) for path in EXPECTED_OUTPUTS if path.exists()
    ]
    before = snapshot_outputs(EXPECTED_OUTPUTS)
    started = datetime.now(timezone.utc)

    paired_rates = build_paired_rates(tables)
    method3 = calculate_method3(paired_rates)
    for atlas in ("Mukin", "TMS"):
        save_magnitude_outputs(magnitude_table(method3, atlas), atlas)

    output_records = validate_refreshed_outputs(EXPECTED_OUTPUTS, before)
    manifest = {
        "master_name": "Objective 3 - death-proliferation dominance",
        "status": "completed",
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "tidyi_input_root": str(input_root),
        "tidyi_version_label": args.tidyi_version_label,
        "tidyi_model_label": args.tidyi_model_label,
        "tidyi_formula_label": args.tidyi_formula_label,
        "rate_inputs": [file_record(path) for path in rate_paths],
        "typical_cv_input": file_record(TYPICAL_CV_PATH),
        "reference_cv_input": file_record(REFERENCE_CV_DETAIL_PATH),
        "legacy_rate_producer": (
            "/scratch/chenxi.sun/thesis_project/TIDYI/TIDYI_pattern/"
            "ty_pipeline_scripts/tidyi_comparison_analysis.py"
        ),
        "legacy_rate_log": "/scratch/chenxi.sun/TIDYI_refine/logs/all_comparisons_refined.log",
        "legacy_annotation_lineage": (
            "tissue_wise_knn; Mukin comparison-specific constrained kNN relabeling"
        ),
        "regulation_magnitude_threshold": REGULATION_MAGNITUDE_THRESHOLD,
        "pre_run_outputs": pre_run_output_records,
        "outputs": output_records,
    }
    (OUTPUT_DIR / "master_death_proliferation_dominance_run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()

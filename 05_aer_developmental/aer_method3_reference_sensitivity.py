#!/usr/bin/env python3
"""Regenerate the retained AER Method 3 Mukin-reference sensitivity outputs.

The current AER helper writes the primary analysis with AER global mean rates
and preserves the Mukin-reference Method 3 weight in
``method3_k_mukin_reference_sensitivity``. This support script converts that
audited primary table into the four historical sensitivity files retained in
output section 07. It is called only by the AER master script.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import PercentFormatter  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402


DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "outputs/07_growth_rate_decomposition"
DEFAULT_PRIMARY_TABLE = (
    DEFAULT_OUTPUT_DIR
    / "aer_method3_aer_global_mean_regulation_magnitude_filter.csv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-table", type=Path, default=DEFAULT_PRIMARY_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def build_sensitivity_table(primary: pd.DataFrame) -> pd.DataFrame:
    required = {
        "abs_delta_death",
        "abs_delta_proliferation",
        "passes_magnitude_filter",
        "method3_k_mukin_reference_sensitivity",
    }
    missing = required.difference(primary.columns)
    if missing:
        raise KeyError(f"Primary Method 3 table is missing columns: {sorted(missing)}")

    sensitivity = primary.copy()
    sensitivity["method3_k"] = sensitivity[
        "method3_k_mukin_reference_sensitivity"
    ]
    denominator = sensitivity["abs_delta_death"] + sensitivity[
        "method3_k"
    ] * sensitivity["abs_delta_proliferation"]
    sensitivity["death_fraction_m3"] = np.where(
        denominator > 0,
        sensitivity["abs_delta_death"] / denominator,
        np.nan,
    )
    sensitivity["death_fraction_m3_magnitude_filtered"] = sensitivity[
        "death_fraction_m3"
    ].where(sensitivity["passes_magnitude_filter"])
    sensitivity["available_method3_fraction"] = sensitivity[
        "death_fraction_m3"
    ].notna()
    sensitivity = sensitivity.drop(
        columns=[
            "method3_k_mukin_reference_sensitivity",
            "global_mean_rate_source",
            "n_aer_rate_rows_for_global_mean",
            "proliferation_global_mean_used",
            "death_global_mean_used",
        ]
    )
    return sensitivity


def plot_sensitivity(results: pd.DataFrame, output_dir: Path) -> None:
    cell_types = results["cell_type"].drop_duplicates().tolist()
    contrast_order = (
        results[["contrast_order", "contrast"]]
        .drop_duplicates()
        .sort_values("contrast_order")["contrast"]
        .tolist()
    )
    original = results.pivot_table(
        index="cell_type",
        columns="contrast",
        values="death_fraction_m3",
        aggfunc="first",
    ).reindex(index=cell_types, columns=contrast_order)
    filtered = results.pivot_table(
        index="cell_type",
        columns="contrast",
        values="death_fraction_m3_magnitude_filtered",
        aggfunc="first",
    ).reindex(index=cell_types, columns=contrast_order)
    availability = original.notna().astype(float)
    figure_height = max(5.2, 0.22 * len(cell_types) + 1.5)
    figure_width = max(6.2, 1.25 * len(contrast_order) + 3.0)
    background_cmap = ListedColormap(["#D9D9D9", "#FFFFFF"])

    def heatmap(values: pd.DataFrame, filtered_plot: bool) -> plt.Figure:
        figure, axis = plt.subplots(
            figsize=(figure_width, figure_height), constrained_layout=True
        )
        sns.heatmap(
            availability,
            ax=axis,
            cmap=background_cmap,
            vmin=0,
            vmax=1,
            cbar=False,
            linewidths=0.35,
            linecolor="white",
        )
        annotations = values.applymap(
            lambda value: "" if pd.isna(value) else f"{value:.2f}"
        )
        sns.heatmap(
            values,
            mask=values.isna(),
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
        threshold = float(results["regulation_threshold"].iloc[0])
        if filtered_plot:
            title = (
                "AER Method 3 death contribution after regulation-magnitude filtering\n"
                f"Colored cells satisfy |Δβ| + |Δα| ≥ {threshold:.4f} h⁻¹"
            )
            handles = [
                Patch(
                    facecolor="white",
                    edgecolor="#BDBDBD",
                    label=f"Filtered: |Δβ| + |Δα| < {threshold:.4f} h⁻¹",
                ),
                Patch(
                    facecolor="#D9D9D9",
                    edgecolor="#BDBDBD",
                    label="Fraction unavailable",
                ),
            ]
        else:
            title = (
                "AER Method 3 death contribution to growth regulation\n"
                "Before regulation-magnitude filtering"
            )
            handles = [
                Patch(
                    facecolor="#D9D9D9",
                    edgecolor="#BDBDBD",
                    label="Fraction unavailable",
                )
            ]
        axis.set_title(title, fontsize=10, pad=28)
        axis.legend(
            handles=handles,
            loc="lower left",
            bbox_to_anchor=(0, 1.0),
            ncol=2 if filtered_plot else 1,
            frameon=False,
            fontsize=7,
        )
        axis.set_xlabel("Contrast: condition A minus condition B", fontsize=8)
        axis.set_ylabel("Cell type", fontsize=8)
        axis.tick_params(axis="x", labelrotation=35, labelsize=8)
        axis.tick_params(axis="y", labelsize=6)
        return figure

    unfiltered_figure = heatmap(original, filtered_plot=False)
    unfiltered_figure.savefig(
        output_dir / "aer_method3_unfiltered_heatmap.png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close(unfiltered_figure)

    filtered_figure = heatmap(filtered, filtered_plot=True)
    filtered_figure.savefig(
        output_dir / "aer_method3_magnitude_filtered_heatmap.png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close(filtered_figure)

    all_fractions = results["death_fraction_m3"].dropna()
    retained_fractions = results["death_fraction_m3_magnitude_filtered"].dropna()
    figure, axes = plt.subplots(
        2, 1, figsize=(3.6, 5.4), sharex=True, sharey=True, constrained_layout=True
    )
    for index, (label, values, color) in enumerate(
        (
            ("Before filter", all_fractions, "#9ECAE1"),
            ("After filter", retained_fractions, "#4C78A8"),
        )
    ):
        axis = axes[index]
        bins = np.linspace(0, 1, 6)
        axis.hist(
            values,
            bins=bins,
            weights=np.ones(len(values)) / len(values),
            color=color,
            edgecolor="white",
            linewidth=0.5,
        )
        axis.set_title(f"{label}; 5 bins (n={len(values)})", fontsize=9)
        axis.set_xlim(0, 1)
        axis.set_xticks(np.linspace(0, 1, 6))
        axis.tick_params(axis="both", labelsize=7)
        axis.yaxis.set_major_formatter(PercentFormatter(1.0))
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_ylabel("Proportion of entries", fontsize=8)
        if index == 1:
            sns.rugplot(x=values, ax=axis, height=0.05, color="black", linewidth=0.7)
    axes[1].set_xlabel("Death contribution fraction (Method 3)", fontsize=8)
    figure.suptitle(
        "AER Method 3 fraction distribution before and after\n"
        "regulation-magnitude filtering",
        fontsize=10,
    )
    figure.savefig(
        output_dir / "aer_method3_magnitude_filter_histogram_5bins.png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if not args.primary_table.exists():
        raise FileNotFoundError(args.primary_table)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sensitivity = build_sensitivity_table(pd.read_csv(args.primary_table))
    sensitivity.to_csv(
        args.output_dir / "aer_method3_regulation_magnitude_filter.csv", index=False
    )
    plot_sensitivity(sensitivity, args.output_dir)


if __name__ == "__main__":
    main()

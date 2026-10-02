from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import LogLocator, MaxNLocator, NullFormatter


ROOT = Path(__file__).resolve().parent
RATES_PATH = ROOT / "outputs/05_tidyi/sample_celltype_rates.csv"
COMPOSITION_PATH = (
    ROOT
    / "inputs/distill_celltype_percentages_by_sample.tsv"
)
OUTPUT_DIR = ROOT / "outputs/pub_ready_fig"

STAGES = ["E9.5", "E10.5", "E11.5"]
GENOTYPES = ["WT", "Dac"]
STAGE_COLORS = {
    "E9.5": "#66C2A5",
    "E10.5": "#FC8D62",
    "E11.5": "#8DA0CB",
}
SAMPLE_COLORS = {
    ("E9.5", "WT"): "#C6DBEF",
    ("E9.5", "Dac"): "#FCBBA1",
    ("E10.5", "WT"): "#6BAED6",
    ("E10.5", "Dac"): "#FB6A4A",
    ("E11.5", "WT"): "#2171B5",
    ("E11.5", "Dac"): "#CB181D",
}
GENOTYPE_MARKERS = {"WT": "o", "Dac": "^"}


def padded_limits(values: pd.Series, fraction: float = 0.09) -> tuple[float, float]:
    lower = float(values.min())
    upper = float(values.max())
    padding = max((upper - lower) * fraction, 1e-5)
    return lower - padding, upper + padding


def save_vector_figure(figure: plt.Figure, stem: str, title: str) -> None:
    common = {"bbox_inches": "tight", "facecolor": "white"}
    figure.savefig(
        OUTPUT_DIR / f"{stem}.pdf",
        metadata={
            "Title": title,
            "Creator": "Matplotlib",
            "Subject": "AER DISTILL publication-ready figure",
        },
        **common,
    )
    figure.savefig(
        OUTPUT_DIR / f"{stem}.svg",
        metadata={
            "Title": title,
            "Creator": "Matplotlib",
            "Description": "AER DISTILL publication-ready figure",
        },
        **common,
    )


def plot_aer_rates(rates: pd.DataFrame) -> plt.Figure:
    required = {
        "stage",
        "genotype",
        "cell_type",
        "n_cells",
        "death_rate",
        "proliferation_rate",
    }
    missing = required.difference(rates.columns)
    if missing:
        raise KeyError(f"Rates table is missing columns: {sorted(missing)}")

    aer = rates.loc[rates["cell_type"].astype(str).eq("AER")].copy()
    aer["stage"] = pd.Categorical(aer["stage"], STAGES, ordered=True)
    aer["genotype"] = pd.Categorical(aer["genotype"], GENOTYPES, ordered=True)
    aer = aer.sort_values(["genotype", "stage"])
    if aer.empty or aer[["stage", "genotype"]].isna().any().any():
        raise AssertionError("AER rates do not match the expected stages/genotypes")
    if aer.duplicated(["stage", "genotype"]).any():
        raise AssertionError("Expected at most one AER rate per stage/genotype")

    fig, ax = plt.subplots(figsize=(4.35, 3.55))
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.15, top=0.88)
    xlim = padded_limits(aer["death_rate"])
    ylim = padded_limits(aer["proliferation_rate"])

    equality_start = max(xlim[0], ylim[0])
    equality_end = min(xlim[1], ylim[1])
    if equality_start < equality_end:
        ax.plot(
            [equality_start, equality_end],
            [equality_start, equality_end],
            color="#B7B7B7",
            linestyle=(0, (3, 2)),
            linewidth=0.8,
            zorder=0,
        )

    for genotype in GENOTYPES:
        trajectory = aer.loc[aer["genotype"].eq(genotype)].sort_values("stage")
        points = trajectory[["death_rate", "proliferation_rate"]].to_numpy()
        for start, end in zip(points[:-1], points[1:]):
            ax.add_patch(
                FancyArrowPatch(
                    start,
                    end,
                    arrowstyle="-|>",
                    mutation_scale=7.5,
                    shrinkA=5.5,
                    shrinkB=7.0,
                    linewidth=0.9,
                    color="#565656",
                    alpha=0.85,
                    zorder=1,
                )
            )

    for _, row in aer.iterrows():
        ax.scatter(
            row["death_rate"],
            row["proliferation_rate"],
            s=68,
            marker=GENOTYPE_MARKERS[str(row["genotype"])],
            facecolor=STAGE_COLORS[str(row["stage"])],
            edgecolor="black",
            linewidth=0.65,
            zorder=3,
        )

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.set_xlabel("Death rate")
    ax.set_ylabel("Proliferation rate")
    ax.set_title("AER cell-state dynamics", fontweight="bold", pad=8)
    ax.grid(color="#D9D9D9", linewidth=0.55, alpha=0.75)
    ax.set_axisbelow(True)

    stage_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=STAGE_COLORS[stage],
            markeredgecolor="black",
            markeredgewidth=0.5,
            markersize=5.7,
            label=stage,
        )
        for stage in STAGES
    ]
    genotype_handles = [
        Line2D(
            [0],
            [0],
            marker=GENOTYPE_MARKERS[genotype],
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=0.65,
            markersize=5.7,
            label=genotype,
        )
        for genotype in GENOTYPES
    ]
    stage_legend = ax.legend(
        handles=stage_handles,
        title="Time point",
        loc="upper right",
        frameon=False,
        borderaxespad=0.35,
        handletextpad=0.45,
        labelspacing=0.35,
    )
    ax.add_artist(stage_legend)
    ax.legend(
        handles=genotype_handles,
        title="Genotype",
        loc="upper right",
        bbox_to_anchor=(1, 0.72),
        frameon=False,
        borderaxespad=0.35,
        handletextpad=0.45,
        labelspacing=0.35,
    )
    return fig


def plot_celltype_composition(table: pd.DataFrame) -> plt.Figure:
    required = {
        "sample",
        "cell_type",
        "n_cells",
        "sample_total_cells",
        "percent_cells",
        "stage",
        "genotype",
    }
    missing = required.difference(table.columns)
    if missing:
        raise KeyError(f"Composition table is missing columns: {sorted(missing)}")
    sums = table.groupby("sample", observed=True)["percent_cells"].sum()
    if not np.allclose(sums.to_numpy(), 100.0):
        raise AssertionError("Cell-type percentages do not sum to 100 per sample")

    cell_types = sorted(table["cell_type"].astype(str).unique())
    sample_meta = table[["sample", "stage", "genotype"]].drop_duplicates().copy()
    sample_meta["stage"] = pd.Categorical(sample_meta["stage"], STAGES, ordered=True)
    sample_meta["genotype"] = pd.Categorical(
        sample_meta["genotype"], GENOTYPES, ordered=True
    )
    sample_meta = sample_meta.sort_values(["stage", "genotype"]).reset_index(drop=True)
    if len(sample_meta) != 6 or sample_meta[["stage", "genotype"]].isna().any().any():
        raise AssertionError("Expected six samples covering three stages x two genotypes")

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.16, top=0.92)
    x = np.arange(len(cell_types), dtype=float)
    group_width = 0.82
    bar_width = group_width / len(sample_meta)
    for sample_index, row in sample_meta.iterrows():
        sample = str(row["sample"])
        stage = str(row["stage"])
        genotype = str(row["genotype"])
        values = (
            table.loc[table["sample"].astype(str).eq(sample)]
            .set_index("cell_type")
            .reindex(cell_types)["percent_cells"]
            .to_numpy()
        )
        offset = (sample_index - (len(sample_meta) - 1) / 2) * bar_width
        ax.bar(
            x + offset,
            np.where(values > 0, values, np.nan),
            width=bar_width * 0.93,
            facecolor=SAMPLE_COLORS[(stage, genotype)],
            edgecolor="#333333",
            linewidth=0.4,
            zorder=2,
        )

    positive = table.loc[table["percent_cells"] > 0, "percent_cells"]
    ax.set_yscale("log", base=10)
    ax.set_ylim(float(positive.min()) * 0.72, float(positive.max()) * 1.35)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=4))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())
    display_labels = {
        "DV": "Dorso-ventral\nectoderm",
        "Vascular cells": "Vascular\ncells",
    }
    ax.set_xticks(x)
    ax.set_xticklabels([display_labels.get(cell_type, cell_type) for cell_type in cell_types])
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Cells within sample (%)\n(log scale)")
    ax.set_title("Cell type composition", fontweight="bold", pad=8)
    ax.grid(axis="y", which="major", color="#CFCFCF", linewidth=0.6, alpha=0.8)
    ax.grid(axis="y", which="minor", color="#E5E5E5", linewidth=0.4, alpha=0.75)
    ax.set_axisbelow(True)

    time_handles = [
        (
            Patch(facecolor=SAMPLE_COLORS[(stage, "WT")], edgecolor="#333333"),
            Patch(facecolor=SAMPLE_COLORS[(stage, "Dac")], edgecolor="#333333"),
        )
        for stage in STAGES
    ]
    time_legend = ax.legend(
        handles=time_handles,
        labels=STAGES,
        title="Time point",
        ncols=1,
        loc="upper right",
        bbox_to_anchor=(0.995, 0.995),
        frameon=False,
        borderaxespad=0.25,
        handlelength=1.3,
        labelspacing=0.3,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)},
    )
    ax.add_artist(time_legend)
    ax.legend(
        handles=[
            Patch(facecolor=SAMPLE_COLORS[("E10.5", "WT")], edgecolor="#333333", label="WT"),
            Patch(facecolor=SAMPLE_COLORS[("E10.5", "Dac")], edgecolor="#333333", label="Dac"),
        ],
        title="Genotype",
        ncols=1,
        loc="upper right",
        bbox_to_anchor=(0.995, 0.69),
        frameon=False,
        borderaxespad=0.25,
        handlelength=1.3,
        labelspacing=0.3,
    )
    return fig


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rates = pd.read_csv(RATES_PATH)
    composition = pd.read_csv(COMPOSITION_PATH, sep="\t")

    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "legend.title_fontsize": 7.5,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    ):
        aer_figure = plot_aer_rates(rates)
        save_vector_figure(
            aer_figure,
            "aer_tidyi_rates_pub_ready",
            "AER cell-state dynamics",
        )
        plt.close(aer_figure)

        composition_figure = plot_celltype_composition(composition)
        save_vector_figure(
            composition_figure,
            "celltype_composition_pub_ready",
            "Cell type composition",
        )
        plt.close(composition_figure)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate Objective 3 reference CVs independently of historical notebooks."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd


RA_ROOT = Path(__file__).resolve().parents[2]
if str(RA_ROOT) not in sys.path:
    sys.path.insert(0, str(RA_ROOT))

from gene_mapping import GENE_MAPPING_FILE, map_mouse_symbols_to_ensembl
from handover_utils import file_record, require_paths


LOGGER = logging.getLogger("reference_cv_generation")
MINIMUM_REFERENCE_CELLS = 100
TYPICAL_CV_FILENAME = "stats_notebook_aging_typical_cv_values.csv"
DETAIL_FILENAME = "stats_notebook_ref_cv_separate_by_dataset_celltype.csv"
OUTPUT_FILENAMES = (
    "reference_sample_rates.csv",
    DETAIL_FILENAME,
    TYPICAL_CV_FILENAME,
    "reference_cv_exclusions.csv",
    "reference_cv_run_manifest.json",
)

# TIDYI_analysis.ipynb cells 13-15 define the historical rate labels.
MUKIN_RATE_LABEL_MAP = {
    "Cd8-effector Tcell": "Cd8 Tcell",
    "Cd8-naive Tcell": "Cd8 Tcell",
    "Cd4-effector Tcell": "Cd4 Tcell",
    "Cd4-naive Tcell": "Cd4 Tcell",
    "Th17-effector Tcell": "Cd4 Tcell",
    "Natural Killer": "ILCs",
    "Myeloid progenitor": "Mye-Neu Progenitor",
    "Neutrophil Progenitor": "Mye-Neu Progenitor",
    "Cd14 Monocyte": "Monocyte",
    "Cd16 Monocyte": "Monocyte",
    "Memory B cell": "B late stage",
    "Naive B cell": "B late stage",
    "Proliferative B cell": "B late stage",
    "Pre-B cell": "B early stage",
    "Immature B cell": "B early stage",
    "HF Bulge cell": "Hair Follicle",
    "HF Quiescennt cell": "Hair Follicle",
    "HF cell": "Hair Follicle",
    "Proliferating KC": "Basal KC",
}
# stats.ipynb cell 7 additionally maps this tissuewise counting label.
MUKIN_COUNT_LABEL_MAP = {
    **MUKIN_RATE_LABEL_MAP,
    "Resident Cd8-effector Tcell": "Cd8 Tcell",
}


REFERENCE_GROUPS = {"Mukin": "Mukin 2m SPF", "TMS": "TMS 1m male"}
RATE_NAMES = {"prolif_rate": "proliferation_rate", "death_rate": "death_rate"}
REFERENCE_KEYS = ["dataset", "sample_unit", "cell_type"]


def build_reference_tables(
    sample_rates: pd.DataFrame, cell_counts: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Apply the historical sample-unit CV and median-CV summary rules.

    Rates are inferred per sample, never from pooled condition rates. As in
    stats.ipynb cell 7, summaries group by reference dataset and cell type,
    without adding tissue to the grouping key. Duplicate sample keys are an
    error rather than silently selecting an arbitrary row.
    """
    for frame, required in (
        (sample_rates, REFERENCE_KEYS + ["tissue", *RATE_NAMES.values()]),
        (cell_counts, REFERENCE_KEYS + ["n_cells"]),
    ):
        missing = set(required).difference(frame.columns)
        if missing:
            raise ValueError(f"Missing reference columns: {sorted(missing)}")
        if frame[REFERENCE_KEYS].isna().any().any():
            raise ValueError("Reference dataset/sample/cell-type keys cannot be missing")
        if frame.duplicated(REFERENCE_KEYS).any():
            raise ValueError("Duplicate reference dataset/sample/cell-type keys")

    for rate in RATE_NAMES.values():
        values = pd.to_numeric(sample_rates[rate], errors="raise")
        if np.isinf(values).any():
            raise ValueError(f"Infinite reference rate: {rate}")
    if not (cell_counts["n_cells"] > 0).all():
        raise ValueError("Reference cell counts must be positive")

    merged = sample_rates.merge(
        cell_counts, on=REFERENCE_KEYS, how="left", validate="one_to_one", indicator=True
    )
    dropped = merged.loc[merged["_merge"].eq("left_only"), sample_rates.columns].copy()
    dropped["exclusion_reason"] = "no_matching_tissue_wise_knn_cell_count"
    matched = merged.loc[merged["_merge"].eq("both")].copy()

    rows = []
    for (dataset, cell_type), group in matched.groupby(["dataset", "cell_type"], sort=True):
        for rate in RATE_NAMES.values():
            values = pd.to_numeric(group[rate], errors="raise").dropna()
            mean_rate = values.mean()
            sd_rate = values.std(ddof=1)
            computed = len(values) >= 2 and pd.notna(mean_rate) and abs(mean_rate) > 1e-9
            rows.append(
                {
                    "dataset": dataset,
                    "cell_type": cell_type,
                    "rate": rate,
                    "n_replicates": len(values),
                    "mean_rate": mean_rate,
                    "sd_rate": sd_rate,
                    "cv": sd_rate / abs(mean_rate) if computed else np.nan,
                    "mean_cell_number": group["n_cells"].mean(),
                    "rate_times_cell_number": mean_rate * group["n_cells"].mean(),
                    "n_cell_count_replicates": group["sample_unit"].nunique(),
                    "cv_status": "computed" if computed else "not_computed",
                    "filter_note": (
                        "CV requires at least 2 sample units and absolute mean > 1e-9; "
                        "unmatched tissue_wise_knn count keys are excluded"
                    ),
                    "tissues": ";".join(sorted(group["tissue"].astype(str).unique())),
                }
            )
    if not rows:
        raise ValueError("No reference rates matched the cell-count keys")
    detail = pd.DataFrame(rows)

    typical_rows = []
    for atlas, reference in REFERENCE_GROUPS.items():
        for output_rate, detail_rate in RATE_NAMES.items():
            valid = detail.loc[
                detail["dataset"].eq(reference)
                & detail["rate"].eq(detail_rate)
                & detail["cv_status"].eq("computed")
            ]
            if valid.empty:
                raise ValueError(f"No computable reference CVs for {reference} / {detail_rate}")
            typical_rows.append(
                {
                    "dataset": atlas,
                    "cv_source_dataset": reference,
                    "rate": output_rate,
                    "cv_rate": detail_rate,
                    "typical_cv": valid["cv"].median(),
                    "n_cell_types": len(valid),
                    "reference_mean_rate": valid["mean_rate"].mean(),
                }
            )
    return detail, pd.DataFrame(typical_rows), dropped


def reference_samples(stats: pd.DataFrame, atlas: str) -> pd.DataFrame:
    """Select the original reference cohorts from the saved QC sample tables."""
    sample_col = "Sample" if atlas == "Mukin" else "sample"
    metadata_col = "model" if atlas == "Mukin" else "mouse.id"
    table = stats[[sample_col, "tissue", "age", metadata_col]].copy()
    if table[sample_col].isna().any() or table.duplicated(sample_col).any():
        raise ValueError(f"Missing or duplicate sample IDs in {atlas} sample metadata")
    if atlas == "Mukin":
        mask = table["age"].eq("2m") & table["model"].eq("SPF")
    else:
        mask = table["age"].eq("1m") & table["mouse.id"].str.split("-").str[1].eq("M")
    selected = table.loc[mask].rename(columns={sample_col: "sample_unit"})
    if selected.empty:
        raise ValueError(f"No reference samples selected for {atlas}")
    LOGGER.info("%s reference sample selection: %s -> %s", atlas, len(table), len(selected))
    return selected


def prepare_reference_cells(adata, samples: pd.DataFrame, atlas: str):
    """Select sample-level groups using historical rate annotations and >=100 cells.

    The saved historical rate-group sets exactly match this cutoff on the
    original H5ADs. The currently saved notebook's default of 10 does not
    reproduce those sets, and the pooled-comparison cutoff of 80 is unrelated.
    """
    sample_col = "Sample" if atlas == "Mukin" else "sample"
    obs = adata.obs[[sample_col, "live-dead_annotation"]].copy()
    mask = obs[sample_col].isin(samples["sample_unit"]) & obs["live-dead_annotation"].notna()
    obs = obs.loc[mask].rename(columns={sample_col: "sample_unit"})
    mapping = MUKIN_RATE_LABEL_MAP if atlas == "Mukin" else {"proerythroblast": "erythroblast"}
    obs["cell_type"] = obs["live-dead_annotation"].astype(str).replace(mapping)
    counts = obs.groupby(["sample_unit", "cell_type"], observed=True).size().rename("n_cells").reset_index()
    excluded = counts.loc[counts["n_cells"] < MINIMUM_REFERENCE_CELLS].copy()
    excluded.insert(0, "dataset", REFERENCE_GROUPS[atlas])
    excluded["exclusion_reason"] = "fewer_than_100_cells_in_rate_annotation"
    eligible = counts.loc[counts["n_cells"] >= MINIMUM_REFERENCE_CELLS]
    keep = pd.MultiIndex.from_frame(obs[["sample_unit", "cell_type"]]).isin(
        pd.MultiIndex.from_frame(eligible[["sample_unit", "cell_type"]])
    )
    obs = obs.loc[keep]
    if obs.empty:
        raise ValueError(f"No {atlas} reference groups meet the 100-cell cutoff")
    selected = adata[obs.index].to_memory()
    selected.obs["sample_unit"] = obs["sample_unit"].to_numpy()
    selected.obs["cell_type"] = obs["cell_type"].to_numpy()
    LOGGER.info(
        "%s reference cells: input=%s retained=%s; groups before=%s after=%s excluded=%s",
        atlas, adata.n_obs, selected.n_obs, len(counts), len(eligible), len(excluded),
    )
    return selected, excluded


def reference_cell_counts(obs: pd.DataFrame, samples: pd.DataFrame, atlas: str) -> pd.DataFrame:
    """Reproduce stats.ipynb's tissuewise count matching, without age repair."""
    sample_col = "Sample" if atlas == "Mukin" else "sample"
    mask = obs[sample_col].isin(samples["sample_unit"]) & obs["tissue_wise_knn"].notna()
    if atlas == "Mukin":
        mask &= obs["age"].eq("2m") & obs["model"].eq("SPF")
    else:
        mask &= obs["age"].eq("1m") & obs["mouse.id"].astype(str).str.split("-").str[1].eq("M")
    selected = obs.loc[mask, [sample_col, "tissue_wise_knn"]].copy()
    selected.columns = ["sample_unit", "cell_type"]
    selected["cell_type"] = selected["cell_type"].astype(str)
    if atlas == "Mukin":
        selected["cell_type"] = selected["cell_type"].replace(MUKIN_COUNT_LABEL_MAP)
    result = selected.groupby(["sample_unit", "cell_type"], observed=True).size().rename("n_cells").reset_index()
    result.insert(0, "dataset", REFERENCE_GROUPS[atlas])
    LOGGER.info("%s count matching: input cells=%s selected=%s count keys=%s", atlas, len(obs), len(selected), len(result))
    return result


def infer_sample_rates(adata, samples: pd.DataFrame, atlas: str, tidyi_module) -> pd.DataFrame:
    """Run the selected TIDYI implementation on each original sample-cell-type group."""
    if "raw" not in adata.layers:
        raise KeyError(f"Missing raw count layer for {atlas} reference inference")
    map_mouse_symbols_to_ensembl(adata, GENE_MAPPING_FILE)
    adata.layers["counts"] = adata.layers["raw"].copy()
    np.random.seed(0)
    model = tidyi_module.tidy_rates_mouse(gene_list=adata.var_names.tolist())
    model.train_classifier()
    tissue = samples.set_index("sample_unit")["tissue"]
    rows = []
    for (sample, cell_type), positions in adata.obs.groupby(
        ["sample_unit", "cell_type"], observed=True, sort=True
    ).indices.items():
        proliferation, death, _, _ = model.apply_classifier(adata[positions].copy())
        rows.append(
            {
                "dataset": REFERENCE_GROUPS[atlas],
                "sample_unit": sample,
                "tissue": tissue.loc[sample],
                "cell_type": cell_type,
                "n_rate_cells": len(positions),
                "proliferation_rate": float(proliferation),
                "death_rate": float(death),
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for atlas in ("mukin", "tms"):
        parser.add_argument(f"--{atlas}-rate-h5ad", type=Path, required=True)
        parser.add_argument(f"--{atlas}-counts-h5ad", type=Path, required=True)
        parser.add_argument(f"--{atlas}-stats-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tidyi-version-label", default="legacy-environment")
    parser.add_argument("--tidyi-model-label", default="legacy-model")
    parser.add_argument("--tidyi-formula-label", default="legacy-formula")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".matplotlib"))
    os.environ.setdefault("NUMBA_CACHE_DIR", str(args.output_dir / ".numba_cache"))
    import anndata as ad
    from tidyi import tidyi_main_mouse

    started = datetime.now(timezone.utc)
    inputs = [GENE_MAPPING_FILE]
    for atlas in ("mukin", "tms"):
        inputs.extend(getattr(args, f"{atlas}_{suffix}") for suffix in ("rate_h5ad", "counts_h5ad", "stats_csv"))
    require_paths(inputs)
    sample_rates, cell_counts, exclusions = [], [], []
    for atlas in REFERENCE_GROUPS:
        prefix = atlas.lower()
        samples = reference_samples(pd.read_csv(getattr(args, f"{prefix}_stats_csv")), atlas)
        original = ad.read_h5ad(getattr(args, f"{prefix}_rate_h5ad"), backed="r")
        try:
            selected, excluded = prepare_reference_cells(original, samples, atlas)
        finally:
            original.file.close()
        sample_rates.append(infer_sample_rates(selected, samples, atlas, tidyi_main_mouse))
        exclusions.append(excluded)
        del selected
        tissuewise = ad.read_h5ad(getattr(args, f"{prefix}_counts_h5ad"), backed="r")
        try:
            cell_counts.append(reference_cell_counts(tissuewise.obs, samples, atlas))
        finally:
            tissuewise.file.close()

    rates = pd.concat(sample_rates, ignore_index=True)
    rates.to_csv(args.output_dir / OUTPUT_FILENAMES[0], index=False)
    detail, typical, unmatched = build_reference_tables(rates, pd.concat(cell_counts, ignore_index=True))
    exclusions.append(unmatched)
    pd.concat(exclusions, ignore_index=True).to_csv(args.output_dir / OUTPUT_FILENAMES[3], index=False)
    detail.to_csv(args.output_dir / DETAIL_FILENAME, index=False)
    typical["cv_summary_source"] = str((args.output_dir / DETAIL_FILENAME).resolve())
    typical.to_csv(args.output_dir / TYPICAL_CV_FILENAME, index=False)
    try:
        version = importlib.metadata.version("tidyi")
    except importlib.metadata.PackageNotFoundError:
        version = "unavailable"
    manifest = {
        "status": "completed",
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "python_executable": sys.executable,
        "tidyi_module": file_record(Path(tidyi_main_mouse.__file__)),
        "tidyi_distribution_version": version,
        "tidyi_version_label": args.tidyi_version_label,
        "tidyi_model_label": args.tidyi_model_label,
        "tidyi_formula_label": args.tidyi_formula_label,
        "inputs": [file_record(path) for path in inputs],
        "configuration": {
            "reference_groups": REFERENCE_GROUPS,
            "rate_annotation": "live-dead_annotation",
            "count_annotation": "tissue_wise_knn",
            "mukin_rate_label_map": MUKIN_RATE_LABEL_MAP,
            "mukin_count_label_map": MUKIN_COUNT_LABEL_MAP,
            "tms_rate_label_map": {"proerythroblast": "erythroblast"},
            "minimum_rate_cells_per_sample_celltype": MINIMUM_REFERENCE_CELLS,
            "cutoff_evidence": "100-cell cutoff reconstructs all 68 Mukin and 455 TMS saved rate groups on original H5ADs",
            "sample_unit": "Sample/sample, not automatically an independent animal",
            "cv_grouping": ["dataset", "cell_type"],
            "cv_ddof": 1,
            "minimum_replicates": 2,
            "minimum_absolute_mean_exclusive": 1e-9,
            "typical_cv": "median across computable cell types",
            "reference_mean": "unweighted mean of computable cell-type mean rates",
            "raw_layer": "raw",
            "numpy_random_seed": 0,
            "age_repair_for_counts": False,
        },
        "n_sample_rate_rows": len(rates),
        "n_unmatched_rate_rows": len(unmatched),
        "n_uncomputable_cv_rows": int(detail["cv_status"].ne("computed").sum()),
        "outputs": [file_record(args.output_dir / name) for name in OUTPUT_FILENAMES[:-1]],
    }
    (args.output_dir / OUTPUT_FILENAMES[-1]).write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Reference CV generation complete: %s", args.output_dir)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run comparison-level TIDYI analysis on ASTRID-annotated inputs."""

from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

RA_ROOT = Path(__file__).resolve().parents[2]
if str(RA_ROOT) not in sys.path:
    sys.path.insert(0, str(RA_ROOT))
from gene_mapping import GENE_MAPPING_FILE, load_mouse_gene_mapping, map_mouse_symbols_to_ensembl

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import anndata as ad
from anndata import ImplicitModificationWarning
import numpy as np
import pandas as pd
import scanpy as sc

from tidyi_comparison_analysis import (
    ANALYSIS_DATASET_COL,
    ANALYSIS_TISSUE_COL,
    ANALYSIS_TISSUE_DELIMITER,
    COMPARISON_GROUP_COL,
    COMPARISON_NAME_COL,
    FINAL_TIDY_COLUMNS,
    PERCENTAGE_COL,
    expand_analysis_sample_column,
    SP_CT_COL,
    require_columns,
)
from tidyi_plot_results import (
    add_derived_metrics,
    ensure_output_dir,
    plot_cross_tissue_heatmaps,
    plot_growth_state_histplot,
    plot_growth_turnover_scatter,
    plot_kinetics_scatter,
)


LOGGER = logging.getLogger("tidyi_astrid_comparison_analysis")

PIPELINE_DIR = Path(__file__).resolve().parent
THESIS_PROJECT_ROOT = PIPELINE_DIR.parents[2]
WORKSPACE_ROOT = THESIS_PROJECT_ROOT.parent
DEFAULT_OUTPUT_ROOT = WORKSPACE_ROOT / "TIDYI_Astrid-anno"

DEFAULT_MUKIN_BM_ROOT = THESIS_PROJECT_ROOT / "mukin_qcdata" / "AllRef_BM_0204"
DEFAULT_MUKIN_SKIN_ROOT = THESIS_PROJECT_ROOT / "mukin_qcdata" / "AllRef_skin_0204"
DEFAULT_TMS_BM_ROOT = THESIS_PROJECT_ROOT / "tms_qcdata" / "AllRef_BM_0204"

MUKIN_RBC_GENES = ("Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Hbb-bh1", "Hbb-bh2")
TMS_REMOVE_CELLTYPES = {
    "fibroblast",
    "endothelial",
    "keratinocyte proliferating",
    "smooth muscle",
    "macrophage alveolar",
    "alveolar epithelial",
    "hair follicle",
}
TMS_AGE_GROUP_MAP = {
    "1m": "young",
    "3m": "young",
    "18m": "mid",
    "21m": "old",
    "24m": "old",
    "30m": "veryold",
}
TMS_ERYTHROCYTE_RULES = {
    "1-M-62": ("1_0", "1_1", "1_2", "1_3"),
    "1-M-63": ("3",),
    "3-F-56": ("0_1", "0_6"),
    "3-F-57": ("2",),
    "18-M-52": ("2",),
    "18-M-53": ("2_2", "2_4"),
    "21-F-54": ("0_9",),
    "21-F-55": ("1_2", "1_4"),
    "24-M-58": ("2", "3"),
    "24-M-59": ("1_1",),
    "24-M-60": ("1", "3"),
    "24-M-61": ("1", "3"),
    "30-M-3": ("2", "3"),
    "30-M-2": ("1", "4"),
    "30-M-5": ("2", "3"),
}
TMS_ERYTHROBLAST_PREFIXES = ("1_2", "1_3", "1_4", "0_6", "2_0", "2_2")
OUTPUT_LAYOUT = {
    "plot_root_dir_name": "plots",
    "ty_plot_dir_name": "ty_comparison",
    "ty_df_dir_name": "ty_df",
    "tidyi_dir_name": "tidyi",
    "tidyi_csv_suffix": "_tidyi.csv",
}


@dataclass(frozen=True)
class ComparisonGroup:
    name: str
    filters: dict[str, tuple[str, ...]]


@dataclass(frozen=True)
class ComparisonDefinition:
    name: str
    dataset: str
    groups: tuple[ComparisonGroup, ...]


COMPARISON_DEFINITIONS = (
    ComparisonDefinition(
        name="tms_young_male_vs_young_female",
        dataset="tms",
        groups=(
            ComparisonGroup("young_male", {"age_group": ("young",), "sex": ("male",)}),
            ComparisonGroup("young_female", {"age_group": ("young",), "sex": ("female",)}),
        ),
    ),
    ComparisonDefinition(
        name="tms_old_male_vs_old_female",
        dataset="tms",
        groups=(
            ComparisonGroup("old_male", {"age_group": ("old",), "sex": ("male",)}),
            ComparisonGroup("old_female", {"age_group": ("old",), "sex": ("female",)}),
        ),
    ),
    ComparisonDefinition(
        name="tms_age_mixed_balanced_young_old",
        dataset="tms",
        groups=(
            ComparisonGroup("young", {"age": ("1m", "3m")}),
            ComparisonGroup("old", {"age": ("21m", "24m")}),
        ),
    ),
    ComparisonDefinition(
        name="mukin_young_spf_vs_young_gf",
        dataset="mukin",
        groups=(
            ComparisonGroup("young_spf", {"model": ("SPF",), "age": ("2m",)}),
            ComparisonGroup("young_gf", {"model": ("GF",), "age": ("2m",)}),
        ),
    ),
    ComparisonDefinition(
        name="mukin_young_spf_vs_old_spf",
        dataset="mukin",
        groups=(
            ComparisonGroup("young_spf", {"model": ("SPF",), "age": ("2m",)}),
            ComparisonGroup("old_spf", {"model": ("SPF",), "age": ("19m",)}),
        ),
    ),
    ComparisonDefinition(
        name="mukin_old_spf_vs_old_gf",
        dataset="mukin",
        groups=(
            ComparisonGroup("old_spf", {"model": ("SPF",), "age": ("19m",)}),
            ComparisonGroup("old_gf", {"model": ("GF",), "age": ("19m",)}),
        ),
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run comparison-level TIDYI analysis on ASTRID-annotated datasets.")
    parser.add_argument("--mukin-bm-root", type=Path, default=DEFAULT_MUKIN_BM_ROOT)
    parser.add_argument("--mukin-skin-root", type=Path, default=DEFAULT_MUKIN_SKIN_ROOT)
    parser.add_argument("--tms-bm-root", type=Path, default=DEFAULT_TMS_BM_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--comparisons",
        nargs="+",
        default=[comparison.name for comparison in COMPARISON_DEFINITIONS],
        help="Subset of supported comparison names to run.",
    )
    parser.add_argument("--raw-layer", default="raw")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return parser.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def require_path(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")


def normalize_string_series(series: pd.Series) -> pd.Series:
    normalized = series.astype("string").str.strip()
    return normalized.replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})


def discover_astrid_h5ad_files(root: Path) -> list[Path]:
    require_path(root, f"ASTRID input root {root}")
    paths = sorted(root.glob("*/*_outASTRID.h5ad"))
    if not paths:
        raise FileNotFoundError(f"No '*_outASTRID.h5ad' files were found under {root}")
    return paths


def load_and_concat_h5ad(paths: list[Path], dataset_name: str) -> sc.AnnData:
    adatas = []
    for path in paths:
        LOGGER.info("Loading %s sample: %s", dataset_name, path)
        adatas.append(sc.read_h5ad(path))
    combined = ad.concat(adatas, join="outer", merge="same", uns_merge="unique", index_unique=None)
    LOGGER.info("%s concatenated shape: %s cells x %s genes", dataset_name, combined.n_obs, combined.n_vars)
    return combined


def require_obs_columns(adata: sc.AnnData, columns: list[str], context: str) -> None:
    missing = [column for column in columns if column not in adata.obs.columns]
    if missing:
        raise KeyError(f"{context} is missing required obs columns: {missing}")


def get_selected_comparisons(selected_names: list[str]) -> list[ComparisonDefinition]:
    selected_set = set(selected_names)
    available = {comparison.name: comparison for comparison in COMPARISON_DEFINITIONS}
    missing = sorted(selected_set - set(available))
    if missing:
        raise ValueError(f"Unknown comparison(s): {missing}. Valid names: {sorted(available)}")
    return [comparison for comparison in COMPARISON_DEFINITIONS if comparison.name in selected_set]


def fill_missing_tms_age(adata: sc.AnnData) -> sc.AnnData:
    require_obs_columns(adata, ["sample", "age"], "TMS age repair input")
    adata = adata.copy()
    adata.obs["age"] = normalize_string_series(adata.obs["age"])
    sample_age_map = (
        adata.obs.dropna(subset=["age"])
        .groupby("sample", observed=True)["age"]
        .first()
        .to_dict()
    )
    missing_mask = adata.obs["age"].isna()
    adata.obs.loc[missing_mask, "age"] = adata.obs.loc[missing_mask, "sample"].map(sample_age_map)
    unresolved = int(adata.obs["age"].isna().sum())
    if unresolved:
        raise ValueError(f"TMS age repair left {unresolved} cells without age")
    return adata


def fill_missing_tms_sex(adata: sc.AnnData) -> sc.AnnData:
    require_obs_columns(adata, ["mouse.id"], "TMS sex repair input")
    adata = adata.copy()
    if "sex" not in adata.obs.columns:
        adata.obs["sex"] = pd.NA
    adata.obs["sex"] = normalize_string_series(adata.obs["sex"])
    missing_mask = adata.obs["sex"].isna()
    derived = adata.obs.loc[missing_mask, "mouse.id"].astype(str).str.extract(r"-([MF])-", expand=False)
    derived = derived.map({"M": "male", "F": "female"})
    adata.obs.loc[missing_mask, "sex"] = derived.values
    unresolved = int(adata.obs["sex"].isna().sum())
    if unresolved:
        raise ValueError(f"TMS sex repair left {unresolved} cells without sex")
    return adata


def add_tms_age_group(adata: sc.AnnData) -> sc.AnnData:
    require_obs_columns(adata, ["age"], "TMS age-group input")
    adata = adata.copy()
    ages = normalize_string_series(adata.obs["age"])
    unknown_ages = sorted(set(ages.dropna().astype(str).unique()) - set(TMS_AGE_GROUP_MAP))
    if unknown_ages:
        raise ValueError(f"Unsupported TMS ages for age-group mapping: {unknown_ages}")
    adata.obs["age"] = ages
    adata.obs["age_group"] = ages.map(TMS_AGE_GROUP_MAP)
    return adata


def remove_blacklisted_tms_celltypes(adata: sc.AnnData) -> sc.AnnData:
    require_obs_columns(adata, ["SingleR_CellType"], "TMS blacklist filtering input")
    labels = normalize_string_series(adata.obs["SingleR_CellType"])
    mask = ~labels.isin(TMS_REMOVE_CELLTYPES)
    filtered = adata[mask.to_numpy()].copy()
    filtered.obs["SingleR_CellType"] = labels.loc[filtered.obs_names].astype(str).values
    LOGGER.info("TMS blacklist filter kept %s / %s cells", filtered.n_obs, adata.n_obs)
    return filtered


def build_tms_celltype_hb(adata: sc.AnnData) -> sc.AnnData:
    require_obs_columns(adata, ["SingleR_CellType", "mouse.id", "final_clustering_level"], "TMS relabeling input")
    adata = adata.copy()
    adata.obs["celltype_hb"] = normalize_string_series(adata.obs["SingleR_CellType"]).astype(str).values

    erythrocyte_mask = pd.Series(False, index=adata.obs_names)
    for mouse_id, prefixes in TMS_ERYTHROCYTE_RULES.items():
        mouse_mask = adata.obs["mouse.id"].astype(str) == mouse_id
        if not mouse_mask.any():
            continue
        cluster_values = adata.obs.loc[mouse_mask, "final_clustering_level"].astype(str)
        prefix_mask = cluster_values.apply(lambda value: any(value.startswith(prefix) for prefix in prefixes))
        erythrocyte_mask.loc[cluster_values.index] = prefix_mask.to_numpy()
    adata.obs.loc[erythrocyte_mask, "celltype_hb"] = "erythrocyte"

    current_labels = adata.obs["celltype_hb"].astype(str)
    erythroblast_mask = current_labels.eq("erythrocyte") & adata.obs["final_clustering_level"].astype(str).apply(
        lambda value: any(value.startswith(prefix) for prefix in TMS_ERYTHROBLAST_PREFIXES)
    )
    adata.obs.loc[erythroblast_mask, "celltype_hb"] = "erythroblast"
    LOGGER.info(
        "TMS erythrocyte relabeling: erythrocyte=%s erythroblast=%s",
        int((adata.obs["celltype_hb"].astype(str) == "erythrocyte").sum()),
        int((adata.obs["celltype_hb"].astype(str) == "erythroblast").sum()),
    )
    return adata


def find_rbc_contaminated_mukin_clusters(combined_csv: Path) -> set[tuple[str, str]]:
    dataframe = pd.read_csv(combined_csv)
    clustering_levels = []
    for _, row in dataframe.iterrows():
        if pd.notna(row.get("clustering_level_6")):
            clustering_levels.append(str(row["clustering_level_6"]))
        elif pd.notna(row.get("clustering_level_5")):
            clustering_levels.append(str(row["clustering_level_5"]))
        else:
            clustering_levels.append(None)
    rbc_mask = dataframe["OddsRatioMarkerGenes"].apply(lambda value: any(gene in str(value) for gene in MUKIN_RBC_GENES))
    flagged = {
        (str(sample_id), str(level))
        for sample_id, level in zip(dataframe.loc[rbc_mask, "sample_id"], clustering_levels)
        if pd.notna(sample_id) and level is not None
    }
    LOGGER.info("Mukin BM RBC contamination detection flagged %s sample/cluster pairs", len(flagged))
    return flagged


def relabel_mukin_rbc_contamination(adata: sc.AnnData, flagged_clusters: set[tuple[str, str]]) -> sc.AnnData:
    require_obs_columns(adata, ["mouseID", "final_clustering_level", "SingleR_CellType"], "Mukin RBC relabeling input")
    adata = adata.copy()
    labels = normalize_string_series(adata.obs["SingleR_CellType"]).astype(str)
    contamination_mask = [
        (str(sample_id), str(cluster_level)) in flagged_clusters
        for sample_id, cluster_level in zip(adata.obs["mouseID"], adata.obs["final_clustering_level"])
    ]
    contamination_mask = pd.Series(contamination_mask, index=adata.obs_names)
    relabeled = "RBC_contam_" + labels.loc[contamination_mask].astype(str)
    labels.loc[contamination_mask] = relabeled.values
    labels = labels.replace({"RBC_contam_erythrocyte": "erythrocyte"})
    adata.obs["SingleR_CellType"] = labels.values
    LOGGER.info("Mukin BM RBC contamination relabeling affected %s cells", int(contamination_mask.sum()))
    return adata


def filter_rare_celltypes(adata: sc.AnnData, column: str, min_cells: int = 100) -> sc.AnnData:
    require_obs_columns(adata, [column], f"Rare-celltype filtering input for {column}")
    labels = normalize_string_series(adata.obs[column])
    counts = labels.value_counts(dropna=True)
    keep = set(counts[counts >= min_cells].index.astype(str).tolist())
    mask = labels.astype(str).isin(keep) & labels.notna()
    filtered = adata[mask.to_numpy()].copy()
    filtered.obs[column] = labels.loc[filtered.obs_names].astype(str).values
    LOGGER.info(
        "Rare-celltype filter on %s kept %s / %s cells and %s labels",
        column,
        filtered.n_obs,
        adata.n_obs,
        len(keep),
    )
    return filtered


def prepare_mukin_dataset(mukin_bm_root: Path, mukin_skin_root: Path) -> sc.AnnData:
    mukin_bm = load_and_concat_h5ad(discover_astrid_h5ad_files(mukin_bm_root), "mukin_bm")
    mukin_skin = load_and_concat_h5ad(discover_astrid_h5ad_files(mukin_skin_root), "mukin_skin")
    flagged_clusters = find_rbc_contaminated_mukin_clusters(mukin_bm_root / "BM_combined_ASTRID_results.csv")
    mukin_bm = relabel_mukin_rbc_contamination(mukin_bm, flagged_clusters)
    mukin_bm = filter_rare_celltypes(mukin_bm, "SingleR_CellType", min_cells=100)
    mukin_skin = filter_rare_celltypes(mukin_skin, "SingleR_CellType", min_cells=100)
    combined = ad.concat([mukin_bm, mukin_skin], join="outer", merge="same", uns_merge="unique", index_unique=None)
    combined.obs["SingleR_CellType"] = normalize_string_series(combined.obs["SingleR_CellType"]).astype(str).values
    return combined


def prepare_tms_dataset(tms_bm_root: Path) -> sc.AnnData:
    tms = load_and_concat_h5ad(discover_astrid_h5ad_files(tms_bm_root), "tms_bm")
    tms = fill_missing_tms_age(tms)
    tms = fill_missing_tms_sex(tms)
    tms = add_tms_age_group(tms)
    tms = remove_blacklisted_tms_celltypes(tms)
    tms = build_tms_celltype_hb(tms)
    tms = filter_rare_celltypes(tms, "celltype_hb", min_cells=100)
    return tms


def ensure_ensembl_ids_with_local_mapping(
    adata: sc.AnnData,
    mapping_csv: Path = GENE_MAPPING_FILE,
) -> sc.AnnData:
    adata = adata.copy()
    map_mouse_symbols_to_ensembl(adata, mapping_csv)
    return adata


def preconvert_dataset_gene_ids(
    datasets: dict[str, sc.AnnData],
    mapping_csv: Path = GENE_MAPPING_FILE,
) -> dict[str, sc.AnnData]:
    converted: dict[str, sc.AnnData] = {}
    for dataset_name, adata in datasets.items():
        LOGGER.info("Pre-converting gene IDs using the local reference: %s", dataset_name)
        converted[dataset_name] = ensure_ensembl_ids_with_local_mapping(adata, mapping_csv)
    return converted


def run_tidyi_analysis(
    adata: sc.AnnData,
    sample_col: str,
    celltype_col: str,
    species: str = "mouse",
    raw_layer: str = "raw",
) -> pd.DataFrame:
    """Retain the original TIDYI grouping/rates, with shared offline mapping."""
    from tidyi import tidyi_main_mouse

    require_columns(adata.obs, [sample_col, celltype_col], "TIDYI input obs")
    if raw_layer not in adata.layers:
        raise KeyError(f"Missing raw layer '{raw_layer}' required for TIDYI analysis.")

    adata = adata.copy()
    adata.obs[SP_CT_COL] = adata.obs[sample_col].astype(str) + "_" + adata.obs[celltype_col].astype(str)
    if species != "mouse":
        raise ValueError("This handover reference supports mouse genes only.")
    map_mouse_symbols_to_ensembl(adata)
    adata.layers["counts"] = adata.layers[raw_layer].copy()

    # TIDYI internally subsets AnnData objects and then writes into .X/.var.
    # Those operations are safe for this workflow but create noisy view warnings in logs.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ImplicitModificationWarning)
        warnings.filterwarnings("ignore", category=UserWarning, message="Received a view of an AnnData.*")
        tidy_model = tidyi_main_mouse.tidy_rates_mouse(gene_list=adata.var_names.tolist())
        tidy_model.train_classifier()

    results: list[dict[str, Any]] = []
    for sample_celltype in adata.obs[SP_CT_COL].astype(str).unique():
        sub_data = adata[adata.obs[SP_CT_COL].astype(str) == sample_celltype].copy()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ImplicitModificationWarning)
            warnings.filterwarnings("ignore", category=UserWarning, message="Received a view of an AnnData.*")
            prolif_rate, death_rate, _, _ = tidy_model.apply_classifier(sub_data)
        sample_name, celltype = sample_celltype.rsplit("_", 1)
        results.append(
            {
                "sample": sample_name,
                "celltype": celltype,
                "prolif_rate": float(prolif_rate),
                "death_rate": float(death_rate),
            }
        )

    return pd.DataFrame(results)


def validate_group_filters(adata: sc.AnnData, group: ComparisonGroup, comparison_name: str) -> None:
    for column_name, allowed_values in group.filters.items():
        require_obs_columns(adata, [column_name], f"{comparison_name} input")
        observed = set(normalize_string_series(adata.obs[column_name]).dropna().astype(str).unique().tolist())
        missing = sorted(set(allowed_values) - observed)
        if missing:
            raise ValueError(
                f"{comparison_name} group '{group.name}' expects values {missing} in column '{column_name}', "
                f"but observed values include {sorted(observed)[:20]}"
            )


def build_comparison_subset(
    adata: sc.AnnData,
    comparison: ComparisonDefinition,
    celltype_col: str,
) -> sc.AnnData:
    require_obs_columns(adata, ["tissue", celltype_col], f"{comparison.name} subset input")
    group_assignments = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
    selected_mask = pd.Series(False, index=adata.obs_names)

    for group in comparison.groups:
        validate_group_filters(adata, group, comparison.name)
        group_mask = pd.Series(True, index=adata.obs_names)
        for column_name, allowed_values in group.filters.items():
            values = normalize_string_series(adata.obs[column_name])
            group_mask &= values.astype(str).isin(allowed_values)
        if not group_mask.any():
            raise ValueError(f"{comparison.name} group '{group.name}' selected zero cells")
        overlap = selected_mask & group_mask
        if overlap.any():
            raise ValueError(f"{comparison.name} produced overlapping comparison groups")
        selected_mask |= group_mask
        group_assignments.loc[group_mask] = group.name

    subset = adata[selected_mask.to_numpy()].copy()
    subset.obs[COMPARISON_GROUP_COL] = group_assignments.loc[subset.obs_names].astype(str).values
    subset.obs[COMPARISON_NAME_COL] = comparison.name
    subset.obs[ANALYSIS_DATASET_COL] = comparison.dataset
    subset.obs[celltype_col] = normalize_string_series(subset.obs[celltype_col]).astype(str).values
    subset.obs["tissue"] = normalize_string_series(subset.obs["tissue"]).astype(str).values
    subset.obs[ANALYSIS_TISSUE_COL] = (
        subset.obs[COMPARISON_GROUP_COL].astype(str) + ANALYSIS_TISSUE_DELIMITER + subset.obs["tissue"].astype(str)
    )
    return subset


def add_percentage(adata: sc.AnnData, celltype_col: str) -> sc.AnnData:
    adata = adata.copy()
    grouped = adata.obs.groupby([ANALYSIS_TISSUE_COL, celltype_col], observed=True).size().astype(float)
    totals = adata.obs.groupby(ANALYSIS_TISSUE_COL, observed=True).size().astype(float)
    mapped_totals = grouped.index.get_level_values(0).map(totals).astype(float)
    percentages = (grouped / mapped_totals) * 100.0
    key_tuples = list(zip(adata.obs[ANALYSIS_TISSUE_COL], adata.obs[celltype_col]))
    adata.obs[PERCENTAGE_COL] = [float(percentages.get(key, 0.0)) for key in key_tuples]
    return adata


def stringify_unique_values(series: pd.Series) -> str | pd.NA:
    values = normalize_string_series(series).dropna().astype(str).unique().tolist()
    if not values:
        return pd.NA
    return "|".join(sorted(values))


def build_final_tidyi_table(
    adata: sc.AnnData,
    tidyi_df: pd.DataFrame,
    dataset_name: str,
    comparison_name: str,
    celltype_col: str,
) -> pd.DataFrame:
    base_keys = [ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, "tissue", celltype_col]
    metadata_frame = (
        adata.obs.groupby(base_keys, observed=True)
        .agg(
            model=("model", stringify_unique_values) if "model" in adata.obs.columns else (celltype_col, lambda _: pd.NA),
            age=("age", stringify_unique_values) if "age" in adata.obs.columns else (celltype_col, lambda _: pd.NA),
            sex=("sex", stringify_unique_values) if "sex" in adata.obs.columns else (celltype_col, lambda _: pd.NA),
            percentage=(PERCENTAGE_COL, "first"),
        )
        .reset_index()
        .rename(columns={celltype_col: "celltype"})
    )
    for optional_column in ["model", "age", "sex"]:
        if optional_column not in metadata_frame.columns:
            metadata_frame[optional_column] = pd.NA

    merged = tidyi_df.merge(
        metadata_frame,
        on=[ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, "tissue", "celltype"],
        how="left",
    )
    merged[ANALYSIS_DATASET_COL] = dataset_name
    merged[COMPARISON_NAME_COL] = comparison_name
    merged = add_derived_metrics(merged, comparison_name)
    return merged[FINAL_TIDY_COLUMNS].sort_values(
        [ANALYSIS_DATASET_COL, COMPARISON_NAME_COL, COMPARISON_GROUP_COL, "tissue", "celltype"]
    ).reset_index(drop=True)


def make_comparison_output_dirs(output_root: Path, comparison_name: str) -> dict[str, Path]:
    comparison_dir = ensure_output_dir(output_root / comparison_name)
    plot_root_dir = ensure_output_dir(output_root / OUTPUT_LAYOUT["plot_root_dir_name"])
    ty_df_dir = ensure_output_dir(output_root / OUTPUT_LAYOUT["ty_df_dir_name"])
    tidyi_dir = ensure_output_dir(comparison_dir / OUTPUT_LAYOUT["tidyi_dir_name"])
    ty_plot_dir = ensure_output_dir(plot_root_dir / OUTPUT_LAYOUT["ty_plot_dir_name"])
    return {
        "comparison_dir": comparison_dir,
        "tidyi_dir": tidyi_dir,
        "ty_df_dir": ty_df_dir,
        "plot_dir": ty_plot_dir,
    }


def save_tidyi_table(table: pd.DataFrame, output_path: Path) -> None:
    ensure_output_dir(output_path.parent)
    table.to_csv(output_path, index=False)
    LOGGER.info("Saved TIDYI table: %s", output_path)


def plot_final_results(dataframes: dict[str, pd.DataFrame], output_root: Path) -> None:
    shared_plot_dir = ensure_output_dir(output_root / OUTPUT_LAYOUT["plot_root_dir_name"] / OUTPUT_LAYOUT["ty_plot_dir_name"])
    comparison_output_dirs = {comparison_name: shared_plot_dir for comparison_name in dataframes}
    plot_kinetics_scatter(dataframes, comparison_output_dirs)
    plot_growth_state_histplot(dataframes, comparison_output_dirs)
    plot_cross_tissue_heatmaps(dataframes, comparison_output_dirs)
    plot_growth_turnover_scatter(dataframes, comparison_output_dirs)


def run_single_comparison(
    adata: sc.AnnData,
    comparison: ComparisonDefinition,
    celltype_col: str,
    output_root: Path,
    raw_layer: str,
) -> pd.DataFrame:
    subset = build_comparison_subset(adata, comparison, celltype_col)
    subset = add_percentage(subset, celltype_col)
    if subset.n_obs == 0:
        raise ValueError(f"{comparison.name} produced an empty subset")

    comparison_dirs = make_comparison_output_dirs(output_root, comparison.name)
    tidyi_df = run_tidyi_analysis(
        adata=subset,
        sample_col=ANALYSIS_TISSUE_COL,
        celltype_col=celltype_col,
        raw_layer=raw_layer,
    )
    tidyi_df = expand_analysis_sample_column(tidyi_df)
    final_table = build_final_tidyi_table(
        adata=subset,
        tidyi_df=tidyi_df,
        dataset_name=comparison.dataset,
        comparison_name=comparison.name,
        celltype_col=celltype_col,
    )
    filename = f"{comparison.name}{OUTPUT_LAYOUT['tidyi_csv_suffix']}"
    save_tidyi_table(final_table, comparison_dirs["tidyi_dir"] / filename)
    save_tidyi_table(final_table, comparison_dirs["ty_df_dir"] / filename)
    return final_table


def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)

    selected_comparisons = get_selected_comparisons(args.comparisons)
    load_mouse_gene_mapping()
    ensure_output_dir(args.output_root)

    datasets = {
        "mukin": prepare_mukin_dataset(args.mukin_bm_root, args.mukin_skin_root),
        "tms": prepare_tms_dataset(args.tms_bm_root),
    }
    datasets = preconvert_dataset_gene_ids(datasets)
    celltype_columns = {
        "mukin": "SingleR_CellType",
        "tms": "celltype_hb",
    }

    final_tables: dict[str, pd.DataFrame] = {}
    for comparison in selected_comparisons:
        LOGGER.info("Running comparison: %s", comparison.name)
        final_table = run_single_comparison(
            adata=datasets[comparison.dataset],
            comparison=comparison,
            celltype_col=celltype_columns[comparison.dataset],
            output_root=args.output_root,
            raw_layer=args.raw_layer,
        )
        final_tables[comparison.name] = final_table

    if final_tables:
        plot_final_results(final_tables, args.output_root)


if __name__ == "__main__":
    main()

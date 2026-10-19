#!/usr/bin/env python3
"""Objective 3 rate producer refactored from TIDYI_analysis_spl-merge.ipynb.

The historical rate calculation uses the handover's frozen offline gene-symbol
mapping instead of a live MyGene query. This producer writes rate CSVs.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault(
    "NUMBA_CACHE_DIR",
    str(Path(__file__).resolve().parents[1] / "outputs/.numba_cache"),
)

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import ImplicitModificationWarning
from scipy import sparse
LOGGER = logging.getLogger("tidyi_comparison_analysis")


RA_ROOT = Path(__file__).resolve().parents[2]
if str(RA_ROOT) not in sys.path:
    sys.path.insert(0, str(RA_ROOT))

from gene_mapping import GENE_MAPPING_FILE, map_mouse_symbols_to_ensembl


PIPELINE_DIR = Path(__file__).resolve().parent
PATTERN_DIR = PIPELINE_DIR.parent
PROJECT_ROOT = PATTERN_DIR.parent.parent


DEFAULT_INPUT_PATHS = {
    "mukin_h5ad": PROJECT_ROOT / "mukin_qcdata" / "mukin_data_tissuewise_knn.h5ad",
    "tms_h5ad": PROJECT_ROOT / "tms_qcdata" / "tms_10x_tissuewise_knn.h5ad",
    "mukin_stats_csv": PATTERN_DIR / "mukin_sample_stas.csv",
    "tms_stats_csv": PATTERN_DIR / "tms_sample_stas.csv",
}


COLUMN_CONFIG = {
    "mukin": {
        "sample": "Sample",
        "sample_stats_sample": "Sample",
        "tissue": "tissue",
        "age": "age",
        "model": "model",
        "base_celltype": "tissue_wise_knn",
        "analysis_celltype": "analysis_celltype",
        "cross_tissue_celltype": "celltye_cross-tissue",
    },
    "tms": {
        "sample": "sample",
        "sample_stats_sample": "sample",
        "sample_stats_mouse_id": "mouse.id",
        "tissue": "tissue",
        "age": "age",
        "sex": "sex",
        "mouse_id": "mouse.id",
        "base_celltype": "tissue_wise_knn",
        "analysis_celltype": "analysis_celltype",
        "cross_tissue_celltype": "celltye_cross-tissue",
    },
}


FILTER_CUTOFFS = {
    "mukin_young_spf_vs_old_spf": 80,
    "mukin_young_spf_vs_young_gf": 80,
    "mukin_old_spf_vs_old_gf": 80,
    "tms_age_mixed": 80,
    "tms_age_mixed_balanced_young_old": 80,
    "tms_young_male_vs_young_female": 80,
    "tms_old_male_vs_old_female": 80,
    "tms_veryold_male_vs_veryold_female": 80,
    "tms_sex_matched_3m_18m_male_vs_female": 80,
}


OUTPUT_LAYOUT = {
    "ty_df_dir_name": "ty_df",
    "tidyi_dir_name": "tidyi",
    "tidyi_csv_suffix": "_tidyi.csv",
}


MUKIN_QUERY_CELLTYPES = (
    "HF Quiescennt cell",
    "Proliferative B cell",
    "Proliferative Tcell",
    "Proliferating KC",
)


MUKIN_REANNOTATION_WHITELISTS = {
    "HF Quiescennt cell": (
        "Basal KC",
        "HF Bulge cell",
        "HF cell",
    ),
    "Proliferating KC": (
        "Basal KC",
        "HF Bulge cell",
        "HF cell",
        "Spinosum KC",
    ),
    "Proliferative B cell": (
        "Immature B cell",
        "Memory B cell",
        "Naive B cell",
        "Pre-B cell",
    ),
    "Proliferative Tcell": (
        "Cd4-effector Tcell",
        "Cd4-naive Tcell",
        "Cd8-naive Tcell",
        "Cd8-effector Tcell",
        "ILCs",
        "Natural Killer",
    ),
}


MUKIN_CELLTYPE_RENAME_MAP = {
    "Pro-monocyte": "promonocyte",
    "Macrophage": "macrophage",
    "Pre-B cell": "precursor B cell",
    "Natural Killer": "NK cell",
}


MUKIN_TISSUE_RENAME_MAP = {
    "Bone-marrow": "Marrow",
    "Epidermis-E": "Epi-e",
}


TMS_AGE_GROUP_MAP = {
    "1m": "young",
    "3m": "young",
    "18m": "mid",
    "21m": "old",
    "24m": "old",
    "30m": "veryold",
}


DISCARD_CROSS_TISSUE_CELLTYPE = "discard"


TIDYI_ALLOWED_TISSUES = {
    "mukin": ("Bone-marrow", "Spleen", "Epidermis-E"),
    "tms": ("Marrow", "Spleen", "Tongue"),
}


MUKIN_CROSS_TISSUE_CELLTYPE_MAPPING = {
    "client progenitor": (
        "Basal KC",
        "Myeloid progenitor",
        "Neutrophil Progenitor",
        "Pro-monocyte",
        "promonocyte",
        "Pre-B cell",
        "precursor B cell",
    ),
    "client": (
        "Spinosum KC",
        "Granulosum KC",
        "HF Bulge cell",
        "HF cell",
        "Dermal Papilla",
        "ILCs",
        "Melanocyte",
        "Marginal B cell",
        "Memory B cell",
        "Cd4-naive Tcell",
        "Cd4-effector Tcell",
    ),
    "Neutrophil": ("Neutrophil",),
    "Immature B cell": ("Immature B cell",),
    "Naive B cell": ("Naive B cell",),
    "monocyte and macrophage": (
        "Cd14 Monocyte",
        "Cd16 Monocyte",
        "Alveolar Macrophage",
        "Macrophage",
        "macrophage",
        "TR-Macrophage",
    ),
    "naive Tcell": ("Cd8-naive Tcell",),
    "effector Tcell": (
        "Epidermal T cell",
        "Cd8-effector Tcell",
        "Resident Cd8-effector Tcell",
        "Th17-effector Tcell",
    ),
    "Dendritic cell": ("Dendritic cell",),
    "NK cell": (
        "NK cell",
        "Natural Killer",
    ),
    "Langerhans cell": ("Langerhans cell",),
    "Plasma Cell": ("Plasma Cell",),
    "Sebocyte": ("Sebocyte",),
    "Age-Associated B cell": ("Age-Associated B cell",),
    "Basophil": ("Basophil",),
    DISCARD_CROSS_TISSUE_CELLTYPE: (
        "Enteroendocrine cell",
        "Mesothelial cell",
        "intestinal stem cell",
    ),
}


TMS_CROSS_TISSUE_CELLTYPE_MAPPING = {
    "client progenitor": (
        "basal cell of epidermis",
        "erythroblast",
        "granulocytopoietic cell",
        "hematopoietic precursor cell",
        "macrophage dendritic cell progenitor",
        "megakaryocyte-erythroid progenitor cell",
        "precursor B cell",
        "proerythroblast",
        "promonocyte",
        "erythroid progenitor",
    ),
    "client": ("keratinocyte",),
    "B cell": ("B cell",),
    "Immature B cell": ("immature B cell",),
    "T cell": ("T cell",),
    "naive T cell": ("naive T cell",),
    "NK cell": (
        "mature NK T cell",
        "NK cell",
        "immature NKT cell",
    ),
    "monocyte and macrophage": (
        "monocyte",
        "macrophage",
    ),
    "Langerhans cell": ("Langerhans cell",),
    "plasma cell": ("plasma cell",),
    "granulocyte": ("granulocyte",),
    "basophil": ("basophil",),
    "late pro-B cell": ("late pro-B cell",),
    "naive B cell": ("naive B cell",),
}


CROSS_TISSUE_CELLTYPE_MAPPINGS = {
    "mukin": MUKIN_CROSS_TISSUE_CELLTYPE_MAPPING,
    "tms": TMS_CROSS_TISSUE_CELLTYPE_MAPPING,
}


ANALYSIS_TISSUE_COL = "analysis_tissue"
COMPARISON_GROUP_COL = "comparison_group"
COMPARISON_NAME_COL = "comparison_name"
ANALYSIS_DATASET_COL = "analysis_dataset"
BASE_CELLTYPE_EXPORT_COL = "base_celltype"
PERCENTAGE_COL = "percentage"
PRE_RELABEL_CELLTYPE_COL = "pre_relabel_celltype"
PRE_RELABEL_PERCENTAGE_COL = "pre_relabel_percentage"
PERCENTAGE_SOURCE_CELLTYPE_COL = "percentage_source_celltype"
PERCENTAGE_COMPONENT_COL = "percentage_component"
SP_CT_COL = "sp_ct"
ANALYSIS_TISSUE_DELIMITER = "||"


FINAL_TIDY_COLUMNS = [
    ANALYSIS_DATASET_COL,
    COMPARISON_NAME_COL,
    COMPARISON_GROUP_COL,
    "tissue",
    "celltype",
    "prolif_rate",
    "death_rate",
    "growth_rate",
    "turn_over",
    PERCENTAGE_COL,
    "model",
    "age",
    "sex",
]


@dataclass(frozen=True)
class ComparisonGroup:
    """One side of a comparison."""

    name: str
    filters: dict[str, tuple[str, ...]]


@dataclass(frozen=True)
class ComparisonDefinition:
    """Comparison configuration near the top of the script."""

    name: str
    dataset: str
    groups: tuple[ComparisonGroup, ...]
    min_cells: int


COMPARISON_DEFINITIONS = (
    ComparisonDefinition(
        name="mukin_young_spf_vs_old_spf",
        dataset="mukin",
        groups=(
            ComparisonGroup("young_spf", {"model": ("SPF",), "age": ("2m",)}),
            ComparisonGroup("old_spf", {"model": ("SPF",), "age": ("19m",)}),
        ),
        min_cells=FILTER_CUTOFFS["mukin_young_spf_vs_old_spf"],
    ),
    ComparisonDefinition(
        name="mukin_young_spf_vs_young_gf",
        dataset="mukin",
        groups=(
            ComparisonGroup("young_spf", {"model": ("SPF",), "age": ("2m",)}),
            ComparisonGroup("young_gf", {"model": ("GF",), "age": ("2m",)}),
        ),
        min_cells=FILTER_CUTOFFS["mukin_young_spf_vs_young_gf"],
    ),
    ComparisonDefinition(
        name="mukin_old_spf_vs_old_gf",
        dataset="mukin",
        groups=(
            ComparisonGroup("old_spf", {"model": ("SPF",), "age": ("19m",)}),
            ComparisonGroup("old_gf", {"model": ("GF",), "age": ("19m",)}),
        ),
        min_cells=FILTER_CUTOFFS["mukin_old_spf_vs_old_gf"],
    ),
    ComparisonDefinition(
        name="tms_age_mixed",
        dataset="tms",
        groups=(
            ComparisonGroup("1m", {"age": ("1m",)}),
            ComparisonGroup("3m", {"age": ("3m",)}),
            ComparisonGroup("18m", {"age": ("18m",)}),
            ComparisonGroup("21m", {"age": ("21m",)}),
            ComparisonGroup("24m", {"age": ("24m",)}),
            ComparisonGroup("30m", {"age": ("30m",)}),
        ),
        min_cells=FILTER_CUTOFFS["tms_age_mixed"],
    ),
    ComparisonDefinition(
        name="tms_age_mixed_balanced_young_old",
        dataset="tms",
        groups=(
            ComparisonGroup("young", {"age": ("1m", "3m")}),
            ComparisonGroup("old", {"age": ("21m", "24m")}),
        ),
        min_cells=FILTER_CUTOFFS["tms_age_mixed_balanced_young_old"],
    ),
    ComparisonDefinition(
        name="tms_young_male_vs_young_female",
        dataset="tms",
        groups=(
            ComparisonGroup("young_male", {"age_group": ("young",), "sex": ("male",)}),
            ComparisonGroup("young_female", {"age_group": ("young",), "sex": ("female",)}),
        ),
        min_cells=FILTER_CUTOFFS["tms_young_male_vs_young_female"],
    ),
    ComparisonDefinition(
        name="tms_old_male_vs_old_female",
        dataset="tms",
        groups=(
            ComparisonGroup("old_male", {"age_group": ("old",), "sex": ("male",)}),
            ComparisonGroup("old_female", {"age_group": ("old",), "sex": ("female",)}),
        ),
        min_cells=FILTER_CUTOFFS["tms_old_male_vs_old_female"],
    ),
    ComparisonDefinition(
        name="tms_veryold_male_vs_veryold_female",
        dataset="tms",
        groups=(
            ComparisonGroup("veryold_male", {"age_group": ("veryold",), "sex": ("male",)}),
        ),
        min_cells=FILTER_CUTOFFS["tms_veryold_male_vs_veryold_female"],
    ),
    ComparisonDefinition(
        name="tms_sex_matched_3m_18m_male_vs_female",
        dataset="tms",
        groups=(
            ComparisonGroup("young_male", {"age": ("3m",), "sex": ("male",)}),
            ComparisonGroup("young_female", {"age": ("3m",), "sex": ("female",)}),
            ComparisonGroup("mid_male", {"age": ("18m",), "sex": ("male",)}),
            ComparisonGroup("mid_female", {"age": ("18m",), "sex": ("female",)}),
        ),
        min_cells=FILTER_CUTOFFS["tms_sex_matched_3m_18m_male_vs_female"],
    ),
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run the comparison-based TIDYI analysis workflow.",
    )
    parser.add_argument(
        "--mukin-h5ad",
        type=Path,
        default=DEFAULT_INPUT_PATHS["mukin_h5ad"],
        help="Path to the Mukin AnnData file.",
    )
    parser.add_argument(
        "--tms-h5ad",
        type=Path,
        default=DEFAULT_INPUT_PATHS["tms_h5ad"],
        help="Path to the TMS AnnData file.",
    )
    parser.add_argument(
        "--mukin-stats-csv",
        type=Path,
        default=DEFAULT_INPUT_PATHS["mukin_stats_csv"],
        help="Path to the Mukin sample statistics table.",
    )
    parser.add_argument(
        "--tms-stats-csv",
        type=Path,
        default=DEFAULT_INPUT_PATHS["tms_stats_csv"],
        help="Path to the TMS sample statistics table.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Base output directory for comparison-specific results.",
    )
    parser.add_argument(
        "--comparisons",
        nargs="+",
        default=["all"],
        help="Comparison names to run, or 'all'.",
    )
    parser.add_argument(
        "--raw-layer",
        default="raw",
        help="Raw counts layer used by TIDYI.",
    )
    parser.add_argument(
        "--species",
        default="mouse",
        help="Species passed to gene ID conversion.",
    )
    parser.add_argument(
        "--entropy-threshold",
        type=float,
        default=0.3,
        help="Maximum normalized entropy kept during Mukin annotation transfer.",
    )
    parser.add_argument(
        "--knn-neighbors",
        type=int,
        default=10,
        help="Number of neighbors for kNN annotation transfer.",
    )
    parser.add_argument(
        "--n-pca",
        type=int,
        default=30,
        help="Number of PCA dimensions for annotation transfer.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def setup_logging(level: str) -> None:
    """Configure application logging."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def require_path(path: Path, description: str) -> None:
    """Raise a clear error if an input path is missing."""
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")


def require_columns(frame: pd.DataFrame, columns: list[str], context: str) -> None:
    """Validate required columns in a DataFrame-like object."""
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise KeyError(f"Missing required columns in {context}: {missing}")


def get_comparison_definitions(selected_names: list[str]) -> list[ComparisonDefinition]:
    """Resolve user-selected comparison names."""
    definitions = {definition.name: definition for definition in COMPARISON_DEFINITIONS}
    if selected_names == ["all"]:
        return list(COMPARISON_DEFINITIONS)

    unknown = [name for name in selected_names if name not in definitions]
    if unknown:
        valid_names = sorted(definitions)
        raise ValueError(f"Unknown comparison(s): {unknown}. Valid names: {valid_names}")

    return [definitions[name] for name in selected_names]


def ensure_output_dir(path: Path) -> Path:
    """Create an output directory if needed."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def require_nonempty_adata(adata: sc.AnnData, context: str) -> None:
    """Raise a clear error if an AnnData object has no cells."""
    if adata.n_obs == 0:
        raise ValueError(f"No cells remain in {context}.")


def require_obs_column(adata: sc.AnnData, column_name: str, context: str) -> None:
    """Require one obs column on AnnData with a clearer helper name."""
    require_columns(adata.obs, [column_name], context)


def validate_celltype_column(
    adata: sc.AnnData,
    celltype_col: str,
    context: str,
) -> None:
    """Check that the configured cell type column exists and contains usable values."""
    require_obs_column(adata, celltype_col, context)
    values = adata.obs[celltype_col]
    non_missing = values.dropna().astype(str)
    if non_missing.empty:
        raise ValueError(f"Cell type column '{celltype_col}' is empty in {context}.")
    invalid_values = {"", "nan", "None"}
    if set(non_missing.unique()).issubset(invalid_values):
        raise ValueError(
            f"Cell type column '{celltype_col}' has no valid annotations in {context}."
        )


def validate_percentage_column(adata: sc.AnnData, context: str) -> None:
    """Validate that percentage information exists after percentage calculation."""
    require_obs_column(adata, PERCENTAGE_COL, context)
    if adata.obs[PERCENTAGE_COL].isna().all():
        raise ValueError(f"Percentage column '{PERCENTAGE_COL}' is entirely missing in {context}.")


def stringify_unique_values(series: pd.Series) -> str | pd.NA:
    """Collapse unique non-missing metadata values into a stable string representation."""
    cleaned = [
        str(value)
        for value in series.dropna().astype(str).tolist()
        if str(value) not in {"", "nan", "None"}
    ]
    unique_values = sorted(set(cleaned))
    if not unique_values:
        return pd.NA
    if len(unique_values) == 1:
        return unique_values[0]
    return "|".join(unique_values)


def normalize_merge_keys(dataframe: pd.DataFrame, key_columns: list[str]) -> pd.DataFrame:
    """Normalize merge keys to stripped strings for stable joins."""
    dataframe = dataframe.copy()
    for column_name in key_columns:
        dataframe[column_name] = dataframe[column_name].astype(str).str.strip()
    return dataframe


def validate_whitelist_targets_exist(
    adata: sc.AnnData,
    celltype_col: str,
) -> None:
    """Log whitelist target labels that are absent from the current annotation space."""
    validate_celltype_column(adata, celltype_col, "whitelist target validation")
    observed_labels = set(adata.obs[celltype_col].dropna().astype(str).unique().tolist())
    missing_targets = sorted(
        {
            target
            for targets in MUKIN_REANNOTATION_WHITELISTS.values()
            for target in targets
            if target not in observed_labels
        }
    )
    if missing_targets:
        LOGGER.warning(
            "Some Mukin whitelist target labels are absent from the current annotation column '%s': %s. "
            "They will not be predicted by kNN for this comparison.",
            celltype_col,
            missing_targets,
        )


def load_inputs(args: argparse.Namespace) -> dict[str, Any]:
    """Load AnnData files and sample statistics tables."""
    require_path(args.mukin_h5ad, "Mukin h5ad")
    require_path(args.tms_h5ad, "TMS h5ad")
    require_path(args.mukin_stats_csv, "Mukin sample statistics CSV")
    require_path(args.tms_stats_csv, "TMS sample statistics CSV")

    LOGGER.info("Loading Mukin AnnData from %s", args.mukin_h5ad)
    mukin = sc.read_h5ad(args.mukin_h5ad)
    LOGGER.info("Loading TMS AnnData from %s", args.tms_h5ad)
    tms = sc.read_h5ad(args.tms_h5ad)

    LOGGER.info("Loading sample statistics tables")
    mukin_stats = pd.read_csv(args.mukin_stats_csv)
    tms_stats = pd.read_csv(args.tms_stats_csv)

    return {
        "mukin": mukin,
        "tms": tms,
        "mukin_stats": mukin_stats,
        "tms_stats": tms_stats,
    }


def filter_to_allowed_samples(
    adata: sc.AnnData,
    stats_df: pd.DataFrame,
    sample_col: str,
    stats_sample_col: str,
    dataset_name: str,
) -> sc.AnnData:
    """Keep only cells whose sample IDs appear in the sample statistics table."""
    require_columns(adata.obs, [sample_col], f"{dataset_name} AnnData obs")
    require_columns(stats_df, [stats_sample_col], f"{dataset_name} sample statistics")

    allowed_samples = set(stats_df[stats_sample_col].astype(str))
    current_samples = adata.obs[sample_col].astype(str)
    mask = current_samples.isin(allowed_samples)

    LOGGER.info(
        "%s sample filter kept %s / %s cells using %s sample IDs",
        dataset_name,
        int(mask.sum()),
        adata.n_obs,
        len(allowed_samples),
    )
    return adata[mask.to_numpy()].copy()


def derive_tms_metadata(
    adata: sc.AnnData,
    stats_df: pd.DataFrame,
    age_col: str,
    mouse_id_col: str,
) -> sc.AnnData:
    """Add TMS-specific metadata columns required by the configured comparisons."""
    require_columns(adata.obs, [age_col], "TMS AnnData obs")
    require_columns(stats_df, ["sample", mouse_id_col], "TMS sample statistics")

    adata = adata.copy()
    if mouse_id_col in adata.obs.columns:
        mouse_id_series = adata.obs[mouse_id_col].astype(str)
    else:
        sample_to_mouse_id = dict(zip(stats_df["sample"].astype(str), stats_df[mouse_id_col].astype(str)))
        mouse_id_series = adata.obs[COLUMN_CONFIG["tms"]["sample"]].astype(str).map(sample_to_mouse_id)
        if mouse_id_series.isna().any():
            missing_samples = sorted(
                adata.obs.loc[mouse_id_series.isna(), COLUMN_CONFIG["tms"]["sample"]].astype(str).unique().tolist()
            )
            raise ValueError(
                "Could not map mouse.id for some TMS samples from the sample statistics table: "
                f"{missing_samples[:10]}"
            )
        adata.obs[mouse_id_col] = mouse_id_series.values

    adata = repair_missing_tms_age_values(
        adata=adata,
        age_col=age_col,
        mouse_id_col=mouse_id_col,
    )

    age_values = adata.obs[age_col]
    non_missing_age_values = age_values.dropna().astype(str)
    missing_ages = sorted(set(non_missing_age_values.unique()) - set(TMS_AGE_GROUP_MAP))
    if missing_ages:
        raise ValueError(
            "Unexpected TMS age values. Update TMS_AGE_GROUP_MAP near the top of the script. "
            f"Observed unsupported ages: {missing_ages}"
        )
    adata.obs["age_group"] = age_values.map(TMS_AGE_GROUP_MAP)

    if "sex" not in adata.obs.columns:
        sex_series = adata.obs[mouse_id_col].astype(str).str.extract(r"-(male|female)-", expand=False)
        if sex_series.isna().any():
            bad_mouse_ids = sorted(adata.obs.loc[sex_series.isna(), mouse_id_col].astype(str).unique().tolist())
            raise ValueError(
                "Could not derive TMS sex from mouse.id. Unexpected values: "
                f"{bad_mouse_ids[:10]}"
            )
        adata.obs["sex"] = sex_series.values

    return adata


def repair_missing_tms_age_values(
    adata: sc.AnnData,
    age_col: str,
    mouse_id_col: str,
) -> sc.AnnData:
    """Repair missing TMS age values by propagating a unique age within each mouse.id."""
    require_columns(adata.obs, [age_col, mouse_id_col], "TMS age repair input obs")
    adata = adata.copy()

    age_series = adata.obs[age_col]
    normalized_age = age_series.replace(r"^\s*$", np.nan, regex=True)
    normalized_age = normalized_age.replace("nan", np.nan)
    normalized_age = normalized_age.replace("None", np.nan)

    missing_before = int(normalized_age.isna().sum())
    LOGGER.info("TMS age repair: %s cells have missing age before repair", missing_before)
    if missing_before == 0:
        adata.obs[age_col] = normalized_age
        return adata

    repaired = 0
    unresolved_conflicts: list[str] = []
    repaired_age = normalized_age.copy()

    for mouse_id, mouse_indices in adata.obs.groupby(mouse_id_col, observed=True).groups.items():
        mouse_age_values = normalized_age.loc[mouse_indices].dropna().astype(str).unique().tolist()
        missing_mask = repaired_age.loc[mouse_indices].isna()
        missing_count = int(missing_mask.sum())
        if missing_count == 0:
            continue
        if len(mouse_age_values) == 1:
            repaired_age.loc[mouse_indices[missing_mask]] = mouse_age_values[0]
            repaired += missing_count
        elif len(mouse_age_values) > 1:
            unresolved_conflicts.append(str(mouse_id))

    unresolved_after = int(repaired_age.isna().sum())
    if unresolved_conflicts:
        LOGGER.warning(
            "TMS age repair found conflicting non-missing ages for %s mouse.id values; "
            "left those cells unchanged. Example mouse.id values: %s",
            len(unresolved_conflicts),
            unresolved_conflicts[:10],
        )

    LOGGER.info(
        "TMS age repair completed: repaired=%s, still_unresolved=%s",
        repaired,
        unresolved_after,
    )
    adata.obs[age_col] = repaired_age
    return adata


def initialize_analysis_celltype_column(
    adata: sc.AnnData,
    base_celltype_col: str,
    analysis_celltype_col: str,
) -> sc.AnnData:
    """Preserve the original annotation and copy it into the analysis cell type column."""
    require_columns(adata.obs, [base_celltype_col], "analysis cell type initialization input obs")
    adata = adata.copy()
    adata.obs[analysis_celltype_col] = adata.obs[base_celltype_col].astype(str)
    LOGGER.info(
        "Initialized analysis cell type column '%s' from base annotation '%s'",
        analysis_celltype_col,
        base_celltype_col,
    )
    return adata


def invert_celltype_mapping(mapping: dict[str, tuple[str, ...]]) -> dict[str, str]:
    """Convert a new-to-old cell type mapping into an old-to-new lookup."""
    inverse: dict[str, str] = {}
    duplicate_old_names: list[tuple[str, str, str]] = []

    for new_celltype, old_celltypes in mapping.items():
        for old_celltype in old_celltypes:
            if old_celltype in inverse and inverse[old_celltype] != new_celltype:
                duplicate_old_names.append((old_celltype, inverse[old_celltype], new_celltype))
            inverse[old_celltype] = new_celltype

    if duplicate_old_names:
        details = "; ".join(
            f"{old}: {first} / {second}" for old, first, second in duplicate_old_names
        )
        raise ValueError(f"Duplicate old cell type names in cross-tissue mapping: {details}")

    return inverse


def add_cross_tissue_celltype_column(
    adata: sc.AnnData,
    dataset_name: str,
    source_celltype_col: str,
    target_celltype_col: str,
) -> sc.AnnData:
    """Add the configured cross-tissue cell type annotation without overwriting the source column."""
    require_columns(adata.obs, [source_celltype_col], f"{dataset_name} cross-tissue mapping input obs")
    mapping = CROSS_TISSUE_CELLTYPE_MAPPINGS[dataset_name]
    inverse_mapping = invert_celltype_mapping(mapping)

    source_values = adata.obs[source_celltype_col].astype(str)
    missing_celltypes = sorted(set(source_values.unique()) - set(inverse_mapping))
    if missing_celltypes:
        raise ValueError(
            f"{dataset_name} has cell types without cross-tissue mapping: {missing_celltypes}. "
            "Update CROSS_TISSUE_CELLTYPE_MAPPINGS near the top of the script."
        )

    adata = adata.copy()
    adata.obs[target_celltype_col] = source_values.map(inverse_mapping).astype(str).values
    discard_mask = adata.obs[target_celltype_col] == DISCARD_CROSS_TISSUE_CELLTYPE
    if discard_mask.any():
        discard_counts = (
            adata.obs.loc[discard_mask, source_celltype_col]
            .astype(str)
            .value_counts()
            .sort_index()
            .to_dict()
        )
        LOGGER.info(
            "%s cross-tissue mapping discarded %s cells by source label: %s",
            dataset_name,
            int(discard_mask.sum()),
            discard_counts,
        )
        adata = adata[~discard_mask.to_numpy()].copy()

    LOGGER.info(
        "%s cross-tissue cell type column '%s' added from '%s': %s source labels -> %s mapped labels",
        dataset_name,
        target_celltype_col,
        source_celltype_col,
        source_values.nunique(),
        adata.obs[target_celltype_col].nunique(),
    )
    return adata


def validate_group_filters(adata: sc.AnnData, group: ComparisonGroup, comparison_name: str) -> None:
    """Validate that a comparison group references existing columns and values."""
    for column_name, allowed_values in group.filters.items():
        require_columns(adata.obs, [column_name], f"{comparison_name} input obs")
        observed_values = set(adata.obs[column_name].astype(str).dropna().unique().tolist())
        missing_values = sorted(set(allowed_values) - observed_values)
        if missing_values:
            raise ValueError(
                f"{comparison_name} group '{group.name}' expects values {missing_values} in column "
                f"'{column_name}', but they are not present. Observed values include: "
                f"{sorted(list(observed_values))[:20]}"
            )


def build_comparison_subset(
    adata: sc.AnnData,
    comparison: ComparisonDefinition,
    tissue_col: str,
) -> sc.AnnData:
    """Subset one dataset into the cells belonging to a configured comparison."""
    group_assignments = pd.Series(index=adata.obs_names, dtype="object")
    selected_mask = pd.Series(False, index=adata.obs_names)

    for group in comparison.groups:
        validate_group_filters(adata, group, comparison.name)
        group_mask = pd.Series(True, index=adata.obs_names)
        for column_name, allowed_values in group.filters.items():
            group_mask &= adata.obs[column_name].astype(str).isin(allowed_values)

        if int(group_mask.sum()) == 0:
            raise ValueError(
                f"Comparison {comparison.name} group '{group.name}' selected zero cells."
            )

        overlap_mask = selected_mask & group_mask
        if overlap_mask.any():
            overlapping_indices = overlap_mask[overlap_mask].index[:10].tolist()
            raise ValueError(
                f"Comparison {comparison.name} produced overlapping groups. Example cell IDs: "
                f"{overlapping_indices}"
            )

        group_assignments.loc[group_mask[group_mask].index] = group.name
        selected_mask |= group_mask

    subset = adata[selected_mask.to_numpy()].copy()
    subset.obs[COMPARISON_GROUP_COL] = group_assignments.loc[subset.obs_names].astype(str).values
    subset.obs[COMPARISON_NAME_COL] = comparison.name
    subset.obs[ANALYSIS_DATASET_COL] = comparison.dataset

    require_columns(subset.obs, [tissue_col], f"{comparison.name} subset obs")
    if subset.obs[tissue_col].astype(str).str.contains(ANALYSIS_TISSUE_DELIMITER, regex=False).any():
        raise ValueError(
            f"Column '{tissue_col}' contains the reserved delimiter '{ANALYSIS_TISSUE_DELIMITER}'."
        )

    subset.obs[ANALYSIS_TISSUE_COL] = (
        subset.obs[COMPARISON_GROUP_COL].astype(str)
        + ANALYSIS_TISSUE_DELIMITER
        + subset.obs[tissue_col].astype(str)
    )
    return subset


def filter_celltypes_by_tissue(
    adata: sc.AnnData,
    tissue_col: str,
    celltype_col: str,
    min_cells: int,
    verbose: bool = False,
) -> sc.AnnData:
    """Remove cell types whose count in a given tissue falls below the minimum threshold."""
    require_columns(adata.obs, [tissue_col, celltype_col], "celltype filter input obs")
    counts = (
        adata.obs.groupby([tissue_col, celltype_col], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    valid = counts[counts["n_cells"] >= min_cells][[tissue_col, celltype_col]]

    if verbose:
        dropped = counts[counts["n_cells"] < min_cells]
        if dropped.empty:
            LOGGER.info("No tissue-celltype groups removed at threshold %s", min_cells)
        else:
            LOGGER.info(
                "Removing %s tissue-celltype groups with < %s cells",
                len(dropped),
                min_cells,
            )
            for _, row in dropped.iterrows():
                LOGGER.info(
                    "Dropped group: %s | %s | %s cells",
                    row[tissue_col],
                    row[celltype_col],
                    row["n_cells"],
                )

    obs_keys = adata.obs[[tissue_col, celltype_col]].reset_index(drop=True)
    valid_index = valid.set_index([tissue_col, celltype_col]).index
    mask = obs_keys.set_index([tissue_col, celltype_col]).index.isin(valid_index)
    return adata[mask].copy()


def export_tissue_celltype_count_summary(
    adata: sc.AnnData,
    comparison_name: str,
    tissue_col: str,
    celltype_col: str,
    output_path: Path,
) -> None:
    """Save one diagnose summary per comparison with tissue and cell type counts."""
    require_columns(
        adata.obs,
        [COMPARISON_GROUP_COL, tissue_col, celltype_col],
        "tissue-celltype count summary input obs",
    )
    ensure_output_dir(output_path.parent)

    summary = (
        adata.obs.groupby([COMPARISON_GROUP_COL, tissue_col, celltype_col], observed=True)
        .size()
        .rename("n_cells")
        .to_dict()
    )

    rows = [
        {
            "comparison": comparison_name,
            COMPARISON_GROUP_COL: comparison_group,
            tissue_col: tissue,
            "celltype": celltype,
            "n_cells": n_cells,
        }
        for (comparison_group, tissue, celltype), n_cells in summary.items()
    ]
    pd.DataFrame(rows).to_csv(output_path, index=False)
    LOGGER.info("Saved tissue-celltype count summary: %s", output_path)


def prepare_percentages_before_relabel(
    adata: sc.AnnData,
    sample_col: str,
    celltype_col: str,
) -> sc.AnnData:
    """Compute percentage before Mukin relabeling and preserve the pre-relabel state."""
    adata = add_celltype_percentage_vectorized(
        adata=adata,
        sample_col=sample_col,
        celltype_col=celltype_col,
    )
    adata.obs[PRE_RELABEL_CELLTYPE_COL] = adata.obs[celltype_col].astype(str)
    adata.obs[PRE_RELABEL_PERCENTAGE_COL] = adata.obs[PERCENTAGE_COL].astype(float)
    adata.obs[PERCENTAGE_SOURCE_CELLTYPE_COL] = adata.obs[PRE_RELABEL_CELLTYPE_COL].astype(str)
    adata.obs[PERCENTAGE_COMPONENT_COL] = adata.obs[PRE_RELABEL_PERCENTAGE_COL].astype(float)
    return adata


def split_reference_and_query_celltypes(
    adata: sc.AnnData,
    celltype_col: str,
) -> tuple[sc.AnnData, sc.AnnData]:
    """Split Mukin cells into reference and query subsets for notebook-style relabeling."""
    require_columns(adata.obs, [celltype_col], "Mukin annotation transfer input obs")
    celltypes = adata.obs[celltype_col].astype(str)
    query_mask = celltypes.isin(MUKIN_QUERY_CELLTYPES)
    query = adata[query_mask.to_numpy(), :].copy()
    ref = adata[(~query_mask).to_numpy(), :].copy()

    if query.n_obs == 0:
        LOGGER.warning("No Mukin query cells found for proliferating/quiescent reannotation")
    if ref.n_obs == 0:
        raise ValueError("No reference cells remain after selecting Mukin query cell types.")
    return ref, query


def ref_to_query_annotation(
    ref: Any,
    query: Any,
    ref_types: pd.Series,
    query_index: pd.Index,
    n_neighbors: int = 10,
    n_pca: int = 30,
) -> pd.DataFrame:
    """Simple and fast kNN-based annotation transfer from the notebook."""
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors

    def entropy(probabilities: np.ndarray) -> float:
        total = 0.0
        for probability in probabilities:
            if probability == 0:
                total += 0.0
            else:
                total += probability * np.log(probability)
        return -total

    def normalize(matrix: Any, scale_factor: float = 1e4) -> Any:
        total_counts = np.ravel(matrix.sum(axis=1))
        scaling_factor = scale_factor / total_counts
        if sparse.issparse(matrix):
            matrix = matrix.multiply(scaling_factor.reshape(-1, 1))
        else:
            matrix = matrix * scaling_factor[:, np.newaxis]
        return np.log1p(matrix)

    if query.shape[0] == 0:
        return pd.DataFrame(
            columns=["cell_type", "entropy", "max_count", "orig_ct"],
            index=query_index,
        )

    median_total_counts = np.median(np.ravel(ref.sum(axis=1)))
    ref = normalize(ref, median_total_counts)
    query = normalize(query, median_total_counts)

    pca = PCA(n_components=n_pca, svd_solver="auto", random_state=0)
    ref_embedding = pca.fit_transform(ref)
    query_embedding = pca.transform(query)

    neighbors = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    neighbors.fit(ref_embedding)
    _, indices_all = neighbors.kneighbors(query_embedding, return_distance=True)

    unique_labels = np.unique(ref_types.astype(str))
    label_to_index = {label: index for index, label in enumerate(unique_labels)}

    label_distribution = np.zeros(len(unique_labels), dtype=float)
    for label in ref_types.astype(str):
        label_distribution[label_to_index[label]] += 1.0
    label_distribution /= label_distribution.sum()
    max_entropy = entropy(label_distribution)

    ref_type_codes = np.array([label_to_index[label] for label in ref_types.astype(str)])
    neighbor_codes = ref_type_codes[indices_all]
    n_labels = len(unique_labels)
    counts = np.zeros((query.shape[0], n_labels), dtype=float)
    for label_index in range(n_labels):
        counts[:, label_index] = (neighbor_codes == label_index).sum(axis=1)
    probabilities = counts / counts.sum(axis=1, keepdims=True)

    norm_entropies = np.array([entropy(probability_row) / max_entropy for probability_row in probabilities])
    predicted_indices = probabilities.argmax(axis=1)
    predicted_labels = unique_labels[predicted_indices]
    confidences = probabilities.max(axis=1)

    return pd.DataFrame(
        {
            "cell_type": predicted_labels,
            "entropy": norm_entropies,
            "max_count": confidences * float(n_neighbors),
        },
        index=query_index,
    )


def apply_constrained_relabel_rules(
    adata: sc.AnnData,
    transferred_df: pd.DataFrame,
    query_adata: sc.AnnData,
    celltype_col: str,
    entropy_threshold: float,
) -> sc.AnnData:
    """Apply whitelist-constrained Mukin relabeling and drop invalid query cells."""
    adata = adata.copy()
    if transferred_df.empty:
        return adata

    transferred_df = transferred_df.copy()
    if ANALYSIS_TISSUE_COL in query_adata.obs.columns:
        relabel_group_col = ANALYSIS_TISSUE_COL
    elif "tissue" in query_adata.obs.columns:
        relabel_group_col = "tissue"
    else:
        raise KeyError(
            "Query AnnData is missing the grouping tissue column required for constrained relabeling. "
            f"Expected '{ANALYSIS_TISSUE_COL}' or 'tissue'."
        )

    # Preserve the tissue/group column on the transfer table so downstream grouping is explicit.
    transferred_df[relabel_group_col] = query_adata.obs[relabel_group_col].astype(str).values
    transferred_df["orig_ct"] = query_adata.obs[celltype_col].astype(str)
    filtered = transferred_df[transferred_df["entropy"] < entropy_threshold].copy()
    filtered["allowed_targets"] = filtered["orig_ct"].map(MUKIN_REANNOTATION_WHITELISTS)
    invalid_query_types = sorted(filtered.loc[filtered["allowed_targets"].isna(), "orig_ct"].astype(str).unique().tolist())
    if invalid_query_types:
        raise ValueError(
            "Mukin constrained relabeling found query cell types without whitelist rules: "
            f"{invalid_query_types}"
        )
    filtered["is_allowed_target"] = filtered.apply(
        lambda row: row["cell_type"] in set(row["allowed_targets"] or ()),
        axis=1,
    )
    valid_assignments = filtered[filtered["is_allowed_target"]].copy()

    # Summarize relabel outcomes by original query cell type to make log inspection easier.
    query_type_counts = query_adata.obs[celltype_col].astype(str).value_counts().to_dict()
    entropy_pass_counts = filtered["orig_ct"].astype(str).value_counts().to_dict()
    kept_counts = valid_assignments["orig_ct"].astype(str).value_counts().to_dict()
    all_query_types = sorted(query_type_counts.keys())
    for query_type in all_query_types:
        total_count = int(query_type_counts.get(query_type, 0))
        entropy_pass_count = int(entropy_pass_counts.get(query_type, 0))
        kept_count = int(kept_counts.get(query_type, 0))
        entropy_drop_count = total_count - entropy_pass_count
        whitelist_drop_count = entropy_pass_count - kept_count
        total_drop_count = total_count - kept_count
        LOGGER.info(
            "Mukin relabel summary for %s: total=%s kept=%s dropped=%s "
            "(entropy_filtered=%s, whitelist_filtered=%s)",
            query_type,
            total_count,
            kept_count,
            total_drop_count,
            entropy_drop_count,
            whitelist_drop_count,
        )

    kept_query_ids = set(valid_assignments.index.tolist())
    all_query_ids = set(query_adata.obs_names.tolist())
    dropped_query_ids = sorted(all_query_ids - kept_query_ids)

    if dropped_query_ids:
        adata = adata[~adata.obs_names.isin(dropped_query_ids)].copy()

    if not valid_assignments.empty:
        percentage_updates: list[dict[str, Any]] = []
        for (
            relabel_group_value,
            orig_ct,
            target_ct,
        ), assignment_df in valid_assignments.groupby([relabel_group_col, "orig_ct", "cell_type"], observed=True):
            source_mask = (
                (query_adata.obs[relabel_group_col].astype(str) == str(relabel_group_value))
                & (query_adata.obs[celltype_col].astype(str) == str(orig_ct))
            )
            source_percentages = query_adata.obs.loc[source_mask, PRE_RELABEL_PERCENTAGE_COL].astype(float).unique()
            if len(source_percentages) != 1:
                raise ValueError(
                    "Expected a unique pre-relabel percentage for each analysis_tissue/original cell type pair. "
                    f"Found {source_percentages.tolist()} for {relabel_group_value} / {orig_ct}"
                )
            source_percentage = float(source_percentages[0])
            valid_for_source = valid_assignments[
                (valid_assignments[relabel_group_col].astype(str) == str(relabel_group_value))
                & (valid_assignments["orig_ct"].astype(str) == str(orig_ct))
            ]
            allocated_percentage = source_percentage * (len(assignment_df) / len(valid_for_source))
            percentage_updates.append(
                {
                    relabel_group_col: relabel_group_value,
                    PRE_RELABEL_CELLTYPE_COL: orig_ct,
                    "target_celltype": target_ct,
                    PERCENTAGE_COMPONENT_COL: allocated_percentage,
                }
            )

        percentage_update_df = pd.DataFrame(percentage_updates)
        for _, row in percentage_update_df.iterrows():
            target_mask = (
                (adata.obs[relabel_group_col].astype(str) == str(row[relabel_group_col]))
                & (adata.obs[PRE_RELABEL_CELLTYPE_COL].astype(str) == str(row[PRE_RELABEL_CELLTYPE_COL]))
                & (adata.obs_names.isin(valid_assignments[
                    (valid_assignments[relabel_group_col].astype(str) == str(row[relabel_group_col]))
                    & (valid_assignments["orig_ct"].astype(str) == str(row[PRE_RELABEL_CELLTYPE_COL]))
                    & (valid_assignments["cell_type"].astype(str) == str(row["target_celltype"]))
                ].index))
            )
            adata.obs.loc[target_mask, PERCENTAGE_COMPONENT_COL] = float(row[PERCENTAGE_COMPONENT_COL])
            adata.obs.loc[target_mask, celltype_col] = str(row["target_celltype"])

    LOGGER.info(
        "Mukin constrained relabeling kept %s query cells and dropped %s query cells at entropy %.3f",
        len(valid_assignments),
        len(dropped_query_ids),
        entropy_threshold,
    )
    return adata


def filter_comparison_groups_separately(
    adata: sc.AnnData,
    comparison: ComparisonDefinition,
    tissue_col: str,
    celltype_col: str,
) -> sc.AnnData:
    """Apply the notebook filter separately to each internal comparison group."""
    retained_obs_names: list[str] = []
    for group in comparison.groups:
        group_adata = adata[adata.obs[COMPARISON_GROUP_COL].astype(str) == group.name].copy()
        if group_adata.n_obs == 0:
            LOGGER.warning("No cells available for filtering in %s / %s", comparison.name, group.name)
            continue
        filtered_group = filter_celltypes_by_tissue(
            adata=group_adata,
            tissue_col=tissue_col,
            celltype_col=celltype_col,
            min_cells=comparison.min_cells,
            verbose=True,
        )
        LOGGER.info(
            "Group-specific cutoff retained %s / %s cells for %s / %s",
            filtered_group.n_obs,
            group_adata.n_obs,
            comparison.name,
            group.name,
        )
        retained_obs_names.extend(filtered_group.obs_names.tolist())

    retained_mask = adata.obs_names.isin(retained_obs_names)
    return adata[retained_mask].copy()


def add_celltype_percentage_vectorized(
    adata: sc.AnnData,
    sample_col: str,
    celltype_col: str,
) -> sc.AnnData:
    """Vectorized notebook-style percentage calculation."""
    require_columns(adata.obs, [sample_col, celltype_col], "percentage input obs")
    adata = adata.copy()

    group_counts = adata.obs.groupby([sample_col, celltype_col]).size()
    sample_totals = adata.obs.groupby(sample_col).size()
    group_counts = group_counts.astype(float)
    sample_totals_mapped = group_counts.index.get_level_values(0).map(sample_totals).astype(float)
    percentages = (group_counts / sample_totals_mapped) * 100.0

    index_tuples = list(zip(adata.obs[sample_col], adata.obs[celltype_col]))
    adata.obs[PERCENTAGE_COL] = [percentages.get(index_tuple, 0.0) for index_tuple in index_tuples]
    return adata


def build_prepared_dataframe(
    adata: sc.AnnData,
    dataset_name: str,
    sample_col: str,
    tissue_col: str,
    base_celltype_col: str,
    analysis_celltype_col: str,
) -> pd.DataFrame:
    """Build the per-comparison prepared dataframe saved before TIDYI output."""
    columns_to_save = [
        sample_col,
        tissue_col,
        base_celltype_col,
        analysis_celltype_col,
        COMPARISON_GROUP_COL,
        COMPARISON_NAME_COL,
        ANALYSIS_DATASET_COL,
        ANALYSIS_TISSUE_COL,
        PERCENTAGE_COL,
    ]
    optional_columns = [
        "age",
        "model",
        "sex",
        "mouse.id",
        PRE_RELABEL_CELLTYPE_COL,
        PRE_RELABEL_PERCENTAGE_COL,
        PERCENTAGE_COMPONENT_COL,
    ]
    existing_optional_columns = [column for column in optional_columns if column in adata.obs.columns]

    require_columns(adata.obs, columns_to_save, f"{dataset_name} prepared dataframe input obs")
    prepared = adata.obs[columns_to_save + existing_optional_columns].copy()
    prepared = prepared.rename(
        columns={
            tissue_col: "tissue",
            sample_col: "sample",
            base_celltype_col: BASE_CELLTYPE_EXPORT_COL,
            analysis_celltype_col: "celltype",
        }
    )
    return prepared


def ensure_ensembl_ids(adata: sc.AnnData, species: str) -> sc.AnnData:
    """Apply the handover's frozen offline mouse symbol-to-Ensembl mapping."""
    if species.lower() != "mouse":
        raise ValueError(
            "Objective 3 rate generation only supports the frozen mouse mapping; "
            f"received species={species!r}"
        )
    adata = adata.copy()
    map_mouse_symbols_to_ensembl(adata, GENE_MAPPING_FILE)
    return adata


def run_tidyi_analysis(
    adata: sc.AnnData,
    sample_col: str,
    celltype_col: str,
    species: str = "mouse",
    raw_layer: str = "raw",
) -> pd.DataFrame:
    """Run TIDYI analysis following the notebook function design."""
    from tidyi import tidyi_main_mouse

    require_columns(adata.obs, [sample_col, celltype_col], "TIDYI input obs")
    if raw_layer not in adata.layers:
        raise KeyError(f"Missing raw layer '{raw_layer}' required for TIDYI analysis.")

    adata = adata.copy()
    adata.obs[SP_CT_COL] = adata.obs[sample_col].astype(str) + "_" + adata.obs[celltype_col].astype(str)
    adata = ensure_ensembl_ids(adata, species)
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


def expand_analysis_sample_column(results_df: pd.DataFrame) -> pd.DataFrame:
    """Split the internal analysis sample column into comparison group and tissue."""
    if "sample" not in results_df.columns:
        raise KeyError("TIDYI results are missing the 'sample' column.")

    results_df = results_df.copy()
    # Use a literal split because the delimiter contains regex metacharacters.
    split_frame = results_df["sample"].astype(str).str.split(
        ANALYSIS_TISSUE_DELIMITER,
        n=1,
        expand=True,
        regex=False,
    )
    if split_frame.shape[1] != 2:
        raise ValueError(
            "Could not split TIDYI sample values into comparison group and tissue. "
            f"Expected delimiter '{ANALYSIS_TISSUE_DELIMITER}'."
        )
    results_df[COMPARISON_GROUP_COL] = split_frame[0].values
    results_df["tissue"] = split_frame[1].values
    results_df = results_df.rename(columns={"sample": ANALYSIS_TISSUE_COL})
    return results_df


def merge_percentage_into_tidyi_results(
    adata: sc.AnnData,
    tidyi_df: pd.DataFrame,
    celltype_col: str,
) -> pd.DataFrame:
    """Merge notebook-style percentage information into the TIDYI result table."""
    if PERCENTAGE_COMPONENT_COL in adata.obs.columns and PERCENTAGE_SOURCE_CELLTYPE_COL in adata.obs.columns:
        require_columns(
            adata.obs,
            [
                ANALYSIS_TISSUE_COL,
                COMPARISON_GROUP_COL,
                "tissue",
                celltype_col,
                PERCENTAGE_SOURCE_CELLTYPE_COL,
                PERCENTAGE_COMPONENT_COL,
            ],
            "percentage merge input obs",
        )
        percentage_info = adata.obs[
            [
                ANALYSIS_TISSUE_COL,
                COMPARISON_GROUP_COL,
                "tissue",
                celltype_col,
                PERCENTAGE_SOURCE_CELLTYPE_COL,
                PERCENTAGE_COMPONENT_COL,
            ]
        ].copy()
        percentage_info = percentage_info.drop_duplicates(
            subset=[
                ANALYSIS_TISSUE_COL,
                COMPARISON_GROUP_COL,
                "tissue",
                celltype_col,
                PERCENTAGE_SOURCE_CELLTYPE_COL,
            ]
        )
        percentage_info = (
            percentage_info.groupby(
                [ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, "tissue", celltype_col],
                observed=True,
            )[PERCENTAGE_COMPONENT_COL]
            .sum()
            .reset_index()
            .rename(columns={celltype_col: "celltype", PERCENTAGE_COMPONENT_COL: PERCENTAGE_COL})
        )
    else:
        require_columns(
            adata.obs,
            [ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, "tissue", celltype_col, PERCENTAGE_COL],
            "percentage merge input obs",
        )
        percentage_info = adata.obs[
            [ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, "tissue", celltype_col, PERCENTAGE_COL]
        ].copy()
        percentage_info = percentage_info.drop_duplicates(
            subset=[ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, "tissue", celltype_col]
        )
        percentage_info = percentage_info.rename(columns={celltype_col: "celltype"})

    merged = tidyi_df.merge(
        percentage_info,
        on=[ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, "tissue", "celltype"],
        how="left",
    )
    return merged


def save_tidyi_table(results_df: pd.DataFrame, output_path: Path) -> None:
    """Save one TIDYI result table."""
    ensure_output_dir(output_path.parent)
    results_df.to_csv(output_path, index=False)
    LOGGER.info("Saved TIDYI results: %s", output_path)


def make_comparison_output_dirs(base_output_dir: Path, comparison_name: str) -> dict[str, Path]:
    """Create the standard output directories for one comparison."""
    comparison_dir = ensure_output_dir(base_output_dir / comparison_name)
    tidyi_dir = ensure_output_dir(comparison_dir / OUTPUT_LAYOUT["tidyi_dir_name"])
    ty_df_dir = ensure_output_dir(base_output_dir / OUTPUT_LAYOUT["ty_df_dir_name"])
    return {
        "comparison_dir": comparison_dir,
        "tidyi_dir": tidyi_dir,
        "ty_df_dir": ty_df_dir,
    }


def normalize_result_labels(results_df: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    """Apply the small label normalizations present in the notebook."""
    results_df = results_df.copy()
    if dataset_name == "mukin":
        results_df["tissue"] = results_df["tissue"].replace(MUKIN_TISSUE_RENAME_MAP)
        results_df["celltype"] = results_df["celltype"].replace(MUKIN_CELLTYPE_RENAME_MAP)
    return results_df


def restrict_to_tidyi_tissues(
    adata: sc.AnnData,
    dataset_name: str,
    tissue_col: str,
) -> sc.AnnData:
    """Restrict the final TIDYI run to the configured tissues for each dataset."""
    require_obs_column(adata, tissue_col, f"{dataset_name} tissue restriction input")
    allowed_tissues = set(TIDYI_ALLOWED_TISSUES[dataset_name])
    mask = adata.obs[tissue_col].astype(str).isin(allowed_tissues)
    restricted = adata[mask.to_numpy()].copy()
    LOGGER.info(
        "Restricted %s to TIDYI tissues %s: kept %s / %s cells",
        dataset_name,
        sorted(allowed_tissues),
        restricted.n_obs,
        adata.n_obs,
    )
    require_nonempty_adata(restricted, f"{dataset_name} after TIDYI tissue restriction")
    return restricted


def require_only_tidyi_tissues(
    adata: sc.AnnData,
    dataset_name: str,
    tissue_col: str,
    context: str,
) -> None:
    """Validate that an intermediate AnnData contains only the selected good tissues."""
    require_obs_column(adata, tissue_col, context)
    allowed_tissues = set(TIDYI_ALLOWED_TISSUES[dataset_name])
    observed_tissues = set(adata.obs[tissue_col].dropna().astype(str).unique().tolist())
    unexpected_tissues = sorted(observed_tissues - allowed_tissues)
    if unexpected_tissues:
        raise ValueError(
            f"{context} contains tissues outside the selected TIDYI tissue set for {dataset_name}: "
            f"{unexpected_tissues}. Allowed tissues: {sorted(allowed_tissues)}"
        )


def summarize_metadata_and_percentage(
    adata: sc.AnnData,
    tissue_col: str,
    celltype_col: str,
) -> pd.DataFrame:
    """Build one clean metadata and percentage table keyed to the final TIDYI output."""
    base_keys = [ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, tissue_col, celltype_col]
    require_columns(
        adata.obs,
        [ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, tissue_col, celltype_col],
        "final metadata summary input obs",
    )

    agg_spec: dict[str, tuple[str, Any]] = {}
    if "model" in adata.obs.columns:
        agg_spec["model"] = ("model", stringify_unique_values)
    if "age" in adata.obs.columns:
        agg_spec["age"] = ("age", stringify_unique_values)
    if "sex" in adata.obs.columns:
        agg_spec["sex"] = ("sex", stringify_unique_values)
    metadata_frame = adata.obs.groupby(base_keys, observed=True).agg(**agg_spec).reset_index()

    if PERCENTAGE_COMPONENT_COL in adata.obs.columns and PERCENTAGE_SOURCE_CELLTYPE_COL in adata.obs.columns:
        percentage_frame = adata.obs[
            base_keys + [PERCENTAGE_SOURCE_CELLTYPE_COL, PERCENTAGE_COMPONENT_COL]
        ].copy()
        percentage_frame = percentage_frame.drop_duplicates(
            subset=base_keys + [PERCENTAGE_SOURCE_CELLTYPE_COL]
        )
        percentage_frame = (
            percentage_frame.groupby(base_keys, observed=True)[PERCENTAGE_COMPONENT_COL]
            .sum()
            .reset_index()
            .rename(columns={PERCENTAGE_COMPONENT_COL: PERCENTAGE_COL})
        )
    else:
        require_obs_column(adata, PERCENTAGE_COL, "final metadata summary input")
        percentage_frame = adata.obs[base_keys + [PERCENTAGE_COL]].copy().drop_duplicates(subset=base_keys)

    summary = metadata_frame.merge(percentage_frame, on=base_keys, how="outer")
    summary = summary.rename(columns={tissue_col: "tissue", celltype_col: "celltype"})
    for optional_column in ["model", "age", "sex"]:
        if optional_column not in summary.columns:
            summary[optional_column] = pd.NA
    return summary


def build_final_tidyi_table(
    adata: sc.AnnData,
    tidyi_df: pd.DataFrame,
    dataset_name: str,
    comparison_name: str,
    tissue_col: str,
    celltype_col: str,
    debug_output_dir: Path,
) -> pd.DataFrame:
    """Merge TIDYI outputs with clean metadata and derived metrics."""
    merge_keys = [ANALYSIS_TISSUE_COL, COMPARISON_GROUP_COL, "tissue", "celltype"]
    metadata_summary = summarize_metadata_and_percentage(
        adata=adata,
        tissue_col=tissue_col,
        celltype_col=celltype_col,
    )
    metadata_summary = normalize_merge_keys(metadata_summary, merge_keys)
    tidyi_df = normalize_merge_keys(tidyi_df, merge_keys)

    # Merge at the same granularity as the final exported tidy table.
    merged = tidyi_df.merge(
        metadata_summary,
        on=merge_keys,
        how="left",
    )
    merged[ANALYSIS_DATASET_COL] = dataset_name
    merged[COMPARISON_NAME_COL] = comparison_name
    merged["growth_rate"] = merged["prolif_rate"] - merged["death_rate"]
    merged["turn_over"] = (np.sqrt(2) / 2) * (merged["prolif_rate"] + merged["death_rate"])
    merged = normalize_result_labels(merged, dataset_name)
    for optional_column in ["model", "age", "sex"]:
        if optional_column not in merged.columns:
            merged[optional_column] = pd.NA
    missing_percentage = int(merged[PERCENTAGE_COL].isna().sum())
    if missing_percentage > 0:
        unmatched_keys = merged.loc[merged[PERCENTAGE_COL].isna(), merge_keys].drop_duplicates()
        unmatched_path = debug_output_dir / f"{comparison_name}_unmatched_percentage_keys.csv"
        ensure_output_dir(debug_output_dir)
        unmatched_keys.to_csv(unmatched_path, index=False)
        LOGGER.warning(
            "Found %s final rows without percentage values for %s. Saved unmatched key combinations to %s",
            missing_percentage,
            comparison_name,
            unmatched_path,
        )
        raise ValueError(
            f"Missing percentage values in final TIDYI table for {dataset_name}: {missing_percentage} rows."
        )
    return merged[FINAL_TIDY_COLUMNS].sort_values(
        [ANALYSIS_DATASET_COL, COMPARISON_NAME_COL, COMPARISON_GROUP_COL, "tissue", "celltype"]
    ).reset_index(drop=True)


def run_single_comparison(
    adata: sc.AnnData,
    comparison: ComparisonDefinition,
    args: argparse.Namespace,
    tissue_col: str,
    base_celltype_col: str,
    analysis_celltype_col: str,
    cross_tissue_celltype_col: str,
    output_base_dir: Path,
) -> pd.DataFrame:
    """Run the full notebook-style workflow for one comparison."""
    LOGGER.info("Running comparison: %s", comparison.name)
    comparison_dirs = make_comparison_output_dirs(output_base_dir, comparison.name)

    comparison_adata = build_comparison_subset(adata, comparison, tissue_col=tissue_col)
    require_nonempty_adata(comparison_adata, f"{comparison.name} after comparison subsetting")
    require_obs_column(comparison_adata, tissue_col, f"{comparison.name} subset validation")
    require_obs_column(comparison_adata, ANALYSIS_TISSUE_COL, f"{comparison.name} subset validation")
    comparison_adata = restrict_to_tidyi_tissues(
        adata=comparison_adata,
        dataset_name=comparison.dataset,
        tissue_col=tissue_col,
    )
    require_only_tidyi_tissues(
        adata=comparison_adata,
        dataset_name=comparison.dataset,
        tissue_col=tissue_col,
        context=f"{comparison.name} after early good-tissue restriction",
    )
    validate_celltype_column(
        comparison_adata,
        analysis_celltype_col,
        f"{comparison.name} comparison subset",
    )

    if comparison.dataset == "mukin":
        # The Mukin transfer step must run on the fine-grained notebook annotation.
        comparison_adata = prepare_percentages_before_relabel(
            adata=comparison_adata,
            sample_col=ANALYSIS_TISSUE_COL,
            celltype_col=analysis_celltype_col,
        )
        validate_percentage_column(comparison_adata, f"{comparison.name} before Mukin relabeling")
        validate_whitelist_targets_exist(comparison_adata, analysis_celltype_col)
        ref_adata, query_adata = split_reference_and_query_celltypes(
            adata=comparison_adata,
            celltype_col=analysis_celltype_col,
        )
        transferred_df = ref_to_query_annotation(
            ref=ref_adata.X,
            query=query_adata.X,
            ref_types=ref_adata.obs[analysis_celltype_col],
            query_index=query_adata.obs_names,
            n_neighbors=args.knn_neighbors,
            n_pca=args.n_pca,
        )
        comparison_adata = apply_constrained_relabel_rules(
            adata=comparison_adata,
            transferred_df=transferred_df,
            query_adata=query_adata,
            celltype_col=analysis_celltype_col,
            entropy_threshold=args.entropy_threshold,
        )
        require_nonempty_adata(comparison_adata, f"{comparison.name} after Mukin constrained relabeling")
        validate_celltype_column(
            comparison_adata,
            analysis_celltype_col,
            f"{comparison.name} after Mukin constrained relabeling",
        )

    # The fixed April rate tables used the fine analysis_celltype labels.
    # Cross-tissue aggregation was added to this script later and is not part of
    # the Objective 3 rate lineage.
    if comparison.dataset != "mukin":
        comparison_adata = prepare_percentages_before_relabel(
            adata=comparison_adata,
            sample_col=ANALYSIS_TISSUE_COL,
            celltype_col=analysis_celltype_col,
        )
    validate_percentage_column(
        comparison_adata,
        f"{comparison.name} before fine-cell-type cutoff filtering",
    )

    comparison_adata = filter_comparison_groups_separately(
        adata=comparison_adata,
        comparison=comparison,
        tissue_col=tissue_col,
        celltype_col=analysis_celltype_col,
    )
    require_nonempty_adata(comparison_adata, f"{comparison.name} after cutoff filtering")
    require_obs_column(comparison_adata, tissue_col, f"{comparison.name} after cutoff filtering")
    require_only_tidyi_tissues(
        adata=comparison_adata,
        dataset_name=comparison.dataset,
        tissue_col=tissue_col,
        context=f"{comparison.name} after cutoff filtering",
    )
    validate_celltype_column(
        comparison_adata,
        analysis_celltype_col,
        f"{comparison.name} after cutoff filtering",
    )

    tidyi_df = run_tidyi_analysis(
        adata=comparison_adata,
        sample_col=ANALYSIS_TISSUE_COL,
        celltype_col=analysis_celltype_col,
        species=args.species,
        raw_layer=args.raw_layer,
    )
    tidyi_df = expand_analysis_sample_column(tidyi_df)
    final_tidyi_df = build_final_tidyi_table(
        adata=comparison_adata,
        tidyi_df=tidyi_df,
        dataset_name=comparison.dataset,
        comparison_name=comparison.name,
        tissue_col=tissue_col,
        celltype_col=analysis_celltype_col,
        debug_output_dir=comparison_dirs["comparison_dir"],
    )

    save_tidyi_table(
        results_df=final_tidyi_df,
        output_path=comparison_dirs["tidyi_dir"] / f"{comparison.name}{OUTPUT_LAYOUT['tidyi_csv_suffix']}",
    )
    save_tidyi_table(
        results_df=final_tidyi_df,
        output_path=comparison_dirs["ty_df_dir"] / f"{comparison.name}{OUTPUT_LAYOUT['tidyi_csv_suffix']}",
    )
    return final_tidyi_df


def main() -> None:
    """Run the terminal analysis workflow."""
    args = parse_args()
    setup_logging(args.log_level)

    selected_comparisons = get_comparison_definitions(args.comparisons)
    inputs = load_inputs(args)

    mukin = filter_to_allowed_samples(
        adata=inputs["mukin"],
        stats_df=inputs["mukin_stats"],
        sample_col=COLUMN_CONFIG["mukin"]["sample"],
        stats_sample_col=COLUMN_CONFIG["mukin"]["sample_stats_sample"],
        dataset_name="Mukin",
    )
    tms = filter_to_allowed_samples(
        adata=inputs["tms"],
        stats_df=inputs["tms_stats"],
        sample_col=COLUMN_CONFIG["tms"]["sample"],
        stats_sample_col=COLUMN_CONFIG["tms"]["sample_stats_sample"],
        dataset_name="TMS",
    )
    tms = derive_tms_metadata(
        adata=tms,
        stats_df=inputs["tms_stats"],
        age_col=COLUMN_CONFIG["tms"]["age"],
        mouse_id_col=COLUMN_CONFIG["tms"]["mouse_id"],
    )
    mukin = initialize_analysis_celltype_column(
        adata=mukin,
        base_celltype_col=COLUMN_CONFIG["mukin"]["base_celltype"],
        analysis_celltype_col=COLUMN_CONFIG["mukin"]["analysis_celltype"],
    )
    tms = initialize_analysis_celltype_column(
        adata=tms,
        base_celltype_col=COLUMN_CONFIG["tms"]["base_celltype"],
        analysis_celltype_col=COLUMN_CONFIG["tms"]["analysis_celltype"],
    )

    datasets = {
        "mukin": mukin,
        "tms": tms,
    }

    ensure_output_dir(args.output_dir)
    for comparison in selected_comparisons:
        dataset_name = comparison.dataset
        dataset_adata = datasets[dataset_name]
        tissue_col = COLUMN_CONFIG[dataset_name]["tissue"]
        base_celltype_col = COLUMN_CONFIG[dataset_name]["base_celltype"]
        analysis_celltype_col = COLUMN_CONFIG[dataset_name]["analysis_celltype"]
        cross_tissue_celltype_col = COLUMN_CONFIG[dataset_name]["cross_tissue_celltype"]

        run_single_comparison(
            adata=dataset_adata,
            comparison=comparison,
            args=args,
            tissue_col=tissue_col,
            base_celltype_col=base_celltype_col,
            analysis_celltype_col=analysis_celltype_col,
            cross_tissue_celltype_col=cross_tissue_celltype_col,
            output_base_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()

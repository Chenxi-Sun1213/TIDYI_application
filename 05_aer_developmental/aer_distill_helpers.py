"""Small, notebook-facing helpers for the AER DISTILL reanalysis.

The module intentionally contains no pipeline runner and performs no work on import.
Each public function corresponds to one inspectable notebook stage.
"""

from __future__ import annotations

import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors


def ensure_output_dirs(output_root: Path) -> dict[str, Path]:
    """Create and return the deliberately small analysis directory layout."""
    paths = {
        "qc": output_root / "01_qc",
        "astrid": output_root / "02_astrid",
        "biolucid": output_root / "03_biolucid",
        "restore": output_root / "04_restore_high_mito",
        "tidyi": output_root / "05_tidyi",
        "final": output_root / "06_final",
        "growth_decomposition": output_root / "07_growth_rate_decomposition",
        "sample_mean_umi": output_root / "08_sample_mean_umi",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def ensure_marker_formula_copy(source_path: Path, project_path: Path) -> Path:
    """Create one editable project-local copy of ASTRID's marker formulas."""
    if not source_path.exists():
        raise FileNotFoundError(f"ASTRID marker formula source not found: {source_path}")
    project_path.parent.mkdir(parents=True, exist_ok=True)
    if not project_path.exists():
        shutil.copyfile(source_path, project_path)

    marker_formulas = pd.read_csv(project_path)
    require_columns(
        marker_formulas,
        ["cell_type", "marker_genes", "expected_markers"],
        "ASTRID marker formula CSV",
    )
    return project_path


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns).difference(frame.columns))
    if missing:
        raise KeyError(f"{label} is missing required columns: {missing}")


def matrix_is_nonnegative_integer_counts(matrix, sample_size: int = 100_000) -> bool:
    """Cheap guard that accepts sparse or dense raw UMI matrices."""
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    if values.size == 0:
        return True
    values = values[:sample_size]
    return bool(np.all(values >= 0) and np.allclose(values, np.round(values)))


def initialize_raw_counts(adata: ad.AnnData) -> ad.AnnData:
    """Validate the reconstructed AER object and preserve its raw X in counts."""
    require_columns(
        adata.obs,
        ["sample", "gsm", "genotype", "stage", "chemistry", "author_cell_annotation"],
        "adata.obs",
    )
    require_columns(adata.var, ["gene_ids", "gene_symbol"], "adata.var")
    if not adata.obs_names.is_unique:
        raise ValueError("AER obs_names are not unique")
    if not adata.var_names.is_unique:
        raise ValueError("AER var_names are not unique")
    if not matrix_is_nonnegative_integer_counts(adata.X):
        raise ValueError("adata.X does not look like non-negative integer UMI counts")

    out = adata.copy()
    out.layers["counts"] = out.X.copy()
    return out


def calculate_qc_metrics(adata: ad.AnnData) -> None:
    """Compute only the metrics needed by DISTILL hard QC and mito splitting."""
    symbols = adata.var["gene_symbol"].astype(str)
    adata.var["mt"] = symbols.str.lower().str.startswith("mt-").to_numpy()
    adata.X = adata.layers["counts"].copy()
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )


def apply_distill_qc(
    adata: ad.AnnData,
    *,
    min_genes_per_cell: int,
    max_genes_per_cell: int,
    min_counts_per_cell: int,
    max_counts_per_cell: int,
    max_pct_mito_for_astrid: float,
) -> tuple[ad.AnnData, ad.AnnData, pd.DataFrame]:
    """Apply DISTILL hard cell QC, then split high-mito cells from ASTRID input.

    No genes are removed. High-mito cells remain in ``full_post_qc`` and are only
    excluded from ``for_astrid``.
    """
    work = adata.copy()
    calculate_qc_metrics(work)

    work.obs["qc_fail_min_genes"] = (
        work.obs["n_genes_by_counts"] < min_genes_per_cell
    )
    work.obs["qc_fail_max_genes"] = (
        work.obs["n_genes_by_counts"] > max_genes_per_cell
    )
    work.obs["qc_fail_min_counts"] = work.obs["total_counts"] < min_counts_per_cell
    work.obs["qc_fail_max_counts"] = work.obs["total_counts"] > max_counts_per_cell
    hard_fail_cols = [
        "qc_fail_min_genes",
        "qc_fail_max_genes",
        "qc_fail_min_counts",
        "qc_fail_max_counts",
    ]
    work.obs["qc_hard_pass"] = ~work.obs[hard_fail_cols].any(axis=1)
    work.obs["qc_high_mito"] = (
        work.obs["pct_counts_mt"] > max_pct_mito_for_astrid
    )

    full_post_qc = work[work.obs["qc_hard_pass"]].copy()
    for_astrid = full_post_qc[~full_post_qc.obs["qc_high_mito"]].copy()

    rows: list[dict[str, Any]] = []
    for sample, sample_obs in work.obs.groupby("sample", observed=True, sort=False):
        hard_pass = sample_obs["qc_hard_pass"]
        rows.append(
            {
                "sample": str(sample),
                "cells_input": int(len(sample_obs)),
                "fail_min_genes": int(sample_obs["qc_fail_min_genes"].sum()),
                "fail_max_genes": int(sample_obs["qc_fail_max_genes"].sum()),
                "fail_min_counts": int(sample_obs["qc_fail_min_counts"].sum()),
                "fail_max_counts": int(sample_obs["qc_fail_max_counts"].sum()),
                "cells_after_hard_qc": int(hard_pass.sum()),
                "high_mito_after_hard_qc": int(
                    (hard_pass & sample_obs["qc_high_mito"]).sum()
                ),
                "cells_for_astrid": int(
                    (hard_pass & ~sample_obs["qc_high_mito"]).sum()
                ),
            }
        )
    summary = pd.DataFrame(rows)

    if full_post_qc.n_vars != adata.n_vars or for_astrid.n_vars != adata.n_vars:
        raise AssertionError("QC unexpectedly changed the gene count")
    if full_post_qc.n_obs != for_astrid.n_obs + int(full_post_qc.obs["qc_high_mito"].sum()):
        raise AssertionError("High-mito split does not reconcile with full_post_qc")
    return full_post_qc, for_astrid, summary


def save_qc_outputs(
    full_post_qc: ad.AnnData,
    for_astrid: ad.AnnData,
    summary: pd.DataFrame,
    qc_dir: Path,
) -> None:
    """Save combined QC objects and the only QC report requested: a table."""
    qc_dir.mkdir(parents=True, exist_ok=True)
    full_post_qc.write_h5ad(qc_dir / "aer_full_post_qc.h5ad", compression="gzip")
    for_astrid.write_h5ad(qc_dir / "aer_for_astrid.h5ad", compression="gzip")
    summary.to_csv(qc_dir / "qc_summary.csv", index=False)


def prepare_astrid_inputs(
    for_astrid: ad.AnnData,
    astrid_dir: Path,
    *,
    author_key: str = "author_cell_annotation",
) -> dict[str, Path]:
    """Write one raw-count ASTRID input per AER sample."""
    require_columns(for_astrid.obs, ["sample", author_key], "for_astrid.obs")
    input_dir = astrid_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    inputs: dict[str, Path] = {}

    for sample in for_astrid.obs["sample"].astype(str).drop_duplicates():
        sample_mask = for_astrid.obs["sample"].astype(str).eq(sample).to_numpy()
        sample_data = for_astrid[sample_mask].copy()
        sample_data.X = sample_data.layers["counts"].copy()
        sample_data.obs["author_validation_label"] = (
            sample_data.obs[author_key].astype(object).fillna("Unannotated")
        )
        # ASTRID expects gene symbols. Repeated symbols retain separate count columns
        # and receive AnnData suffixes such as "-1" to keep var_names unique.
        sample_data.var_names = pd.Index(
            sample_data.var["gene_symbol"].astype(str).to_numpy(), name=None
        )
        sample_data.var_names_make_unique()
        path = input_dir / f"{sample}_astrid_input.h5ad"
        sample_data.write_h5ad(path, compression="gzip")
        inputs[sample] = path
    return inputs


def build_astrid_commands(
    astrid_inputs: dict[str, Path],
    astrid_dir: Path,
    *,
    astrid_script: Path,
    astrid_reference: Path,
    marker_formula: Path,
    cutoff_level: int | None,
    python_executable: str = sys.executable,
) -> dict[str, list[str]]:
    """Build, but do not run, the per-sample commands used in the notebook."""
    commands: dict[str, list[str]] = {}
    for sample, input_path in astrid_inputs.items():
        sample_dir = astrid_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        command = [
            "nice",
            "-n",
            "19",
            python_executable,
            str(astrid_script),
            "--clustering",
            "--annotation",
            "--validation",
            "--input_file",
            str(input_path),
            "--input_prefix",
            sample,
            "--output_file",
            str(sample_dir / f"{sample}_outASTRID.h5ad"),
            "--output_clustering_results",
            str(sample_dir / f"{sample}_ASTRID_results.csv"),
            "--author_type",
            "author_validation_label",
            "--out_dir",
            str(sample_dir),
            "--reference",
            str(astrid_reference),
            "--marker_formula",
            str(marker_formula),
            "--species",
            "mouse",
        ]
        if cutoff_level is not None:
            command.extend(["--cutoff_level", str(cutoff_level)])
        commands[sample] = command
    return commands


def astrid_command_table(commands: dict[str, list[str]]) -> pd.DataFrame:
    return pd.DataFrame(
        [{"sample": sample, "command": shlex.join(command)} for sample, command in commands.items()]
    )


def run_astrid_commands(
    commands: dict[str, list[str]], astrid_dir: Path
) -> pd.DataFrame:
    """Run ASTRID sequentially and retain stdout/stderr in one log per sample."""
    rows: list[dict[str, Any]] = []
    for sample, command in commands.items():
        sample_dir = astrid_dir / sample
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        log_path = sample_dir / "astrid.log"
        log_path.write_text(
            "COMMAND:\n"
            + shlex.join(command)
            + f"\nRETURN_CODE: {completed.returncode}\n\nSTDOUT:\n"
            + completed.stdout
            + "\nSTDERR:\n"
            + completed.stderr,
            encoding="utf-8",
        )
        rows.append(
            {
                "sample": sample,
                "return_code": completed.returncode,
                "log": str(log_path),
            }
        )
        if completed.returncode != 0:
            raise RuntimeError(f"ASTRID failed for {sample}; inspect {log_path}")
    return pd.DataFrame(rows)


def load_astrid_outputs(
    astrid_dir: Path,
    samples: Iterable[str],
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Load uncurated ASTRID outputs and restore raw counts to X."""
    datasets: list[ad.AnnData] = []
    rows: list[dict[str, Any]] = []
    for sample in samples:
        path = astrid_dir / sample / f"{sample}_outASTRID.h5ad"
        if not path.exists():
            raise FileNotFoundError(f"Missing ASTRID output: {path}")
        sample_data = sc.read_h5ad(path)
        require_columns(
            sample_data.obs,
            [
                "SingleR_CellType",
                "SingleR_Pruned_CellType",
                "SingleR_All_deltanext",
                "final_clustering_level",
            ],
            f"{sample} ASTRID obs",
        )
        if "raw" not in sample_data.layers:
            raise KeyError(f"{sample} ASTRID output has no layers['raw']")
        if "sample" not in sample_data.obs:
            raise KeyError(f"{sample} ASTRID output has no obs['sample']")
        observed_samples = set(sample_data.obs["sample"].astype(str).unique())
        if observed_samples != {str(sample)}:
            raise ValueError(
                f"{sample} ASTRID output has unexpected sample values: {sorted(observed_samples)}"
            )
        sample_data.obs["astrid_suggested_cell_type"] = sample_data.obs[
            "SingleR_CellType"
        ].astype("string")
        sample_data.obs["astrid_pruned_cell_type"] = sample_data.obs[
            "SingleR_Pruned_CellType"
        ].astype("string")
        sample_data.layers["counts"] = sample_data.layers["raw"].copy()
        sample_data.X = sample_data.layers["counts"].copy()
        rows.append(
            {
                "sample": sample,
                "n_cells": int(sample_data.n_obs),
                "n_genes": int(sample_data.n_vars),
                "n_clusters": int(sample_data.obs["final_clustering_level"].nunique()),
                "n_suggested_cell_types": int(
                    sample_data.obs["astrid_suggested_cell_type"].nunique()
                ),
                "n_missing_suggested_labels": int(
                    sample_data.obs["astrid_suggested_cell_type"].isna().sum()
                ),
                "n_missing_pruned_labels": int(
                    sample_data.obs["astrid_pruned_cell_type"].isna().sum()
                ),
                "n_unassigned_clusters": int(
                    sample_data.obs["final_clustering_level"].astype(str).eq("unassigned").sum()
                ),
            }
        )
        datasets.append(sample_data)

    combined = ad.concat(datasets, axis=0, join="inner", merge="same", uns_merge="same")
    if not combined.obs_names.is_unique:
        raise ValueError("Combined ASTRID output has duplicate obs_names")
    return combined, pd.DataFrame(rows)


def build_astrid_cluster_review(
    astrid_unreviewed: ad.AnnData,
    astrid_dir: Path,
    samples: Iterable[str],
    *,
    cluster_key: str = "final_clustering_level",
) -> pd.DataFrame:
    """Combine per-sample ASTRID validation CSVs into a manual-review table."""
    require_columns(astrid_unreviewed.obs, ["sample", cluster_key], "ASTRID obs")
    rows: list[pd.DataFrame] = []
    for sample in samples:
        result_path = astrid_dir / sample / f"{sample}_ASTRID_results.csv"
        if not result_path.exists():
            raise FileNotFoundError(f"Missing ASTRID validation table: {result_path}")
        result = pd.read_csv(result_path, dtype={cluster_key: "string"})
        require_columns(result, [cluster_key, "SingleR_CellType"], f"{sample} ASTRID CSV")
        result[cluster_key] = result[cluster_key].astype("string")
        if result[cluster_key].isna().any() or result[cluster_key].duplicated().any():
            raise ValueError(f"{sample} ASTRID CSV must contain one non-empty row per cluster")

        expected = set(
            astrid_unreviewed.obs.loc[
                astrid_unreviewed.obs["sample"].astype(str).eq(str(sample)), cluster_key
            ].astype(str)
        )
        observed = set(result[cluster_key].astype(str))
        if expected != observed:
            raise ValueError(
                f"{sample} ASTRID CSV clusters do not match h5ad: "
                f"missing={sorted(expected - observed)}, extra={sorted(observed - expected)}"
            )
        result.insert(0, "sample", str(sample))
        rows.append(result)

    review = pd.concat(rows, ignore_index=True)
    review.insert(2, "cell_type", pd.Series(pd.NA, index=review.index, dtype="string"))
    review.insert(3, "review_notes", pd.Series(pd.NA, index=review.index, dtype="string"))
    return review


def apply_curated_astrid_annotations(
    astrid_unreviewed: ad.AnnData,
    curated_csv: Path,
    *,
    sample_key: str = "sample",
    cluster_key: str = "final_clustering_level",
    label_key: str = "cell_type",
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Attach any currently reviewed labels and report incomplete/ambiguous rows."""
    if not curated_csv.exists():
        raise FileNotFoundError(f"Curated ASTRID annotation CSV not found: {curated_csv}")
    require_columns(astrid_unreviewed.obs, [sample_key, cluster_key], "ASTRID obs")
    mapping = pd.read_csv(curated_csv, dtype={sample_key: "string", cluster_key: "string"})
    require_columns(mapping, [sample_key, cluster_key, label_key], "curated ASTRID CSV")

    key_columns = [sample_key, cluster_key]
    mapping[label_key] = mapping[label_key].astype("string").str.strip()
    invalid_labels = mapping[label_key].isna() | mapping[label_key].str.lower().isin(
        ["", "nan", "none", "null"]
    )
    mapping.loc[invalid_labels, label_key] = pd.NA
    duplicate_mask = mapping.duplicated(key_columns, keep=False)
    duplicate_key_count = int(
        mapping.loc[duplicate_mask, key_columns].drop_duplicates().shape[0]
    )
    # Do not block exploratory review. If a key occurs more than once, the last
    # non-empty label is used provisionally and the ambiguity remains in the audit.
    mapping_for_use = (
        mapping[key_columns + [label_key]]
        .groupby(key_columns, as_index=False, sort=False, dropna=False)[label_key]
        .last()
    )

    expected_frame = pd.DataFrame(
        {
            sample_key: astrid_unreviewed.obs[sample_key].astype("string").to_numpy(),
            cluster_key: astrid_unreviewed.obs[cluster_key].astype("string").to_numpy(),
        }
    ).drop_duplicates()
    expected_keys = set(map(tuple, expected_frame[key_columns].to_numpy()))
    mapping_keys = set(map(tuple, mapping_for_use[key_columns].to_numpy()))
    missing_keys = expected_keys - mapping_keys
    extra_keys = mapping_keys - expected_keys

    lookup = mapping_for_use.set_index(key_columns)[label_key]
    cell_keys = pd.MultiIndex.from_arrays(
        [
            astrid_unreviewed.obs[sample_key].astype("string"),
            astrid_unreviewed.obs[cluster_key].astype("string"),
        ],
        names=key_columns,
    )
    curated = astrid_unreviewed.copy()
    curated.obs[label_key] = lookup.reindex(cell_keys).to_numpy(dtype=object)
    curated.obs["cell_type_source"] = np.where(
        curated.obs[label_key].notna(),
        "manual_astrid_cluster_review",
        "unreviewed_astrid_cluster",
    )

    rows: list[dict[str, Any]] = []
    for sample in curated.obs[sample_key].astype(str).drop_duplicates():
        sample_obs = curated.obs[curated.obs[sample_key].astype(str).eq(sample)]
        sample_expected = {
            (sample, str(cluster)) for cluster in sample_obs[cluster_key].astype(str).unique()
        }
        sample_mapping = mapping_for_use[
            mapping_for_use[sample_key].astype(str).eq(sample)
        ]
        sample_labeled_keys = set(
            map(
                tuple,
                sample_mapping.loc[
                    sample_mapping[label_key].notna(), key_columns
                ].to_numpy(),
            )
        ).intersection(sample_expected)
        rows.append(
            {
                sample_key: sample,
                "clusters_total": len(sample_expected),
                "clusters_with_curated_label": len(sample_labeled_keys),
                "clusters_without_curated_label": len(sample_expected - sample_labeled_keys),
                "cells_total": int(len(sample_obs)),
                "cells_with_curated_label": int(sample_obs[label_key].notna().sum()),
                "cells_without_curated_label": int(sample_obs[label_key].isna().sum()),
                "mapping_missing_keys_total": len(missing_keys),
                "mapping_extra_keys_total": len(extra_keys),
                "mapping_duplicate_keys_total": duplicate_key_count,
            }
        )
    summary = pd.DataFrame(rows)
    return curated, summary


def preview_biolucid_groups(
    adata: ad.AnnData,
    *,
    batch_key: str,
    cell_type_key: str,
    min_cells: int,
) -> tuple[pd.DataFrame, list[str]]:
    """Show exactly which cell types satisfy BioLUCID's all-batches rule."""
    counts = pd.crosstab(adata.obs[batch_key], adata.obs[cell_type_key])
    valid = counts.columns[(counts >= min_cells).all(axis=0)].astype(str).tolist()
    return counts, valid


def plot_biolucid_scores(results_df: pd.DataFrame, output_path: Path):
    """Plot saved BioLUCID sample scores without relying on ``plt.show()``."""
    import seaborn as sns

    plot_data = results_df.copy().reset_index()
    plot_data = plot_data.rename(columns={plot_data.columns[0]: "sample"})
    plot_data["sample"] = plot_data["sample"].astype(str)
    plot_data["genotype"] = plot_data["sample"].str.extract(
        r"_(WT|Dac)_", expand=False
    )
    plot_data["stage"] = (
        plot_data["sample"]
        .str.extract(r"_(E\d+)$", expand=False)
        .replace({"E95": "E9.5", "E105": "E10.5", "E115": "E11.5"})
    )
    plot_data["sample_label"] = plot_data["genotype"] + " " + plot_data["stage"]

    figure, axis = plt.subplots(figsize=(7, 5.5))
    sns.scatterplot(
        data=plot_data,
        x="q_sp_score_per_batch",
        y="q_sh_score_per_batch",
        hue="genotype",
        style="stage",
        s=110,
        alpha=0.85,
        edgecolor="black",
        linewidth=0.4,
        ax=axis,
    )

    q_max = float(plot_data["q_sh_score_per_batch"].max()) * 1.1
    q = np.linspace(0.01, q_max, 100)
    axis.plot(q * np.sqrt(1 / 0.8 - 1), q, "r--", alpha=0.7, label="b=80%")
    axis.plot(q * np.sqrt(1 / 0.6 - 1), q, "b--", alpha=0.7, label="b=60%")

    label_offsets = {
        ("WT", "E9.5"): (6, 8),
        ("Dac", "E9.5"): (-6, 8),
        ("WT", "E10.5"): (6, 8),
        ("Dac", "E10.5"): (6, -14),
        ("WT", "E11.5"): (6, 8),
        ("Dac", "E11.5"): (6, 8),
    }
    for row in plot_data.itertuples(index=False):
        offset = label_offsets[(row.genotype, row.stage)]
        axis.annotate(
            row.sample_label,
            (row.q_sp_score_per_batch, row.q_sh_score_per_batch),
            xytext=offset,
            textcoords="offset points",
            fontsize=8,
            alpha=0.8,
            ha="right" if offset[0] < 0 else "left",
        )

    axis.set_xlabel("q$_{sp}$")
    axis.set_ylabel("q$_{sh}$")
    axis.set_title("BioLUCID sample-level batch effect scores")
    axis.grid(True, alpha=0.25)
    axis.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    return figure


def run_biolucid(
    adata: ad.AnnData,
    output_dir: Path,
    *,
    batch_key: str,
    cell_type_key: str,
    min_cells: int,
    abundant_gene_threshold: float,
    min_abundant_genes: int,
):
    """Run BioLUCID and save its tables and sample-level score plot."""
    import biolucid

    output_dir.mkdir(parents=True, exist_ok=True)
    bio_data = adata.copy()
    bio_data.X = bio_data.layers["counts"].copy()
    bio_data.obs[cell_type_key] = bio_data.obs[cell_type_key].astype("category")
    params = {
        "batch_key": batch_key,
        "celltype_key": cell_type_key,
        "min_cells": int(min_cells),
        "abundant_gene_threshold": float(abundant_gene_threshold),
        "min_abundant_genes": int(min_abundant_genes),
    }
    analyzer = biolucid.core.BatchEffectAnalyzer(bio_data, params=params)
    analyzer.run_analysis()
    global_df, per_sample_df = biolucid.visualization.results_to_df(analyzer.results)
    global_df.to_csv(output_dir / "global_results.csv", index=False)
    per_sample_df.to_csv(output_dir / "per_sample_results.csv", index=True)

    figure = plot_biolucid_scores(
        per_sample_df,
        output_dir / "biolucid_scatter.png",
    )
    return analyzer, global_df, per_sample_df, figure


def _normalize_to_reference_median(matrix, target_sum: float):
    totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    if np.any(totals <= 0):
        raise ValueError("kNN input contains a zero-count cell")
    scale = target_sum / totals
    if sparse.issparse(matrix):
        return matrix.multiply(scale[:, None]).tocsr().log1p()
    return np.log1p(np.asarray(matrix, dtype=float) * scale[:, None])


def _entropy(probabilities: np.ndarray) -> float:
    nonzero = probabilities > 0
    return float(-np.sum(probabilities[nonzero] * np.log(probabilities[nonzero])))


def transfer_labels_knn(
    ref_matrix,
    query_matrix,
    ref_labels: np.ndarray,
    *,
    n_neighbors: int,
    n_pca: int,
    random_state: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """DISTILL-compatible median-depth/log1p/PCA/kNN label transfer."""
    effective_neighbors = max(1, min(n_neighbors, ref_matrix.shape[0]))
    effective_pca = min(n_pca, ref_matrix.shape[0] - 1, ref_matrix.shape[1] - 1)
    if effective_pca < 1:
        raise ValueError("Not enough reference cells or genes for label transfer")
    target_sum = float(np.median(np.asarray(ref_matrix.sum(axis=1)).ravel()))
    ref_norm = _normalize_to_reference_median(ref_matrix, target_sum)
    query_norm = _normalize_to_reference_median(query_matrix, target_sum)

    pca = PCA(n_components=effective_pca, svd_solver="arpack", random_state=random_state)
    ref_pca = pca.fit_transform(ref_norm)
    query_pca = pca.transform(query_norm)
    neighbors = NearestNeighbors(n_neighbors=effective_neighbors, metric="euclidean")
    neighbors.fit(ref_pca)
    _, indices = neighbors.kneighbors(query_pca)

    labels = np.unique(ref_labels)
    label_to_index = {label: index for index, label in enumerate(labels)}
    ref_codes = np.array([label_to_index[label] for label in ref_labels])
    neighbor_codes = ref_codes[indices]
    votes = np.column_stack(
        [(neighbor_codes == label_index).sum(axis=1) for label_index in range(len(labels))]
    ).astype(float)
    probabilities = votes / votes.sum(axis=1, keepdims=True)
    predicted = labels[probabilities.argmax(axis=1)]
    confidence = probabilities.max(axis=1)

    reference_probabilities = np.array([(ref_labels == label).mean() for label in labels])
    max_entropy = _entropy(reference_probabilities)
    entropy = (
        np.array([_entropy(row) / max_entropy for row in probabilities])
        if max_entropy > 0
        else np.zeros(query_matrix.shape[0])
    )
    return predicted, entropy, confidence


def restore_high_mito_labels(
    full_post_qc: ad.AnnData,
    astrid_annotated: ad.AnnData,
    *,
    cell_type_key: str,
    n_neighbors: int,
    n_pca: int,
    random_state: int,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Restore only QC-flagged high-mito cells; never restore unrelated exclusions."""
    require_columns(full_post_qc.obs, ["sample", "qc_high_mito"], "full_post_qc.obs")
    require_columns(astrid_annotated.obs, ["sample", cell_type_key], "ASTRID obs")

    expected_reference = full_post_qc.obs_names[~full_post_qc.obs["qc_high_mito"]]
    missing_astrid = expected_reference.difference(astrid_annotated.obs_names)
    extra_astrid = astrid_annotated.obs_names.difference(expected_reference)
    if len(missing_astrid) or len(extra_astrid):
        raise AssertionError(
            "ASTRID/reference cells do not match the low-mito QC set: "
            f"missing={len(missing_astrid)}, extra={len(extra_astrid)}"
        )

    combined = full_post_qc.copy()
    combined.X = combined.layers["counts"].copy()
    combined.obs[cell_type_key] = pd.Series(index=combined.obs_names, dtype="object")
    combined.obs["cell_type_source"] = pd.Series(index=combined.obs_names, dtype="object")
    combined.obs["cell_type_transfer_confidence"] = np.nan
    combined.obs["cell_type_transfer_entropy"] = np.nan
    combined.obs.loc[expected_reference, cell_type_key] = (
        astrid_annotated.obs.loc[expected_reference, cell_type_key].astype(str)
    )
    if "cell_type_source" in astrid_annotated.obs:
        combined.obs.loc[expected_reference, "cell_type_source"] = (
            astrid_annotated.obs.loc[expected_reference, "cell_type_source"].astype(str)
        )
    else:
        combined.obs.loc[expected_reference, "cell_type_source"] = (
            "manual_astrid_cluster_review"
        )

    rows: list[dict[str, Any]] = []
    for sample in combined.obs["sample"].astype(str).drop_duplicates():
        in_sample = combined.obs["sample"].astype(str).eq(sample)
        query_mask = in_sample & combined.obs["qc_high_mito"]
        ref_mask = in_sample & ~combined.obs["qc_high_mito"]
        query_names = combined.obs_names[query_mask]
        ref_names = combined.obs_names[ref_mask]
        if len(query_names) == 0:
            rows.append({"sample": sample, "n_reference": len(ref_names), "n_restored": 0})
            continue

        expression = combined.layers["counts"]
        ref_positions = combined.obs_names.get_indexer(ref_names)
        query_positions = combined.obs_names.get_indexer(query_names)
        predicted, entropy, confidence = transfer_labels_knn(
            expression[ref_positions, :],
            expression[query_positions, :],
            combined.obs.loc[ref_names, cell_type_key].astype(str).to_numpy(),
            n_neighbors=n_neighbors,
            n_pca=n_pca,
            random_state=random_state,
        )
        combined.obs.loc[query_names, cell_type_key] = predicted
        combined.obs.loc[query_names, "cell_type_source"] = "high_mito_knn"
        combined.obs.loc[query_names, "cell_type_transfer_confidence"] = confidence
        combined.obs.loc[query_names, "cell_type_transfer_entropy"] = entropy
        rows.append(
            {
                "sample": sample,
                "n_reference": len(ref_names),
                "n_restored": len(query_names),
                "median_confidence": float(np.median(confidence)),
                "median_entropy": float(np.median(entropy)),
            }
        )

    if combined.obs[cell_type_key].isna().any():
        raise AssertionError("Restored object still contains missing cell-type labels")
    return combined, pd.DataFrame(rows)


def prepare_mouse_tidyi(
    adata: ad.AnnData,
    *,
    min_pct_genes_mapped: float = 70.0,
    fail_pct_genes_mapped_lt: float = 30.0,
    allow_low_mapping: bool = False,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Keep valid mouse Ensembl genes and preserve symbols for final restoration."""
    require_columns(adata.var, ["gene_ids", "gene_symbol"], "adata.var")
    gene_ids = adata.var["gene_ids"].astype(str)
    valid = gene_ids.str.startswith("ENSMUSG").to_numpy()
    prepared = adata[:, valid].copy()
    prepared.var["original_var_name"] = prepared.var_names.astype(str)
    prepared.var_names = pd.Index(prepared.var["gene_ids"].astype(str))
    if not prepared.var_names.is_unique:
        raise ValueError("Mouse Ensembl IDs are not unique after TIDYI preparation")
    prepared.X = prepared.layers["counts"].copy()
    prepared.layers["counts"] = prepared.X.copy()
    pct_mapped = 100.0 * prepared.n_vars / max(1, adata.n_vars)
    report = pd.DataFrame(
        [
            {
                "genes_input": int(adata.n_vars),
                "genes_mapped_ensmusg": int(prepared.n_vars),
                "pct_genes_mapped": pct_mapped,
                "genes_removed": int(adata.n_vars - prepared.n_vars),
                "min_pct_genes_mapped": float(min_pct_genes_mapped),
                "fail_pct_genes_mapped_lt": float(fail_pct_genes_mapped_lt),
                "allow_low_mapping": bool(allow_low_mapping),
                "mapping_warning": pct_mapped < min_pct_genes_mapped,
            }
        ]
    )
    if pct_mapped < min_pct_genes_mapped:
        import warnings

        warnings.warn(
            f"Only {pct_mapped:.2f}% of genes have ENSMUSG IDs "
            f"(<{min_pct_genes_mapped}%)",
            stacklevel=2,
        )
    if pct_mapped < fail_pct_genes_mapped_lt and not allow_low_mapping:
        raise RuntimeError(
            f"Only {pct_mapped:.2f}% of genes have ENSMUSG IDs; "
            f"DISTILL failure threshold is {fail_pct_genes_mapped_lt}%"
        )
    return prepared, report


def train_mouse_tidyi(gene_names: Iterable[str], *, random_state: int = 0):
    """Construct and train the mouse TIDYI model (an intentionally explicit step)."""
    from tidyi import tidyi_main_mouse

    np.random.seed(random_state)
    model = tidyi_main_mouse.tidy_rates_mouse(gene_list=list(gene_names))
    model.train_classifier()
    return model


def _positive_score(values) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim > 1 and array.shape[1] > 1:
        return np.asarray(array[:, 1], dtype=float)
    return np.asarray(array, dtype=float).ravel()


def run_tidyi_by_sample_celltype(
    adata: ad.AnnData,
    model,
    *,
    sample_key: str,
    cell_type_key: str,
    min_cells_exclusive: int,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """Apply TIDYI to eligible sample x cell-type groups and retain cell scores."""
    require_columns(adata.obs, [sample_key, cell_type_key, "genotype", "stage"], "TIDYI obs")
    group_counts = (
        adata.obs.groupby(
            [sample_key, "genotype", "stage", cell_type_key], observed=True
        )
        .size()
        .reset_index(name="n_cells")
    )
    group_counts["keep_for_tidyi"] = group_counts["n_cells"] > min_cells_exclusive
    scored = adata.copy()
    scored.obs["division_score"] = np.nan
    scored.obs["death_score"] = np.nan
    rows: list[dict[str, Any]] = []

    # DISTILL-compatible part: every hard-QC cell receives a score. Rates from
    # these six whole-sample calls are not used for the cell-type comparison.
    for sample in scored.obs[sample_key].astype(str).drop_duplicates():
        sample_mask = scored.obs[sample_key].astype(str).eq(sample)
        names = scored.obs_names[sample_mask]
        subset = scored[names].copy()
        _, _, death_prediction, division_prediction = model.apply_classifier(subset)
        scored.obs.loc[names, "death_score"] = _positive_score(death_prediction)
        scored.obs.loc[names, "division_score"] = _positive_score(division_prediction)

    # Earlier AER reporting rule: only n > threshold groups receive a rate row.
    for group in group_counts.loc[group_counts["keep_for_tidyi"]].itertuples(index=False):
        sample = str(getattr(group, sample_key))
        cell_type = str(getattr(group, cell_type_key))
        mask = (
            scored.obs[sample_key].astype(str).eq(sample)
            & scored.obs[cell_type_key].astype(str).eq(cell_type)
        )
        names = scored.obs_names[mask]
        subset = scored[names].copy()
        division_rate, death_rate, _, _ = model.apply_classifier(subset)
        rows.append(
            {
                sample_key: sample,
                "genotype": str(group.genotype),
                "stage": str(group.stage),
                cell_type_key: cell_type,
                "n_cells": int(group.n_cells),
                "proliferation_rate": float(division_rate),
                "death_rate": float(death_rate),
            }
        )
    rates = pd.DataFrame(rows)
    return scored, group_counts, rates


def restore_gene_symbols(adata: ad.AnnData) -> ad.AnnData:
    """Restore final var_names while keeping TIDYI Ensembl IDs in var metadata."""
    out = adata.copy()
    out.var["tidyi_gene_id"] = out.var_names.astype(str)
    out.var_names = pd.Index(out.var["original_var_name"].astype(str))
    out.var_names_make_unique()
    out.var_names.name = None
    return out


def plot_tidyi_rates(rates: pd.DataFrame, output_path: Path):
    """Plot rates in the same form as ``TIDYI_cellrangerAER.ipynb``."""
    import seaborn as sns

    stage_order = ["E9.5", "E10.5", "E11.5"]
    cell_types = sorted(rates["cell_type"].astype(str).unique())
    markers = ["o", "s", "^", "D", "P", "X", "v", "<", ">", "*"]
    if len(cell_types) > len(markers):
        raise ValueError(f"Too many cell types for the marker map: {len(cell_types)}")
    marker_map = dict(zip(cell_types, markers[: len(cell_types)]))
    stage_palette = dict(
        zip(stage_order, sns.color_palette("Set2", n_colors=len(stage_order)))
    )
    x_max = max(float(rates["death_rate"].max()) * 1.1, 1e-6)
    y_max = max(float(rates["proliferation_rate"].max()) * 1.1, 1e-6)

    figure, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
    for axis, genotype in zip(axes, ["WT", "Dac"]):
        subset = rates[rates["genotype"].astype(str).eq(genotype)]
        sns.scatterplot(
            data=subset,
            x="death_rate",
            y="proliferation_rate",
            hue="stage",
            hue_order=stage_order,
            palette=stage_palette,
            style="cell_type",
            style_order=cell_types,
            markers=marker_map,
            s=90,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.4,
            ax=axis,
        )
        axis.plot(
            [0, min(x_max, y_max)],
            [0, min(x_max, y_max)],
            color="gray",
            linestyle="--",
            linewidth=1,
            alpha=0.5,
        )
        axis.set_xlim(0, x_max)
        axis.set_ylim(0, y_max)
        axis.set_title(genotype)
        axis.set_xlabel("Death rate")
        axis.set_ylabel("Proliferation rate")
        axis.grid(alpha=0.25)
        if axis.legend_ is not None:
            axis.legend_.remove()
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, bbox_to_anchor=(1.01, 0.5), loc="center left", frameon=False)
    figure.suptitle(
        "TIDYI: ASTRID-reannotated reference + KNN-predicted cells with "
        "entropy < 0.3; sample × cell type >50",
        y=1.02,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    return figure


def plot_celltype_percentages_by_sample(
    input_path: Path,
    output_dir: Path,
) -> tuple[pd.DataFrame, plt.Figure]:
    """Plot final DISTILL cell-type percentages with the legacy AER settings."""
    stage_order = ["E9.5", "E10.5", "E11.5"]
    genotype_order = ["WT", "Dac"]
    sample_colors = {
        ("E9.5", "WT"): "#c6dbef",
        ("E9.5", "Dac"): "#fcbba1",
        ("E10.5", "WT"): "#6baed6",
        ("E10.5", "Dac"): "#fb6a4a",
        ("E11.5", "WT"): "#2171b5",
        ("E11.5", "Dac"): "#cb181d",
    }

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_path, backed="r")
    required = {"sample", "stage", "genotype", "cell_type"}
    missing = required.difference(adata.obs.columns)
    if missing:
        adata.file.close()
        raise KeyError(f"Missing required obs columns: {sorted(missing)}")
    if adata.obs["cell_type"].isna().any():
        adata.file.close()
        raise AssertionError("Final DISTILL cell_type contains missing values")

    obs = adata.obs[["sample", "stage", "genotype", "cell_type"]].copy()
    adata.file.close()
    sample_meta = obs[["sample", "stage", "genotype"]].astype(str).drop_duplicates()
    if sample_meta["sample"].duplicated().any():
        raise AssertionError("A sample maps to multiple stage/genotype combinations")
    sample_meta["stage"] = pd.Categorical(
        sample_meta["stage"], stage_order, ordered=True
    )
    sample_meta["genotype"] = pd.Categorical(
        sample_meta["genotype"], genotype_order, ordered=True
    )
    sample_meta = sample_meta.sort_values(["stage", "genotype"]).reset_index(drop=True)
    if len(sample_meta) != 6 or sample_meta[["stage", "genotype"]].isna().any().any():
        raise AssertionError("Expected six samples covering three stages x two genotypes")
    sample_meta["display_label"] = (
        sample_meta["stage"].astype(str) + " " + sample_meta["genotype"].astype(str)
    )

    cell_types = sorted(obs["cell_type"].astype(str).unique())
    counts = (
        obs.groupby(["sample", "cell_type"], observed=True)
        .size()
        .reindex(
            pd.MultiIndex.from_product(
                [sample_meta["sample"].astype(str), cell_types],
                names=["sample", "cell_type"],
            ),
            fill_value=0,
        )
        .rename("n_cells")
        .reset_index()
    )
    totals = (
        obs.groupby("sample", observed=True)
        .size()
        .rename("sample_total_cells")
        .reset_index()
    )
    table = counts.merge(totals, on="sample", validate="many_to_one")
    table["percent_cells"] = 100 * table["n_cells"] / table["sample_total_cells"]
    table = table.merge(sample_meta, on="sample", validate="many_to_one")
    percentage_sums = table.groupby("sample", observed=True)["percent_cells"].sum()
    if not np.allclose(percentage_sums.to_numpy(), 100.0):
        raise AssertionError("DISTILL sample percentages do not sum to 100")
    table.to_csv(
        output_dir / "distill_celltype_percentages_by_sample.tsv",
        sep="\t",
        index=False,
    )

    x = np.arange(len(cell_types), dtype=float)
    n_samples = len(sample_meta)
    group_width = 0.84
    bar_width = group_width / n_samples
    figure, axis = plt.subplots(figsize=(15, 6.2), constrained_layout=True)
    for sample_index, sample_row in sample_meta.iterrows():
        sample = str(sample_row["sample"])
        stage = str(sample_row["stage"])
        genotype = str(sample_row["genotype"])
        values = (
            table.loc[table["sample"].astype(str).eq(sample)]
            .set_index("cell_type")
            .reindex(cell_types)["percent_cells"]
            .to_numpy()
        )
        # A true zero has no finite position on a log axis, so its bar is absent.
        plot_values = np.where(values > 0, values, np.nan)
        offset = (sample_index - (n_samples - 1) / 2) * bar_width
        axis.bar(
            x + offset,
            plot_values,
            width=bar_width * 0.94,
            color=sample_colors[(stage, genotype)],
            edgecolor="black",
            linewidth=0.35,
            label=str(sample_row["display_label"]),
        )

    positive = table.loc[table["percent_cells"] > 0, "percent_cells"]
    axis.set_yscale("log", base=10)
    axis.set_ylim(float(positive.min()) * 0.7, float(positive.max()) * 1.35)
    axis.set_xticks(x)
    display_cell_types = [
        "Dorso-ventral ectoderm" if cell_type == "DV" else cell_type
        for cell_type in cell_types
    ]
    axis.set_xticklabels(display_cell_types, rotation=25, ha="right")
    axis.set_xlabel("Cell type")
    axis.set_ylabel("Cells within sample (%) — log10 scale")
    axis.grid(axis="y", which="both", alpha=0.25)
    axis.set_axisbelow(True)
    axis.set_title(
        "Cell-type composition: ASTRID-reannotated + KNN-predicted cells "
        "with entropy < 0.3"
    )
    axis.legend(
        title="Stage / genotype",
        ncols=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0),
        frameon=False,
    )
    figure.savefig(
        output_dir / "distill_entropy_lt_0_3_celltype_percent_log10_barplot.png",
        dpi=300,
        bbox_inches="tight",
    )
    return table, figure


def run_growth_rate_decomposition_method3(
    rates_path: Path,
    typical_cv_path: Path,
    reference_cv_detail_path: Path,
    output_dir: Path,
    contrasts: list[dict[str, str]],
    *,
    cv_dataset: str = "Mukin",
    regulation_magnitude_threshold: float = 2 * np.log(2) / (24 * 3),
) -> tuple[pd.DataFrame, plt.Figure, plt.Figure, plt.Figure]:
    """Apply the legacy Method 3 decomposition to AER TIDYI rates."""
    import seaborn as sns
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    from matplotlib.ticker import PercentFormatter

    rates = pd.read_csv(rates_path)
    require_columns(
        rates,
        [
            "sample",
            "genotype",
            "stage",
            "cell_type",
            "n_cells",
            "proliferation_rate",
            "death_rate",
        ],
        "AER TIDYI rate table",
    )
    if rates.duplicated(["genotype", "stage", "cell_type"]).any():
        raise ValueError("AER TIDYI rates contain duplicate genotype/stage/cell-type rows")

    cv_values = pd.read_csv(typical_cv_path)
    reference_cv = pd.read_csv(reference_cv_detail_path)
    require_columns(
        cv_values,
        ["dataset", "cv_source_dataset", "rate", "cv_rate", "typical_cv"],
        "typical CV table",
    )
    require_columns(
        reference_cv,
        ["dataset", "rate", "mean_rate", "cv_status"],
        "reference CV detail table",
    )

    cv_subset = cv_values[cv_values["dataset"].eq(cv_dataset)].copy()
    if set(cv_subset["rate"]) != {"prolif_rate", "death_rate"}:
        raise ValueError(f"Expected proliferation and death CV rows for {cv_dataset}")
    cv_subset = cv_subset.set_index("rate")

    reference_global_means: dict[str, float] = {}
    source_datasets: dict[str, str] = {}
    source_rates: dict[str, str] = {}
    for rate in ["prolif_rate", "death_rate"]:
        source_dataset = str(cv_subset.at[rate, "cv_source_dataset"])
        source_rate = str(cv_subset.at[rate, "cv_rate"])
        source_rows = reference_cv[
            reference_cv["dataset"].eq(source_dataset)
            & reference_cv["rate"].eq(source_rate)
            & reference_cv["cv_status"].eq("computed")
        ]
        if source_rows.empty:
            raise ValueError(
                f"No computed reference rows for {source_dataset} {source_rate}"
            )
        reference_global_means[rate] = float(source_rows["mean_rate"].mean())
        source_datasets[rate] = source_dataset
        source_rates[rate] = source_rate

    proliferation_cv = float(cv_subset.at["prolif_rate", "typical_cv"])
    death_cv = float(cv_subset.at["death_rate", "typical_cv"])
    aer_global_means = {
        "prolif_rate": float(rates["proliferation_rate"].mean()),
        "death_rate": float(rates["death_rate"].mean()),
    }
    method3_k = (
        death_cv * aer_global_means["death_rate"]
        / (proliferation_cv * aer_global_means["prolif_rate"])
    )
    method3_k_mukin_reference = (
        death_cv * reference_global_means["death_rate"]
        / (proliferation_cv * reference_global_means["prolif_rate"])
    )

    required_contrast_keys = {
        "contrast",
        "genotype_a",
        "stage_a",
        "genotype_b",
        "stage_b",
    }
    for contrast in contrasts:
        missing = required_contrast_keys.difference(contrast)
        if missing:
            raise KeyError(f"Contrast is missing keys: {sorted(missing)}")

    cell_types = sorted(rates["cell_type"].astype(str).unique())
    rows: list[dict[str, Any]] = []
    for contrast_index, contrast in enumerate(contrasts):
        condition_a = rates[
            rates["genotype"].astype(str).eq(contrast["genotype_a"])
            & rates["stage"].astype(str).eq(contrast["stage_a"])
        ].set_index("cell_type")
        condition_b = rates[
            rates["genotype"].astype(str).eq(contrast["genotype_b"])
            & rates["stage"].astype(str).eq(contrast["stage_b"])
        ].set_index("cell_type")

        for cell_type in cell_types:
            has_a = cell_type in condition_a.index
            has_b = cell_type in condition_b.index
            if has_a and has_b:
                proliferation_a = float(
                    condition_a.at[cell_type, "proliferation_rate"]
                )
                proliferation_b = float(
                    condition_b.at[cell_type, "proliferation_rate"]
                )
                death_a = float(condition_a.at[cell_type, "death_rate"])
                death_b = float(condition_b.at[cell_type, "death_rate"])
                delta_proliferation = proliferation_a - proliferation_b
                delta_death = death_a - death_b
                denominator = abs(delta_death) + method3_k * abs(delta_proliferation)
                death_fraction = abs(delta_death) / denominator if denominator > 0 else np.nan
                regulation_magnitude = abs(delta_proliferation) + abs(delta_death)
            else:
                proliferation_a = np.nan
                proliferation_b = np.nan
                death_a = np.nan
                death_b = np.nan
                delta_proliferation = np.nan
                delta_death = np.nan
                death_fraction = np.nan
                regulation_magnitude = np.nan

            passes_filter = bool(
                pd.notna(death_fraction)
                and regulation_magnitude >= regulation_magnitude_threshold
            )
            missing_conditions = []
            if not has_a:
                missing_conditions.append("condition_a")
            if not has_b:
                missing_conditions.append("condition_b")
            rows.append(
                {
                    "contrast_order": contrast_index,
                    "contrast": contrast["contrast"],
                    "condition_a": f'{contrast["genotype_a"]} {contrast["stage_a"]}',
                    "condition_b": f'{contrast["genotype_b"]} {contrast["stage_b"]}',
                    "cell_type": cell_type,
                    "n_cells_a": (
                        int(condition_a.at[cell_type, "n_cells"]) if has_a else np.nan
                    ),
                    "n_cells_b": (
                        int(condition_b.at[cell_type, "n_cells"]) if has_b else np.nan
                    ),
                    "proliferation_rate_a": proliferation_a,
                    "proliferation_rate_b": proliferation_b,
                    "death_rate_a": death_a,
                    "death_rate_b": death_b,
                    "delta_proliferation": delta_proliferation,
                    "delta_death": delta_death,
                    "abs_delta_proliferation": abs(delta_proliferation),
                    "abs_delta_death": abs(delta_death),
                    "regulation_magnitude": regulation_magnitude,
                    "regulation_threshold": regulation_magnitude_threshold,
                    "passes_magnitude_filter": passes_filter,
                    "method3_k": method3_k,
                    "method3_k_mukin_reference_sensitivity": (
                        method3_k_mukin_reference
                    ),
                    "death_fraction_m3": death_fraction,
                    "death_fraction_m3_magnitude_filtered": (
                        death_fraction if passes_filter else np.nan
                    ),
                    "available_method3_fraction": bool(pd.notna(death_fraction)),
                    "unavailable_reason": (
                        ";".join(missing_conditions) if missing_conditions else pd.NA
                    ),
                    "cv_dataset": cv_dataset,
                    "proliferation_typical_cv": proliferation_cv,
                    "death_typical_cv": death_cv,
                    "global_mean_rate_source": (
                        "AER sample_celltype_rates unweighted"
                    ),
                    "n_aer_rate_rows_for_global_mean": len(rates),
                    "proliferation_global_mean_used": aer_global_means[
                        "prolif_rate"
                    ],
                    "death_global_mean_used": aer_global_means["death_rate"],
                    "proliferation_global_mean_reference": reference_global_means[
                        "prolif_rate"
                    ],
                    "death_global_mean_reference": reference_global_means[
                        "death_rate"
                    ],
                    "proliferation_cv_source_dataset": source_datasets[
                        "prolif_rate"
                    ],
                    "death_cv_source_dataset": source_datasets["death_rate"],
                    "proliferation_cv_source_rate": source_rates["prolif_rate"],
                    "death_cv_source_rate": source_rates["death_rate"],
                }
            )

    results = pd.DataFrame(rows)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(
        output_dir
        / "aer_method3_aer_global_mean_regulation_magnitude_filter.csv",
        index=False,
    )

    contrast_order = [contrast["contrast"] for contrast in contrasts]
    original_heat = results.pivot_table(
        index="cell_type",
        columns="contrast",
        values="death_fraction_m3",
        aggfunc="first",
    ).reindex(index=cell_types, columns=contrast_order)
    filtered_heat = results.pivot_table(
        index="cell_type",
        columns="contrast",
        values="death_fraction_m3_magnitude_filtered",
        aggfunc="first",
    ).reindex(index=cell_types, columns=contrast_order)
    availability = original_heat.notna().astype(float)

    figure_height = max(5.2, 0.22 * len(cell_types) + 1.5)
    figure_width = max(6.2, 1.25 * len(contrast_order) + 3.0)
    background_cmap = ListedColormap(["#D9D9D9", "#FFFFFF"])

    unfiltered_figure, unfiltered_axis = plt.subplots(
        figsize=(figure_width, figure_height), constrained_layout=True
    )
    sns.heatmap(
        availability,
        ax=unfiltered_axis,
        cmap=background_cmap,
        vmin=0,
        vmax=1,
        cbar=False,
        linewidths=0.35,
        linecolor="white",
    )
    unfiltered_annotations = original_heat.map(
        lambda value: "" if pd.isna(value) else f"{value:.2f}"
    )
    sns.heatmap(
        original_heat,
        mask=original_heat.isna(),
        ax=unfiltered_axis,
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
        linewidths=0,
        annot=unfiltered_annotations,
        fmt="",
        annot_kws={"fontsize": 5.8},
        cbar_kws={"label": "Death contribution fraction"},
    )
    unfiltered_axis.set_title(
        "AER Method 3: Mukin CV + AER global mean rates\n"
        "Before regulation-magnitude filtering",
        fontsize=10,
        pad=28,
    )
    unfiltered_axis.legend(
        handles=[
            Patch(
                facecolor="#D9D9D9",
                edgecolor="#BDBDBD",
                label="Fraction unavailable",
            )
        ],
        loc="lower left",
        bbox_to_anchor=(0, 1.0),
        frameon=False,
        fontsize=7,
    )
    unfiltered_axis.set_xlabel("Contrast: condition A minus condition B", fontsize=8)
    unfiltered_axis.set_ylabel("Cell type", fontsize=8)
    unfiltered_axis.tick_params(axis="x", labelrotation=35, labelsize=8)
    unfiltered_axis.tick_params(axis="y", labelsize=6)
    unfiltered_figure.savefig(
        output_dir / "aer_method3_aer_global_mean_unfiltered_heatmap.png",
        bbox_inches="tight",
        dpi=300,
    )

    filtered_figure, filtered_axis = plt.subplots(
        figsize=(figure_width, figure_height), constrained_layout=True
    )
    sns.heatmap(
        availability,
        ax=filtered_axis,
        cmap=background_cmap,
        vmin=0,
        vmax=1,
        cbar=False,
        linewidths=0.35,
        linecolor="white",
    )
    filtered_annotations = filtered_heat.map(
        lambda value: "" if pd.isna(value) else f"{value:.2f}"
    )
    sns.heatmap(
        filtered_heat,
        mask=filtered_heat.isna(),
        ax=filtered_axis,
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
        linewidths=0,
        annot=filtered_annotations,
        fmt="",
        annot_kws={"fontsize": 5.8},
        cbar_kws={"label": "Death contribution fraction"},
    )
    filtered_axis.set_title(
        "AER Method 3: Mukin CV + AER global mean rates\n"
        "After regulation-magnitude filtering; "
        f"Colored cells satisfy |Δβ| + |Δα| ≥ "
        f"{regulation_magnitude_threshold:.4f} h⁻¹",
        fontsize=10,
        pad=28,
    )
    filtered_axis.legend(
        handles=[
            Patch(
                facecolor="white",
                edgecolor="#BDBDBD",
                label=(
                    f"Filtered: |Δβ| + |Δα| < "
                    f"{regulation_magnitude_threshold:.4f} h⁻¹"
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
    filtered_axis.set_xlabel("Contrast: condition A minus condition B", fontsize=8)
    filtered_axis.set_ylabel("Cell type", fontsize=8)
    filtered_axis.tick_params(axis="x", labelrotation=35, labelsize=8)
    filtered_axis.tick_params(axis="y", labelsize=6)
    filtered_figure.savefig(
        output_dir
        / "aer_method3_aer_global_mean_magnitude_filtered_heatmap.png",
        bbox_inches="tight",
        dpi=300,
    )

    all_fractions = results["death_fraction_m3"].dropna()
    filtered_fractions = results[
        "death_fraction_m3_magnitude_filtered"
    ].dropna()
    histogram_rows = [
        ("Before filter", all_fractions, "#9ECAE1"),
        ("After filter", filtered_fractions, "#4C78A8"),
    ]
    histogram_figure, histogram_axes = plt.subplots(
        2,
        1,
        figsize=(3.6, 5.4),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    for row_index, (filter_label, values, color) in enumerate(histogram_rows):
        axis = histogram_axes[row_index]
        bins = np.linspace(0, 1, 6)
        axis.hist(
            values,
            bins=bins,
            weights=np.ones(len(values)) / len(values),
            color=color,
            edgecolor="white",
            linewidth=0.5,
        )
        axis.set_title(
            f"{filter_label}; 5 bins (n={len(values)})",
            fontsize=9,
        )
        axis.set_xlim(0, 1)
        axis.set_xticks(np.linspace(0, 1, 6))
        axis.tick_params(axis="both", labelsize=7)
        axis.yaxis.set_major_formatter(PercentFormatter(1.0))
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_ylabel("Proportion of entries", fontsize=8)
        if row_index == 1:
            sns.rugplot(
                x=values,
                ax=axis,
                height=0.05,
                color="black",
                linewidth=0.7,
            )
    histogram_axes[1].set_xlabel(
        "Death contribution fraction (Method 3)", fontsize=8
    )
    histogram_figure.suptitle(
        "AER Method 3: Mukin CV + AER global mean rates\n"
        "Before and after regulation-magnitude filtering",
        fontsize=10,
    )
    histogram_figure.savefig(
        output_dir
        / "aer_method3_aer_global_mean_magnitude_filter_histogram_5bins.png",
        bbox_inches="tight",
        dpi=300,
    )
    return results, unfiltered_figure, filtered_figure, histogram_figure


def plot_sample_mean_umi(
    adata_path: Path,
    sample_col: str,
    output_path: Path,
    title: str,
    *,
    include_samples: Iterable[str] | None = None,
    color: str = "#4C78A8",
) -> tuple[pd.DataFrame, plt.Figure]:
    """Plot mean raw UMI counts per retained cell for every selected sample."""
    from matplotlib.ticker import FuncFormatter

    adata = ad.read_h5ad(adata_path, backed="r")
    try:
        require_columns(adata.obs, [sample_col, "total_counts"], str(adata_path))
        obs = adata.obs[[sample_col, "total_counts"]].copy()
    finally:
        adata.file.close()

    obs[sample_col] = obs[sample_col].astype(str)
    if include_samples is None:
        sample_order = sorted(obs[sample_col].unique())
    else:
        sample_order = list(dict.fromkeys(str(sample) for sample in include_samples))
        missing_samples = sorted(set(sample_order).difference(obs[sample_col]))
        if missing_samples:
            raise ValueError(
                f"{adata_path} is missing selected samples: {missing_samples}"
            )
        obs = obs[obs[sample_col].isin(sample_order)].copy()

    if obs["total_counts"].isna().any() or (obs["total_counts"] < 0).any():
        raise ValueError(f"Invalid total_counts values in {adata_path}")

    overall_mean_umi = float(obs["total_counts"].mean())
    summary = (
        obs.groupby(sample_col, observed=True)["total_counts"]
        .agg(n_cells="size", mean_umi="mean")
        .reindex(sample_order)
        .reset_index()
        .rename(columns={sample_col: "sample"})
    )

    figure_width = min(20.0, max(6.0, 2.6 + 0.24 * len(summary)))
    figure, axis = plt.subplots(
        figsize=(figure_width, 5.2), constrained_layout=True
    )
    x_positions = np.arange(len(summary))
    bars = axis.bar(
        x_positions,
        summary["mean_umi"],
        color=color,
        edgecolor="white",
        linewidth=0.4,
    )
    axis.axhline(
        overall_mean_umi,
        color="#D62728",
        linestyle="--",
        linewidth=1.0,
        label=f"Overall mean = {overall_mean_umi:,.0f}",
    )
    axis.set_title(title, fontsize=10)
    axis.set_xlabel("Sample", fontsize=8)
    axis.set_ylabel("Mean UMI counts per cell", fontsize=8)
    axis.set_xticks(x_positions)
    axis.set_xticklabels(
        summary["sample"], rotation=90, ha="center", fontsize=6
    )
    axis.tick_params(axis="y", labelsize=7)
    axis.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value:,.0f}")
    )
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(loc="upper right", frameon=False, fontsize=7)
    if len(summary) <= 12:
        axis.bar_label(
            bars,
            labels=[f"{value:,.0f}" for value in summary["mean_umi"]],
            padding=2,
            fontsize=7,
        )
        axis.margins(y=0.12)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    return summary, figure


def final_audit(
    raw: ad.AnnData,
    full_post_qc: ad.AnnData,
    for_astrid: ad.AnnData,
    restored: ad.AnnData,
    final: ad.AnnData,
) -> pd.DataFrame:
    """Create a compact, per-sample reconciliation table."""
    rows = []
    samples = raw.obs["sample"].astype(str).drop_duplicates()
    for sample in samples:
        rows.append(
            {
                "sample": sample,
                "input_cells": int(raw.obs["sample"].astype(str).eq(sample).sum()),
                "hard_qc_cells": int(full_post_qc.obs["sample"].astype(str).eq(sample).sum()),
                "astrid_cells": int(for_astrid.obs["sample"].astype(str).eq(sample).sum()),
                "restored_cells": int(restored.obs["sample"].astype(str).eq(sample).sum()),
                "final_cells": int(final.obs["sample"].astype(str).eq(sample).sum()),
                "final_genes": int(final.n_vars),
            }
        )
    return pd.DataFrame(rows)

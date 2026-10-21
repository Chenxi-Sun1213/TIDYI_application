#!/usr/bin/env python3
"""Generate ASTRID abundance and versioned TIDYI rates; never run model fitting."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
RA_ROOT = ROOT.parent
THESIS_ROOT = Path("/scratch/chenxi.sun/thesis_project")
LEGACY_PIPELINE = THESIS_ROOT / "TIDYI/TIDYI_pattern/ty_pipeline_scripts"
COMPARISONS = (
    "mukin_young_spf_vs_old_spf", "mukin_young_spf_vs_young_gf",
    "mukin_old_spf_vs_old_gf", "tms_young_male_vs_young_female",
    "tms_old_male_vs_old_female",
)


def file_record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "size": path.stat().st_size, "sha256": digest.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/input_preparation")
    parser.add_argument("--tidyi-source-root", type=Path)
    parser.add_argument("--tidyi-version-label", default="installed-baseline")
    args = parser.parse_args()
    output = args.output_dir.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Use an empty output directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    for name, child in [("NUMBA_CACHE_DIR", "numba"), ("MPLCONFIGDIR", "matplotlib"), ("IPYTHONDIR", "ipython")]:
        os.environ.setdefault(name, str(output / ".runtime_cache" / child))
    os.environ.setdefault("MPLBACKEND", "Agg")
    # The historical plotting helper also provides the rate-table derived metrics.
    for directory in (LEGACY_PIPELINE, RA_ROOT / "03_death_proliferation/support", ROOT / "support", RA_ROOT):
        sys.path.insert(0, str(directory))
    if args.tidyi_source_root:
        sys.path.insert(0, str(args.tidyi_source_root.expanduser().resolve()))

    import logging
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import tidyi
    from tidyi import tidyi_main_mouse
    from tidyi.data_utils import check_data_exists
    from gene_mapping import GENE_MAPPING_FILE
    import tidyi_astrid_comparison_analysis as pipeline
    import tidyi_comparison_analysis as comparison_helper
    import tidyi_plot_results as plot_helper
    from system_biology_inputs import attach_astrid_abundance, pool_astrid_abundance

    pipeline.setup_logging("INFO")
    logger = logging.getLogger("prepare_systems_biology_inputs")
    if args.tidyi_source_root:
        Path(tidyi.__file__).resolve().relative_to(args.tidyi_source_root.resolve())
    if not check_data_exists():
        raise FileNotFoundError("Local TIDYI classifier data are required; no download is performed")
    manifest = {
        "started_utc": datetime.now(timezone.utc).isoformat(), "status": "running",
        "python": sys.executable, "tidyi_version_label": args.tidyi_version_label,
        "installed_distribution_version": importlib.metadata.version("tidyi"),
        "tidyi_module": file_record(tidyi_main_mouse.__file__),
        "gene_mapping": file_record(GENE_MAPPING_FILE),
        "source_files": [file_record(p) for p in (
            __file__, ROOT / "system_biology.ipynb", pipeline.__file__,
            comparison_helper.__file__, plot_helper.__file__, ROOT / "support/system_biology_inputs.py",
        )],
        "comparisons": list(COMPARISONS), "model_fitting_executed": False,
    }
    manifest_path = output / "input_preparation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    roots = {"mukin_bm": THESIS_ROOT / "mukin_qcdata/AllRef_BM_0204",
             "mukin_skin": THESIS_ROOT / "mukin_qcdata/AllRef_skin_0204",
             "tms_bm": THESIS_ROOT / "tms_qcdata/AllRef_BM_0204"}
    input_paths = [p for root in roots.values() for p in pipeline.discover_astrid_h5ad_files(root)]
    input_paths.append(roots["mukin_bm"] / "BM_combined_ASTRID_results.csv")
    manifest["astrid_inputs"] = [file_record(p) for p in input_paths]
    # Preserve the historical Mukin BM+skin gene universe during TIDYI training.
    datasets = {
        "mukin": pipeline.prepare_mukin_dataset(roots["mukin_bm"], roots["mukin_skin"]),
        "tms": pipeline.prepare_tms_dataset(roots["tms_bm"]),
    }
    manifest["prepared_dataset_shapes"] = {k: list(v.shape) for k, v in datasets.items()}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    notebook = json.loads((ROOT / "system_biology.ipynb").read_text())
    abundance_source = "".join(next(c for c in notebook["cells"] if c.get("id") == "5250c396")["source"])
    abundance_source = abundance_source[abundance_source.index("mukin_obs = "):]
    abundance_source = abundance_source.replace(str(ROOT / "outputs"), str(output))
    abundance_namespace = {
        "mukin_astrid": datasets["mukin"][datasets["mukin"].obs["tissue"].eq("Bone-marrow")],
        "tms_astrid": datasets["tms"], "Path": Path, "pd": pd,
        "pool_astrid_abundance": pool_astrid_abundance,
    }
    exec(compile(abundance_source, "notebook-abundance-only", "exec"), abundance_namespace)
    pooled = abundance_namespace["astrid_condition_abundance"]
    logger.info("Generated fresh ASTRID abundance; starting TIDYI inference only")
    datasets = pipeline.preconvert_dataset_gene_ids(datasets)
    model_inputs = []
    selected_groups = {
        COMPARISONS[0]: ["young_spf", "old_spf"], COMPARISONS[1]: ["young_gf"],
        COMPARISONS[2]: ["old_gf"], COMPARISONS[3]: ["young_male", "young_female"],
        COMPARISONS[4]: ["old_male", "old_female"],
    }
    for comparison in pipeline.get_selected_comparisons(list(COMPARISONS)):
        logger.info("Recomputing %s", comparison.name)
        celltype_col = "SingleR_CellType" if comparison.dataset == "mukin" else "celltype_hb"
        subset = pipeline.build_comparison_subset(datasets[comparison.dataset], comparison, celltype_col)
        subset = pipeline.add_percentage(subset, celltype_col)
        rates = pipeline.run_tidyi_analysis(subset, pipeline.ANALYSIS_TISSUE_COL, celltype_col, raw_layer="raw")
        generated = pipeline.build_final_tidyi_table(
            subset, pipeline.expand_analysis_sample_column(rates),
            comparison.dataset, comparison.name, celltype_col,
        )
        pipeline.save_tidyi_table(generated, output / "ty_df" / f"{comparison.name}_tidyi.csv")
        plt.close("all")
        selected = generated[generated.tissue.isin(["Bone-marrow", "Marrow"]) & generated.comparison_group.isin(selected_groups[comparison.name])].copy()
        selected["dataset"] = "Mukin" if comparison.dataset == "mukin" else "TMS"
        model_inputs.append(attach_astrid_abundance(selected, pooled))
    model_input = pd.concat(model_inputs, ignore_index=True)
    model_input["celltype_key"] = model_input.celltype.astype("string").str.strip().str.lower()
    model_input["gamma"] = model_input.prolif_rate - model_input.death_rate
    model_input["abundance_fraction"] = model_input.percentage / 100
    model_input = model_input[model_input.celltype_key.isin(["b pre", "dendritic", "monocyte pro", "myeloid progenitor", "myelomonocyte", "t cd8"])]
    assert len(model_input) == 48 and np.isfinite(model_input[["gamma", "abundance_fraction"]]).all().all()
    model_input.to_csv(output / "six_celltype_model_inputs.csv", index=False)
    manifest.update(status="complete", completed_utc=datetime.now(timezone.utc).isoformat(),
                    six_celltype_input_rows=len(model_input), outputs=[file_record(p) for p in sorted(output.rglob("*.csv"))])
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    logger.info("Input preparation complete. No systems-biology model fitting executed.")


if __name__ == "__main__":
    main()

"""Keep observed ASTRID abundance separate from versioned TIDYI rates."""

import numpy as np
import pandas as pd


CONDITION_GROUPS = {
    "Mukin SPF young": "young_spf", "Mukin SPF old": "old_spf",
    "Mukin GF young": "young_gf", "Mukin GF old": "old_gf",
    "TMS male young": "young_male", "TMS female young": "young_female",
    "TMS male old": "old_male", "TMS female old": "old_female",
}


def pool_astrid_abundance(mouse_abundance):
    """Pool cell counts across mice, retaining the full filtered-cell denominator."""
    groups = ["dataset", "tissue", "condition"]
    mouse_keys = groups + ["mouse_id"]
    if mouse_abundance.duplicated(mouse_keys + ["cell_type"]).any():
        raise ValueError("Duplicate mouse/cell-type abundance rows")
    mouse_totals = mouse_abundance[mouse_keys + ["n_total"]].drop_duplicates()
    if mouse_totals.duplicated(mouse_keys).any():
        raise ValueError("Conflicting cell-count denominators for a mouse")
    totals = mouse_totals.groupby(groups, observed=True, as_index=False)["n_total"].sum()
    counts = mouse_abundance.groupby(
        groups + ["cell_type"], observed=True, as_index=False,
    )["n_celltype"].sum()
    pooled = counts.merge(totals, on=groups, validate="many_to_one")
    pooled = pooled.rename(columns={"cell_type": "celltype"})
    pooled["comparison_group"] = pooled["condition"].map(CONDITION_GROUPS)
    if pooled["comparison_group"].isna().any():
        raise ValueError("ASTRID abundance contains an unsupported condition")
    if (pooled["n_total"] <= 0).any() or (pooled["n_celltype"] > pooled["n_total"]).any():
        raise ValueError("Invalid ASTRID cell counts")
    pooled["percentage"] = 100 * pooled["n_celltype"] / pooled["n_total"]
    pooled["abundance_fraction"] = pooled["n_celltype"] / pooled["n_total"]
    return pooled.sort_values(groups + ["celltype"]).reset_index(drop=True)


def attach_astrid_abundance(rates, pooled_abundance):
    """Join selected TIDYI rates to observed abundance without using CSV x values."""
    rates = rates.copy()
    pooled = pooled_abundance.copy()
    for frame in (rates, pooled):
        frame["_tissue_key"] = frame["tissue"].replace({"Bone-marrow": "Marrow"})
        frame["_celltype_key"] = frame["celltype"].astype("string").str.strip().str.lower()
    keys = ["dataset", "comparison_group", "_tissue_key", "_celltype_key"]
    if "percentage" in rates:
        rates = rates.rename(columns={"percentage": "percentage_tidyi_csv"})
    pooled = pooled[keys + ["n_celltype", "n_total", "percentage"]]
    merged = rates.merge(pooled, on=keys, how="left", validate="many_to_one", indicator=True)
    if merged["_merge"].ne("both").any():
        missing = merged.loc[merged["_merge"].ne("both"), keys].to_dict("records")
        raise ValueError(f"TIDYI rate rows have no matching ASTRID abundance: {missing}")
    if "percentage_tidyi_csv" in merged and not np.allclose(
        merged["percentage_tidyi_csv"], merged["percentage"], rtol=0, atol=1e-8,
    ):
        raise ValueError(
            "TIDYI CSV percentages disagree with the ASTRID counts. "
            "Recompute rates using the same ASTRID cells, labels and condition groups."
        )
    return merged.drop(columns=["_tissue_key", "_celltype_key", "_merge"])

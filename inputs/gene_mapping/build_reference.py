"""Rebuild the frozen reference from existing local sources, without networking."""

from collections import defaultdict
import hashlib
import json
from pathlib import Path

import h5py
import pandas as pd


OUTPUT_DIR = Path(__file__).resolve().parent
LEGACY_CACHE = Path(
    "/scratch/chenxi.sun/TIDYI_thesis_code_organized/data/reference/"
    "gene_cache/mouse_symbol_to_ensembl.csv"
)
AER_H5AD = Path(
    "/srv/mfs/hausserlab/data/AER2025_Juliane/aer_six_samples_raw_counts.h5ad"
)
BIOMART = OUTPUT_DIR.parent / "classifier_data/mart_export_gene_names.txt"


def read_h5_strings(node):
    if isinstance(node, h5py.Group):
        categories = node["categories"].asstr()[:]
        return [categories[code] if code >= 0 else "" for code in node["codes"][:]]
    return node.asstr()[:]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    legacy = pd.read_csv(LEGACY_CACHE, keep_default_na=False).rename(
        columns={"ensembl_gene_clean": "ensembl_gene"}
    )
    with h5py.File(AER_H5AD, "r") as handle:
        aer = pd.DataFrame({
            "symbol": read_h5_strings(handle["var/gene_symbol"]),
            "ensembl_gene": read_h5_strings(handle["var/gene_ids"]),
        })
    biomart = pd.read_csv(
        BIOMART, sep="\t", keep_default_na=False,
        usecols=["Mouse gene name", "Mouse gene stable ID"],
    ).rename(columns={"Mouse gene name": "symbol", "Mouse gene stable ID": "ensembl_gene"})
    sources = {"legacy_cache": legacy, "aer_raw_features": aer, "biomart_orthologues": biomart}
    candidates = defaultdict(lambda: defaultdict(set))
    for source, frame in sources.items():
        for symbol, gene_id in frame[["symbol", "ensembl_gene"]].itertuples(index=False, name=None):
            symbol, gene_id = symbol.strip(), gene_id.strip().split(".")[0]
            if symbol:
                candidates[symbol]  # Keep documented legacy misses as unresolved rows.
                if gene_id.startswith("ENSMUSG"):
                    candidates[symbol][source].add(gene_id)

    rows = []
    for symbol, by_source in sorted(candidates.items()):
        all_ids = sorted(set().union(*by_source.values())) if by_source else []
        legacy_ids = by_source.get("legacy_cache", set())
        selected = ""
        if len(legacy_ids) == 1:
            selected = next(iter(legacy_ids))
            status = "resolved_legacy_cache"
        elif len(all_ids) == 1:
            selected = all_ids[0]
            status = "resolved_local_supplement"
        else:
            status = "ambiguous" if all_ids else "unmapped"
        rows.append({
            "symbol": symbol, "ensembl_gene": selected, "status": status,
            "sources": ";".join(by_source), "candidate_ensembl_ids": ";".join(all_ids),
        })
    table = pd.DataFrame(rows)
    destination = OUTPUT_DIR / "mouse_symbol_to_ensembl.csv"
    table.to_csv(destination, index=False)
    paths = {"legacy_cache": LEGACY_CACHE, "aer_raw_features": AER_H5AD, "biomart_orthologues": BIOMART}
    provenance = {
        "species": "Mus musculus", "taxon_id": 10090,
        "policy": "Unique legacy-cache result wins; otherwise only a single candidate ID is resolved.",
        "annotation_release": "Not recorded in the source files; this is a frozen mixed-source reference.",
        "sources": {name: {"path": str(path), "sha256": sha256(path)} for name, path in paths.items()},
        "mapping_sha256": sha256(destination),
        "symbol_count": len(table), "status_counts": table.status.value_counts().to_dict(),
        "symbols_with_multiple_candidate_ids": int(table.candidate_ensembl_ids.str.contains(";").sum()),
    }
    (OUTPUT_DIR / "mapping_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({key: provenance[key] for key in ("symbol_count", "status_counts")}, indent=2))


if __name__ == "__main__":
    main()

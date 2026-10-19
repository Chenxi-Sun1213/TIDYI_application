"""Apply the frozen replay3 mouse cache without supplementation or networking."""

from __future__ import annotations

import csv
from functools import lru_cache
import hashlib
from pathlib import Path
import re


GENE_MAPPING_FILE = (
    Path(__file__).resolve().parent / "inputs/gene_mapping/mouse_symbol_to_ensembl.csv"
)
MOUSE_ENSEMBL_PATTERN = re.compile(r"^ENSMUSG\d+(?:\.\d+)?$")


@lru_cache(maxsize=None)
def load_mouse_gene_mapping(mapping_path: Path = GENE_MAPPING_FILE) -> dict[str, str]:
    """Read valid mouse IDs; empty or non-Ensembl archived values stay unresolved."""
    mapping: dict[str, str] = {}
    seen: set[str] = set()
    with Path(mapping_path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"symbol", "ensembl_gene_clean"}.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing symbol/ensembl_gene_clean columns: {mapping_path}")
        for row in reader:
            symbol = row["symbol"].strip()
            gene_id = row["ensembl_gene_clean"].strip()
            if not symbol or symbol in seen:
                raise ValueError(f"Empty or duplicate symbol in {mapping_path}: {symbol!r}")
            seen.add(symbol)
            if MOUSE_ENSEMBL_PATTERN.fullmatch(gene_id):
                mapping[symbol] = gene_id.split(".")[0]
    return mapping


def map_mouse_symbols_to_ensembl(adata, mapping_path: Path = GENE_MAPPING_FILE) -> None:
    """Rename in place, retaining unmapped names and uniquifying duplicate IDs."""
    mapping_path = Path(mapping_path).resolve()
    mapping = load_mouse_gene_mapping(mapping_path)
    original = list(map(str, adata.var_names))
    mapped = [
        name if MOUSE_ENSEMBL_PATTERN.fullmatch(name)
        else mapping.get(name, name)
        for name in original
    ]
    unresolved_set = {
        name for name in original
        if not MOUSE_ENSEMBL_PATTERN.fullmatch(name) and name not in mapping
    }
    adata.var_names = mapped
    if adata.var_names.duplicated().any():
        adata.var_names_make_unique()
    report = {
        "mapping_file": str(mapping_path),
        "mapping_sha256": hashlib.sha256(mapping_path.read_bytes()).hexdigest(),
        "genes_input": len(original),
        "symbols_mapped": sum(name in mapping for name in original),
        "ensembl_ids_retained": sum(bool(MOUSE_ENSEMBL_PATTERN.fullmatch(name)) for name in original),
        "unmapped_symbols": sorted(unresolved_set),
        "unmapped_gene_count": sum(name in unresolved_set for name in original),
    }
    adata.uns["local_gene_mapping"] = report
    print(
        f"Local gene mapping: {report['symbols_mapped']} symbols mapped, "
        f"{report['ensembl_ids_retained']} Ensembl IDs retained, "
        f"{report['unmapped_gene_count']} names unresolved; {mapping_path}"
    )

"""Freeze the replay3 archived mouse mapping cache without supplementation."""
import hashlib
import json
from pathlib import Path
import shutil

import pandas as pd

OUTPUT_DIR = Path(__file__).resolve().parent
LEGACY_CACHE = Path('/scratch/chenxi.sun/.tidyi_gene_cache/mouse_symbol_to_ensembl.csv')


def main():
    table = pd.read_csv(LEGACY_CACHE, keep_default_na=False)
    assert table.symbol.is_unique
    valid = table.ensembl_gene_clean.str.fullmatch(r'ENSMUSG\d+(?:\.\d+)?')
    destination = OUTPUT_DIR / 'mouse_symbol_to_ensembl.csv'
    shutil.copyfile(LEGACY_CACHE, destination)
    checksum = hashlib.sha256(destination.read_bytes()).hexdigest()
    provenance = {
        'species': 'Mus musculus', 'taxon_id': 10090,
        'policy': 'Replay3: archived cache only; valid mouse IDs mapped; unresolved symbols retained; duplicate names made unique.',
        'annotation_release': 'Not recorded in the archived cache.',
        'source': {'path': str(LEGACY_CACHE), 'sha256': checksum},
        'mapping_sha256': checksum, 'symbol_count': len(table),
        'resolved_mouse_ids': int(valid.sum()),
        'empty_entries': int(table.ensembl_gene_clean.eq('').sum()),
        'invalid_nonempty_entries': table.loc[~valid & table.ensembl_gene_clean.ne('')].to_dict('records'),
        'supplementation': 'None; no AER, BioMart, or online additions.',
    }
    (OUTPUT_DIR / 'mapping_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    main()

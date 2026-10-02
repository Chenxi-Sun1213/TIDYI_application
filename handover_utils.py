"""Small shared utilities for the RA handover master scripts.

This module deliberately contains orchestration and provenance only. Scientific
calculations remain in the active analysis notebook or in the corresponding
master script.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Any, Iterable, Mapping, Sequence


TIDYI_RATE_COLUMNS = (
    "analysis_dataset",
    "comparison_name",
    "comparison_group",
    "tissue",
    "celltype",
    "prolif_rate",
    "death_rate",
    "growth_rate",
    "turn_over",
    "percentage",
    "model",
    "age",
    "sex",
)


@dataclass(frozen=True)
class SourceReplacement:
    """An exact, audited source replacement for one notebook cell."""

    old: str
    new: str
    expected_count: int = 1


def require_paths(paths: Iterable[Path]) -> None:
    """Fail before analysis if an authoritative input is missing."""

    missing = [str(Path(path)) for path in paths if not Path(path).exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n- " + "\n- ".join(missing))


def configure_tidyi_source(source_root: Path | None) -> None:
    """Make an optional source checkout importable by the notebook kernel.

    ``source_root`` must be the directory that contains the importable ``tidyi``
    package/module. When omitted, the kernel's existing environment is used.
    """

    if source_root is None:
        return
    source_root = source_root.expanduser().resolve()
    if not source_root.exists():
        raise FileNotFoundError(f"TIDYI source root does not exist: {source_root}")
    current = os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONPATH"] = (
        str(source_root) if not current else f"{source_root}{os.pathsep}{current}"
    )


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return a SHA-256 checksum without loading the whole file into memory."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, checksum_limit_bytes: int = 100 * 1024 * 1024) -> dict[str, Any]:
    """Record lightweight provenance; hash only files at or below 100 MB."""

    path = Path(path).resolve()
    stat = path.stat()
    if path.is_dir():
        return {
            "path": str(path),
            "type": "directory",
            "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
            "sha256": None,
            "sha256_status": "not_applicable_directory",
        }
    return {
        "path": str(path),
        "type": "file",
        "size_bytes": stat.st_size,
        "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        "sha256": sha256_file(path) if stat.st_size <= checksum_limit_bytes else None,
        "sha256_status": (
            "computed" if stat.st_size <= checksum_limit_bytes else "skipped_file_over_100_mb"
        ),
    }


def snapshot_outputs(paths: Sequence[Path]) -> dict[str, int | None]:
    """Capture output modification times before a master run."""

    return {
        str(Path(path).resolve()): (
            Path(path).stat().st_mtime_ns if Path(path).exists() else None
        )
        for path in paths
    }


def validate_refreshed_outputs(
    paths: Sequence[Path], before: Mapping[str, int | None]
) -> list[dict[str, Any]]:
    """Require every declared output to exist, be non-empty, and be refreshed."""

    problems: list[str] = []
    records: list[dict[str, Any]] = []
    for path_value in paths:
        path = Path(path_value).resolve()
        if not path.exists():
            problems.append(f"missing: {path}")
            continue
        if path.stat().st_size == 0:
            problems.append(f"empty: {path}")
            continue
        previous_mtime = before.get(str(path))
        if previous_mtime is not None and path.stat().st_mtime_ns <= previous_mtime:
            problems.append(f"not refreshed by this run: {path}")
            continue
        records.append(file_record(path))
    if problems:
        raise RuntimeError("Expected output validation failed:\n- " + "\n- ".join(problems))
    return records


def _replace_cell_source(
    source: str,
    replacements: Sequence[SourceReplacement],
    cell_number: int,
) -> str:
    for replacement in replacements:
        count = source.count(replacement.old)
        if count != replacement.expected_count:
            raise ValueError(
                f"Cell {cell_number}: expected {replacement.expected_count} occurrence(s) "
                f"of {replacement.old!r}, found {count}. The source notebook changed."
            )
        source = source.replace(replacement.old, replacement.new)
    return source


def tidyi_runtime_capture_cell(output_path: Path) -> str:
    """Return a notebook cell that records the imported TIDYI implementation."""

    return f'''\
import importlib
import importlib.metadata
import json
from pathlib import Path

_handover_tidyi = importlib.import_module("tidyi")
try:
    _handover_tidyi_distribution_version = importlib.metadata.version("tidyi")
except importlib.metadata.PackageNotFoundError:
    _handover_tidyi_distribution_version = None
_handover_tidyi_provenance = {{
    "module_file": str(Path(_handover_tidyi.__file__).resolve()),
    "module_version_attribute": getattr(_handover_tidyi, "__version__", None),
    "distribution_version": _handover_tidyi_distribution_version,
}}
Path({str(output_path)!r}).write_text(
    json.dumps(_handover_tidyi_provenance, indent=2) + "\\n"
)
'''


def run_notebook_master(
    *,
    master_name: str,
    source_notebook: Path,
    cell_numbers: Sequence[int],
    expected_outputs: Sequence[Path],
    authoritative_inputs: Sequence[Path],
    manifest_path: Path,
    kernel_name: str,
    source_replacements: Mapping[int, Sequence[SourceReplacement]] | None = None,
    prepend_cells: Sequence[str] = (),
    append_cells: Sequence[str] = (),
    post_commands: Sequence[Sequence[str]] = (),
    run_metadata: Mapping[str, Any] | None = None,
) -> None:
    """Execute an explicit subset of notebook cells in a fresh kernel.

    Cell numbers are one-based positions in the source notebook. Selecting the
    cells here makes the formerly hidden notebook execution order inspectable.
    """

    import nbformat
    from nbclient import NotebookClient

    source_notebook = Path(source_notebook).resolve()
    manifest_path = Path(manifest_path).resolve()
    replacements = source_replacements or {}
    require_paths([source_notebook, *authoritative_inputs])

    source = nbformat.read(source_notebook, as_version=4)
    selected_cells = []
    for cell_number in cell_numbers:
        if cell_number < 1 or cell_number > len(source.cells):
            raise IndexError(f"Cell {cell_number} is outside {source_notebook}")
        original = source.cells[cell_number - 1]
        if original.cell_type != "code":
            raise TypeError(f"Cell {cell_number} is not a code cell")
        cell_source = _replace_cell_source(
            original.source,
            replacements.get(cell_number, ()),
            cell_number,
        )
        selected_cells.append(nbformat.v4.new_code_cell(cell_source))

    notebook = nbformat.v4.new_notebook(
        cells=[
            *[nbformat.v4.new_code_cell(value) for value in prepend_cells],
            *selected_cells,
            *[nbformat.v4.new_code_cell(value) for value in append_cells],
        ],
        metadata={
            "kernelspec": {
                "name": kernel_name,
                "display_name": kernel_name,
                "language": "python",
            }
        },
    )

    pre_run_output_records = [
        file_record(path) for path in expected_outputs if Path(path).exists()
    ]
    before = snapshot_outputs(expected_outputs)
    started = datetime.now(timezone.utc)
    client = NotebookClient(
        notebook,
        kernel_name=kernel_name,
        timeout=None,
        allow_errors=False,
        resources={"metadata": {"path": str(source_notebook.parent)}},
    )
    client.execute()

    for command in post_commands:
        subprocess.run(
            list(command),
            cwd=source_notebook.parent,
            check=True,
            env=os.environ.copy(),
        )

    output_records = validate_refreshed_outputs(expected_outputs, before)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "master_name": master_name,
        "status": "completed",
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "source_notebook": file_record(source_notebook),
        "selected_one_based_cell_numbers": list(cell_numbers),
        "authoritative_inputs": [file_record(path) for path in authoritative_inputs],
        "pre_run_outputs": pre_run_output_records,
        "outputs": output_records,
        "run_metadata": dict(run_metadata or {}),
    }
    manifest_path.write_text(json.dumps(payload, indent=2) + "\n")

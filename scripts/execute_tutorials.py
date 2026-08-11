#!/usr/bin/env python3
"""Execute tutorial notebooks in place so RTD can render stored outputs.

Read the Docs uses ``nb_execution_mode = "off"`` and displays *stored* notebook
outputs. Re-run this script after changing tutorial code, then commit the
``.ipynb`` files so https://scatrans.readthedocs.io shows figures again.

Prefer the same Python as the Jupyter ``python3`` kernel (here often the
``scatrans`` conda env)::

    /path/to/scatrans/bin/python scripts/execute_tutorials.py
    /path/to/scatrans/bin/python scripts/execute_tutorials.py t_synthetic_visualization.ipynb

Requires data files at the repo root for real-data notebooks, plus optional
extras (``pseudobulk``, ``memento``, ``gsea``) for the full suite.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
TUT = ROOT / "docs" / "tutorials"

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("PYTHONWARNINGS", "ignore")
# Ensure local package is importable even if not installed editable.
sys.path.insert(0, str(ROOT / "src"))


# Patterns from third-party packages that drown notebook outputs (not scatrans).
_SPAM_SUBSTRINGS = (
    "where' used without",
    'where" used without',
    "unitialized memory",
    "uninitialized memory",
    "hypothesis_test.py",
    "condition = np.less_equal",
)


def _strip_third_party_warning_spam(nb) -> None:
    """Drop known third-party warning spam from stream outputs.

    Mutates notebook nodes in place (keeps nbformat NotebookNode types).
    """
    for cell in nb.cells:
        if cell.get("cell_type") != "code":
            continue
        cleaned_outputs = []
        for out in cell.get("outputs") or []:
            if out.get("output_type") == "stream" and "text" in out:
                text = out["text"]
                lines = text if isinstance(text, list) else text.splitlines(keepends=True)
                kept = [ln for ln in lines if not any(s in ln for s in _SPAM_SUBSTRINGS)]
                if not kept:
                    continue
                out["text"] = kept
            cleaned_outputs.append(out)
        cell["outputs"] = cleaned_outputs


DEFAULT = [
    "t_synthetic_visualization.ipynb",
    "t_ec_active_transcription.ipynb",
    "t_ec_gene_upset.ipynb",
    "t_ec_standalone_de_enrichment.ipynb",
    "t_ga_active_transcription.ipynb",
    "t_gse226488_partition_mechanism.ipynb",
]


def run_one(name: str, timeout: int = 7200) -> None:
    path = TUT / name
    if not path.is_file():
        raise SystemExit(f"missing notebook: {path}")
    print(f"\n==== executing {name} (timeout={timeout}s) ====", flush=True)
    nb = nbformat.read(path, as_version=4)
    for cell in nb.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    client = NotebookClient(
        nb,
        timeout=timeout,
        kernel_name="python3",
        resources={"metadata": {"path": str(TUT)}},
    )
    client.execute()
    _strip_third_party_warning_spam(nb)
    nbformat.write(nb, path)
    n_img = n_out = 0
    for cell in nb.cells:
        for o in cell.get("outputs") or []:
            n_out += 1
            data = o.get("data") or {}
            if any(str(k).startswith("image/") for k in data):
                n_img += 1
    print(f"OK {name}: outputs={n_out} images={n_img}", flush=True)


def main() -> None:
    targets = sys.argv[1:] or DEFAULT
    for name in targets:
        run_one(name)
    print("\nAll requested notebooks executed.", flush=True)


if __name__ == "__main__":
    main()

"""Package layout guards — prevent pre-split shadow modules from shipping again.

After the tl/enrich package split, top-level ``scatrans/tl.py`` and
``scatrans/enrich.py`` must never coexist with ``scatrans/tl/`` and
``scatrans/enrich/``: Python prefers the package, leaving the flat modules
as unreachable dead code that diverges from real fixes.
"""

from __future__ import annotations

from pathlib import Path

import scatrans
import scatrans.enrich
import scatrans.tl


def test_tl_and_enrich_resolve_to_packages_not_flat_modules():
    """Imports must bind to package __init__.py under tl/ and enrich/."""
    assert scatrans.tl.__file__ is not None
    assert scatrans.enrich.__file__ is not None
    assert (
        scatrans.tl.__file__.endswith(f"{Path('tl') / '__init__.py'}")
        or Path(scatrans.tl.__file__).name == "__init__.py"
    )
    assert Path(scatrans.tl.__file__).parent.name == "tl"
    assert Path(scatrans.enrich.__file__).parent.name == "enrich"


def test_no_shadow_flat_modules_beside_packages():
    """Source (or installed) tree must not ship scatrans/tl.py or enrich.py."""
    root = Path(scatrans.__file__).resolve().parent
    assert (root / "tl").is_dir(), f"expected package dir {root / 'tl'}"
    assert (root / "enrich").is_dir(), f"expected package dir {root / 'enrich'}"
    assert not (root / "tl.py").exists(), (
        f"dead shadow module {root / 'tl.py'} would be packaged but never imported "
        f"(package {root / 'tl'} takes precedence). Delete it."
    )
    assert not (root / "enrich.py").exists(), (
        f"dead shadow module {root / 'enrich.py'} would be packaged but never imported "
        f"(package {root / 'enrich'} takes precedence). Delete it."
    )


def test_no_zip_strict_keyword_in_package_source():
    """``zip(..., strict=True)`` is Python 3.10+; package supports 3.9+ (CI matrix)."""
    import re

    root = Path(scatrans.__file__).resolve().parent
    pat = re.compile(r"\bzip\s*\([^)]*\bstrict\s*=")
    offenders: list[str] = []
    for path in root.rglob("*.py"):
        text = path.read_text(encoding="utf-8", errors="replace")
        for i, line in enumerate(text.splitlines(), 1):
            code = line.split("#", 1)[0]
            if pat.search(code):
                offenders.append(f"{path.relative_to(root)}:{i}:{line.strip()}")
    assert not offenders, (
        "zip(..., strict=) is not valid on Python 3.9; found:\n  " + "\n  ".join(offenders)
    )


def test_gsea_mapping_rate_warns_and_empty_on_zero_overlap():
    """GSEA must warn on low/zero symbol overlap (same class of check as ORA)."""
    import importlib.util

    import pandas as pd
    import pytest

    import scatrans as scat

    if importlib.util.find_spec("gseapy") is None:
        pytest.skip("gseapy not installed")

    # Mixed-case mouse-like symbols vs uppercase gene-set members (Enrichr style)
    ranked = pd.Series({f"Gene{i}": float(10 - i) for i in range(20)})  # Gene0..Gene19
    gene_sets = {
        "PATH_A": [f"GENE{i}" for i in range(30)],  # all UPPERCASE → zero overlap
        "PATH_B": [f"GENE{i}" for i in range(10, 40)],
    }
    with pytest.warns(UserWarning, match="mapping rate|overlap|gene_case|Low mapping"):
        res = scat.run_gsea(
            ranked,
            gene_sets=gene_sets,
            gene_case=None,
            nperm=10,
            min_size=1,
            max_size=50,
            verbose=False,
        )
    assert res.empty
    assert res.attrs.get("reason") == "no_ranked_genes_mapped"
    mapping = (res.attrs.get("gene_set_info") or {}).get("mapping") or {}
    assert mapping.get("n_mapped") == 0

    # Same symbols with gene_case='upper' should map (overlap > 0)
    res_ok = scat.run_gsea(
        ranked,
        gene_sets=gene_sets,
        gene_case="upper",
        nperm=10,
        min_size=1,
        max_size=50,
        verbose=False,
    )
    assert res_ok.attrs.get("reason") != "no_ranked_genes_mapped"
    m = (res_ok.attrs.get("gsea_info") or {}).get("n_genes_overlap")
    if m is None:
        m = ((res_ok.attrs.get("gene_set_info") or {}).get("mapping") or {}).get("n_mapped")
    assert m is not None and m > 0


def test_gsea_does_not_auto_pick_active_score():
    """GSEA auto-rank must prefer signed logFC over non-negative active_score."""
    import numpy as np
    import pandas as pd
    import pytest

    from scatrans.enrich.gsea import _coerce_ranked_genes_dataframe, _pick_gsea_score_column

    df = pd.DataFrame(
        {
            "active_score": np.linspace(90, 10, 20),
            "logFC": np.linspace(2, -2, 20),
            "p_adj": np.full(20, 0.01),
        },
        index=[f"G{i}" for i in range(20)],
    )
    assert _pick_gsea_score_column(df, prefer=None) == "logFC"

    # No logFC → must not silently choose active_score when a residual/p column exists;
    # last-resort may still pick a numeric column, but active_score alone warns.
    no_lfc = df.drop(columns=["logFC"])
    # With only unsigned columns, pick may fall back — coercing must warn.
    with pytest.warns(UserWarning, match="one-sided|non-negative|signed"):
        _coerce_ranked_genes_dataframe(no_lfc[["active_score"]], score_column="active_score")

    # Explicit active_score also warns.
    with pytest.warns(UserWarning, match="one-sided|non-negative|signed"):
        series = _coerce_ranked_genes_dataframe(df, score_column="active_score")
    assert series is not None
    assert (series >= 0).all()

"""Tests for abundance-/length-normalized unspliced-excess residual (tl.bias)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import scatrans as scat
from scatrans.tl.bias import RESID_COL, add_abundance_normalized_residual


def _make_results(seed: int = 0):
    """Induced genes + a MALAT1-like abundance artifact + a long-gene artifact."""
    rng = np.random.default_rng(seed)
    n_ind, n_bg = 20, 78
    rows = []

    for i in range(n_ind):  # induced: moderate everything
        rows.append(
            {
                "gene": f"IND{i}",
                "unspliced_excess_delta": rng.uniform(4, 6),
                "logFC": 2.0,
                "p_adj": 1e-6,
                "total_us_counts": rng.uniform(800, 1500),
                "gene_length": rng.uniform(3000, 8000),
                "intron_number": int(rng.integers(5, 15)),
            }
        )
    for i in range(n_bg):  # background: small delta, not induced
        rows.append(
            {
                "gene": f"BG{i}",
                "unspliced_excess_delta": rng.normal(0, 1),
                "logFC": rng.normal(0, 0.2),
                "p_adj": rng.uniform(0.2, 0.9),
                "total_us_counts": rng.uniform(100, 900),
                "gene_length": rng.uniform(1000, 6000),
                "intron_number": int(rng.integers(0, 10)),
            }
        )
    # MALAT1-like: huge abundance + huge delta, NOT induced
    rows.append(
        {
            "gene": "MALAT1",
            "unspliced_excess_delta": 100.0,
            "logFC": 0.05,
            "p_adj": 0.9,
            "total_us_counts": 200000.0,
            "gene_length": 150000.0,
            "intron_number": 7,
        }
    )
    # long-gene artifact: very long, moderate-high delta, NOT induced
    rows.append(
        {
            "gene": "LONGGENE",
            "unspliced_excess_delta": 14.0,
            "logFC": -0.1,
            "p_adj": 0.8,
            "total_us_counts": 3000.0,
            "gene_length": 450000.0,
            "intron_number": 60,
        }
    )

    df = pd.DataFrame(rows).set_index("gene")
    df["unspliced_excess_residual"] = df["unspliced_excess_delta"]
    df["valid_expr"] = True
    return df


def _rank(df, col, gene):
    return int(df[col].rank(ascending=False)[gene])


def test_abundance_demotes_malat1():
    df = _make_results()
    assert _rank(df, "unspliced_excess_residual", "MALAT1") == 1  # #1 before
    out, diag = add_abundance_normalized_residual(df, method="abundance")
    assert diag["method"] == "abundance"
    assert out[RESID_COL].loc[out["valid_expr"]].notna().all()
    # MALAT1 pushed well out of the top; an induced gene now outranks it
    assert _rank(out, RESID_COL, "MALAT1") > 10
    ind_best = min(_rank(out, RESID_COL, f"IND{i}") for i in range(20))
    assert ind_best < _rank(out, RESID_COL, "MALAT1")


def test_abundance_length_also_demotes_long_gene():
    df = _make_results()
    out_ab, _ = add_abundance_normalized_residual(df, method="abundance")
    out_al, diag = add_abundance_normalized_residual(df, method="abundance_length")
    assert diag["method"] == "abundance_length"
    # the length residualization drops the long-gene artifact further than abundance-only
    assert _rank(out_al, RESID_COL, "LONGGENE") >= _rank(out_ab, RESID_COL, "LONGGENE")
    assert _rank(out_al, RESID_COL, "MALAT1") > 10


def test_floor_quantile_controls_deflation():
    df = _make_results()
    strong, _ = add_abundance_normalized_residual(df, method="abundance", floor_quantile=0.5)
    gentle, _ = add_abundance_normalized_residual(df, method="abundance", floor_quantile=0.99)
    # a lower floor deflates the extreme-abundance MALAT1 harder (worse rank)
    assert _rank(strong, RESID_COL, "MALAT1") >= _rank(gentle, RESID_COL, "MALAT1")


def test_validation_errors():
    df = _make_results()
    with pytest.raises(ValueError):
        add_abundance_normalized_residual(df, method="bogus")
    with pytest.raises(KeyError):
        add_abundance_normalized_residual(pd.DataFrame({"valid_expr": [True]}))


def test_public_api():
    assert hasattr(scat, "add_abundance_normalized_residual")
    out, _ = scat.add_abundance_normalized_residual(_make_results(), method="abundance")
    assert RESID_COL in out.columns

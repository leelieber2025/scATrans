"""Tests for reliability-adaptive weighting of the nascent leg (tl.adaptive)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

import scatrans as scat
from scatrans.tl.adaptive import (
    _auc,
    adaptive_weight,
    add_adaptive_score,
    labeling_anchor,
)


def _make_results(reliable: bool, n_induced: int = 40, n_bg: int = 160, seed: int = 0):
    """Synthetic all_results where the nascent residual either agrees (reliable)
    or disagrees (anti-correlated) with the DE-induced genes."""
    rng = np.random.default_rng(seed)
    n = n_induced + n_bg
    logFC = np.concatenate([rng.normal(2.5, 0.3, n_induced), rng.normal(0.0, 0.2, n_bg)])
    p_adj = np.concatenate([rng.uniform(1e-8, 1e-4, n_induced), rng.uniform(0.2, 0.9, n_bg)])
    induced = np.arange(n) < n_induced
    # residual: high on induced (reliable) or low on induced (anti-correlated)
    base = rng.normal(0, 1, n)
    resid = base + (3.0 if reliable else -3.0) * induced
    df = pd.DataFrame(
        {
            "active_score": rng.uniform(0, 100, n),
            "unspliced_excess_residual": resid,
            "logFC": logFC,
            "p_adj": p_adj,
            "valid_expr": True,
        },
        index=[f"g{i}" for i in range(n)],
    )
    return df, induced


@pytest.mark.parametrize(
    "reliability,expected",
    [(0.5, 0.0), (0.4, 0.0), (0.75, 1.0), (1.0, 2.0), (float("nan"), 0.0)],
)
def test_adaptive_weight_map(reliability, expected):
    assert adaptive_weight(reliability, k=4.0, w_max=2.0) == pytest.approx(expected)


def test_reliable_proxy_is_upweighted():
    df, _ = _make_results(reliable=True)
    out, diag = add_adaptive_score(df)
    assert diag["reliability_auc"] > 0.9
    assert diag["w_proxy"] > 1.0  # strong proxy allowed to lead
    assert "up-weighted" in diag["verdict"] or "leading" in diag["verdict"]
    assert out["adaptive_score"].notna().all()
    assert out["adaptive_score"].between(0, 100).all()


def test_anticorrelated_proxy_is_disabled():
    df, induced = _make_results(reliable=False)
    out, diag = add_adaptive_score(df)
    assert diag["reliability_auc"] < 0.5
    assert diag["w_proxy"] == 0.0
    assert "DISABLED" in diag["verdict"]
    # with the proxy off, ranking follows the DE signal (induced genes on top)
    # and is NOT driven by the anti-correlated residual.
    assert _auc(out["adaptive_score"].to_numpy(), induced.astype(int)) > 0.95
    assert abs(spearmanr(out["adaptive_score"], out["unspliced_excess_residual"]).correlation) < 0.5


def test_missing_columns_raises():
    df = pd.DataFrame({"logFC": [1.0], "valid_expr": [True]})
    with pytest.raises(KeyError):
        add_adaptive_score(df)


def test_no_anchor_does_not_crash():
    # no DE-induced genes -> reliability NaN -> weight 0, no exception
    df, _ = _make_results(reliable=True)
    df["logFC"] = 0.0
    df["p_adj"] = 0.9
    out, diag = add_adaptive_score(df)
    assert diag["w_proxy"] == 0.0
    assert out["adaptive_score"].notna().all()


def _make_anchor_mismatch(seed: int = 1):
    """A frame where the proxy tracks the LABELING truth but NOT the DE anchor.

    Mirrors the scNT-seq finding: DE-induced genes (fast IEGs) have depleted
    unspliced signal, so the DE anchor mis-scores reliability, while a broader
    labeling-induced set the proxy actually follows says the proxy is reliable.
    """
    rng = np.random.default_rng(seed)
    n_de, n_lab_only, n_bg = 20, 60, 160
    n = n_de + n_lab_only + n_bg
    idx = np.arange(n)
    de_induced = idx < n_de
    lab_induced = idx < (n_de + n_lab_only)  # DE genes are a subset of labeling genes
    logFC = np.where(de_induced, rng.normal(2.5, 0.3, n), rng.normal(0.0, 0.2, n))
    p_adj = np.where(de_induced, rng.uniform(1e-8, 1e-4, n), rng.uniform(0.2, 0.9, n))
    # new-RNA log2FC: high on ALL labeling-induced genes (the ground truth)
    new_log2fc = np.where(lab_induced, rng.uniform(1.5, 4.0, n), rng.uniform(-0.3, 0.5, n))
    # residual (proxy): HIGH on the labeling-only genes (still transcribing),
    # LOW on the fast DE/IEG genes (unspliced depleted) -> tracks labeling, not DE.
    resid = rng.normal(0, 1, n)
    resid = np.where(lab_induced & ~de_induced, resid + 3.0, resid)
    resid = np.where(de_induced, resid - 2.0, resid)
    df = pd.DataFrame(
        {
            "active_score": rng.uniform(0, 100, n),
            "unspliced_excess_residual": resid,
            "logFC": logFC,
            "p_adj": p_adj,
            "new_log2fc": new_log2fc,
            "valid_expr": True,
        },
        index=[f"g{i}" for i in range(n)],
    )
    return df


def test_labeling_anchor_rescues_proxy_de_anchor_disables():
    df = _make_anchor_mismatch()
    _, de_diag = add_adaptive_score(df, anchor="de")
    _, lab_diag = add_adaptive_score(df, anchor=labeling_anchor())
    # DE anchor mis-scores reliability and disables the proxy...
    assert de_diag["reliability_auc"] < 0.5
    assert de_diag["w_proxy"] == 0.0
    # ...while the labeling anchor recovers it and keeps the proxy on.
    assert lab_diag["reliability_auc"] > 0.7
    assert lab_diag["w_proxy"] > 0.0
    assert lab_diag["anchor"].startswith("labeling_anchor")
    assert lab_diag["n_anchor_induced"] > de_diag["n_anchor_induced"]


def test_anchor_accepts_callable_and_array():
    df, induced = _make_results(reliable=True)
    # callable
    _, d1 = add_adaptive_score(df, anchor=lambda e: (e["logFC"] >= 1.0).to_numpy(int))
    assert d1["anchor"] == "<lambda>"
    # plain array aligned to valid_expr rows
    _, d2 = add_adaptive_score(df, anchor=induced.astype(int))
    assert d2["anchor"] == "array"
    assert d1["w_proxy"] > 0 and d2["w_proxy"] > 0


def test_labeling_anchor_missing_column_raises():
    df, _ = _make_results(reliable=True)  # no new_log2fc column
    with pytest.raises(KeyError):
        add_adaptive_score(df, anchor=labeling_anchor())


def test_unknown_string_anchor_raises():
    df, _ = _make_results(reliable=True)
    with pytest.raises(ValueError):
        add_adaptive_score(df, anchor="labeling")


def test_public_api_and_end_to_end(adata_basic):
    assert hasattr(scat, "adaptive_active_score")
    ar, diag = scat.adaptive_active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        organism="mouse",
    )
    assert "adaptive_score" in ar.columns
    assert {"reliability_auc", "w_proxy", "verdict"} <= set(diag)

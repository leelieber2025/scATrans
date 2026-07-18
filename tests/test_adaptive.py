"""Tests for reliability-adaptive weighting of the nascent leg (tl.adaptive)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

import scatrans as scat
from scatrans.tl.adaptive import _auc, adaptive_weight, add_adaptive_score


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

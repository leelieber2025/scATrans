"""Tests for qc.regime_diagnosis — proxy-reliability pre-flight from the global
unspliced fraction, and its wiring into run_default_pipeline."""

from __future__ import annotations

import anndata as ad
import pytest
import scipy.sparse as sp

import scatrans as scat
from scatrans.qc import _reliability_from_unspliced_fraction as rel
from scatrans.qc import regime_diagnosis


def _adata_with_unspliced_fraction(frac: float, n: int = 200, g: int = 50, seed: int = 0):
    S = (sp.random(n, g, density=0.3, random_state=seed) * 10).tocsr()
    U = sp.random(n, g, density=0.3, random_state=seed + 1)
    scale = S.sum() * frac / (1 - frac) / max(U.sum(), 1e-9)
    U = (U * scale).tocsr()
    a = ad.AnnData(X=S.copy())
    a.layers["spliced"] = S
    a.layers["unspliced"] = U
    return a


def test_reliability_mapping_u_shaped():
    assert rel(0.0) == 0.0
    assert rel(0.30) == 1.0
    assert rel(0.90) == 0.0
    # monotone increasing on the low ramp, decreasing on the high ramp
    assert rel(0.05) < rel(0.10)
    assert rel(0.55) > rel(0.65)
    assert 0.0 <= rel(0.55) <= 1.0


def test_reliability_invalid_thresholds_raise():
    with pytest.raises(ValueError):
        rel(0.3, low_ok=0.5, high_ok=0.4)  # low_ok > high_ok


@pytest.mark.parametrize(
    "frac,regime", [(0.03, "low_unspliced"), (0.30, "ok"), (0.68, "high_unspliced")]
)
def test_regime_labels(frac, regime):
    d = regime_diagnosis(_adata_with_unspliced_fraction(frac))
    assert d["regime"] == regime
    assert 0.0 <= d["reliability"] <= 1.0
    assert d["basis"] == "unspliced_fraction"
    assert abs(d["unspliced_fraction"] - frac) < 0.02


def test_regime_ok_is_full_reliability_high_is_low():
    ok = regime_diagnosis(_adata_with_unspliced_fraction(0.30))
    hi = regime_diagnosis(_adata_with_unspliced_fraction(0.68))
    assert ok["reliability"] == 1.0
    assert hi["reliability"] < 0.3


def test_regime_missing_layers_raises():
    a = ad.AnnData(X=sp.random(10, 5, density=0.5).tocsr())
    with pytest.raises(ValueError, match="not found"):
        regime_diagnosis(a)


def test_pipeline_reports_regime_and_scales_confidence(adata_basic):
    res = scat.run_default_pipeline(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        run_go_enrichment=False,
        show_plot=False,
        annotate_mechanism=True,
    )
    assert "regime" in res.meta
    assert 0.0 <= res.meta["regime"]["reliability"] <= 1.0
    # annotation reliability comes from the regime pre-flight
    assert res.meta["mechanism"]["reliability"] == res.meta["regime"]["reliability"]


def test_pipeline_regime_present_even_without_mechanism(adata_basic):
    res = scat.run_default_pipeline(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        run_go_enrichment=False,
        show_plot=False,
    )
    assert "regime" in res.meta

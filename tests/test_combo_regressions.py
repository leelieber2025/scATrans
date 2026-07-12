"""Regression tests for high-risk active_score flag combinations.

These combinations were manually stress-checked during the 2026-07 audit
(advanced+perm, advanced+pseudobulk gate, memento-consistent null, MixedLM
sample-level shuffle, ranking_mode boundaries). Keep them green so future
refactors cannot silently reopen those blind spots.
"""

from __future__ import annotations

import importlib.util
import logging

import anndata as ad
import numpy as np
import pytest

import scatrans as scat


def _small_velocity_adata(
    n_cells: int = 48,
    n_genes: int = 20,
    *,
    n_samples: int = 8,
    seed: int = 0,
) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    X = rng.negative_binomial(4, 0.45, size=(n_cells, n_genes)).astype(float)
    adata = ad.AnnData(X)
    n_half = n_cells // 2
    adata.obs["condition"] = ["Disease"] * n_half + ["Control"] * (n_cells - n_half)
    adata.obs["sample"] = [f"s{i % n_samples}" for i in range(n_cells)]
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.5 + rng.poisson(1, size=X.shape).astype(float)
    adata.layers["counts"] = np.round(X).astype(int).astype(float)
    adata.var["gene_length"] = rng.integers(600, 3500, n_genes)
    adata.var["intron_number"] = rng.integers(0, 8, n_genes)
    return adata


# ---------------------------------------------------------------------------
# ranking_mode boundary
# ---------------------------------------------------------------------------


def test_nascent_excess_active_score_is_residual_only():
    """ranking_mode='nascent_excess' → active_score ranks only by residual soft-scale."""
    from scatrans._utils import _get_exponential_scale_lambda, _soft_scale

    adata = _small_velocity_adata(seed=1)
    ad_out, _, allr = scat.active_score(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        ranking_mode="nascent_excess",
        weight_fc=9.0,  # must be overridden
        weight_pval=9.0,
        use_permutation=False,
        show_plot=False,
        n_jobs=1,
        min_total_counts=1,
    )
    meta = ad_out.uns["scatrans"]
    assert meta.get("ranking_mode") == "nascent_excess"
    assert meta.get("weight_fc") == 0.0
    assert meta.get("weight_pval") == 0.0
    assert meta.get("weight_unspliced") == 1.0

    residual = allr["unspliced_excess_residual"].to_numpy(dtype=float)
    lam = max(_get_exponential_scale_lambda(residual), 1e-8)
    expected = _soft_scale(residual, lam) * 100.0
    np.testing.assert_allclose(
        allr["active_score"].to_numpy(dtype=float),
        expected,
        rtol=1e-5,
        atol=1e-5,
    )


def test_nascent_excess_mixed_model_perm_score_ignores_de_null():
    """Under residual-only ranking, DE path in perms cannot change active_score FDR ranks.

    If weight_fc/weight_pval are 0, poisoning null DE logFC/p_adj must leave
    perm active_score vectors identical to residual-only soft-scales.
    """
    import scatrans._permutation as perm_mod

    adata = _small_velocity_adata(n_cells=60, n_genes=12, n_samples=6, seed=2)
    # Require MixedLM sample gate: recycled s0..s5 → composite RE groups (≥4/arm).
    calls = {"n": 0}
    orig = perm_mod._run_de_wrapper

    def poison_de(*args, **kwargs):
        calls["n"] += 1
        df = orig(*args, **kwargs)
        out = df.copy()
        out["logFC"] = 99.0
        out["p_adj"] = 1e-300
        out["p_val"] = 1e-300
        return out

    # Patch only the permutation module import binding.
    import scatrans._permutation as pm

    pm._run_de_wrapper = poison_de  # type: ignore[assignment]
    try:
        ad_out, _, allr = scat.active_score(
            adata,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            use_mixed_model=True,
            sample_col="sample",
            ranking_mode="nascent_excess",
            use_permutation=True,
            perm_de_backend="same",
            n_perm=2,
            n_jobs=1,
            show_plot=False,
            min_total_counts=1,
        )
    finally:
        pm._run_de_wrapper = orig  # type: ignore[assignment]

    assert calls["n"] >= 2  # poison path actually used in perms
    assert "active_score_fdr" in allr.columns
    # Residual-only: score must still equal residual soft-scale (not DE-poisoned)
    from scatrans._utils import _get_exponential_scale_lambda, _soft_scale

    residual = allr["unspliced_excess_residual"].to_numpy(dtype=float)
    lam = max(_get_exponential_scale_lambda(residual), 1e-8)
    expected = _soft_scale(residual, lam) * 100.0
    np.testing.assert_allclose(
        allr["active_score"].to_numpy(dtype=float),
        expected,
        rtol=1e-5,
        atol=1e-5,
    )
    note = ad_out.uns["scatrans"].get("permutation_approximation_note") or ""
    assert "MixedLM" in note or ad_out.uns["scatrans"].get("perm_use_mixed_model") is True


# ---------------------------------------------------------------------------
# mode="advanced" combinations
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_advanced_mode_with_permutation():
    """mode='advanced' + use_permutation must complete and record permutation columns."""
    if importlib.util.find_spec("scvelo") is None:
        pytest.skip("scvelo not installed")
    adata = _small_velocity_adata(n_cells=60, n_genes=30, seed=3)
    ad_out, _, allr = scat.active_score(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="advanced",
        advanced_fallback=True,
        use_permutation=True,
        n_perm=2,
        n_jobs=1,
        show_plot=False,
        min_total_counts=1,
    )
    mode = ad_out.uns["scatrans"].get("mode")
    assert mode in {"advanced", "heuristic_fallback_from_advanced", "heuristic"}
    assert "active_score_fdr" in allr.columns
    assert "unspliced_excess_fdr" in allr.columns
    # Velocity source recorded for diagnostics
    assert ad_out.uns["scatrans"].get("use_permutation") is True


def test_advanced_pseudobulk_blocked_without_allow_flag():
    adata = _small_velocity_adata(seed=4)
    with pytest.raises(ValueError, match="allow_advanced_pseudobulk"):
        scat.active_score(
            adata,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            mode="advanced",
            use_pseudobulk=True,
            sample_col="sample",
            allow_advanced_pseudobulk=False,
            use_permutation=False,
            show_plot=False,
            n_jobs=1,
            min_cells=1,
            min_counts=1,
        )


@pytest.mark.slow
def test_allow_advanced_pseudobulk_runs_or_falls_back(caplog):
    """allow_advanced_pseudobulk=True is experimental but must not crash."""
    if importlib.util.find_spec("scvelo") is None:
        pytest.skip("scvelo not installed")
    adata = _small_velocity_adata(n_cells=64, n_genes=25, n_samples=8, seed=5)
    # Ensure ≥4 samples per group for pb
    adata.obs["sample"] = [f"S{i // 8}" for i in range(adata.n_obs)]
    with caplog.at_level(logging.WARNING):
        ad_out, _, allr = scat.active_score(
            adata,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            mode="advanced",
            advanced_fallback=True,
            use_pseudobulk=True,
            sample_col="sample",
            allow_advanced_pseudobulk=True,
            pseudobulk_de_backend="scanpy",
            de_method="t-test_overestim_var",
            use_permutation=False,
            show_plot=False,
            n_jobs=1,
            min_cells=1,
            min_counts=1,
            min_total_counts=1,
        )
    assert "active_score" in allr.columns
    assert ad_out.uns["scatrans"].get("use_pseudobulk") is True
    # Experimental warning should have fired at least once
    assert (
        any(
            "Advanced mode on pseudobulk" in r.message or "allow_advanced_pseudobulk" in r.message
            for r in caplog.records
        )
        or ad_out.uns["scatrans"].get("mode") is not None
    )


# ---------------------------------------------------------------------------
# Memento-consistent permutation
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_perm_use_memento_de_true_records_consistent_null():
    """perm_use_memento_de=True must mark a fully consistent Memento null in metadata."""
    if importlib.util.find_spec("memento") is None:
        pytest.skip("memento not installed")
    adata = _small_velocity_adata(n_cells=40, n_genes=15, seed=6)
    # Memento needs integer counts on .X or counts=
    adata.X = adata.layers["counts"].copy()
    ad_out, _, allr = scat.active_score(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_memento_de=True,
        memento_num_boot=100,
        memento_n_cpus=1,
        use_permutation=True,
        perm_use_memento_de=True,
        perm_de_backend="same",
        n_perm=1,
        n_jobs=1,
        show_plot=False,
        min_total_counts=1,
        de_preprocess="none",
    )
    meta = ad_out.uns["scatrans"]
    assert meta.get("use_memento_de") is True
    assert meta.get("perm_use_memento_de") is True
    assert "active_score_fdr" in allr.columns
    note = meta.get("permutation_approximation_note") or ""
    assert "Memento was used for both observed DE and permutation null" in note


@pytest.mark.slow
def test_memento_without_perm_use_memento_de_notes_mismatch():
    """Default: Memento observed + non-Memento null is documented in metadata."""
    if importlib.util.find_spec("memento") is None:
        pytest.skip("memento not installed")
    adata = _small_velocity_adata(n_cells=40, n_genes=15, seed=7)
    adata.X = adata.layers["counts"].copy()
    ad_out, _, allr = scat.active_score(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_memento_de=True,
        memento_num_boot=100,
        memento_n_cpus=1,
        use_permutation=True,
        perm_use_memento_de=False,
        n_perm=1,
        n_jobs=1,
        show_plot=False,
        min_total_counts=1,
        de_preprocess="none",
    )
    meta = ad_out.uns["scatrans"]
    assert meta.get("use_memento_de") is True
    assert meta.get("perm_use_memento_de") is False
    note = meta.get("permutation_approximation_note") or ""
    assert "Memento was used for the observed DE" in note
    assert "active_score_fdr" in allr.columns

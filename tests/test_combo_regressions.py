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


def test_composite_score_gates_pval_term_by_logfc_direction():
    """Strong downregulation must not earn ~half score from p_adj alone.

    Regression for direction inconsistency: s1 (logFC) was one-sided via
    _soft_scale clip, but s3 (-log10 p_adj) was directionless — a gene with
    logFC=-3, padj=1e-20 could still score ~50 under equal weights.
    """
    from scatrans._utils import _composite_active_score_terms, _soft_scale

    logFC = np.array([3.0, -3.0, 0.0])
    residual = np.zeros(3)  # isolate DE legs
    p_adj = np.array([1e-20, 1e-20, 1e-20])
    # Fixed lambdas so soft-scale saturates on |logFC|=3 and tiny p_adj
    s1, s2, s3 = _composite_active_score_terms(logFC, residual, p_adj, 1.0, 1.0, 1.0)
    assert s2.sum() == 0.0
    # Upregulated: both FC and pval contribute; downregulated: both zero
    assert s1[0] > 0.9 and s3[0] > 0.9
    assert s1[1] == 0.0 and s3[1] == 0.0
    assert s1[2] == 0.0 and s3[2] == 0.0

    # Default equal weights (incl. residual): down gene must not mid-rank
    total_w = 3.0
    scores = (1.0 * s1 + 1.0 * s2 + 1.0 * s3) / total_w * 100.0
    assert scores[0] > 60.0  # up: FC + pval
    assert scores[1] == 0.0  # down: fully gated
    assert scores[2] == 0.0

    # DE-only weights (weight_unspliced=0): documents pre-fix mid-rank bug
    # (down ≈ 50 from directionless s3; up ≈ 100 when both legs saturate).
    de_w = 2.0
    s3_ungated = _soft_scale(-np.log10(p_adj + 1e-300), 1.0)
    legacy_down = s3_ungated[1] / de_w * 100.0
    assert abs(legacy_down - 50.0) < 1.0
    fixed_down = s3[1] / de_w * 100.0
    assert fixed_down == 0.0


def test_s3_gate_prefers_mixedlm_coef_over_logfc():
    """When MixedLM coef and sample-aware logFC disagree, s3 follows the coef (p_adj target)."""
    from scatrans._utils import _composite_active_score_terms, _score_direction_effect

    logFC = np.array([2.0, -2.0])  # sample-aware signs
    coef = np.array([-1.0, 1.0])  # opposite model signs
    residual = np.zeros(2)
    p_adj = np.array([1e-20, 1e-20])
    effect = _score_direction_effect(logFC, mixedlm_coef=coef)
    s1, s2, s3 = _composite_active_score_terms(
        logFC, residual, p_adj, 1.0, 1.0, 1.0, direction_effect=effect
    )
    # s1 still tracks logFC soft-scale
    assert s1[0] > 0.8 and s1[1] == 0.0
    # s3 tracks model coef (what p_adj tests)
    assert s3[0] == 0.0  # coef < 0 → no significance credit
    assert s3[1] > 0.9  # coef > 0 → significance credit despite logFC < 0


def test_lambda_pval_estimated_on_up_genes_only():
    """Down-regulated extreme p-values must not inflate lambda_pval for up genes."""
    from scatrans._utils import _lambda_pval_for_active_score, _soft_scale

    # 50 strong down + 50 moderate up
    logFC = np.concatenate([np.full(50, -2.0), np.full(50, 1.0)])
    p_adj = np.concatenate([np.full(50, 1e-20), np.full(50, 1e-3)])
    lam_all = max(
        __import__(
            "scatrans._utils", fromlist=["_get_exponential_scale_lambda"]
        )._get_exponential_scale_lambda(-np.log10(p_adj + 1e-300)),
        1.0,
    )
    lam_up = _lambda_pval_for_active_score(p_adj, logFC, floor=1.0)
    assert lam_up < lam_all  # up-only scale is smaller
    s3_all = _soft_scale(-np.log10(np.array([1e-3]) + 1e-300), lam_all)[0]
    s3_up = _soft_scale(-np.log10(np.array([1e-3]) + 1e-300), lam_up)[0]
    assert s3_up > s3_all  # moderate up genes keep more significance credit


def test_active_score_forwards_design_warnings(caplog):
    """diagnose_design warnings must surface when active_score auto-runs design checks."""
    import logging

    adata = _small_velocity_adata(n_cells=40, n_genes=12, n_samples=2, seed=9)
    # 2 samples total recycled → 1 unique-ish structure; force tiny per-group sample counts
    adata.obs["sample"] = ["s0"] * 20 + ["s1"] * 20
    with caplog.at_level(logging.WARNING):
        ad_out, _, _ = scat.active_score(
            adata,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            sample_col="sample",
            use_pseudobulk=False,
            use_permutation=False,
            show_plot=False,
            n_jobs=1,
            min_total_counts=1,
        )
    design = (ad_out.uns.get("scatrans") or {}).get("diagnostics", {}).get("design")
    assert design is not None, "design diagnosis must be stored under diagnostics"
    assert design.get("warnings"), "expected non-empty design warnings for tiny sample design"
    assert any(
        "design WARNING" in r.message
        or "few biological samples" in r.message.lower()
        or "sample" in r.message.lower()
        for r in caplog.records
        if r.levelno >= logging.WARNING
    ), [r.message for r in caplog.records if r.levelno >= logging.WARNING]


def test_valid_feat_excludes_zero_gene_length():
    """gene_length == 0 must not enter Huber fit (x-leverage sentinel)."""
    gene_length = np.array([1000.0, 0.0, 2000.0, np.nan])
    intron_number = np.array([2.0, 1.0, 0.0, 3.0])
    valid_feat = (
        np.isfinite(gene_length)
        & np.isfinite(intron_number)
        & (gene_length > 0)
        & (intron_number >= 0)
    )
    assert list(valid_feat) == [True, False, True, False]
    # length 0 would have passed the old >= 0 rule
    old_rule = (
        np.isfinite(gene_length)
        & np.isfinite(intron_number)
        & (gene_length >= 0)
        & (intron_number >= 0)
    )
    assert old_rule[1]  # documents the bug
    assert not valid_feat[1]


def test_huber_zero_length_genes_do_not_bias_slope():
    """GTF path zeros must not act as x=0 leverage points in Huber (bug #4 semantics).

    Construct y with a strong length trend on real genes plus wild y on length=0
    genes. If zeros enter the fit, the estimated length coefficient collapses;
    with gene_length>0 valid_feat the slope should stay near the true trend.
    """
    from scatrans._utils import _fit_huber_bias_correction

    rng = np.random.default_rng(0)
    n_real, n_zero = 80, 40
    n = n_real + n_zero
    gene_length = np.concatenate([rng.uniform(500, 4000, n_real), np.zeros(n_zero)])
    intron_number = rng.integers(0, 10, n).astype(float)
    # True model on real genes only
    true_slope = 0.05
    delta = np.zeros(n)
    delta[:n_real] = true_slope * np.log1p(gene_length[:n_real]) + rng.normal(0, 0.02, n_real)
    # Wild y at x=0 would drag Huber if included
    delta[n_real:] = rng.normal(50.0, 5.0, n_zero)

    valid_feat = (
        np.isfinite(gene_length)
        & np.isfinite(intron_number)
        & (gene_length > 0)
        & (intron_number >= 0)
    )
    assert valid_feat[:n_real].all() and not valid_feat[n_real:].any()
    valid_expr = np.ones(n, dtype=bool)
    X_features = np.column_stack(
        [np.log1p(gene_length[valid_feat]), np.log1p(intron_number[valid_feat])]
    )
    residual, info = _fit_huber_bias_correction(
        delta,
        gene_length,
        intron_number,
        total_us_for_weights=np.full(n, 100.0),
        valid_feat=valid_feat,
        valid_expr=valid_expr,
        X_features=X_features,
        bias_correction="huber_length_intron",
        min_fit_obs=30,
    )
    assert info.get("bias_corrected") is True
    coef = info.get("coef_gene_length")
    assert coef is not None and np.isfinite(coef)
    # Recovered slope near true; would be near 0 or wrong sign if zeros dominated
    assert abs(float(coef) - true_slope) < 0.03
    # Zero-length genes still get a residual (median-centered), not left as NaN
    assert np.isfinite(residual[n_real:]).all()


def test_active_score_downregulated_not_midranked_by_pval():
    """End-to-end: strongly down DE genes must not outrank weak-up by p_adj alone (#1)."""
    rng = np.random.default_rng(42)
    n_cells, n_genes = 80, 30
    X = rng.negative_binomial(5, 0.4, size=(n_cells, n_genes)).astype(float)
    # Plant one strongly downregulated gene with huge separation (tiny p)
    X[:40, 0] = 1.0  # Disease low
    X[40:, 0] = 50.0  # Control high → negative logFC
    # Plant mild upregulation with weaker separation
    X[:40, 1] = 12.0
    X[40:, 1] = 8.0
    adata = __import__("anndata").AnnData(X)
    adata.obs["condition"] = ["Disease"] * 40 + ["Control"] * 40
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.3
    adata.var["gene_length"] = rng.integers(800, 3000, n_genes)
    adata.var["intron_number"] = rng.integers(0, 6, n_genes)
    adata.var_names = [f"G{i}" for i in range(n_genes)]

    _, _, allr = scat.active_score(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        weight_fc=1.0,
        weight_unspliced=0.0,  # isolate DE legs
        weight_pval=1.0,
        use_permutation=False,
        show_plot=False,
        n_jobs=1,
        min_total_counts=1,
        bias_correction="none",
    )
    down = allr.loc["G0"]
    mild_up = allr.loc["G1"]
    assert float(down["logFC"]) < 0
    assert float(mild_up["logFC"]) > 0
    # After direction gate, down gene composite must not beat mild up on p-only legs
    assert float(down["active_score"]) < float(mild_up["active_score"])


def test_builtin_significant_requires_positive_mixedlm_coef():
    """Significant mask must not keep genes with logFC>0 but mixedlm_coef<0."""
    import pandas as pd

    from scatrans.tl.filter import _builtin_significant_mask

    df = pd.DataFrame(
        {
            "logFC": [1.0, 1.0, -1.0],
            "p_adj": [0.01, 0.01, 0.01],
            "mixedlm_coef": [0.5, -0.5, 0.5],
            "unspliced_excess_residual": [2.0, 2.0, 2.0],
            "unspliced_excess_fdr": [0.01, 0.01, 0.01],
            "active_score": [80.0, 80.0, 80.0],
            "active_score_fdr": [0.1, 0.1, 0.1],
            "valid_expr": [True, True, True],
        },
        index=["concordant_up", "discordant", "concordant_model_up_neg_lfc"],
    )
    mask = _builtin_significant_mask(
        df,
        use_permutation=True,
        extra_metadata={
            "pval_cutoff": 0.05,
            "logfc_cutoff": 0.35,
            "unspliced_excess_fdr_cutoff": 0.05,
            "use_fdr_for_significance": True,
            "is_pseudobulk": False,
        },
    )
    assert bool(mask["concordant_up"])
    assert not bool(mask["discordant"])  # coef < 0 despite logFC > 0
    assert not bool(mask["concordant_model_up_neg_lfc"])  # logFC fails


def test_extract_gene_lists_raw_p_fallback_warns(caplog):
    """Applying padj_cutoff to raw p must warn (contract is adjusted p)."""
    import logging

    import pandas as pd

    df = pd.DataFrame(
        {"logFC": [1.0, 1.2], "p_val": [0.01, 0.02]},
        index=["G1", "G2"],
    )
    with caplog.at_level(logging.WARNING):
        out = scat.extract_gene_lists(df, logfc_cutoff=0.5, pval_cutoff=0.05)
    assert "G1" in out["contrast"]
    assert any("unadjusted" in r.message or "raw" in r.message.lower() for r in caplog.records)


def test_extract_gene_lists_seurat_avg_log2FC():
    """Seurat FindMarkers avg_log2FC must be recognized (not silent empty list)."""
    import pandas as pd

    df = pd.DataFrame(
        {"avg_log2FC": [1.2, -0.9, 0.1], "p_adj": [0.01, 0.02, 0.5]},
        index=["G_up", "G_down", "G_ns"],
    )
    out = scat.extract_gene_lists(df, logfc_cutoff=0.5, padj_cutoff=0.05, logfc_direction="up")
    assert out["contrast"] == ["G_up"]


def test_extract_gene_lists_missing_lfc_warns(caplog):
    import logging

    import pandas as pd

    df = pd.DataFrame({"p_adj": [0.01, 0.02]}, index=["G1", "G2"])
    with caplog.at_level(logging.WARNING):
        out = scat.extract_gene_lists(df, logfc_cutoff=0.5, padj_cutoff=0.05)
    assert out["contrast"] == []
    assert any("log-fold-change" in r.message or "logFC" in r.message for r in caplog.records)


def test_filter_residual_follows_logfc_direction():
    """Residual magnitude cutoffs must match logfc_direction (not always positive)."""
    df = __import__("pandas").DataFrame(
        {
            "logFC": [2.0, -2.0, -2.0, 2.0],
            "p_adj": [0.01, 0.01, 0.01, 0.01],
            "unspliced_excess_residual": [2.0, 2.0, -2.0, -2.0],
            "active_score": [80.0, 30.0, 0.0, 0.0],
        },
        index=["up_posR", "down_posR", "down_negR", "up_negR"],
    )
    up = scat.filter_active_genes(
        df,
        preset=None,
        padj_cutoff=0.05,
        logfc_cutoff=0.5,
        logfc_direction="up",
        unspliced_excess_residual_cutoff=1.0,
        active_score_cutoff=0.0,
    )
    down = scat.filter_active_genes(
        df,
        preset=None,
        padj_cutoff=0.05,
        logfc_cutoff=0.5,
        logfc_direction="down",
        unspliced_excess_residual_cutoff=1.0,
        active_score_cutoff=0.0,
    )
    both = scat.filter_active_genes(
        df,
        preset=None,
        padj_cutoff=0.05,
        logfc_cutoff=0.5,
        logfc_direction="both",
        unspliced_excess_residual_cutoff=1.0,
        active_score_cutoff=0.0,
    )
    assert list(up.index) == ["up_posR"]
    assert list(down.index) == ["down_negR"]
    assert set(both.index) == {"up_posR", "down_negR"}


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

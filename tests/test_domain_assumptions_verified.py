"""Empirical checks that documented domain assumptions match runtime behavior.

These are the rules users most often misunderstand (semantic, not crash bugs).
Keep green so docs and code cannot drift apart.
"""

from __future__ import annotations

import logging
import warnings

import anndata as ad
import numpy as np
import pandas as pd

import scatrans as scat
from scatrans._utils import (
    _composite_active_score_terms,
    _fit_huber_bias_correction,
    _get_exponential_scale_lambda,
    _lambda_pval_for_active_score,
    _score_direction_effect,
    _soft_scale,
)
from scatrans.enrich.compare import extract_gene_lists
from scatrans.enrich.gsea import _coerce_ranked_genes_dataframe, _pick_gsea_score_column
from scatrans.tl.filter import _builtin_significant_mask


def test_assumption_lambda_data_adaptive_rescales_same_raw_value():
    core = np.array([1.0, 0.5, 0.2])
    lam_small = _get_exponential_scale_lambda(np.concatenate([core, np.full(30, 0.1)]))
    lam_large = _get_exponential_scale_lambda(np.concatenate([core, np.full(30, 4.0)]))
    assert lam_large > lam_small * 2
    s_small = float(_soft_scale(np.array([1.0]), lam_small)[0])
    s_large = float(_soft_scale(np.array([1.0]), lam_large)[0])
    assert s_small > s_large + 0.1


def test_assumption_s1_s3_gated_s2_independent_of_de_sign():
    logFC = np.array([2.0, -2.0])
    residual = np.array([0.0, 3.0])
    p_adj = np.array([1e-20, 1e-20])
    s1, s2, s3 = _composite_active_score_terms(logFC, residual, p_adj, 1.0, 1.0, 1.0)
    assert s1[0] > 0.8 and s1[1] == 0.0
    assert s3[0] > 0.9 and s3[1] == 0.0
    assert s2[1] > 0.9
    # Down gene score entirely from s2 under equal weights
    score_down = (s1[1] + s2[1] + s3[1]) / 3.0 * 100.0
    assert 20.0 < score_down < 40.0


def test_assumption_mixedlm_s3_follows_coef_not_logfc():
    effect = _score_direction_effect(np.array([2.0, -2.0]), mixedlm_coef=np.array([-1.0, 1.0]))
    _s1, _s2, s3 = _composite_active_score_terms(
        np.array([2.0, -2.0]),
        np.zeros(2),
        np.array([1e-20, 1e-20]),
        1.0,
        1.0,
        1.0,
        direction_effect=effect,
    )
    assert s3[0] == 0.0 and s3[1] > 0.9


def test_assumption_lambda_pval_estimated_on_up_genes():
    logfc = np.concatenate([np.full(40, -2.0), np.full(40, 1.0)])
    padj = np.concatenate([np.full(40, 1e-20), np.full(40, 1e-3)])
    lam_all = max(_get_exponential_scale_lambda(-np.log10(padj + 1e-300)), 1.0)
    lam_up = _lambda_pval_for_active_score(padj, logfc, floor=1.0)
    assert lam_up < lam_all


def test_assumption_huber_excludes_nonpositive_length():
    n_good, n_bad = 40, 20
    gl = np.concatenate([np.full(n_good, 2000.0), np.zeros(n_bad)])
    intr = np.ones(n_good + n_bad)
    delta = np.concatenate([np.full(n_good, 0.1), np.full(n_bad, 50.0)])
    vf = np.isfinite(gl) & (gl > 0)
    Xfeat = np.column_stack([np.log1p(gl[vf]), np.log1p(intr[vf])])
    _, info = _fit_huber_bias_correction(
        delta,
        gl,
        intr,
        np.full(len(gl), 100.0),
        vf,
        np.ones(len(gl), dtype=bool),
        Xfeat,
        bias_correction="huber_length_intron",
        min_fit_obs=20,
    )
    assert info.get("bias_corrected")
    assert int(info["n_genes_used_for_fit"]) == n_good


def test_assumption_significant_empty_without_permutation():
    rng = np.random.default_rng(0)
    n_c, n_g = 40, 20
    X = rng.negative_binomial(4, 0.45, size=(n_c, n_g)).astype(float)
    adata = ad.AnnData(X)
    adata.obs["condition"] = ["D"] * 20 + ["C"] * 20
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.4
    adata.var["gene_length"] = rng.integers(800, 3000, n_g)
    adata.var["intron_number"] = rng.integers(0, 5, n_g)
    _ad, sig, _allr = scat.active_score(
        adata,
        groupby="condition",
        target_group="D",
        reference_group="C",
        use_permutation=False,
        show_plot=False,
        n_jobs=1,
        min_total_counts=1,
        bias_correction="none",
    )
    assert len(sig) == 0
    scoring = _ad.uns["scatrans"]["diagnostics"]["scoring"]
    assert scoring["lambda_fc"] > 0 and "scale_note" in scoring


def test_assumption_filter_residual_follows_direction_and_padj():
    df = pd.DataFrame(
        {
            "logFC": [1.0, -1.0, -1.5],
            "p_adj": [0.01, 0.01, 0.2],
            "unspliced_excess_residual": [2.0, -2.0, -3.0],
            "active_score": [80.0, 10.0, 5.0],
        },
        index=["up", "down_res", "ns"],
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
    assert list(up.index) == ["up"]
    assert list(down.index) == ["down_res"]


def test_assumption_extract_warns_on_raw_p_and_accepts_avg_log2FC(caplog):
    df_raw = pd.DataFrame({"logFC": [1.0], "p_val": [0.01]}, index=["g"])
    with caplog.at_level(logging.WARNING):
        out = extract_gene_lists(df_raw, logfc_cutoff=0.5, pval_cutoff=0.05)
    assert out["contrast"] == ["g"]
    assert any("unadjusted" in r.message or "raw" in r.message.lower() for r in caplog.records)

    df_s = pd.DataFrame({"avg_log2FC": [1.2, -0.9], "p_adj": [0.01, 0.01]}, index=["U", "D"])
    out_s = extract_gene_lists(df_s, logfc_cutoff=0.5, padj_cutoff=0.05, logfc_direction="up")
    assert out_s["contrast"] == ["U"]


def test_assumption_gsea_prefers_signed_and_warns_unsigned():
    df = pd.DataFrame(
        {
            "active_score": np.linspace(90, 10, 20),
            "logFC": np.linspace(2, -2, 20),
        },
        index=[f"G{i}" for i in range(20)],
    )
    assert _pick_gsea_score_column(df, prefer=None) == "logFC"
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _coerce_ranked_genes_dataframe(df[["active_score"]], score_column="active_score")
    assert any(
        issubclass(x.category, UserWarning)
        and any(k in str(x.message).lower() for k in ("one-sided", "non-negative", "signed"))
        for x in w
    )


def test_assumption_significant_requires_positive_mixedlm_coef():
    df = pd.DataFrame(
        {
            "logFC": [1.0, 1.0],
            "p_adj": [0.01, 0.01],
            "mixedlm_coef": [0.5, -0.5],
            "unspliced_excess_residual": [2.0, 2.0],
            "unspliced_excess_fdr": [0.01, 0.01],
            "active_score": [80.0, 80.0],
            "active_score_fdr": [0.1, 0.1],
            "valid_expr": [True, True],
        },
        index=["ok", "discordant"],
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
    assert bool(mask["ok"]) and not bool(mask["discordant"])

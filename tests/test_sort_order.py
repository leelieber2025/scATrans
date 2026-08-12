"""Regression tests for the DE-first default row ordering of ``active_score``
outputs, and the ``preset="significant"`` order contract.

Covers three points:
  1. The paired sort-key builder ``_de_first_sort_keys`` — in particular the
     edge where ``p_adj`` is absent must NOT flip ``logFC`` to ascending.
  2. ``active_score`` / ``active_score_simple`` return tables sorted DE-first
     (``p_adj`` ascending) — the new default after demoting the composite score.
  3. ``filter_active_genes(preset="significant")`` reproduces the built-in
     ``significant`` list in the SAME order (no row-order drift).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

import scatrans as scat
from scatrans._utils import UNSPLICED_EXCESS_FDR_COL, UNSPLICED_EXCESS_RESIDUAL_COL
from scatrans.tl._common import _de_first_sort_keys
from scatrans.tl.filter import filter_active_genes


# --- 1. unit: paired sort-key builder -------------------------------------
def test_de_first_sort_keys_both_columns():
    assert _de_first_sort_keys(["active_score", "logFC", "p_adj", "x"]) == (
        ["p_adj", "logFC"],
        [True, False],
    )


def test_de_first_sort_keys_logfc_only_stays_descending():
    # The bug this guards: with p_adj absent, positional pairing would sort
    # logFC ASCENDING. It must remain descending.
    assert _de_first_sort_keys(["logFC", "active_score"]) == (["logFC"], [False])


def test_de_first_sort_keys_padj_only():
    assert _de_first_sort_keys(["p_adj", "active_score"]) == (["p_adj"], [True])


def test_de_first_sort_keys_active_score_fallback():
    assert _de_first_sort_keys(["active_score", "x"]) == (["active_score"], [False])


def test_de_first_sort_keys_none_present():
    assert _de_first_sort_keys(["x", "y"]) == ([], [])


# --- integration fixture with planted induced genes -----------------------
def _adata_with_induced():
    rng = np.random.default_rng(0)
    n_cells, n_genes = 240, 100
    base = rng.negative_binomial(6, 0.4, size=(n_cells, n_genes)).astype(float)
    cond = np.array(["Disease"] * 120 + ["Control"] * 120)
    disease = np.where(cond == "Disease")[0]
    induced = np.arange(15)
    base[np.ix_(disease, induced)] *= 3.0  # up in Disease -> positive logFC
    ad = sc.AnnData(base)
    ad.obs["condition"] = cond
    ad.layers["spliced"] = base.copy()
    uns = base * 0.4
    # Strong extra nascent (unspliced) for induced genes so the nascent excess is
    # clearly POSITIVE (U_target >> gamma_ref * S_target), passing the strict gate.
    uns[np.ix_(disease, induced)] *= 8.0
    ad.layers["unspliced"] = uns
    ad.var["gene_length"] = rng.integers(700, 4500, n_genes)
    ad.var["intron_number"] = rng.integers(0, 12, n_genes)
    return ad


# --- 2. new default: active_score outputs are DE-first ordered -------------
def test_active_score_all_results_is_de_first_sorted():
    _, _significant, all_results = scat.active_score(
        _adata_with_induced(),
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        show_plot=False,
    )
    pa = all_results["p_adj"].to_numpy(dtype=float)
    pa = pa[~np.isnan(pa)]
    # p_adj non-decreasing top-to-bottom (ties broken by logFC desc), i.e. NOT
    # ordered by the (legacy) composite active_score anymore.
    assert np.all(pa[:-1] <= pa[1:] + 1e-12), "all_results must be sorted p_adj ascending"
    assert not all_results["active_score"].is_monotonic_decreasing or len(all_results) <= 1


# --- 3. preset='significant' order contract (no drift vs built-in) ---------
def test_preset_significant_orders_de_first_not_by_active_score():
    # Three genes that pass every strict gate. Their active_score order (B, C, A)
    # differs from the DE-first order (A, B, C) — so the result distinguishes the
    # two. The built-in ``significant`` list is DE-first, and preset='significant'
    # must reproduce that SAME order (this test guards the row-order drift).
    df = pd.DataFrame(
        {
            "logFC": [3.0, 2.5, 2.0],
            "p_adj": [1e-9, 1e-6, 1e-3],
            "active_score": [80.0, 99.0, 90.0],
            "active_score_fdr": [0.01, 0.01, 0.01],
            UNSPLICED_EXCESS_RESIDUAL_COL: [5.0, 5.0, 5.0],
            UNSPLICED_EXCESS_FDR_COL: [0.01, 0.01, 0.01],
            "valid_expr": [True, True, True],
        },
        index=["A", "B", "C"],
    )
    df.attrs["scatrans_filter_context"] = {
        "use_permutation": True,
        "use_fdr_for_significance": True,
        "is_pseudobulk": False,
    }
    out = filter_active_genes(df, preset="significant")
    assert list(out.index) == ["A", "B", "C"]  # DE-first, not active_score order (B, C, A)

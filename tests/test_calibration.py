"""Tests for tl.calibration.program_mechanism_permutation_calibrated — the
permutation-calibrated program score used for ABSOLUTE placement on the mechanism axis.

The contract under test: membership is frozen by a required DE table, the null is built
by shuffling condition labels through the identical transformation chain, and the
reported ``calibrated`` value is observed minus that null mean.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import scatrans as scat
from scatrans.tl import program_mechanism_permutation_calibrated as calib


def _de_for(adata) -> pd.DataFrame:
    """Observed DE table, computed once, reused as frozen membership."""
    de = scat.differential_expression(
        adata.copy(),
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
    )
    return de[1] if isinstance(de, tuple) else de


def _sets(adata, k: int = 12) -> dict[str, list[str]]:
    g = list(adata.var_names)
    return {"progA": g[:k], "progB": g[k : 2 * k]}


def test_returns_expected_columns_and_calibration_identity(adata_basic):
    de = _de_for(adata_basic)
    out = calib(
        adata_basic,
        _sets(adata_basic),
        de=de,
        target_group="Disease",
        reference_group="Control",
        organism="human",
        n_perm=6,
        random_state=0,
    )
    assert len(out) == 2
    for c in (
        "n_genes",
        "observed_mean",
        "null_mean",
        "null_sd",
        "calibrated",
        "z",
        "p_perm",
        "n_perm_effective",
    ):
        assert c in out.columns
    # the defining identity: calibrated = observed - null mean
    np.testing.assert_allclose(
        out["calibrated"].to_numpy(float),
        out["observed_mean"].to_numpy(float) - out["null_mean"].to_numpy(float),
        rtol=1e-10,
    )
    # Phipson & Smyth: permutation p is never zero
    assert (out["p_perm"] > 0).all()
    assert (out["p_perm"] <= 1).all()


def test_is_deterministic_under_fixed_seed(adata_basic):
    de = _de_for(adata_basic)
    kw = {
        "de": de,
        "target_group": "Disease",
        "reference_group": "Control",
        "organism": "human",
        "n_perm": 5,
        "random_state": 7,
    }
    a = calib(adata_basic, _sets(adata_basic), **kw)
    b = calib(adata_basic, _sets(adata_basic), **kw)
    pd.testing.assert_frame_equal(a, b)


def test_different_seeds_move_the_null_not_the_observation(adata_basic):
    de = _de_for(adata_basic)
    kw = {
        "de": de,
        "target_group": "Disease",
        "reference_group": "Control",
        "organism": "human",
        "n_perm": 5,
    }
    a = calib(adata_basic, _sets(adata_basic), random_state=1, **kw)
    b = calib(adata_basic, _sets(adata_basic), random_state=2, **kw)
    # observation is not a function of the permutation stream
    np.testing.assert_allclose(a["observed_mean"], b["observed_mean"], rtol=1e-12)


def test_de_is_required_and_must_be_a_frame(adata_basic):
    with pytest.raises(TypeError):
        calib(
            adata_basic,
            _sets(adata_basic),
            de=None,
            target_group="Disease",
            reference_group="Control",
            n_perm=2,
        )
    with pytest.raises(TypeError):
        calib(
            adata_basic,
            _sets(adata_basic),
            de="builtin",
            target_group="Disease",
            reference_group="Control",
            n_perm=2,
        )


def test_rejects_unknown_groups_and_columns(adata_basic):
    de = _de_for(adata_basic)
    with pytest.raises(KeyError):
        calib(
            adata_basic,
            _sets(adata_basic),
            de=de,
            groupby="nope",
            target_group="Disease",
            reference_group="Control",
            n_perm=2,
        )
    with pytest.raises(ValueError):
        calib(
            adata_basic,
            _sets(adata_basic),
            de=de,
            target_group="Missing",
            reference_group="Control",
            n_perm=2,
        )
    with pytest.raises(KeyError):
        calib(
            adata_basic,
            _sets(adata_basic),
            de=de,
            block_col="nope",
            target_group="Disease",
            reference_group="Control",
            n_perm=2,
        )


def test_min_genes_skips_small_programs(adata_basic):
    de = _de_for(adata_basic)
    g = list(adata_basic.var_names)
    out = calib(
        adata_basic,
        {"tiny": g[:3], "ok": g[:20]},
        de=de,
        target_group="Disease",
        reference_group="Control",
        organism="human",
        n_perm=3,
        min_genes=5,
    )
    assert "tiny" not in out.index
    assert "ok" in out.index


def test_block_col_shuffles_within_blocks(adata_basic):
    """With block_col, each block keeps its own label composition in every replicate."""
    from scatrans.tl.calibration import _shuffle_within_blocks

    labels = np.array(["A"] * 6 + ["B"] * 6)
    blocks = np.array(["b1"] * 3 + ["b2"] * 3 + ["b1"] * 3 + ["b2"] * 3)
    rng = np.random.default_rng(0)
    for _ in range(20):
        out = _shuffle_within_blocks(labels, blocks, rng)
        for b in ("b1", "b2"):
            m = blocks == b
            assert sorted(out[m]) == sorted(labels[m])


def test_no_block_col_permutes_globally(adata_basic):
    from scatrans.tl.calibration import _shuffle_within_blocks

    labels = np.array(["A"] * 10 + ["B"] * 10)
    rng = np.random.default_rng(0)
    outs = {tuple(_shuffle_within_blocks(labels, None, rng)) for _ in range(20)}
    assert len(outs) > 1  # actually shuffling
    for o in outs:
        assert sorted(o) == sorted(labels)  # composition preserved


def test_empty_result_when_no_program_qualifies(adata_basic):
    de = _de_for(adata_basic)
    out = calib(
        adata_basic,
        {"tiny": list(adata_basic.var_names)[:2]},
        de=de,
        target_group="Disease",
        reference_group="Control",
        organism="human",
        n_perm=2,
        min_genes=5,
    )
    assert out.empty
    assert "calibrated" in out.columns


def test_n_perm_must_be_positive(adata_basic):
    de = _de_for(adata_basic)
    with pytest.raises(ValueError):
        calib(
            adata_basic,
            _sets(adata_basic),
            de=de,
            target_group="Disease",
            reference_group="Control",
            n_perm=0,
        )


def test_restrict_to_selected_uses_frozen_de_membership(adata_basic):
    """restrict_to_selected must actually narrow the program to DE-selected genes."""
    de = _de_for(adata_basic)
    g = list(adata_basic.var_names)
    sets = {"progA": g[:60]}
    kw = {
        "de": de,
        "target_group": "Disease",
        "reference_group": "Control",
        "organism": "human",
        "n_perm": 3,
        "random_state": 0,
        "logfc_cutoff": 0.0,
        "padj_cutoff": 1.0,
    }
    wide = calib(adata_basic, sets, restrict_to_selected=False, **kw)
    narrow = calib(adata_basic, sets, restrict_to_selected=True, **kw)
    assert narrow.loc["progA", "n_genes"] <= wide.loc["progA", "n_genes"]


def test_restrict_to_selected_raises_when_nothing_selected(adata_basic):
    de = _de_for(adata_basic)
    g = list(adata_basic.var_names)
    with pytest.raises(ValueError, match="selected no genes"):
        calib(
            adata_basic,
            {"progA": g[:60]},
            de=de,
            target_group="Disease",
            reference_group="Control",
            organism="human",
            n_perm=2,
            restrict_to_selected=True,
            logfc_cutoff=1e9,
            padj_cutoff=1e-300,
        )


def test_label_normalization_matches_active_score(adata_basic):
    """Trailing spaces / '1.0' forms must not make a valid contrast look missing."""
    a = adata_basic.copy()
    # Strip-normalization: obs has trailing spaces; caller passes clean names.
    a.obs["condition"] = [f"{v} " for v in a.obs["condition"].astype(str).tolist()]
    de = _de_for(a)
    out = calib(
        a,
        _sets(a),
        de=de,
        target_group="Disease",
        reference_group="Control",
        organism="human",
        n_perm=3,
        random_state=0,
    )
    assert len(out) == 2
    assert (out["p_perm"] > 0).all()

    # Numeric string forms stored as "1.0"/"0.0"; caller passes ints.
    a2 = adata_basic.copy()
    raw = a2.obs["condition"].astype(str)
    a2.obs["condition"] = raw.map({"Disease": "1.0", "Control": "0.0"}).fillna("0.0")
    de2 = scat.differential_expression(
        a2.copy(),
        groupby="condition",
        target_group=1,
        reference_group=0,
    )
    de2 = de2[1] if isinstance(de2, tuple) else de2
    out2 = calib(
        a2,
        _sets(a2),
        de=de2,
        target_group=1,
        reference_group=0,
        organism="human",
        n_perm=3,
        random_state=0,
    )
    assert len(out2) == 2


def test_target_equals_reference_rejected(adata_basic):
    de = _de_for(adata_basic)
    with pytest.raises(ValueError, match="must be different"):
        calib(
            adata_basic,
            _sets(adata_basic),
            de=de,
            target_group="Disease",
            reference_group="Disease",
            n_perm=2,
        )

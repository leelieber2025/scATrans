"""Tests for tl.nascent.nascent_activity_score and its additive, DECOUPLED
integration into partition_de_by_mechanism(add_nascent_score=True).

Covers: column contract, nascent-detection discrimination, DE reproducibility flag,
error paths, and the partition wiring — crucially that the nascent detection columns
are additive and do NOT drive the transcription-vs-stabilization mechanism call
(which stays on the induction-normalized residual).
"""

from __future__ import annotations

import numpy as np
import pytest
import scanpy as sc

import scatrans as scat
from scatrans.tl import nascent_activity_score
from scatrans.tl.mechanism import SUPPORT_COL
from scatrans.tl.nascent import NASCENT_Z_COL, REPRO_COL, REPRO_FRAC_COL


@pytest.fixture(scope="module")
def adata_mech():
    """Known-truth contrast: g0-19 transcription-driven (unspliced+spliced up),
    g20-39 stabilization-driven (spliced up, unspliced flat), rest background."""
    rng = np.random.default_rng(0)
    ng, nc = 200, 400
    cond = np.array(["ctrl"] * (nc // 2) + ["trt"] * (nc // 2))
    S = rng.poisson(5, size=(nc, ng)).astype(float)
    U = rng.poisson(2, size=(nc, ng)).astype(float)
    trt = cond == "trt"
    U[np.ix_(trt, np.arange(20))] += 8
    S[np.ix_(trt, np.arange(20))] += 6
    S[np.ix_(trt, np.arange(20, 40))] += 6
    a = sc.AnnData(S.copy())
    a.obs["condition"] = cond
    a.obs["sample"] = ["s" + str(i % 6) for i in range(nc)]
    a.layers["spliced"] = S
    a.layers["unspliced"] = U
    a.var_names = [f"g{i}" for i in range(ng)]
    a.var["gene_length"] = rng.integers(700, 40000, ng)
    a.var["intron_number"] = rng.integers(0, 30, ng)
    return a


def test_columns_and_index(adata_mech):
    nz = nascent_activity_score(adata_mech, "condition", "trt", "ctrl")
    assert list(nz.columns) == [
        NASCENT_Z_COL,
        "dlog_unspliced",
        "dlog_spliced",
        REPRO_COL,
        REPRO_FRAC_COL,
    ]
    assert list(nz.index) == list(adata_mech.var_names)
    assert nz[REPRO_FRAC_COL].between(0, 1).all()
    assert nz[REPRO_COL].dtype == bool


def test_discriminates_transcription_from_stabilization(adata_mech):
    nz = nascent_activity_score(adata_mech, "condition", "trt", "ctrl")
    txn = nz.iloc[:20][NASCENT_Z_COL].median()
    stab = nz.iloc[20:40][NASCENT_Z_COL].median()
    bg = nz.iloc[40:][NASCENT_Z_COL].median()
    # transcription-driven (nascent gain) >> stabilization (spliced-only) ~ background
    assert txn > stab + 5
    assert abs(stab - bg) < abs(txn - bg)


def test_requires_explicit_groups(adata_mech):
    # treatment direction is never guessed
    with pytest.raises(ValueError, match="required"):
        nascent_activity_score(adata_mech, "condition")


def test_missing_layer_raises(adata_mech):
    a = adata_mech.copy()
    del a.layers["unspliced"]
    with pytest.raises(KeyError):
        nascent_activity_score(a, "condition", "trt", "ctrl")


def test_too_few_cells_raises(adata_mech):
    a = adata_mech[adata_mech.obs.condition == "trt"].copy()
    a.obs["condition"] = ["trt"] * (a.n_obs - 1) + ["ctrl"]  # 1 ref cell
    with pytest.raises(ValueError):
        nascent_activity_score(a, "condition", "trt", "ctrl")


def test_sample_col_used_for_reproducibility(adata_mech):
    # with a sample column the folds are the samples (proper cross-replicate check)
    nz = nascent_activity_score(adata_mech, "condition", "trt", "ctrl", sample_col="sample")
    assert nz[REPRO_COL].any()


def test_bad_sample_col_raises(adata_mech):
    # an explicit but missing sample_col must not silently fall back to random folds
    with pytest.raises(KeyError, match="sample_col"):
        nascent_activity_score(adata_mech, "condition", "trt", "ctrl", sample_col="nope")


def test_unknown_group_label_lists_available(adata_mech):
    with pytest.raises(ValueError, match="not found"):
        nascent_activity_score(adata_mech, "condition", "Treated", "ctrl")


def test_n_splits_validation(adata_mech):
    with pytest.raises(ValueError, match="n_splits"):
        nascent_activity_score(adata_mech, "condition", "trt", "ctrl", n_splits=0)


def test_flat_gene_not_reproducible():
    # genes with no spliced change (overall sign == 0) must not be flagged reproducible
    ng, nc = 10, 40
    cond = np.array(["ctrl"] * 20 + ["trt"] * 20)
    S = np.ones((nc, ng))  # identical spliced everywhere -> dlog_spliced == 0
    U = np.ones((nc, ng))
    U[cond == "trt", 0] += 5  # a nascent gain, but spliced is flat
    a = sc.AnnData(S.copy())
    a.obs["condition"] = cond
    a.layers["spliced"] = S
    a.layers["unspliced"] = U
    a.var_names = [f"g{i}" for i in range(ng)]
    nz = nascent_activity_score(a, "condition", "trt", "ctrl")
    assert (nz["dlog_spliced"].abs() < 1e-9).all()
    assert not nz[REPRO_COL].any()


def test_mature_nascent_layers_resolved(adata_mech):
    # kb_python-style mature/nascent layers resolve for standalone use
    a = adata_mech.copy()
    a.layers["mature"] = a.layers.pop("spliced")
    a.layers["nascent"] = a.layers.pop("unspliced")
    nz = nascent_activity_score(a, "condition", "trt", "ctrl")
    assert nz.iloc[:20][NASCENT_Z_COL].median() > nz.iloc[40:][NASCENT_Z_COL].median()


# --------------------------- partition integration --------------------------
def _partition(a, **kw):
    return scat.partition_de_by_mechanism(
        a,
        groupby="condition",
        target_group="trt",
        reference_group="ctrl",
        organism="human",
        **kw,
    )


def test_partition_add_nascent_score_columns(adata_mech):
    r = _partition(adata_mech, add_nascent_score=True)
    assert r.meta["nascent_score"]["enabled"] is True
    assert r.meta["nascent_score"]["status"] == "ok"
    # detection columns land on the full gene table
    for c in (NASCENT_Z_COL, REPRO_COL):
        assert c in r.gene_table.columns
    # mechanism support is STILL present (from the residual, not the nascent z)
    assert SUPPORT_COL in r.gene_table.columns


def test_nascent_score_does_not_drive_mechanism(adata_mech):
    # DECOUPLING: the mechanism partition (support + program directions) must be
    # identical with and without the nascent detection columns.
    gs = {"txn": [f"g{i}" for i in range(20)], "stab": [f"g{i}" for i in range(20, 40)]}
    r0 = _partition(adata_mech, gene_sets=gs)
    r1 = _partition(adata_mech, add_nascent_score=True, gene_sets=gs)
    assert (
        r0.programs.set_index("program")["direction"].to_dict()
        == r1.programs.set_index("program")["direction"].to_dict()
    )
    s0 = r0.gene_table[SUPPORT_COL].to_numpy()
    s1 = r1.gene_table[SUPPORT_COL].to_numpy()
    assert np.allclose(s0, s1, equal_nan=True)


def test_partition_default_no_nascent(adata_mech):
    r = _partition(adata_mech)  # default add_nascent_score=False
    assert r.meta["nascent_score"]["enabled"] is False
    assert NASCENT_Z_COL not in r.gene_table.columns
    assert SUPPORT_COL in r.gene_table.columns  # annotated from the residual

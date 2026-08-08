"""Tests for the scat.datasets toy-data loader."""

import numpy as np
import pytest

import scatrans as scat


def test_load_toy_default_shape():
    adata = scat.datasets.load_toy()
    assert adata.shape == (240, 200)
    assert set(adata.obs["condition"]) == {"Control", "Disease"}
    assert "spliced" in adata.layers and "unspliced" in adata.layers
    assert "ground_truth_mechanism" in adata.var
    assert adata.X.shape == adata.layers["spliced"].shape


def test_load_toy_deterministic():
    a = scat.datasets.load_toy(seed=7)
    b = scat.datasets.load_toy(seed=7)
    assert np.allclose(a.X, b.X)
    assert list(a.var_names) == list(b.var_names)


def test_load_toy_different_seed_differs():
    a = scat.datasets.load_toy(seed=0)
    b = scat.datasets.load_toy(seed=1)
    assert not np.allclose(a.X, b.X)


def test_load_toy_real_gene_symbols_no_mapping_warning(recwarn):
    adata = scat.datasets.load_toy(n_genes=50, seed=0)
    scat.add_gene_features(adata, organism="mouse")
    assert not adata.var["gene_length"].isna().all()
    assert not any("mapping rate" in str(w.message) for w in recwarn.list)


def test_load_toy_human_organism():
    adata = scat.datasets.load_toy(organism="human", n_genes=50, seed=0)
    assert adata.n_vars == 50
    assert len(set(adata.var_names)) == 50


def test_load_toy_sample_col_structure():
    adata = scat.datasets.load_toy(n_samples=4, seed=0)
    counts = adata.obs.groupby(["condition", "sample"]).size()
    assert len(counts) == 4
    for cond in ("Control", "Disease"):
        samples_in_cond = adata.obs.loc[adata.obs["condition"] == cond, "sample"].unique()
        assert len(samples_in_cond) == 2


@pytest.mark.parametrize("n_samples", [1, 3, 5])
def test_load_toy_rejects_odd_or_too_small_n_samples(n_samples):
    with pytest.raises(ValueError, match="n_samples"):
        scat.datasets.load_toy(n_samples=n_samples)


def test_load_toy_rejects_tiny_n_genes():
    with pytest.raises(ValueError, match="n_genes"):
        scat.datasets.load_toy(n_genes=2)


def test_load_toy_partition_recovers_ground_truth():
    """The whole point of the toy set: partition_de_by_mechanism should
    recover the simulated transcription/stabilization labels cleanly."""
    adata = scat.datasets.load_toy(seed=0)
    result = scat.partition_de_by_mechanism(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        organism="mouse",
        de="builtin",
    )
    assert result.regime["regime"] == "ok"
    assert len(result.selected) > 0

    truth = adata.var.loc[result.selected.index, "ground_truth_mechanism"]
    calls = result.selected["mechanism_class"]

    txn_calls = calls[truth == "transcription_up"]
    stab_calls = calls[truth == "stabilization_up"]
    assert len(txn_calls) > 0
    assert len(stab_calls) > 0
    assert (txn_calls == "transcription-driven").mean() > 0.7
    assert (stab_calls == "stabilization-driven").mean() > 0.7

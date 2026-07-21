"""Pytest plugin: shared AnnData fixtures for the scATrans test suite.

Loaded via ``pytest_plugins`` from root ``conftest.py`` (and ``tests/conftest.py``)
so fixtures are registered even when nested conftest discovery fails in CI.
"""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc


@pytest.fixture(scope="module")
def adata_basic():
    """AnnData with spliced/unspliced layers + gene features."""
    np.random.seed(42)
    n_cells, n_genes = 120, 250
    X = np.random.negative_binomial(4, 0.45, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["Disease"] * 60 + ["Control"] * 60
    ad.obs["sample"] = ["s" + str(i % 8) for i in range(n_cells)]
    ad.layers["spliced"] = X.copy()
    ad.layers["unspliced"] = X * 0.55
    ad.var["gene_length"] = np.random.randint(700, 4500, n_genes)
    ad.var["intron_number"] = np.random.randint(0, 12, n_genes)
    return ad


@pytest.fixture(scope="module")
def adata_de_only():
    """Count AnnData without velocity layers (pure DE path)."""
    np.random.seed(99)
    n_cells, n_genes = 80, 120
    X = np.random.negative_binomial(5, 0.4, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["Disease"] * 40 + ["Control"] * 40
    ad.obs["sample"] = ["s" + str(i % 6) for i in range(n_cells)]
    return ad


@pytest.fixture(scope="module")
def adata_small_reference():
    """Synthetic data with a small reference group (empirical Bayes gamma)."""
    np.random.seed(123)
    n_ref, n_tgt, n_genes = 18, 80, 200
    n_cells = n_ref + n_tgt
    X = np.random.negative_binomial(4, 0.45, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["Control"] * n_ref + ["Disease"] * n_tgt
    ad.layers["spliced"] = X.copy()
    ad.layers["unspliced"] = X * 0.5
    ad.var["gene_length"] = np.random.randint(700, 4500, n_genes)
    ad.var["intron_number"] = np.random.randint(0, 12, n_genes)
    return ad


@pytest.fixture(scope="module")
def adata_high_unspliced():
    """Very high unspliced fraction for qc warning path."""
    np.random.seed(7)
    n_cells, n_genes = 40, 50
    X = np.random.negative_binomial(3, 0.5, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["A"] * 20 + ["B"] * 20
    ad.layers["spliced"] = X * 0.2
    ad.layers["unspliced"] = X * 2.0
    return ad


@pytest.fixture(scope="module")
def adata_mature_nascent():
    """kb_python style layer names."""
    np.random.seed(123)
    n_cells, n_genes = 80, 180
    X = np.random.negative_binomial(3, 0.5, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["GA"] * 40 + ["Ctrl"] * 40
    ad.layers["mature"] = X.copy()
    ad.layers["nascent"] = X * 0.5
    ad.var["gene_length"] = np.random.randint(800, 4000, n_genes)
    ad.var["intron_number"] = np.random.randint(1, 9, n_genes)
    return ad


@pytest.fixture(scope="module")
def adata_mixed_small():
    """Small fixture for mixed-model + filter_active_genes tests."""
    np.random.seed(42)
    n_cells, n_genes = 60, 70
    X = np.random.negative_binomial(3, 0.5, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["Disease"] * 30 + ["Control"] * 30
    ad.obs["sample"] = ["s" + str(i % 6) for i in range(n_cells)]
    ad.layers["spliced"] = X.copy()
    ad.layers["unspliced"] = X * 0.45
    ad.var["gene_length"] = np.random.randint(600, 3500, n_genes)
    ad.var["intron_number"] = np.random.randint(0, 8, n_genes)
    return ad


@pytest.fixture(scope="module")
def adata_pb():
    """Pseudobulk-oriented fixture (4 samples per group)."""
    np.random.seed(11)
    n_cells, n_genes = 96, 80
    X = np.random.negative_binomial(6, 0.35, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["Disease"] * 48 + ["Control"] * 48
    ad.obs["sample"] = [f"S{i // 12}" for i in range(n_cells)]
    ad.layers["spliced"] = X.copy()
    ad.layers["unspliced"] = X * 0.5
    return ad


@pytest.fixture(scope="module")
def results_df(adata_basic):
    """Shared active_score table for plotting tests (one run per module)."""
    import scatrans as scat

    _, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_permutation=False,
        show_plot=False,
    )
    return allr


@pytest.fixture(scope="module")
def enrich_df():
    """Minimal enrichment result for plotting smoke tests."""
    return pd.DataFrame(
        {
            "Term": ["T1", "T2", "T3"],
            "Description": ["desc1", "desc2", "desc3"],
            "Count": [5, 3, 8],
            "GeneRatio": [0.1, 0.05, 0.2],
            "FoldEnrichment": [2.0, 1.5, 3.0],
            "pvalue": [0.001, 0.01, 0.0001],
            "p.adjust": [0.01, 0.05, 0.001],
            "neg_log10_padj": [2.0, 1.3, 3.0],
            "Genes": ["G1;G2", "G2;G3", "G4;G5;G6"],
        }
    )

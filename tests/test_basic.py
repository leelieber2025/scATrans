import pytest
import numpy as np
import scanpy as sc
import scatrans as scat

def make_test_adata():
    np.random.seed(42)
    n_cells, n_genes = 150, 300
    X = np.random.negative_binomial(4, 0.4, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs['condition'] = ['GA']*75 + ['Ctrl']*75
    ad.layers['spliced'] = X.copy()
    ad.layers['unspliced'] = X * 0.6
    ad.var['gene_length'] = np.random.randint(800, 4000, n_genes)
    ad.var['intron_number'] = np.random.randint(1, 8, n_genes)
    return ad

def test_heuristic_runs():
    ad = make_test_adata()
    res, sig, allr = scat.active_score(ad, mode='heuristic', show_plot=False, use_permutation=False)
    assert 'active_score' in res.var.columns

def test_advanced_runs_or_skips():
    ad = make_test_adata()
    try:
        res, sig, allr = scat.active_score(ad, mode='advanced', advanced_fallback=True, show_plot=False, use_permutation=False)
        assert res.uns['scatrans']['mode'] in ['advanced', 'heuristic_fallback_from_advanced']
    except ImportError:
        pytest.skip("scvelo not available")

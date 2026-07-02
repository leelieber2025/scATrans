"""Additional enrichment tests (run_go ALL, compare errors)."""

import pandas as pd
import pytest

import scatrans as scat


@pytest.mark.slow
def test_run_go_all_smoke():
    res = scat.run_go(
        ["ZZZ_not_a_gene"],
        ontology="ALL",
        organism="mouse",
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        adjust_across_all=True,
        verbose=False,
    )
    assert isinstance(res, pd.DataFrame)
    assert "per_ontology_attrs" in res.attrs


def test_compare_enrichment_raise_on_error():
    def _boom(*_args, **_kwargs):
        raise ValueError("simulated enrichment failure")

    with pytest.raises(ValueError, match="simulated enrichment failure"):
        scat.compare_enrichment(
            {"C1": ["GeneA"]},
            gene_sets={"T1": ["GeneA", "GeneB"]},
            fun=_boom,
            raise_on_error=True,
            verbose=False,
        )


def test_extract_gene_lists_empty_df():
    empty = pd.DataFrame(columns=["logFC", "p_adj"])
    out = scat.extract_gene_lists(empty, pval_cutoff=0.05)
    assert isinstance(out, dict)

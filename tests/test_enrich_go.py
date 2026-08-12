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


def test_run_go_adjust_across_all_honors_p_adjust_method(monkeypatch):
    """Cross-ontology re-adjust must use p_adjust_method (not hardcode BH)."""
    import numpy as np

    from scatrans.enrich import ora as ora_mod
    from scatrans.enrich._data import _apply_p_adjust

    def _fake_run_enrichment(*_args, **kwargs):
        ont = kwargs.get("gene_sets", "")
        # Distinct raw pvalues so BH vs Bonferroni differ
        if "Biological" in str(ont) or ont.endswith("BP") or "Process" in str(ont):
            pvals = [0.01, 0.02]
            terms = ["BP_t1", "BP_t2"]
            ont_tag = "BP"
        elif "Cellular" in str(ont) or "Component" in str(ont):
            pvals = [0.03]
            terms = ["CC_t1"]
            ont_tag = "CC"
        else:
            pvals = [0.04]
            terms = ["MF_t1"]
            ont_tag = "MF"
        df = pd.DataFrame(
            {
                "Term": terms,
                "pvalue": pvals,
                "p.adjust": list(pvals),  # dummy within-ontology
                "Count": [1] * len(pvals),
            }
        )
        df.attrs["ontology_tag"] = ont_tag
        return df

    monkeypatch.setattr(ora_mod, "run_enrichment", _fake_run_enrichment)
    res = scat.run_go(
        ["A", "B"],
        ontology="ALL",
        organism="mouse",
        return_all=True,
        adjust_across_all=True,
        p_adjust_method="bonferroni",
        verbose=False,
    )
    assert not res.empty
    assert res.attrs.get("p_adjust_method") == "bonferroni"
    assert "p.adjust.within_ontology" in res.columns
    expected = _apply_p_adjust(res["pvalue"].to_numpy(dtype=float), method="bonferroni")
    np.testing.assert_allclose(
        res["p.adjust"].to_numpy(dtype=float), expected, rtol=1e-10, atol=1e-12
    )


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

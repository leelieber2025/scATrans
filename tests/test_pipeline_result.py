"""Tests for PipelineResult and package layout after tl/enrich split."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
import pytest

import scatrans as scat


def _empty_result(**overrides):
    base = {
        "adata": None,
        "significant": pd.DataFrame(),
        "all_results": pd.DataFrame(),
        "candidates": pd.DataFrame({"g": [1]}),
        "enrichment": None,
        "filter_preset": "heuristic",
        "backend": {"use_pseudobulk": False},
        "meta": {"scatrans_version": "test"},
    }
    base.update(overrides)
    return scat.PipelineResult(**base)


def test_run_default_pipeline_returns_pipeline_result(adata_basic):
    result = scat.run_default_pipeline(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        run_go_enrichment=False,
        show_plot=False,
    )
    assert isinstance(result, scat.PipelineResult)
    assert isinstance(result, dict)
    assert isinstance(result, Mapping)
    assert "candidates" in result
    assert result["filter_preset"] in ("heuristic", "pseudobulk")
    assert isinstance(result.candidates, pd.DataFrame)
    assert isinstance(result.all_results, pd.DataFrame)
    assert result.enrichment is None
    summary = result.summary()
    assert summary["n_all_results"] == len(result.all_results)
    assert summary["filter_preset"] == result.filter_preset
    d = result.to_dict()
    assert set(d) >= {"adata", "candidates", "backend", "meta"}
    assert "scatrans_version" in result.meta
    assert result.meta["organism"] == "mouse"
    # Diagnostics promised by docstring / CHANGELOG must be folded into meta
    scatrans_uns = result.adata.uns.get("scatrans", {})
    assert isinstance(scatrans_uns, dict) and scatrans_uns
    assert "diagnostics" in result.meta
    assert result.meta["diagnostics"] is scatrans_uns["diagnostics"]
    assert "use_permutation" in result.meta
    assert result.meta["use_permutation"] == scatrans_uns.get("use_permutation")
    assert "gamma_method" in result.meta
    assert result.meta["gamma_method"] == scatrans_uns.get("gamma_method")


def test_pipeline_result_keyerror():
    pr = _empty_result(candidates=pd.DataFrame())
    with pytest.raises(KeyError):
        _ = pr["not_a_field"]
    assert pr.get("missing", 123) == 123


def test_pipeline_result_mapping_and_dict_protocol():
    """dict subclass: isinstance(dict), iteration, items, **unpacking."""
    pr = _empty_result()
    assert isinstance(pr, dict)
    assert isinstance(pr, Mapping)
    assert len(pr) == 8
    assert list(pr) == list(pr.keys())
    assert ("candidates", pr.candidates) in list(pr.items())
    assert dict(pr)["filter_preset"] == "heuristic"
    assert {**pr}["backend"] == {"use_pseudobulk": False}


def test_pipeline_result_is_readonly():
    pr = _empty_result()
    with pytest.raises(TypeError, match="read-only"):
        pr["filter_preset"] = "other"
    with pytest.raises(TypeError, match="read-only"):
        pr.candidates = pd.DataFrame()
    with pytest.raises(TypeError, match="read-only"):
        pr.update({"filter_preset": "x"})
    with pytest.raises(TypeError, match="read-only"):
        pr.pop("backend")
    with pytest.raises(TypeError, match="read-only"):
        pr.clear()
    # Nested plain dicts remain mutable (shallow freeze only)
    pr.backend["extra"] = 1
    assert pr.backend["extra"] == 1
    # Mutable shallow copy
    d = pr.to_dict()
    d["filter_preset"] = "changed"
    assert pr["filter_preset"] == "heuristic"


def test_pipeline_result_copy_and_deepcopy():
    """copy / deepcopy must not hit frozen __setitem__ during rebuild."""
    import copy
    import pickle

    pr = _empty_result(
        candidates=pd.DataFrame({"g": [1, 2]}),
        backend={"use_pseudobulk": False, "nested": {"a": 1}},
        meta={"scatrans_version": "test", "tags": ["x"]},
    )
    shallow = copy.copy(pr)
    assert isinstance(shallow, scat.PipelineResult)
    assert isinstance(shallow, dict)
    assert shallow is not pr
    assert shallow["filter_preset"] == pr["filter_preset"]
    # Shallow: nested dicts shared
    assert shallow.backend is pr.backend
    with pytest.raises(TypeError, match="read-only"):
        shallow["filter_preset"] = "nope"

    deep = copy.deepcopy(pr)
    assert isinstance(deep, scat.PipelineResult)
    assert deep is not pr
    assert deep.backend is not pr.backend
    assert deep.backend == pr.backend
    deep.backend["nested"]["a"] = 99
    assert pr.backend["nested"]["a"] == 1
    with pytest.raises(TypeError, match="read-only"):
        deep.update({"x": 1})

    # method .copy() is mutable dict (alias of to_dict); copy.copy is PipelineResult
    method_copy = pr.copy()
    assert type(method_copy) is dict
    assert not isinstance(method_copy, scat.PipelineResult)
    method_copy["filter_preset"] = "ok"
    assert pr["filter_preset"] == "heuristic"

    # Explicit reduce path used by pickle (and joblib-style dumps)
    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        roundtrip = pickle.loads(pickle.dumps(pr, protocol=protocol))
        assert isinstance(roundtrip, scat.PipelineResult), protocol
        assert roundtrip["filter_preset"] == "heuristic"
        assert list(roundtrip.candidates["g"]) == [1, 2]
        with pytest.raises(TypeError, match="read-only"):
            roundtrip["filter_preset"] = "x"

    # __reduce__ owns serialization; do not reintroduce dead __getstate__/__setstate__
    assert "__getstate__" not in scat.PipelineResult.__dict__
    assert "__setstate__" not in scat.PipelineResult.__dict__
    assert "__reduce__" in scat.PipelineResult.__dict__
    assert "__deepcopy__" in scat.PipelineResult.__dict__


def test_pipeline_result_ior_raises():
    """``|=`` must not silently mutate via C-level dict.__ior__."""
    pr = _empty_result()
    original_adata = pr["adata"]
    with pytest.raises(TypeError, match="read-only"):
        pr |= {"adata": "OVERWRITTEN"}
    assert pr["adata"] is original_adata
    assert pr["adata"] != "OVERWRITTEN"


def test_pipeline_result_or_returns_mutable_dict():
    """``|`` returns a new plain dict; left operand stays read-only."""
    pr = _empty_result()
    merged = pr | {"extra": 1, "filter_preset": "other"}
    assert type(merged) is dict
    assert not isinstance(merged, scat.PipelineResult)
    assert merged["extra"] == 1
    assert merged["filter_preset"] == "other"
    assert pr["filter_preset"] == "heuristic"
    with pytest.raises(TypeError, match="read-only"):
        pr["filter_preset"] = "nope"

    # reverse or
    merged2 = {"extra": 2} | pr
    assert type(merged2) is dict
    assert merged2["extra"] == 2
    assert merged2["filter_preset"] == "heuristic"


def test_tl_enrich_public_all_no_privates():
    import scatrans.enrich as enrich
    import scatrans.tl as tl

    for name in tl.__all__:
        assert not name.startswith("_"), name
    for name in enrich.__all__:
        assert not name.startswith("_"), name
    # dir() follows __all__
    assert "_materialize_if_view" not in dir(tl)
    assert "_clean_gene_list" not in dir(enrich)
    # Still reachable from implementation modules for advanced/tests use
    from scatrans.enrich._data import _clean_gene_list
    from scatrans.tl._common import _materialize_if_view

    assert callable(_materialize_if_view)
    assert callable(_clean_gene_list)


def test_tl_enrich_are_packages():
    import scatrans.enrich as enrich
    import scatrans.tl as tl

    assert hasattr(tl, "active_score")
    assert hasattr(tl, "PipelineResult")
    assert hasattr(enrich, "run_enrichment")
    import scatrans.enrich.compare as _c  # noqa: F401
    import scatrans.enrich.gsea as _g  # noqa: F401
    import scatrans.enrich.ora as _o  # noqa: F401
    import scatrans.tl.active as _a  # noqa: F401
    import scatrans.tl.de as _d  # noqa: F401
    import scatrans.tl.design as _des  # noqa: F401
    import scatrans.tl.filter as _f  # noqa: F401
    import scatrans.tl.pipeline as _p  # noqa: F401

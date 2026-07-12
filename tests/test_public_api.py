"""Smoke tests: every symbol in scatrans.__all__ is importable and minimally callable."""

import scatrans as scat


def test_all_exports_importable():
    missing = [name for name in scat.__all__ if not hasattr(scat, name)]
    assert not missing, f"Missing exports: {missing}"


def test_version_string():
    assert isinstance(scat.__version__, str)
    assert len(scat.__version__) > 0


def test_submodules_exposed():
    assert hasattr(scat, "pl")
    assert hasattr(scat, "qc")


def test_generate_gene_features_cli_entry_importable():
    """CLI entry lives on the submodule, not the top-level public surface."""
    from scatrans.generate_gene_features import main

    assert callable(main)


def test_pipeline_result_exported():
    assert hasattr(scat, "PipelineResult")
    assert "PipelineResult" in scat.__all__


def test_pl_public_dir_no_typing_leaks():
    names = set(dir(scat.pl))
    for leak in ("Any", "Optional", "Union", "Iterable", "Mapping", "Normalize"):
        assert leak not in names, f"pl dir() leaks {leak}"


def test_qc_public_dir_minimal():
    names = set(dir(scat.qc))
    assert "unspliced_global" in names
    assert "logger" not in names
    assert "logging" not in names


def test_list_available_gene_features():
    feats = scat.list_available_gene_features(verbose=False)
    assert isinstance(feats, list)


def test_compare_enrichment_callable():
    assert callable(scat.compare_enrichment)


def test_extract_gene_lists_callable():
    assert callable(scat.extract_gene_lists)


def test_concat_compare_results_callable():
    assert callable(scat.concat_compare_results)


def test_public_callables_have_docstrings():
    """Exported callables (not modules) should document user-facing behavior."""
    skip = {"pl", "qc", "__version__", "PipelineResult", "WORKFLOW_PRESETS"}
    for name in scat.__all__:
        if name in skip:
            continue
        obj = getattr(scat, name)
        if callable(obj) and not isinstance(obj, type):
            assert obj.__doc__, f"{name} missing docstring"

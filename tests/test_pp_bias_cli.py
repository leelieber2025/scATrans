"""Tests for gene feature attachment and CLI entry point."""

import sys

import pytest

import scatrans as scat


def test_add_gene_features_custom_path(adata_basic, tmp_path):
    import pandas as pd

    feats = pd.DataFrame(
        {
            "gene_name": adata_basic.var_names[:10],
            "gene_length": range(1000, 1010),
            "intron_number": range(1, 11),
        }
    )
    path = tmp_path / "custom_feats.parquet"
    feats.to_parquet(path)
    out = scat.add_gene_features(adata_basic.copy(), gene_features_path=str(path))
    assert out.var["gene_length"].notna().sum() >= 1


def test_add_gene_features_human_organism(adata_basic):
    out = scat.add_gene_features(adata_basic.copy(), organism="human")
    assert "gene_length" in out.var.columns


def test_generate_gene_features_cli_help():
    from scatrans.generate_gene_features import main

    with pytest.raises(SystemExit) as exc:
        sys.argv = ["generate-gene-features", "--help"]
        main()
    assert exc.value.code == 0

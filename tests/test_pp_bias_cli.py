"""Tests for gene feature attachment and CLI entry point."""

import sys

import pytest

import scatrans as scat


def test_generate_gene_features_no_exon_rows(tmp_path):
    pytest.importorskip("gtfparse")
    from scatrans.pp_bias import generate_gene_features_from_gtf

    gtf = tmp_path / "no_exon.gtf"
    rows = [
        'chr1\t.\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "OnlyGene"; gene_type "protein_coding";',
    ]
    gtf.write_text("\n".join(rows) + "\n")
    out = tmp_path / "features.parquet"
    gene_df = generate_gene_features_from_gtf(str(gtf), str(out))
    assert len(gene_df) == 1
    assert gene_df.iloc[0]["gene_length"] == 0


def test_generate_gene_features_dedup_gene_names(tmp_path):
    pytest.importorskip("gtfparse")
    from scatrans.pp_bias import generate_gene_features_from_gtf

    gtf = tmp_path / "dup.gtf"
    rows = [
        'chr1\t.\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "DupGene"; gene_type "protein_coding";',
        'chr1\t.\tgene\t300\t400\t.\t+\t.\tgene_id "G2"; gene_name "DupGene"; gene_type "protein_coding";',
        'chr1\t.\texon\t100\t150\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
        'chr1\t.\texon\t300\t350\t.\t+\t.\tgene_id "G2"; transcript_id "T2";',
    ]
    gtf.write_text("\n".join(rows) + "\n")
    out = tmp_path / "features.parquet"
    gene_df = generate_gene_features_from_gtf(str(gtf), str(out))
    assert len(gene_df) == 1
    assert gene_df.iloc[0]["gene_name"] == "DupGene"


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

"""
scATrans internal differential expression helpers.

Contains the wrapper that supports both scanpy rank_genes_groups and PyDESeq2
for pseudobulk. Extracted so tl.py stays small.
"""

from __future__ import annotations

import logging
import warnings
from typing import Any, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from ._utils import (
    _is_integer_counts_like,
    _warn_if_low_counts_matrix,
)

logger = logging.getLogger(__name__)


def _run_de_wrapper(
    adata: ad.AnnData,
    groupby: str,
    target_group: str,
    reference_group: str,
    de_method: str = "t-test_overestim_var",
    is_pseudobulk: bool = False,
    pb_backend: str = "pydeseq2",
    n_jobs: int = 1,
    labels: Optional[Any] = None,
    strict_pydeseq2_counts: bool = True,
) -> pd.DataFrame:
    """Run DE and return a minimal DataFrame with logFC, p_val, p_adj indexed by var_names."""
    ad_temp = adata.copy() if labels is not None else adata
    use_groupby = groupby

    if labels is not None:
        use_groupby = "_de_temp_group"
        ad_temp.obs[use_groupby] = pd.Categorical(
            np.asarray(labels).astype(str), categories=[reference_group, target_group]
        )

    if is_pseudobulk and pb_backend == "pydeseq2":
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
        except ImportError as e:
            raise ImportError(
                "pydeseq2 is required when pseudobulk_de_backend='pydeseq2'. "
                "Install with: pip install pydeseq2 or 'scatrans[pseudobulk]'"
            ) from e

        n_t = (ad_temp.obs[use_groupby] == target_group).sum()
        n_r = (ad_temp.obs[use_groupby] == reference_group).sum()
        if n_t < 2 or n_r < 2:
            raise ValueError(
                f"PyDESeq2 requires >=2 replicates per group. Found {n_t} target, {n_r} ref."
            )

        is_count_like = _is_integer_counts_like(ad_temp.X)

        if not is_count_like:
            msg = (
                "Data passed to PyDESeq2 does not look like raw non-negative integer counts. "
                "PyDESeq2 requires unnormalized integer counts in adata.X. "
                "If you intentionally want to allow rounding, set strict_pydeseq2_counts=False."
            )
            if strict_pydeseq2_counts:
                raise ValueError(msg)
            logger.warning(msg)
        else:
            _warn_if_low_counts_matrix(ad_temp.X)

        if sparse.issparse(ad_temp.X):
            gene_sums = np.asarray(ad_temp.X.sum(axis=0)).ravel()
            gene_keep = gene_sums >= 10
            if gene_keep.sum() == 0:
                raise ValueError("No genes passed the DESeq2 count filter (sum(counts) >= 10).")
            X_filtered = ad_temp.X[:, gene_keep].toarray()
            X_filtered = np.clip(np.round(np.nan_to_num(X_filtered)), 0, None).astype(int)
            counts_use = pd.DataFrame(
                X_filtered, index=ad_temp.obs_names, columns=ad_temp.var_names[gene_keep]
            )
        else:
            X = np.asarray(ad_temp.X)
            X = np.clip(np.round(np.nan_to_num(X)), 0, None).astype(int)
            counts_df = pd.DataFrame(X, index=ad_temp.obs_names, columns=ad_temp.var_names)
            gene_keep = counts_df.sum(axis=0) >= 10
            counts_use = counts_df.loc[:, gene_keep].copy()

        if counts_use.shape[1] == 0:
            raise ValueError("No genes passed the DESeq2 count filter (sum(counts) >= 10).")

        condition = ad_temp.obs[use_groupby].astype(str).values
        metadata = pd.DataFrame(
            {use_groupby: pd.Categorical(condition, categories=[reference_group, target_group])},
            index=counts_use.index,
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                dds = DeseqDataSet(
                    counts=counts_use,
                    metadata=metadata,
                    design_factors=use_groupby,
                    ref_level=[use_groupby, reference_group],
                    quiet=True,
                    n_cpus=n_jobs,
                )
            except TypeError:
                dds = DeseqDataSet(
                    counts=counts_use,
                    metadata=metadata,
                    design=f"~{use_groupby}",
                    refit_cooks=True,
                    quiet=True,
                    n_cpus=n_jobs,
                )
            dds.deseq2()

            try:
                stat_res = DeseqStats(
                    dds,
                    contrast=[use_groupby, target_group, reference_group],
                    quiet=True,
                    n_cpus=n_jobs,
                )
            except TypeError:
                stat_res = DeseqStats(
                    dds,
                    contrast=[use_groupby, target_group, reference_group],
                    n_cpus=n_jobs,
                )
            stat_res.summary()

        res2 = stat_res.results_df.copy().reindex(ad_temp.var_names)
        de_df = pd.DataFrame(index=ad_temp.var_names)
        de_df["logFC"] = res2["log2FoldChange"].fillna(0.0)
        de_df["p_val"] = res2.get("pvalue", pd.Series(1.0, index=res2.index)).fillna(1.0)
        de_df["p_adj"] = res2.get("padj", pd.Series(1.0, index=res2.index)).fillna(1.0)
        return de_df

    else:
        # Standard scanpy path (works for both regular and pseudobulk when not using pydeseq2)
        rank_key = "_scatrans_rank_genes_groups"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.rank_genes_groups(
                ad_temp,
                groupby=use_groupby,
                groups=[target_group],
                reference=reference_group,
                method=de_method,
                key_added=rank_key,
            )
        de_raw = sc.get.rank_genes_groups_df(ad_temp, group=target_group, key=rank_key).set_index(
            "names"
        )
        de_df = pd.DataFrame(index=ad_temp.var_names)
        de_df["logFC"] = de_raw["logfoldchanges"].reindex(ad_temp.var_names).fillna(0.0)
        de_df["p_val"] = de_raw["pvals"].reindex(ad_temp.var_names).fillna(1.0)
        de_df["p_adj"] = de_raw["pvals_adj"].reindex(ad_temp.var_names).fillna(1.0)
        return de_df

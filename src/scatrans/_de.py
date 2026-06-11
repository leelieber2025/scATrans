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
from joblib import Parallel, delayed
from scipy import sparse
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

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
    use_mixed_model: bool = False,
    sample_col: Optional[str] = None,
    mixed_model_pval: str = "wald",
    # Memento (Cell 2024 method-of-moments) as independent cell-level DE backend
    use_memento_de: bool = False,
    memento_capture_rate: float = 0.07,
    memento_num_boot: int = 5000,
    memento_n_cpus: int = -1,
) -> pd.DataFrame:
    """Run DE and return a DataFrame with logFC, p_val, p_adj (and optionally delta_variance, delta_var_pval when mixed).

    When use_memento_de=True, Memento is used for the primary DE statistics (logFC/p_adj).
    This is treated as a third parallel cell-level backend (alongside scanpy-style and mixed-model).
    """
    if use_mixed_model:
        if sample_col is None:
            raise ValueError("sample_col must be provided when use_mixed_model=True")
        return _run_mixedlm_de(
            adata,
            groupby=groupby,
            target_group=target_group,
            reference_group=reference_group,
            sample_col=sample_col,
            n_jobs=n_jobs,
            labels=labels,
            mixed_model_pval=mixed_model_pval,
        )

    if use_memento_de:
        if is_pseudobulk:
            raise ValueError(
                "use_memento_de=True is not supported with use_pseudobulk=True "
                "(Memento is a cell-level method-of-moments estimator; use PyDESeq2 for pseudobulk)."
            )
        return _run_memento_de(
            adata,
            groupby=groupby,
            target_group=target_group,
            reference_group=reference_group,
            labels=labels,
            capture_rate=memento_capture_rate,
            num_boot=memento_num_boot,
            n_cpus=memento_n_cpus,
        )

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


def _run_mixedlm_de(
    adata: ad.AnnData,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str,
    n_jobs: int = 1,
    labels: Optional[Any] = None,
    mixed_model_pval: str = "wald",
) -> pd.DataFrame:
    """
    Mixed linear model (LMM) DE + Delta Variance using statsmodels mixedlm.

    Models: y_log ~ C(condition) + (1 | sample)
    - logFC: coefficient for the target condition (on log1p scale)
    - p_val / p_adj: Wald p for the condition fixed effect (BH adj across genes)
    - delta_variance: fraction of total variance (var_fe + re_var + resid) attributable to the fixed condition effect.
    - delta_var_pval: LRT p-value for the contribution of condition (full vs reduced ~1 + (1|sample))

    This provides a lightweight Python analogue to variancePartition/dream (fraction of variation explained)
    + LMM DE, suitable for cell-level data with sample-level random effects (addresses pseudoreplication).
    For full voom + precision weights + dreampy/dreamlet on pseudobulk, or NEBULA NB-GLMM, see external packages.
    """
    try:
        import statsmodels.formula.api as smf
    except ImportError as e:
        raise ImportError(
            "statsmodels is required for use_mixed_model=True (it is a core dependency of scatrans)."
        ) from e

    ad_temp = adata.copy() if labels is not None else adata
    use_groupby = groupby
    if labels is not None:
        use_groupby = "_de_temp_group"
        ad_temp.obs[use_groupby] = pd.Categorical(
            np.asarray(labels).astype(str), categories=[reference_group, target_group]
        )

    if sample_col not in ad_temp.obs.columns:
        raise ValueError(f"sample_col='{sample_col}' not found in adata.obs")

    # Prepare expression for LMM: always use log1p library-size normalized (LMM assumes approx Gaussian on log scale)
    # Work on a temp copy to avoid mutating caller adata state for this auxiliary norm
    ad_expr = ad_temp.copy()
    # If very sparse or raw counts, normalize; safe for already-log too (will just re-log)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.normalize_total(ad_expr, target_sum=1e4)
            sc.pp.log1p(ad_expr)
    except Exception:
        # Fallback: manual log1p on X if pp fails (e.g. already transformed or edge data)
        if sparse.issparse(ad_expr.X):
            Xn = np.log1p(ad_expr.X.toarray())
        else:
            Xn = np.log1p(np.asarray(ad_expr.X))
        ad_expr.X = Xn

    if sparse.issparse(ad_expr.X):
        expr_mat = ad_expr.X.toarray()
    else:
        expr_mat = np.asarray(ad_expr.X)

    obs = ad_temp.obs
    condition = obs[use_groupby].astype(str).values
    samples = obs[sample_col].astype(str).values

    n_genes = expr_mat.shape[1]
    var_names = ad_temp.var_names

    # Per-gene worker (returns idx, logfc, wald_p, lrt_p, delta_var)
    def _fit_gene_mixed(idx: int):
        y = expr_mat[:, idx].astype(float)
        # guard against all-zero / constant (mixedlm will be singular)
        if np.allclose(y, y[0]):
            return idx, 0.0, 1.0, 1.0, 0.0
        df = pd.DataFrame({"y": y, "condition": condition, "sample": samples})
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # Full model
                md_full = smf.mixedlm("y ~ C(condition)", df, groups=df["sample"])
                m_full = md_full.fit(reml=False, maxiter=200, disp=False)

                # Reduced (null) for LRT on condition contribution
                md_null = smf.mixedlm("y ~ 1", df, groups=df["sample"])
                m_null = md_null.fit(reml=False, maxiter=200, disp=False)

            # LRT statistic and p (chi2 df=1 for the added fixed effect term(s))
            lrt_stat = -2.0 * (m_null.llf - m_full.llf)
            lrt_p = float(chi2.sf(max(lrt_stat, 0.0), 1))

            # Extract condition coef (target vs ref)
            # The param name is typically "C(condition)[T.<target>]" or similar
            coef_name = None
            for pname in m_full.params.index:
                if "condition" in str(pname) and target_group in str(pname):
                    coef_name = pname
                    break
            if coef_name is None:
                # fallback: take the second coef if intercept + one more
                if len(m_full.params) >= 2:
                    coef_name = m_full.params.index[1]
                else:
                    coef_name = m_full.params.index[0]
            logfc = float(m_full.params.get(coef_name, 0.0))
            p_wald = float(m_full.pvalues.get(coef_name, 1.0))

            # Delta variance: var attributable to fixed effects / total modeled var
            exog = m_full.model.exog
            beta = np.asarray(m_full.fe_params)
            fe_contrib = exog @ beta
            var_fe = float(np.var(fe_contrib))
            re_var = 0.0
            try:
                if hasattr(m_full, "cov_re") and m_full.cov_re is not None and len(m_full.cov_re) > 0:
                    re_var = float(np.diag(m_full.cov_re)[0])  # first (only) RE variance
            except Exception:
                re_var = 0.0
            resid_var = float(getattr(m_full, "scale", 0.0))
            total_v = var_fe + max(re_var, 0.0) + max(resid_var, 0.0)
            delta_var = var_fe / total_v if total_v > 1e-12 else 0.0

            return idx, logfc, p_wald, lrt_p, float(np.clip(delta_var, 0.0, 1.0))
        except Exception:
            # Degenerate fit (few samples per group, collinear, etc.) -> non-informative
            return idx, 0.0, 1.0, 1.0, 0.0

    # Parallel execution (loky or threading; mixedlm releases GIL-ish via numpy)
    effective_jobs = max(1, n_jobs) if n_jobs and n_jobs > 0 else 1

    results = Parallel(n_jobs=effective_jobs, backend="loky")(
        delayed(_fit_gene_mixed)(i) for i in range(n_genes)
    )

    # Assemble
    results = sorted(results, key=lambda t: t[0])
    logfcs = np.array([r[1] for r in results], dtype=float)
    p_walds = np.array([r[2] for r in results], dtype=float)
    p_lrts = np.array([r[3] for r in results], dtype=float)
    dvars = np.array([r[4] for r in results], dtype=float)

    # Choose which p-value to expose as the main "p_val" for active_score weighting and default filtering.
    # "wald": the coefficient test (standard for logFC-like effect)
    # "lrt": the likelihood ratio test for the condition term contribution (ties directly to delta_variance)
    if mixed_model_pval == "lrt":
        main_pvals = p_lrts
    else:
        if mixed_model_pval != "wald":
            logger.warning("mixed_model_pval must be 'wald' or 'lrt'; falling back to 'wald'.")
        main_pvals = p_walds

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p_adjs = multipletests(main_pvals, method="fdr_bh")[1]

    de_df = pd.DataFrame(index=var_names)
    de_df["logFC"] = pd.Series(logfcs, index=var_names)
    de_df["p_val"] = pd.Series(main_pvals, index=var_names)
    de_df["p_adj"] = pd.Series(p_adjs, index=var_names)
    de_df["delta_variance"] = pd.Series(dvars, index=var_names)
    de_df["delta_var_pval"] = pd.Series(p_lrts, index=var_names)
    return de_df


def _run_memento_de(
    adata: ad.AnnData,
    groupby: str,
    target_group: str,
    reference_group: str,
    labels: Optional[Any] = None,
    capture_rate: float = 0.07,
    num_boot: int = 5000,
    n_cpus: int = -1,
) -> pd.DataFrame:
    """Memento (method of moments) cell-level DE backend.

    Returns a DataFrame with at minimum 'logFC', 'p_val', 'p_adj' (plus optional
    memento_de_* and memento_dv_* columns for advanced inspection).
    This replaces the scanpy rank_genes_groups path when use_memento_de=True.
    """
    try:
        import memento
    except ImportError as e:
        raise ImportError(
            "memento-de is required when use_memento_de=True. "
            'Install with: pip install "scatrans[memento]" (or pip install memento-de)'
        ) from e

    ad_temp = adata.copy() if labels is not None else adata
    use_groupby = groupby

    if labels is not None:
        use_groupby = "_memento_temp_group"
        ad_temp.obs[use_groupby] = pd.Categorical(
            np.asarray(labels).astype(str), categories=[reference_group, target_group]
        )

    # Restrict to the two groups being compared
    keep = ad_temp.obs[use_groupby].astype(str).isin([target_group, reference_group])
    ad_temp = ad_temp[keep].copy()

    # Binary treatment column expected by memento.binary_test_1d
    ad_temp.obs["stim"] = (ad_temp.obs[use_groupby].astype(str) == target_group).astype(int)

    # Memento expects (or strongly prefers) raw non-negative integer counts in .X
    if not _is_integer_counts_like(ad_temp.X):
        if "counts" in getattr(ad_temp, "layers", {}):
            ad_temp.X = ad_temp.layers["counts"].copy()
            logger.info("Memento: using 'counts' layer as raw counts input.")
        else:
            logger.warning(
                "Memento input does not look like raw integer counts and no 'counts' layer was found. "
                "Memento (method-of-moments + hypergeometric model) performs best on raw counts; "
                "results may be unreliable. Provide a 'counts' layer or ensure adata.X contains unnormalized counts."
            )

    # Memento (via its wrappers) requires scipy CSR for .X
    from scipy.sparse import csr_matrix, issparse
    if not (issparse(ad_temp.X) and type(ad_temp.X) == csr_matrix):
        ad_temp.X = csr_matrix(ad_temp.X)
        logger.debug("Memento: converted .X to CSR matrix.")

    # Effective cpus
    effective_cpus = n_cpus if n_cpus and n_cpus > 0 else -1

    result = memento.binary_test_1d(
        adata=ad_temp,
        treatment_col="stim",
        capture_rate=capture_rate,
        num_boot=num_boot,
        num_cpus=effective_cpus,
    )

    # result may be indexed by gene or have a 'gene' column (handle both)
    if isinstance(result, pd.DataFrame):
        if "gene" in result.columns:
            result = result.set_index("gene")
        res_index = result.index
    else:
        # Fallback (should not happen)
        res_index = ad_temp.var_names
        result = pd.DataFrame(index=res_index)

    de_df = pd.DataFrame(index=res_index)
    de_df["logFC"] = pd.to_numeric(result.get("de_coef", 0.0), errors="coerce").reindex(res_index).fillna(0.0)
    pvals = pd.to_numeric(result.get("de_pval", 1.0), errors="coerce").reindex(res_index).fillna(1.0)
    de_df["p_val"] = pvals

    # BH adjustment (Memento may return raw p; we make p_adj consistent with other backends)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        de_df["p_adj"] = multipletests(pvals.values, method="fdr_bh")[1]

    # Expose Memento's native columns for users who want mean + variability signals
    for src, dst in [
        ("de_se", "memento_de_se"),
        ("dv_coef", "memento_dv_coef"),
        ("dv_se", "memento_dv_se"),
        ("dv_pval", "memento_dv_pval"),
    ]:
        if src in result.columns:
            de_df[dst] = pd.to_numeric(result[src], errors="coerce").reindex(res_index)

    # Re-align to the var_names of the adata object that was passed into the wrapper
    # (important for the labels= permutation case and any internal subsetting)
    de_df = de_df.reindex(adata.var_names).fillna({"logFC": 0.0, "p_val": 1.0, "p_adj": 1.0})

    return de_df

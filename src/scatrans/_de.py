"""
scATrans internal differential expression helpers.

Contains the wrapper that supports both scanpy rank_genes_groups and PyDESeq2
for pseudobulk. Extracted so tl.py stays small.
"""

from __future__ import annotations

import logging
import warnings
from contextlib import contextmanager
from importlib.metadata import version
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

from ._utils import (
    _is_integer_counts_like,
    _prepare_log_normalized_expression,
    _warn_if_low_counts_matrix,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

_DE_REQUIRED_COLS = frozenset({"logFC", "p_val", "p_adj"})

# Suppress noisy deprecation/future warnings from scanpy/anndata/PyDESeq2 during DE,
# but keep convergence/runtime/user warnings visible to callers.
_DE_SUPPRESS_WARNING_CATEGORIES = (
    DeprecationWarning,
    FutureWarning,
    PendingDeprecationWarning,
)


@contextmanager
def _de_warning_context():
    with warnings.catch_warnings():
        for cat in _DE_SUPPRESS_WARNING_CATEGORIES:
            warnings.simplefilter("ignore", category=cat)
        yield


def _validate_de_result(de_df: pd.DataFrame, *, backend: str) -> pd.DataFrame:
    """Assert all DE backends return the minimal schema expected downstream."""
    missing = _DE_REQUIRED_COLS - set(de_df.columns)
    if missing:
        raise RuntimeError(
            f"DE backend '{backend}' returned incomplete results; "
            f"missing columns: {sorted(missing)}"
        )
    if len(de_df) == 0:
        return de_df
    for col in _DE_REQUIRED_COLS:
        vals = pd.to_numeric(de_df[col], errors="coerce")
        if not np.isfinite(vals.to_numpy()).any():
            raise RuntimeError(
                f"DE backend '{backend}' returned no finite values in column {col!r}."
            )
    return de_df


def _pydeseq2_uses_design_factors() -> bool:
    """Detect whether the installed PyDESeq2 supports the modern design_factors= API.

    Uses version parsing instead of blind except TypeError to avoid swallowing
    unrelated future errors.
    """
    try:
        vstr = version("pydeseq2")
    except Exception:
        # If we cannot determine (e.g. not installed or metadata issue), prefer the
        # modern path and let the try/except in caller be a last resort.
        return True
    # pydeseq2 0.4+ standardized on design_factors; earlier used design=
    # Be defensive: parse major.minor
    try:
        parts = vstr.split(".")
        major = int(parts[0]) if parts else 0
        minor = int(parts[1]) if len(parts) > 1 else 0
        return (major, minor) >= (0, 4)
    except Exception:
        return True


def _run_de_wrapper(
    adata: ad.AnnData,
    groupby: str,
    target_group: str,
    reference_group: str,
    de_method: str = "t-test_overestim_var",
    is_pseudobulk: bool = False,
    pb_backend: str = "pydeseq2",
    n_jobs: int = 1,
    labels: Any | None = None,
    strict_pydeseq2_counts: bool = True,
    use_mixed_model: bool = False,
    sample_col: str | None = None,
    mixed_model_pval: str = "wald",
    # Memento (Cell 2024 method-of-moments) as independent cell-level DE backend
    use_memento_de: bool = False,
    memento_capture_rate: float = 0.07,
    memento_num_boot: int = 5000,
    memento_n_cpus: int = -1,
    # Allow providing raw counts separately (common when adata.X is already HVG + log1p)
    counts: str | np.ndarray | sparse.spmatrix | pd.DataFrame | ad.AnnData | None = None,
    # Minimum total counts (across observations) to keep a gene for PyDESeq2.
    # The original hard-coded 10 can be too aggressive for very small pseudobulk
    # (few samples) or low-depth data, and too lenient for huge datasets.
    min_counts_per_gene: int = 10,
) -> pd.DataFrame:
    """Run DE and return a DataFrame with logFC, p_val, p_adj (and optionally delta_variance, delta_var_pval when mixed).

    When use_memento_de=True, Memento is used for the primary DE statistics (logFC/p_adj).
    This is treated as a third parallel cell-level backend (alongside scanpy-style and mixed-model).

    logFC is normalized toward log2 scale for cross-backend comparability of logfc_cutoff:
      - PyDESeq2 + scanpy backends (wilcoxon, t-test, logreg, etc.): native log2
      - mixedlm + memento: converted from natural/log1p scale (/ log(2))

    Internal function: type hints strengthened for mypy/pyright.
    """
    if use_mixed_model:
        if sample_col is None:
            raise ValueError("sample_col must be provided when use_mixed_model=True")
        return _validate_de_result(
            _run_mixedlm_de(
                adata,
                groupby=groupby,
                target_group=target_group,
                reference_group=reference_group,
                sample_col=sample_col,
                n_jobs=n_jobs,
                labels=labels,
                mixed_model_pval=mixed_model_pval,
            ),
            backend="mixedlm",
        )

    if use_memento_de:
        if is_pseudobulk:
            raise ValueError(
                "use_memento_de=True is not supported with use_pseudobulk=True "
                "(Memento is a cell-level method-of-moments estimator; use PyDESeq2 for pseudobulk)."
            )
        return _validate_de_result(
            _run_memento_de(
                adata,
                groupby=groupby,
                target_group=target_group,
                reference_group=reference_group,
                labels=labels,
                capture_rate=memento_capture_rate,
                num_boot=memento_num_boot,
                n_cpus=memento_n_cpus,
                counts=counts,
            ),
            backend="memento",
        )

    ad_temp = adata.copy() if labels is not None else adata
    use_groupby = groupby

    if labels is not None:
        use_groupby = "_de_temp_group"
        # Quiet anndata "storing ... as categorical" note during internal label injection (perm + shuffle)
        _ann_log = logging.getLogger("anndata")
        _prev_ann = _ann_log.level
        _ann_log.setLevel(logging.WARNING)
        try:
            ad_temp.obs[use_groupby] = pd.Categorical(
                np.asarray(labels).astype(str), categories=[reference_group, target_group]
            )
        finally:
            _ann_log.setLevel(_prev_ann)

    if is_pseudobulk and pb_backend == "pydeseq2":
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
        except ImportError as e:
            raise ImportError(
                "PyDESeq2 backend requested but 'pydeseq2' is not installed.\n"
                "Install with:\n"
                '    pip install "scatrans[pseudobulk]"\n'
                "or\n"
                "    pip install pydeseq2"
            ) from e

        n_t = (ad_temp.obs[use_groupby] == target_group).sum()
        n_r = (ad_temp.obs[use_groupby] == reference_group).sum()
        if n_t < 2 or n_r < 2:
            raise ValueError(
                f"PyDESeq2 requires >=2 replicates per group. Found {n_t} target, {n_r} ref."
            )

        # For pseudobulk the X is an aggregation (sum) of counts; values may arrive as float64
        # but are integer-valued. Use a tolerant check + pre-coercion so users don't have to set
        # strict_pydeseq2_counts=False for normal aggregated data.
        if is_pseudobulk:
            # Coerce a tolerant view for the "looks like counts" decision
            try:
                if sparse.issparse(ad_temp.X):
                    X_check_arr = np.asarray(ad_temp.X.todense())
                else:
                    X_check_arr = np.asarray(ad_temp.X)
                X_check_arr = np.clip(np.round(np.nan_to_num(X_check_arr)), 0, None)
                is_count_like = _is_integer_counts_like(X_check_arr)
            except Exception as _e:
                logger.debug("Count-like check fallback: %s", _e)
                is_count_like = _is_integer_counts_like(ad_temp.X)
        else:
            is_count_like = _is_integer_counts_like(ad_temp.X)

        if not is_count_like:
            msg = (
                "Data passed to PyDESeq2 does not look like raw non-negative integer counts. "
                "PyDESeq2 requires unnormalized integer counts in adata.X. "
                "For pseudobulk data we automatically round to integer counts; "
                "set strict_pydeseq2_counts=False to allow the (rounded) data anyway."
            )
            if strict_pydeseq2_counts:
                raise ValueError(msg)
            logger.warning(msg)
        else:
            _warn_if_low_counts_matrix(ad_temp.X)

        if sparse.issparse(ad_temp.X):
            gene_sums = np.asarray(ad_temp.X.sum(axis=0)).ravel()
            gene_keep = gene_sums >= min_counts_per_gene
            if gene_keep.sum() == 0:
                raise ValueError(
                    f"No genes passed the DESeq2 count filter (sum(counts) >= {min_counts_per_gene})."
                )
            X_filtered = ad_temp.X[:, gene_keep].toarray()
            X_filtered = np.clip(np.round(np.nan_to_num(X_filtered)), 0, None).astype(int)
            counts_use = pd.DataFrame(
                X_filtered, index=ad_temp.obs_names, columns=ad_temp.var_names[gene_keep]
            )
        else:
            X = np.asarray(ad_temp.X)
            X = np.clip(np.round(np.nan_to_num(X)), 0, None).astype(int)
            counts_df = pd.DataFrame(X, index=ad_temp.obs_names, columns=ad_temp.var_names)
            gene_keep = counts_df.sum(axis=0) >= min_counts_per_gene
            counts_use = counts_df.loc[:, gene_keep].copy()

        if counts_use.shape[1] == 0:
            raise ValueError(
                f"No genes passed the DESeq2 count filter (sum(counts) >= {min_counts_per_gene})."
            )

        condition = ad_temp.obs[use_groupby].astype(str).values
        metadata = pd.DataFrame(
            {use_groupby: pd.Categorical(condition, categories=[reference_group, target_group])},
            index=counts_use.index,
        )

        with _de_warning_context():
            if _pydeseq2_uses_design_factors():
                dds = DeseqDataSet(
                    counts=counts_use,
                    metadata=metadata,
                    design_factors=use_groupby,
                    ref_level=[use_groupby, reference_group],
                    quiet=True,
                    n_cpus=n_jobs,
                )
            else:
                dds = DeseqDataSet(
                    counts=counts_use,
                    metadata=metadata,
                    design=f"~{use_groupby}",
                    refit_cooks=True,
                    quiet=True,
                    n_cpus=n_jobs,
                )
            dds.deseq2()

            if _pydeseq2_uses_design_factors():
                stat_res = DeseqStats(
                    dds,
                    contrast=[use_groupby, target_group, reference_group],
                    quiet=True,
                    n_cpus=n_jobs,
                )
            else:
                stat_res = DeseqStats(
                    dds,
                    contrast=[use_groupby, target_group, reference_group],
                    n_cpus=n_jobs,
                )
            stat_res.summary()

        res_df = stat_res.results_df.copy()
        n_genes_filtered_low_count = int((~gene_keep).sum())
        n_genes_nan_from_deseq2 = (
            int(res_df["padj"].isna().sum()) if "padj" in res_df.columns else 0
        )
        res2 = res_df.reindex(ad_temp.var_names)
        de_df = pd.DataFrame(index=ad_temp.var_names)
        de_df["logFC"] = res2["log2FoldChange"].reindex(ad_temp.var_names)
        de_df["p_val"] = res2.get("pvalue", pd.Series(np.nan, index=res2.index)).reindex(
            ad_temp.var_names
        )
        de_df["p_adj"] = res2.get("padj", pd.Series(np.nan, index=res2.index)).reindex(
            ad_temp.var_names
        )
        _validate_de_result(de_df, backend="pydeseq2")
        de_df["logFC"] = de_df["logFC"].fillna(0.0)
        de_df["p_val"] = de_df["p_val"].fillna(1.0)
        de_df["p_adj"] = de_df["p_adj"].fillna(1.0)
        n_total = len(de_df)
        de_df.attrs["n_genes_filtered_low_count"] = n_genes_filtered_low_count
        de_df.attrs["n_genes_nan_from_deseq2"] = n_genes_nan_from_deseq2
        de_df.attrs["pydeseq2_neutral_fill"] = True
        if n_total > 0 and (n_genes_filtered_low_count > 0 or n_genes_nan_from_deseq2 > 0):
            logger.warning(
                "PyDESeq2: %d/%d genes (%.1f%%) skipped by min_counts_per_gene filter; "
                "%d/%d genes (%.1f%%) had NaN padj from DESeq2 independent filtering/outliers. "
                "These genes appear as neutral values (logFC=0, p_adj=1) and are NOT "
                "'tested and non-significant'. "
                "See de_df.attrs['n_genes_filtered_low_count'] and ['n_genes_nan_from_deseq2'].",
                n_genes_filtered_low_count,
                n_total,
                100.0 * n_genes_filtered_low_count / n_total,
                n_genes_nan_from_deseq2,
                n_total,
                100.0 * n_genes_nan_from_deseq2 / n_total,
            )
        return de_df

    else:
        # Standard scanpy path (works for both regular and pseudobulk when not using pydeseq2)
        rank_key = "_scatrans_rank_genes_groups"
        with _de_warning_context():
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
        # scanpy rank_genes_groups always returns logfoldchanges on log2 scale:
        # log2( (expm1(mean_t) + eps) / (expm1(mean_r) + eps) ), independent of the
        # statistical method (wilcoxon / t-test / etc.). No secondary conversion needed.
        de_df["logFC"] = de_raw["logfoldchanges"].reindex(ad_temp.var_names)
        de_df["p_val"] = de_raw["pvals"].reindex(ad_temp.var_names)
        de_df["p_adj"] = de_raw["pvals_adj"].reindex(ad_temp.var_names)
        _validate_de_result(de_df, backend=f"scanpy:{de_method}")
        de_df["logFC"] = de_df["logFC"].fillna(0.0)
        de_df["p_val"] = de_df["p_val"].fillna(1.0)
        de_df["p_adj"] = de_df["p_adj"].fillna(1.0)
        return de_df


def _run_mixedlm_de(
    adata: ad.AnnData,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str,
    n_jobs: int = 1,
    labels: Any | None = None,
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

    Performance note: fits a MixedLM (full + null) per gene via LRT (or Wald). This is
    O(n_genes) non-linear optimizations; use n_jobs for parallelism, or prefer faster
    backends (scanpy rank_genes or pydeseq2) for very large gene sets unless you need
    the per-gene delta_variance and sample-aware modeling.
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

    # Prepare expression for LMM: log1p-normalized, without double-transforming already-log data
    expr_mat = _prepare_log_normalized_expression(ad_temp)

    obs = ad_temp.obs
    condition = pd.Categorical(
        obs[use_groupby].astype(str),
        categories=[reference_group, target_group],
    )
    samples = obs[sample_col].astype(str).values

    n_genes = expr_mat.shape[1]
    var_names = ad_temp.var_names

    # Per-gene worker (returns idx, logfc, wald_p, lrt_p, delta_var, failed_fit_flag)
    # The explicit flag avoids counting truly neutral (logFC~0, p~1, dvar~0) biological genes as "failed".
    def _fit_gene_mixed(idx: int):
        y = expr_mat[:, idx].astype(float)
        # guard against near-constant expression (mixedlm will be singular)
        if float(np.nanvar(y)) < 1e-12:
            return idx, 0.0, 1.0, 1.0, 0.0, True
        df = pd.DataFrame({"y": y, "condition": condition, "sample": samples})
        try:
            with _de_warning_context():
                # Full model
                md_full = smf.mixedlm("y ~ C(condition)", df, groups=df["sample"])
                m_full = md_full.fit(reml=False, maxiter=200, disp=False)

                # Reduced (null) for LRT on condition contribution
                md_null = smf.mixedlm("y ~ 1", df, groups=df["sample"])
                m_null = md_null.fit(reml=False, maxiter=200, disp=False)

            if not getattr(m_full, "converged", True) or not getattr(m_null, "converged", True):
                return idx, 0.0, 1.0, 1.0, 0.0, True

            # LRT statistic and p (chi2 df=1 for the added fixed effect term(s))
            lrt_stat = -2.0 * (m_null.llf - m_full.llf)
            lrt_p = float(chi2.sf(max(lrt_stat, 0.0), 1))

            # Extract condition coef (target vs ref) — exact statsmodels name only
            expected_coef = f"C(condition)[T.{target_group}]"
            coef_name = expected_coef if expected_coef in m_full.params.index else None
            if coef_name is None:
                return idx, 0.0, 1.0, 1.0, 0.0, True
            logfc = float(m_full.params.get(coef_name, 0.0))
            p_wald = float(m_full.pvalues.get(coef_name, 1.0))

            # Delta variance: var attributable to fixed effects / total modeled var
            exog = m_full.model.exog
            beta = np.asarray(m_full.fe_params)
            fe_contrib = exog @ beta
            var_fe = float(np.var(fe_contrib))
            re_var = 0.0
            try:
                if (
                    hasattr(m_full, "cov_re")
                    and m_full.cov_re is not None
                    and len(m_full.cov_re) > 0
                ):
                    re_var = float(np.diag(m_full.cov_re)[0])  # first (only) RE variance
            except Exception:
                re_var = 0.0
            resid_var = float(getattr(m_full, "scale", 0.0))
            total_v = var_fe + max(re_var, 0.0) + max(resid_var, 0.0)
            delta_var = var_fe / total_v if total_v > 1e-12 else 0.0

            return idx, logfc, p_wald, lrt_p, float(np.clip(delta_var, 0.0, 1.0)), False
        except np.linalg.LinAlgError as e:
            logger.debug("MixedLM singular matrix for gene %d: %s", idx, e)
            return idx, 0.0, 1.0, 1.0, 0.0, True
        except Exception:
            # Degenerate fit (few samples per group, collinear, etc.) -> non-informative
            return idx, 0.0, 1.0, 1.0, 0.0, True

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
    failed_flags = np.array([r[5] for r in results], dtype=bool)

    # Choose which p-value to expose as the main "p_val" for active_score weighting and default filtering.
    # "wald": the coefficient test (standard for logFC-like effect)
    # "lrt": the likelihood ratio test for the condition term contribution (ties directly to delta_variance)
    if mixed_model_pval == "lrt":
        main_pvals = p_lrts
    else:
        if mixed_model_pval != "wald":
            logger.warning("mixed_model_pval must be 'wald' or 'lrt'; falling back to 'wald'.")
        main_pvals = p_walds

    with _de_warning_context():
        p_adjs = multipletests(main_pvals, method="fdr_bh")[1]

    de_df = pd.DataFrame(index=var_names)
    # mixedlm coefficient is on log1p scale (natural log). Convert to log2 for comparability.
    de_df["logFC"] = pd.Series(logfcs, index=var_names) / np.log(2)
    de_df["p_val"] = pd.Series(main_pvals, index=var_names)
    de_df["p_adj"] = pd.Series(p_adjs, index=var_names)
    de_df["delta_variance"] = pd.Series(dvars, index=var_names)
    de_df["delta_var_pval"] = pd.Series(p_lrts, index=var_names)

    # Use explicit failure flag returned from workers (avoids counting true biological
    # neutral genes (logFC~0, p~1, dvar~0) as "failed fits").
    n_failed = int(np.sum(failed_flags))
    n_total = len(de_df)
    failed_rate = (n_failed / n_total) if n_total else 0.0
    de_df.attrs["n_genes_failed_fit"] = n_failed
    de_df.attrs["failed_fit_rate"] = failed_rate
    if n_failed > 0:
        logger.warning(
            "MixedLM: %d/%d genes (%.1f%%) had degenerate or non-convergent fits "
            "(near-constant expression, singular Hessian, missing condition coefficient, etc.) "
            "and received neutral values (logFC=0, p_val=1, delta_variance=0). "
            "See diagnostics['mixed_model']['n_genes_failed_fit'].",
            n_failed,
            n_total,
            100.0 * failed_rate,
        )
    return de_df


def _run_memento_de(
    adata: ad.AnnData,
    groupby: str,
    target_group: str,
    reference_group: str,
    labels: Any | None = None,
    capture_rate: float = 0.07,
    num_boot: int = 5000,
    n_cpus: int = -1,
    counts: str | np.ndarray | sparse.spmatrix | pd.DataFrame | ad.AnnData | None = None,
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
            "Memento backend requested but 'memento-de' is not installed.\n"
            "Install with:\n"
            '    pip install "scatrans[memento]"\n'
            "or\n"
            "    pip install memento-de"
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

    # --- Resolve raw counts for Memento ---
    # Priority:
    # 1. Explicit `counts` argument (most flexible)
    # 2. adata.layers["counts"]
    # 3. adata.raw (if it has the counts)
    # 4. Current .X if it already looks like raw counts

    def _to_csr(x):
        from scipy.sparse import csr_matrix, issparse

        if issparse(x):
            return x.tocsr()
        return csr_matrix(np.asarray(x))

    raw_counts = None

    if counts is not None:
        if isinstance(counts, str):
            if counts in ad_temp.layers:
                raw_counts = ad_temp.layers[counts]
            else:
                raise ValueError(f"counts='{counts}' layer not found in adata.layers")
        elif isinstance(counts, ad.AnnData):
            if not np.array_equal(counts.obs_names, ad_temp.obs_names):
                missing = ad_temp.obs_names.difference(counts.obs_names)
                if len(missing):
                    raise ValueError(
                        f"counts AnnData is missing {len(missing)} cell(s) required by the "
                        f"comparison subset (first missing: {missing[0]!r}). "
                        "Pass counts aligned to adata.obs_names or use a matrix layer."
                    )
                counts = counts[ad_temp.obs_names]
            raw_counts = counts.X
            if counts.var_names.tolist() != ad_temp.var_names.tolist():
                common = counts.var_names.intersection(ad_temp.var_names)
                if len(common) == 0:
                    raise ValueError(
                        "No overlapping genes between provided counts AnnData and current adata"
                    )
                raw_counts = counts[:, common].X
                ad_temp = ad_temp[:, common].copy()
        else:
            raw_counts = counts

        if raw_counts is not None:
            raw_counts = _to_csr(raw_counts)
            logger.info("Memento: using explicitly provided counts.")

    if raw_counts is None and "counts" in getattr(ad_temp, "layers", {}):
        raw_counts = ad_temp.layers["counts"]
        logger.info("Memento: using 'counts' layer.")
        raw_counts = _to_csr(raw_counts)

    if (
        raw_counts is None
        and hasattr(ad_temp, "raw")
        and ad_temp.raw is not None
        and getattr(ad_temp.raw, "shape", (0, 0))[1] == ad_temp.n_vars
        and hasattr(ad_temp.raw, "var_names")
        and np.array_equal(ad_temp.raw.var_names, ad_temp.var_names)
    ):
        raw_counts = ad_temp.raw.X
        logger.info("Memento: using counts from adata.raw (exact match).")
        raw_counts = _to_csr(raw_counts)

    if raw_counts is None and _is_integer_counts_like(ad_temp.X):
        raw_counts = ad_temp.X
        logger.info("Memento: using current .X (looks like raw counts).")
        raw_counts = _to_csr(raw_counts)

    if raw_counts is None:
        logger.warning(
            "Could not obtain raw counts for Memento. "
            "Memento works best with raw integer UMI counts. "
            "Please call scat.store_raw_counts(adata) early (before HVG + log), "
            "or provide via the `counts` parameter, or ensure adata.raw / layers['counts'] has raw counts."
        )
        raw_counts = _to_csr(ad_temp.X)

    ad_temp.X = raw_counts

    from scipy.sparse import issparse

    if not (issparse(ad_temp.X) and isinstance(ad_temp.X, csr_matrix)):
        ad_temp.X = csr_matrix(ad_temp.X)

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

    # Schema guard: memento.binary_test_1d contract is not version-pinned strictly.
    # Without this, missing 'de_coef'/'de_pval' leads to silent/wrong or crashing fallback.
    if isinstance(result, pd.DataFrame):
        expected_cols = {"de_coef", "de_pval"}
        missing = expected_cols - set(result.columns)
        if missing:
            logger.warning(
                "Memento result missing required columns %s (possible memento-de API/version drift). "
                "logFC/p_val will fall back to neutral values (0/1). Results unreliable for use_memento_de=True. "
                "Pin compatible version e.g. 'memento-de>=0.1.0,<0.3.0' if needed.",
                sorted(missing),
            )

    de_df = pd.DataFrame(index=res_index)
    # Safe access to avoid scalar .reindex crash or silent zeroing on column absence
    if isinstance(result, pd.DataFrame) and "de_coef" in result.columns:
        m_lfc_raw = result["de_coef"]
    else:
        m_lfc_raw = pd.Series(0.0, index=res_index)
    m_lfc = pd.to_numeric(m_lfc_raw, errors="coerce").reindex(res_index).fillna(0.0)
    # memento de_coef is typically on natural log scale; convert to log2
    de_df["logFC"] = m_lfc / np.log(2)

    if isinstance(result, pd.DataFrame) and "de_pval" in result.columns:
        pval_raw = result["de_pval"]
    else:
        pval_raw = pd.Series(1.0, index=res_index)
    pvals_for_correction = pd.to_numeric(pval_raw, errors="coerce").reindex(res_index)
    valid_pval_mask = pvals_for_correction.notna()
    de_df["p_val"] = pvals_for_correction.fillna(1.0)

    # Preserve Memento-native adjusted p if present (audit trail); package p_adj uses BH for consistency.
    for native_col in ("de_padj", "de_padjs", "padj", "p_adj", "de_fdr"):
        if isinstance(result, pd.DataFrame) and native_col in result.columns:
            de_df["memento_p_adj_native"] = pd.to_numeric(
                result[native_col], errors="coerce"
            ).reindex(res_index)
            break

    de_df["p_adj"] = 1.0
    if valid_pval_mask.sum() > 0:
        with _de_warning_context():
            de_df.loc[valid_pval_mask, "p_adj"] = multipletests(
                pvals_for_correction[valid_pval_mask].values, method="fdr_bh"
            )[1]

    # Memento typically drops low-coverage genes entirely (not NaN rows). Genes absent
    # from the memento result are reindexed below with neutral fill — same UX issue as
    # PyDESeq2 filtered genes, but via row omission rather than NaN padj.
    n_genes_not_returned = int(len(adata.var_names.difference(res_index)))
    de_df.attrs["n_genes_not_returned_by_memento"] = n_genes_not_returned
    de_df.attrs["n_genes_missing_pval"] = int((~valid_pval_mask).sum())
    if n_genes_not_returned > 0:
        logger.warning(
            "Memento: %d/%d genes were not returned by binary_test_1d (e.g. internal "
            "min_perc_group filtering) and will appear as neutral values (logFC=0, p_adj=1) "
            "after reindexing. These were not tested. "
            "See de_df.attrs['n_genes_not_returned_by_memento'].",
            n_genes_not_returned,
            len(adata.var_names),
        )

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

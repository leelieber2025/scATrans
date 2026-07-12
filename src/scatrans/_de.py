"""
scATrans internal differential expression helpers.

Contains the wrapper that supports both scanpy rank_genes_groups and PyDESeq2
for pseudobulk. Extracted so tl.py stays small.
"""

from __future__ import annotations

import inspect
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
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from ._utils import (
    _as_contrast_categorical,
    _dense_expression_matrix,
    _is_integer_counts_like,
    _matrix_sum_axis0,
    _normalize_label_array,
    _prepare_log_normalized_expression,
    _require_matrix,
    _resolve_aligned_raw_counts,
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
def _de_warning_context(extra_categories: tuple = ()):
    with warnings.catch_warnings():
        for cat in (*_DE_SUPPRESS_WARNING_CATEGORIES, *extra_categories):
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


def _pydeseq2_filter_init_kwargs(cls: type, kwargs: dict[str, Any]) -> dict[str, Any]:
    """Drop kwargs not accepted by ``cls.__init__`` (API varies across pydeseq2 versions).

    Some releases accept ``n_cpus`` / ``quiet`` on ``DeseqDataSet`` but not on
    ``DeseqStats`` (or the reverse). Filtering by signature avoids
    ``TypeError: unexpected keyword argument`` on CI / older pins.
    """
    try:
        params = inspect.signature(cls.__init__).parameters
    except (TypeError, ValueError):
        return kwargs
    if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in params.values()):
        return kwargs
    return {k: v for k, v in kwargs.items() if k in params}


def _pydeseq2_preagg_count_like_verdict(
    ad_temp: ad.AnnData,
    *,
    count_source_name: str,
    X_for_deseq: Any,
    is_pseudobulk: bool,
) -> bool | None:
    """Return pre-aggregation count-likeness from ``uns`` when applicable, else None.

    Pseudobulk aggregation always rounds, so post-hoc ``_is_integer_counts_like``
    on the matrix that will feed PyDESeq2 is not a reliable safety net. Prefer:

    - ``pb_x_is_count_like`` when the matrix is ``.X``
    - ``pb_counts_is_count_like`` when the matrix is the aggregated ``layers['counts']``
      (including when passed as an already-resolved ``counts matrix`` that is the
      same object as ``layers['counts']``)
    """
    if not is_pseudobulk:
        return None

    if count_source_name == "adata.X" and "pb_x_is_count_like" in ad_temp.uns:
        return bool(ad_temp.uns["pb_x_is_count_like"])

    # Aggregated counts layer: string path, or resolved matrix that still *is* the layer.
    using_counts_layer = count_source_name == "layers['counts']"
    if not using_counts_layer and "counts" in getattr(ad_temp, "layers", {}):
        try:
            using_counts_layer = X_for_deseq is ad_temp.layers["counts"]
        except Exception:
            using_counts_layer = False

    if using_counts_layer:
        if "pb_counts_is_count_like" in ad_temp.uns:
            return bool(ad_temp.uns["pb_counts_is_count_like"])
        layer_map = ad_temp.uns.get("pb_layer_is_count_like")
        if isinstance(layer_map, dict) and "counts" in layer_map:
            return bool(layer_map["counts"])

    return None


def _coerce_pydeseq2_counts_matrix(
    counts: str | np.ndarray | sparse.spmatrix | pd.DataFrame | ad.AnnData,
    ad_temp: ad.AnnData,
) -> tuple[Any, str] | None:
    """Resolve ``counts=`` to a matrix aligned with ``ad_temp`` (n_obs × n_vars).

    Returns ``(matrix, source_name)`` or ``None`` if counts cannot be aligned safely.
    """
    if isinstance(counts, str):
        if counts not in ad_temp.layers:
            logger.warning(
                "PyDESeq2 counts='%s' not found in adata.layers; falling back to .X.",
                counts,
            )
            return None
        mat = ad_temp.layers[counts]
        source = f"layers['{counts}']"
    elif isinstance(counts, ad.AnnData):
        counts_ad: ad.AnnData = counts
        # Always align by names when they differ (same shape ≠ same order).
        if not np.array_equal(counts_ad.obs_names, ad_temp.obs_names):
            missing = ad_temp.obs_names.difference(counts_ad.obs_names)
            if len(missing):
                logger.warning(
                    "PyDESeq2 counts AnnData missing %d obs; falling back to .X.",
                    len(missing),
                )
                return None
            counts_ad = counts_ad[ad_temp.obs_names]
        if not np.array_equal(counts_ad.var_names, ad_temp.var_names):
            common = counts_ad.var_names.intersection(ad_temp.var_names)
            if len(common) != ad_temp.n_vars:
                logger.warning(
                    "PyDESeq2 counts AnnData gene set does not match adata; falling back to .X."
                )
                return None
            counts_ad = counts_ad[:, ad_temp.var_names]
        mat = counts_ad.X
        source = "counts AnnData"
    elif isinstance(counts, pd.DataFrame):
        if list(counts.index) != list(ad_temp.obs_names) or list(counts.columns) != list(
            ad_temp.var_names
        ):
            try:
                reindexed = counts.reindex(index=ad_temp.obs_names, columns=ad_temp.var_names)
            except Exception:
                logger.warning(
                    "PyDESeq2 counts DataFrame could not be reindexed to adata; falling back to .X."
                )
                return None
            if reindexed.shape != (ad_temp.n_obs, ad_temp.n_vars) or reindexed.isna().any().any():
                logger.warning(
                    "PyDESeq2 counts DataFrame reindex introduced missing values or wrong shape; "
                    "falling back to .X."
                )
                return None
            mat = reindexed.to_numpy()
            source = "counts DataFrame (reindexed)"
        else:
            mat = counts.to_numpy()
            source = "counts DataFrame"
    else:
        mat = counts
        source = "counts matrix"

    try:
        n0, n1 = mat.shape[0], mat.shape[1]
    except Exception:
        logger.warning("PyDESeq2 counts object has no shape; falling back to .X.")
        return None
    if n0 != ad_temp.n_obs or n1 != ad_temp.n_vars:
        logger.warning(
            "PyDESeq2 counts shape (%s, %s) != adata (%s, %s); falling back to .X.",
            n0,
            n1,
            ad_temp.n_obs,
            ad_temp.n_vars,
        )
        return None
    return mat, source


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
    paired_replicates: bool = False,
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
      - PyDESeq2 + scanpy backends (wilcoxon, t-test, t-test_overestim_var): native log2
      - mixedlm + memento: converted from natural/log1p scale (/ log(2))

    Internal function: type hints strengthened for mypy/pyright.
    """
    if use_mixed_model and use_memento_de:
        raise ValueError(
            "use_mixed_model=True and use_memento_de=True are incompatible. "
            "Choose one cell-level DE backend (MixedLM with sample_col, or Memento)."
        )

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
                paired_replicates=paired_replicates,
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
            ad_temp.obs[use_groupby] = _as_contrast_categorical(
                labels, reference_group, target_group
            )
        finally:
            _ann_log.setLevel(_prev_ann)

    if is_pseudobulk and pb_backend == "pydeseq2":
        # Validate design before importing the optional dependency so missing
        # pydeseq2 does not mask a clearer design-error ValueError, and so
        # base CI (.[dev] only) can exercise this branch without extras.
        n_t = int((ad_temp.obs[use_groupby] == target_group).sum())
        n_r = int((ad_temp.obs[use_groupby] == reference_group).sum())
        if n_t < 2 or n_r < 2:
            raise ValueError(
                f"PyDESeq2 requires >=2 replicates per group. Found {n_t} target, {n_r} ref."
            )

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

        # Prefer explicit counts= when aligned to current obs×var (layer / matrix /
        # AnnData). Callers often keep log1p in .X and raw counts elsewhere.
        count_source_name = "adata.X"
        X_for_deseq: Any = ad_temp.X
        if counts is not None:
            resolved = _coerce_pydeseq2_counts_matrix(counts, ad_temp)
            if resolved is not None:
                X_for_deseq, count_source_name = resolved
                logger.info("PyDESeq2: using %s for count matrix.", count_source_name)
            else:
                logger.warning(
                    "PyDESeq2: provided counts= could not be aligned to current "
                    "adata shape %s; falling back to .X.",
                    ad_temp.shape,
                )

        # Count-likeness of the matrix that will actually be used for DESeq2.
        # After _pseudobulk_with_layers, every aggregated matrix is np.round()'d and
        # therefore always looks integer — so prefer pre-aggregation verdicts in uns
        # for both .X (pb_x_is_count_like) and layers['counts'] (pb_counts_is_count_like).
        is_count_like = _pydeseq2_preagg_count_like_verdict(
            ad_temp,
            count_source_name=count_source_name,
            X_for_deseq=X_for_deseq,
            is_pseudobulk=is_pseudobulk,
        )
        if is_count_like is None:
            try:
                X_check_arr = _dense_expression_matrix(X_for_deseq)
                X_check_arr = np.clip(np.nan_to_num(X_check_arr), 0, None)
                is_count_like = _is_integer_counts_like(X_check_arr)
            except Exception as _e:
                logger.debug("Count-like check fallback: %s", _e)
                is_count_like = _is_integer_counts_like(X_for_deseq)

        if not is_count_like:
            msg = (
                f"Data passed to PyDESeq2 ({count_source_name}) does not look like raw "
                "non-negative integer counts. PyDESeq2 requires unnormalized integer counts. "
                "For pseudobulk data we automatically round to integer counts; "
                "set strict_pydeseq2_counts=False to allow the (rounded) data anyway. "
                "If raw counts live in layers['counts'], pass counts= or ensure "
                "pseudobulk aggregates that layer into .X."
            )
            if strict_pydeseq2_counts:
                raise ValueError(msg)
            logger.warning(msg)
        else:
            _warn_if_low_counts_matrix(X_for_deseq)

        # Narrow matrix type once, then densify / filter genes.
        X_mat = _require_matrix(X_for_deseq, name=count_source_name)
        if sparse.issparse(X_mat):
            gene_sums = _matrix_sum_axis0(X_mat)
            gene_keep = gene_sums >= min_counts_per_gene
            if gene_keep.sum() == 0:
                raise ValueError(
                    f"No genes passed the DESeq2 count filter (sum(counts) >= {min_counts_per_gene})."
                )
            X_filtered = _dense_expression_matrix(X_mat[:, gene_keep])
            X_filtered = np.clip(np.round(np.nan_to_num(X_filtered)), 0, None).astype(int)
            counts_use = pd.DataFrame(
                X_filtered, index=ad_temp.obs_names, columns=ad_temp.var_names[gene_keep]
            )
        else:
            X = _dense_expression_matrix(X_mat)
            X = np.clip(np.round(np.nan_to_num(X)), 0, None).astype(int)
            counts_df = pd.DataFrame(X, index=ad_temp.obs_names, columns=ad_temp.var_names)
            gene_keep = counts_df.sum(axis=0) >= min_counts_per_gene
            counts_use = counts_df.loc[:, gene_keep].copy()

        if counts_use.shape[1] == 0:
            raise ValueError(
                f"No genes passed the DESeq2 count filter (sum(counts) >= {min_counts_per_gene})."
            )

        condition = _normalize_label_array(ad_temp.obs[use_groupby])
        metadata = pd.DataFrame(
            {use_groupby: pd.Categorical(condition, categories=[reference_group, target_group])},
            index=counts_use.index,
        )

        with _de_warning_context():
            if _pydeseq2_uses_design_factors():
                dds_kw = _pydeseq2_filter_init_kwargs(
                    DeseqDataSet,
                    {
                        "counts": counts_use,
                        "metadata": metadata,
                        "design_factors": use_groupby,
                        "ref_level": [use_groupby, reference_group],
                        "quiet": True,
                        "n_cpus": n_jobs,
                    },
                )
                dds = DeseqDataSet(**dds_kw)
            else:
                dds_kw = _pydeseq2_filter_init_kwargs(
                    DeseqDataSet,
                    {
                        "counts": counts_use,
                        "metadata": metadata,
                        "design": f"~{use_groupby}",
                        "refit_cooks": True,
                        "quiet": True,
                        "n_cpus": n_jobs,
                    },
                )
                dds = DeseqDataSet(**dds_kw)
            dds.deseq2()

            # DeseqStats: some pydeseq2 builds reject n_cpus and/or quiet.
            stats_kw = _pydeseq2_filter_init_kwargs(
                DeseqStats,
                {
                    "contrast": [use_groupby, target_group, reference_group],
                    "quiet": True,
                    "n_cpus": n_jobs,
                },
            )
            stat_res = DeseqStats(dds, **stats_kw)
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
        if de_method == "logreg":
            raise ValueError(
                "de_method='logreg' is not supported: scanpy's logreg ranks genes by logistic "
                "regression scores and does not produce logFC, p-values, or adjusted p-values. "
                "Use 'wilcoxon', 't-test', or 't-test_overestim_var' instead."
            )
        labels = _normalize_label_array(ad_temp.obs[use_groupby])
        n_target = int((labels == target_group).sum())
        n_reference = int((labels == reference_group).sum())
        if n_target < 2 or n_reference < 2:
            raise ValueError(
                f"scanpy DE (method='{de_method}') requires at least 2 cells per group. "
                f"Found {n_target} in target '{target_group}' and {n_reference} in reference "
                f"'{reference_group}'. With a single cell per group, no valid statistics can be "
                f"computed. Consider pseudobulk aggregation (use_pseudobulk=True with sample_col) "
                f"if you have multiple biological replicates per condition."
            )
        rank_key = "_scatrans_rank_genes_groups"
        with _de_warning_context():
            try:
                sc.tl.rank_genes_groups(
                    ad_temp,
                    groupby=use_groupby,
                    groups=[target_group],
                    reference=reference_group,
                    method=de_method,
                    key_added=rank_key,
                )
            except Exception as exc:
                msg = str(exc)
                if "only contain one sample" in msg or "one sample" in msg.lower():
                    raise ValueError(
                        f"scanpy could not compute DE statistics: each group needs at least 2 "
                        f"cells (found {n_target} target, {n_reference} reference for "
                        f"'{use_groupby}'). For single-cell-per-group designs, use pseudobulk "
                        f"with biological replicates (use_pseudobulk=True, sample_col=...) or "
                        f"add more cells per group."
                    ) from exc
                raise
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
        # Neutral-fill missing/non-finite values (scanpy can emit ±inf logFC when
        # rank_genes_groups is run on raw counts and expm1 overflows).
        logfc = pd.to_numeric(de_df["logFC"], errors="coerce")
        de_df["logFC"] = logfc.where(np.isfinite(logfc), 0.0).fillna(0.0)
        pval = pd.to_numeric(de_df["p_val"], errors="coerce")
        de_df["p_val"] = pval.where(np.isfinite(pval), 1.0).fillna(1.0)
        padj = pd.to_numeric(de_df["p_adj"], errors="coerce")
        de_df["p_adj"] = padj.where(np.isfinite(padj), 1.0).fillna(1.0)
        return de_df


def _resolve_mixedlm_random_groups(
    obs: pd.DataFrame,
    groupby: str,
    sample_col: str,
    *,
    paired_replicates: bool = False,
    quiet: bool = False,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Build per-cell random-effect group IDs for MixedLM.

    When the same ``sample_col`` string appears under multiple ``groupby`` levels
    (e.g. both conditions use ``rep1``/``rep2`` for *different* individuals),
    pooling them in one random intercept is wrong. Default: composite
    ``{condition}::{sample}`` groups. Set ``paired_replicates=True`` when the same
    ID intentionally denotes the *same* biological replicate in every condition
    (paired / repeated-measures design).

    Parameters
    ----------
    quiet
        If True, suppress informational/warning logs (used inside tight loops
        such as sample-level permutation shuffles).
    """
    cond = obs[groupby].astype(str)
    raw = obs[sample_col].astype(str)
    cross = (
        pd.DataFrame({"condition": cond, "sample": raw})
        .groupby("sample", observed=True)["condition"]
        .nunique()
    )
    overlapping = cross[cross > 1].index.astype(str).tolist()

    if overlapping and not paired_replicates:
        groups = (cond + "::" + raw).to_numpy()
        grouping = "condition_sample_composite"
        if not quiet:
            logger.warning(
                "sample_col=%r has %d label(s) reused across %r levels (e.g. %s). "
                "MixedLM random effects use composite 'condition::sample' IDs so "
                "distinct biological replicates are not pooled. For paired designs "
                "where the same ID is the same individual in every condition, pass "
                "paired_replicates=True.",
                sample_col,
                len(overlapping),
                groupby,
                overlapping[:5],
            )
    else:
        groups = raw.to_numpy()
        grouping = "sample_col_raw"
        if not quiet:
            if overlapping and paired_replicates:
                logger.info(
                    "sample_col=%r labels overlap across %r levels; paired_replicates=True "
                    "so identical sample IDs share one random effect.",
                    sample_col,
                    groupby,
                )
            elif paired_replicates and not overlapping:
                logger.warning(
                    "paired_replicates=True but no sample_col labels are shared across %r "
                    "levels — pairing has no effect (each sample ID appears in only one "
                    "condition). Check that paired replicates use the same sample ID string "
                    "in every condition.",
                    groupby,
                )

    meta = {
        "grouping": grouping,
        "paired_replicates": bool(paired_replicates),
        "overlapping_sample_labels": overlapping,
        "n_random_groups": int(pd.Series(groups).nunique()),
        "n_random_groups_raw_sample_col": int(raw.nunique()),
    }
    return groups, meta


def _mixedlm_sample_aware_logfc(
    expr_mat: np.ndarray,
    *,
    condition: Any,
    samples: Any,
    target_group: str,
    reference_group: str,
    eps: float = 1e-9,
) -> np.ndarray:
    """Scanpy-style log2FC from **per-sample means**, then mean of sample means.

    Gives each biological sample equal weight so a large-N sample cannot dominate
    the effect size (aligned with the LMM's anti-pseudoreplication intent). Still
    on the log1p → expm1 → log2 scale used by scanpy rank_genes_groups.
    """
    expr = np.asarray(expr_mat, dtype=float)
    n_genes = expr.shape[1]
    cond = np.asarray(condition)
    samp = np.asarray(samples)

    def _mean_of_sample_means(mask: np.ndarray) -> np.ndarray | None:
        if not np.any(mask):
            return None
        sub_x = expr[mask]
        sub_s = samp[mask]
        # pd.unique preserves order; works for str/categorical codes
        uniq = pd.unique(sub_s)
        means: list[np.ndarray] = []
        for s in uniq:
            m = sub_s == s
            if np.any(m):
                means.append(sub_x[m].mean(axis=0))
        if not means:
            return None
        return np.mean(np.stack(means, axis=0), axis=0)

    mean_t = _mean_of_sample_means(cond == target_group)
    mean_r = _mean_of_sample_means(cond == reference_group)
    if mean_t is None or mean_r is None:
        return np.zeros(n_genes, dtype=float)
    logfc = np.log2((np.expm1(mean_t) + eps) / (np.expm1(mean_r) + eps))
    return np.where(np.isfinite(logfc), logfc, 0.0)


def _run_mixedlm_de(
    adata: ad.AnnData,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str,
    n_jobs: int = 1,
    labels: Any | None = None,
    mixed_model_pval: str = "wald",
    paired_replicates: bool = False,
) -> pd.DataFrame:
    """
    Mixed linear model (LMM) DE + Delta Variance using statsmodels mixedlm.

    Models: y_log ~ C(condition) + (1 | sample)
    - logFC: **sample-aware** scanpy-style log2 fold-change — mean of per-sample
      means within each condition, then log2((expm1(mean_t)+eps)/(expm1(mean_r)+eps)).
      Equal weight per biological sample (not dominated by high-cell samples);
      magnitude stays comparable to wilcoxon / PyDESeq2 cutoffs.
    - mixedlm_coef: fixed-effect coefficient on log1p-expression scale (natural log)
    - p_val / p_adj: Wald (or LRT) p for the condition fixed effect (BH adj across genes)
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
        ad_temp.obs[use_groupby] = _as_contrast_categorical(labels, reference_group, target_group)

    if sample_col not in ad_temp.obs.columns:
        raise ValueError(f"sample_col='{sample_col}' not found in adata.obs")

    # Keep in sync with scatrans.tl.MIXED_MODEL_MIN_SAMPLES_PER_GROUP / _TOTAL_SAMPLES.
    _min_per_group = 4
    _min_total = 6
    obs = ad_temp.obs
    cond_str = _normalize_label_array(obs[use_groupby])
    group_ids, group_meta = _resolve_mixedlm_random_groups(
        obs,
        use_groupby,
        sample_col,
        paired_replicates=paired_replicates,
    )
    group_s = pd.Series(group_ids, index=obs.index)
    samples_per_group = [
        int(group_s.loc[cond_str == g].nunique()) for g in (reference_group, target_group)
    ]
    total_samples = int(group_s.nunique())
    min_per_group = min(samples_per_group) if samples_per_group else 0
    if min_per_group < _min_per_group or total_samples < _min_total:
        raise ValueError(
            f"Mixed linear model requires >= {_min_per_group} biological samples per group "
            f"(found min={min_per_group}) and >= {_min_total} total random-effect groups "
            f"(found {total_samples}). With few replicates, use use_pseudobulk=True with "
            "pseudobulk_de_backend='pydeseq2' instead of use_mixed_model=True."
        )

    # Prepare expression for LMM: log1p-normalized, without double-transforming already-log data
    expr_mat = _prepare_log_normalized_expression(ad_temp)

    condition = pd.Categorical(
        cond_str,
        categories=[reference_group, target_group],
    )
    samples = group_ids

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
            # ConvergenceWarning is suppressed here: on small-sample/genome-wide runs,
            # many genes hit near-singular covariance or non-convergence, which floods
            # stdout with tens of thousands of per-gene statsmodels warnings. This is
            # already tracked and surfaced to the caller in aggregate via the
            # n_genes_failed_fit / failed_fit_rate metadata, so the raw per-gene
            # warning is redundant noise rather than new information.
            with _de_warning_context(extra_categories=(ConvergenceWarning,)):
                # Full model
                md_full = smf.mixedlm("y ~ C(condition)", df, groups=df["sample"])
                m_full = md_full.fit(reml=False, maxiter=200, disp=False)

                # Reduced (null) for LRT on condition contribution
                md_null = smf.mixedlm("y ~ 1", df, groups=df["sample"])
                m_null = md_null.fit(reml=False, maxiter=200, disp=False)

            # Extract condition coef (target vs ref). Categories are
            # [reference_group, target_group] so the name is always T.<target>.
            expected_coef = f"C(condition)[T.{target_group}]"
            if expected_coef not in m_full.params.index:
                return idx, 0.0, 1.0, 1.0, 0.0, True
            logfc = float(m_full.params.get(expected_coef, np.nan))
            p_wald = float(m_full.pvalues.get(expected_coef, np.nan))
            # statsmodels often reports converged=False when RE variance is on the
            # boundary (singular cov_re) even though fixed-effect coefs are usable.
            # Only discard when the coefficient itself is missing/non-finite.
            if not (np.isfinite(logfc) and np.isfinite(p_wald)):
                return idx, 0.0, 1.0, 1.0, 0.0, True

            # LRT: require finite llf on both models; else fall back to Wald p
            llf_full = float(getattr(m_full, "llf", np.nan))
            llf_null = float(getattr(m_null, "llf", np.nan))
            if np.isfinite(llf_full) and np.isfinite(llf_null):
                lrt_stat = -2.0 * (llf_null - llf_full)
                lrt_p = float(chi2.sf(max(lrt_stat, 0.0), 1))
                if not np.isfinite(lrt_p):
                    lrt_p = p_wald
            else:
                lrt_p = p_wald

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

    # Neutral-fill non-finite values BEFORE multipletests. statsmodels multipletests
    # is not NaN-safe: a single NaN p-value makes the entire corrected vector NaN,
    # which would then be filled to 1.0 for every gene (silent total loss of DE signal).
    logfcs = np.where(np.isfinite(logfcs), logfcs, 0.0)
    main_pvals = np.where(np.isfinite(main_pvals), main_pvals, 1.0)
    dvars = np.where(np.isfinite(dvars), dvars, 0.0)
    p_lrts = np.where(np.isfinite(p_lrts), p_lrts, 1.0)

    with _de_warning_context():
        p_adjs = multipletests(main_pvals, method="fdr_bh")[1]
    p_adjs = np.where(np.isfinite(p_adjs), p_adjs, 1.0)

    de_df = pd.DataFrame(index=var_names)
    # LMM coefficient on log1p-expression scale (natural-log mean difference).
    # Exposed as mixedlm_coef; logFC is a sample-aware scanpy-style log2FC so
    # shared cutoffs (e.g. 0.35) remain usable while matching anti-pseudoreplication.
    de_df["mixedlm_coef"] = pd.Series(logfcs, index=var_names)
    logfc_means = _mixedlm_sample_aware_logfc(
        expr_mat,
        condition=condition,
        samples=samples,
        target_group=target_group,
        reference_group=reference_group,
    )
    de_df["logFC"] = pd.Series(logfc_means, index=var_names)
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
    de_df.attrs["mixedlm_grouping"] = group_meta
    de_df.attrs["logFC_method"] = "sample_mean_of_means_log2"
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
    # Sign discordance is rare after sample-aware logFC but still possible because
    # p-values test the LMM fixed effect (mixedlm_coef), not the mean log2FC.
    finite_both = np.isfinite(logfc_means) & np.isfinite(logfcs)
    nonzero = finite_both & (np.abs(logfc_means) > 1e-12) & (np.abs(logfcs) > 1e-12)
    n_sign_flip = int(np.sum(nonzero & (np.sign(logfc_means) != np.sign(logfcs))))
    de_df.attrs["n_genes_logFC_mixedlm_sign_discordant"] = n_sign_flip
    if n_sign_flip > 0:
        logger.info(
            "MixedLM: %d/%d genes have opposite signs for sample-aware logFC vs "
            "mixedlm_coef (p_val tests the LMM coefficient). Prefer mixedlm_coef "
            "when the model direction matters.",
            n_sign_flip,
            n_total,
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
    # Fail fast on non-integer count layers before the optional import so
    # design/data errors are not masked by missing memento-de on base installs.
    if isinstance(counts, str):
        if counts not in adata.layers:
            raise ValueError(f"counts='{counts}' layer not found in adata.layers")
        if not _is_integer_counts_like(adata.layers[counts]):
            raise ValueError(
                f"Memento counts layer {counts!r} does not look like raw integer counts. "
                "Call store_raw_counts() before normalize/log1p, or pass a true count matrix."
            )

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
        ad_temp.obs[use_groupby] = _as_contrast_categorical(labels, reference_group, target_group)

    # Restrict to the two groups being compared (normalized labels)
    keep = pd.Series(
        _normalize_label_array(ad_temp.obs[use_groupby]), index=ad_temp.obs_names
    ).isin([target_group, reference_group])
    ad_temp = ad_temp[keep].copy()

    # Binary treatment column expected by memento.binary_test_1d
    ad_temp.obs["stim"] = (_normalize_label_array(ad_temp.obs[use_groupby]) == target_group).astype(
        int
    )

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
            if counts not in ad_temp.layers:
                raise ValueError(f"counts='{counts}' layer not found in adata.layers")
            raw_counts = ad_temp.layers[counts]
            # Integer check already applied on full adata above; re-check after subset.
            if not _is_integer_counts_like(raw_counts):
                raise ValueError(
                    f"Memento counts layer {counts!r} does not look like raw integer counts. "
                    "Call store_raw_counts() before normalize/log1p, or pass a true count matrix."
                )
            logger.info("Memento: using explicitly provided counts layer %r.", counts)
        else:
            # AnnData / ndarray / sparse / DataFrame — shape + name alignment
            coerced = _coerce_pydeseq2_counts_matrix(counts, ad_temp)
            if coerced is None:
                raise ValueError(
                    "Memento counts= is missing cells/genes or could not be aligned to the "
                    f"comparison subset shape {ad_temp.shape}. Pass a matrix/AnnData/layer "
                    "with matching obs×var order (or use counts='counts' after store_raw_counts)."
                )
            raw_counts, src = coerced
            if not _is_integer_counts_like(raw_counts):
                raise ValueError(f"Memento counts= ({src}) does not look like raw integer counts.")
            logger.info("Memento: using %s.", src)
        raw_counts = _to_csr(raw_counts)

    # Prefer shared resolve path (same integer + alignment guards as PyDESeq2).
    if raw_counts is None:
        resolved = _resolve_aligned_raw_counts(ad_temp, layer="counts", require_integer=True)
        if resolved is not None:
            raw_counts = _to_csr(resolved)
            logger.info("Memento: using aligned integer counts from layers['counts']/raw.")

    if raw_counts is None and _is_integer_counts_like(ad_temp.X):
        raw_counts = ad_temp.X
        logger.info("Memento: using current .X (looks like raw counts).")
        raw_counts = _to_csr(raw_counts)

    if raw_counts is None:
        raise ValueError(
            "Could not obtain raw integer counts for Memento. "
            "Call scat.store_raw_counts(adata) early (before HVG + log), "
            "or provide counts= / layers['counts'] with unnormalized integer UMIs. "
            "Refusing to run Memento on log-normalized or non-count matrices."
        )

    # Final shape guard (layers / .X paths)
    try:
        rc_shape = raw_counts.shape
    except Exception:
        rc_shape = None
    if rc_shape is not None and (rc_shape[0] != ad_temp.n_obs or rc_shape[1] != ad_temp.n_vars):
        raise ValueError(
            f"Memento counts shape {rc_shape} does not match AnnData subset "
            f"{ad_temp.shape}. Align genes/cells before calling."
        )

    ad_temp.X = raw_counts

    from scipy.sparse import issparse

    if not (issparse(ad_temp.X) and isinstance(ad_temp.X, csr_matrix)):
        ad_temp.X = csr_matrix(ad_temp.X)

    # Effective cpus
    effective_cpus = n_cpus if n_cpus and n_cpus > 0 else -1

    # memento-de's bootstrap internals emit a large volume of third-party
    # UserWarning/FutureWarning noise (unrelated to anything the caller can act
    # on); suppress it here so it doesn't flood stdout on genome-wide runs.
    with _de_warning_context(extra_categories=(UserWarning,)):
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
    # (important for the labels= permutation case and any internal subsetting).
    # DataFrame.attrs propagation across reindex/fillna is experimental in pandas
    # and not guaranteed across versions — capture and restore explicitly so
    # n_genes_not_returned_by_memento / n_genes_missing_pval cannot silently drop.
    saved_attrs = dict(de_df.attrs)
    de_df = de_df.reindex(adata.var_names).fillna({"logFC": 0.0, "p_val": 1.0, "p_adj": 1.0})
    de_df.attrs.update(saved_attrs)
    de_df.attrs["n_genes_not_returned_by_memento"] = n_genes_not_returned
    de_df.attrs["n_genes_missing_pval"] = int((~valid_pval_mask).sum())

    return de_df

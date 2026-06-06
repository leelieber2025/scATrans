"""
scATrans tl.py - Final Complete Dual-Track Version

This is the complete, ready-to-integrate tl.py with all review fixes applied.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse
from sklearn.linear_model import HuberRegressor
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
import warnings
import logging
from math import comb
import joblib
from joblib import Parallel, delayed
from typing import Optional, Union, List, Tuple, Dict, Any

try:
    from . import _version
    VERSION = _version.version
except (ImportError, AttributeError):
    VERSION = "0.7.0-dev"

logger = logging.getLogger(__name__)


# ==================== ORIGINAL HELPER FUNCTIONS (COMPLETE) ====================

def _is_integer_counts_like(X, max_check=100000):
    if sparse.issparse(X):
        data = X.data
        if data.size == 0:
            return True
        if not np.all(np.isfinite(data)):
            return False
        vals = data
    else:
        arr = np.asarray(X)
        if not np.all(np.isfinite(arr)):
            return False
        vals = arr.ravel()

    if vals.size == 0:
        return True

    if vals.size > max_check:
        rng = np.random.default_rng(0)
        vals = rng.choice(vals, size=max_check, replace=False)

    return np.all(vals >= 0) and np.allclose(vals, np.round(vals))


def _warn_if_not_integer_counts_matrix(X, max_check=100000):
    if not _is_integer_counts_like(X, max_check=max_check):
        logger.warning(
            "Data passed to PyDESeq2 may not be raw non-negative integer counts. "
            "Please ensure the input contains unnormalized counts."
        )


def _warn_if_low_counts_matrix(X, max_check=100000):
    if sparse.issparse(X):
        vals = X.data
    else:
        vals = np.asarray(X).ravel()

    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return

    if vals.size > max_check:
        rng = np.random.default_rng(0)
        vals = rng.choice(vals, size=max_check, replace=False)

    if vals.max() < 30:
        logger.warning(
            "Maximum count passed to PyDESeq2 is <30. This may be valid for small datasets, "
            "but please verify that the matrix contains raw counts, not normalized/log-transformed values."
        )


def _safe_add_matrices(a, b):
    if sparse.issparse(a) or sparse.issparse(b):
        return sparse.csr_matrix(a) + sparse.csr_matrix(b)
    return np.asarray(a) + np.asarray(b)


def _normalize_velocity_layers_by_size_factor(uns_layer, spl_layer, target_sum=None):
    total_layer = _safe_add_matrices(uns_layer, spl_layer)
    row_totals = np.asarray(total_layer.sum(axis=1)).ravel()
    positive = row_totals > 0

    if positive.sum() == 0:
        return uns_layer, spl_layer, row_totals, np.ones_like(row_totals, dtype=float)

    if target_sum is None:
        target_sum = np.median(row_totals[positive])

    factors = target_sum / np.maximum(row_totals, 1e-8)

    if sparse.issparse(uns_layer) or sparse.issparse(spl_layer):
        uns_layer = sparse.csr_matrix(uns_layer)
        spl_layer = sparse.csr_matrix(spl_layer)
        D = sparse.diags(factors)
        return D @ uns_layer, D @ spl_layer, row_totals, factors

    return (
        np.asarray(uns_layer) * factors[:, None],
        np.asarray(spl_layer) * factors[:, None],
        row_totals,
        factors,
    )


def _get_group_mean(matrix, mask):
    if np.sum(mask) == 0:
        raise ValueError("Cannot compute group mean for an empty group.")
    sub = matrix[mask]
    if sparse.issparse(sub):
        return np.asarray(sub.mean(axis=0)).ravel()
    return np.asarray(sub.mean(axis=0)).ravel()


def _pseudobulk_with_layers(
    adata, sample_col, groupby, layers=("spliced", "unspliced"),
    x_layer=None, use_total_for_x=False, min_cells=10, min_counts=1000,
):
    """Internal helper. After layer normalization in active_score(), layers are always 'spliced'/'unspliced'."""
    if sample_col not in adata.obs.columns:
        raise ValueError(f"sample_col='{sample_col}' not found.")
    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby='{groupby}' not found.")
    for layer in layers:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' not found in adata.layers")

    if use_total_for_x:
        X_source = _safe_add_matrices(adata.layers["spliced"], adata.layers["unspliced"])
        x_source_name = "spliced + unspliced"
    else:
        if x_layer is not None and x_layer not in adata.layers:
            raise ValueError(f"x_layer '{x_layer}' not found in adata.layers")
        X_source = adata.X if x_layer is None else adata.layers[x_layer]
        x_source_name = "adata.X" if x_layer is None else f"layer '{x_layer}'"

    group_df = adata.obs[[sample_col, groupby]].copy()
    group_df[sample_col] = group_df[sample_col].astype(str)
    group_df[groupby] = group_df[groupby].astype(str)
    pb_key = group_df[sample_col] + "||" + group_df[groupby]
    unique_keys = pd.Index(pb_key.unique())

    X_rows, obs_rows = [], []
    layer_rows = {layer: [] for layer in layers}

    for key in unique_keys:
        mask = pb_key.values == key
        n_cells = int(mask.sum())
        if n_cells < min_cells:
            continue
        x_sum = np.nan_to_num(np.asarray(X_source[mask].sum(axis=0)).ravel())
        if float(x_sum.sum()) < min_counts:
            continue

        sample_id, group_value = key.split("||", 1)
        X_rows.append(sparse.csr_matrix(x_sum.reshape(1, -1)))
        obs_rows.append({
            sample_col: sample_id, groupby: group_value,
            "n_cells": n_cells, "total_counts": float(x_sum.sum()),
            "pb_x_source": x_source_name,
        })
        for layer in layers:
            l_sum = np.nan_to_num(np.asarray(adata.layers[layer][mask].sum(axis=0)).ravel())
            layer_rows[layer].append(sparse.csr_matrix(l_sum.reshape(1, -1)))

    if not X_rows:
        raise ValueError("No samples remained after pseudobulk filtering.")

    adata_pb = ad.AnnData(X=sparse.vstack(X_rows).tocsr(),
                          obs=pd.DataFrame(obs_rows), var=adata.var.copy())
    adata_pb.obs.index = adata_pb.obs[sample_col].astype(str) + "_" + adata_pb.obs[groupby].astype(str)
    for layer in layers:
        adata_pb.layers[layer] = sparse.vstack(layer_rows[layer]).tocsr()
    adata_pb.obs_names_make_unique()
    return adata_pb


def _run_de_wrapper(adata, groupby, target_group, reference_group,
                    de_method="t-test_overestim_var", is_pseudobulk=False,
                    pb_backend="pydeseq2", n_jobs=1, labels=None,
                    strict_pydeseq2_counts=True):
    ad_temp = adata.copy() if labels is not None else adata
    use_groupby = groupby

    if labels is not None:
        use_groupby = "_de_temp_group"
        ad_temp.obs[use_groupby] = pd.Categorical(
            np.asarray(labels).astype(str), categories=[reference_group, target_group])

    if is_pseudobulk and pb_backend == "pydeseq2":
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
        except ImportError as e:
            raise ImportError(
                "pydeseq2 is required when pseudobulk_de_backend='pydeseq2'. "
                "Install with: pip install pydeseq2"
            ) from e

        n_t = (ad_temp.obs[use_groupby] == target_group).sum()
        n_r = (ad_temp.obs[use_groupby] == reference_group).sum()
        if n_t < 2 or n_r < 2:
            raise ValueError(f"PyDESeq2 requires >=2 replicates per group. Found {n_t} target, {n_r} ref.")

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
            counts_use = pd.DataFrame(X_filtered, index=ad_temp.obs_names, columns=ad_temp.var_names[gene_keep])
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
            index=counts_use.index
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                dds = DeseqDataSet(
                    counts=counts_use, metadata=metadata, design_factors=use_groupby,
                    ref_level=[use_groupby, reference_group], quiet=True, n_cpus=n_jobs
                )
            except TypeError:
                dds = DeseqDataSet(
                    counts=counts_use, metadata=metadata, design=f"~{use_groupby}",
                    refit_cooks=True, quiet=True, n_cpus=n_jobs
                )
            dds.deseq2()

            try:
                stat_res = DeseqStats(
                    dds, contrast=[use_groupby, target_group, reference_group],
                    quiet=True, n_cpus=n_jobs
                )
            except TypeError:
                stat_res = DeseqStats(
                    dds, contrast=[use_groupby, target_group, reference_group], n_cpus=n_jobs
                )
            stat_res.summary()

        res2 = stat_res.results_df.copy().reindex(ad_temp.var_names)
        de_df = pd.DataFrame(index=ad_temp.var_names)
        de_df["logFC"] = res2["log2FoldChange"].fillna(0.0)
        de_df["p_val"] = res2.get("pvalue", pd.Series(1.0, index=res2.index)).fillna(1.0)
        de_df["p_adj"] = res2.get("padj", pd.Series(1.0, index=res2.index)).fillna(1.0)
        return de_df

    else:
        rank_key = "_scatrans_rank_genes_groups"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.rank_genes_groups(
                ad_temp, groupby=use_groupby, groups=[target_group],
                reference=reference_group, method=de_method, key_added=rank_key
            )
        de_raw = sc.get.rank_genes_groups_df(ad_temp, group=target_group, key=rank_key).set_index("names")
        de_df = pd.DataFrame(index=ad_temp.var_names)
        de_df["logFC"] = de_raw["logfoldchanges"].reindex(ad_temp.var_names).fillna(0.0)
        de_df["p_val"] = de_raw["pvals"].reindex(ad_temp.var_names).fillna(1.0)
        de_df["p_adj"] = de_raw["pvals_adj"].reindex(ad_temp.var_names).fillna(1.0)
        return de_df


def _compute_velocity_delta(uns_layer, spl_layer, t_mask, r_mask, prior_weight=5.0):
    U_t = _get_group_mean(uns_layer, t_mask)
    S_t = _get_group_mean(spl_layer, t_mask)
    U_r = _get_group_mean(uns_layer, r_mask)
    S_r = _get_group_mean(spl_layer, r_mask)

    global_gamma = (np.sum(U_r) + 1e-8) / (np.sum(S_r) + 1e-8)
    beta = prior_weight
    alpha = global_gamma * beta
    gamma_ref = (U_r + alpha) / (S_r + beta)
    delta_velocity = U_t - (gamma_ref * S_t)

    total_uns = np.asarray(uns_layer.sum(axis=0)).ravel()
    total_spl = np.asarray(spl_layer.sum(axis=0)).ravel()
    total_us = total_uns + total_spl
    return np.nan_to_num(delta_velocity), np.nan_to_num(total_us)


def _get_exponential_scale_lambda(x):
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    x_pos = np.clip(x, 0.0, None)
    nonzero_x = x_pos[x_pos > 0]
    if len(nonzero_x) < 2:
        return 1e-8
    med = np.median(nonzero_x)
    return med / np.log(2.0) if med > 0 else 1e-8


def _soft_scale(x, lam):
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    x_pos = np.clip(x, 0.0, None)
    if lam <= 1e-8:
        return np.zeros_like(x)
    return 1.0 - np.exp(-x_pos / lam)


def _single_permutation_task(seed, original_labels, target_group, reference_group,
                             adata_subset, X_features, valid_feat, uns_layer, spl_layer,
                             total_us_for_filter, min_total_counts,
                             weight_fc, weight_unspliced, weight_pval,
                             lambda_fc, lambda_res, lambda_pval, is_pseudobulk, pb_backend,
                             de_method, prior_weight, de_preprocess, strict_pydeseq2_counts):
    rng = np.random.default_rng(seed)
    for _ in range(50):
        shuffled_labels = rng.permutation(original_labels)
        if not np.array_equal(shuffled_labels, original_labels):
            break
    else:
        logger.warning("Failed to generate a different permutation after 50 attempts.")

    ad_temp = adata_subset.copy()

    if de_preprocess == "normalize_log1p":
        sc.pp.normalize_total(ad_temp, target_sum=1e4)
        sc.pp.log1p(ad_temp)
    elif de_preprocess == "auto" and not (is_pseudobulk and pb_backend == "pydeseq2"):
        if "log1p" not in ad_temp.uns:
            sc.pp.normalize_total(ad_temp, target_sum=1e4)
            sc.pp.log1p(ad_temp)
    elif de_preprocess == "none":
        pass

    perm_de_df = _run_de_wrapper(
        ad_temp, groupby="_unused_when_labels_provided", target_group=target_group,
        reference_group=reference_group, de_method=de_method,
        is_pseudobulk=is_pseudobulk, pb_backend=pb_backend,
        n_jobs=1, labels=shuffled_labels,
        strict_pydeseq2_counts=strict_pydeseq2_counts
    )

    t_mask = shuffled_labels == target_group
    r_mask = shuffled_labels == reference_group
    delta_velocity, total_us_velocity = _compute_velocity_delta(
        uns_layer, spl_layer, t_mask, r_mask, prior_weight
    )

    total_us_for_filter = np.asarray(total_us_for_filter)
    valid_expr = total_us_for_filter >= min_total_counts
    residual = np.zeros(adata_subset.n_vars, dtype=float)

    if X_features is not None and (valid_feat & valid_expr).sum() >= 30:
        try:
            fit_mask = valid_feat & valid_expr
            X_fit = np.column_stack([
                np.log1p(adata_subset.var["gene_length"].values[fit_mask]),
                np.log1p(adata_subset.var["intron_number"].values[fit_mask])
            ])
            weights = np.clip(
                total_us_for_filter[fit_mask],
                a_min=None,
                a_max=np.percentile(total_us_for_filter[fit_mask], 95)
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = HuberRegressor(epsilon=1.35, max_iter=100).fit(
                    X_fit, delta_velocity[fit_mask], sample_weight=weights)
            pred = model.predict(X_features)
            residual[valid_feat] = delta_velocity[valid_feat] - pred
        except Exception as e:
            logger.debug(f"Bias correction failed in permutation: {e}")

    if np.all(residual == 0) and valid_expr.sum() > 0:
        residual[valid_expr] = delta_velocity[valid_expr] - np.nanmedian(delta_velocity[valid_expr])
    residual[~valid_expr] = 0.0

    s1 = _soft_scale(perm_de_df["logFC"].values, lambda_fc)
    s2 = _soft_scale(residual, lambda_res)
    s3 = _soft_scale(-np.log10(perm_de_df["p_adj"].values + 1e-300), lambda_pval)

    total_w = weight_fc + weight_unspliced + weight_pval
    return (weight_fc * s1 + weight_unspliced * s2 + weight_pval * s3) / total_w * 100.0


# ==================== IMPROVED ADVANCED HELPER ====================

def _compute_moments_velocity_delta(
    adata_comp: ad.AnnData,
    t_mask: np.ndarray,
    r_mask: np.ndarray,
    prior_weight: float = 5.0,
    n_neighbors: int = 30,
    n_pcs: int = 30,
    use_precomputed: bool = False,
    recompute_neighbors: bool = True,
    random_state: int = 42,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    try:
        import scvelo as scv
    except ImportError as e:
        raise ImportError(
            "mode='advanced' requires the 'scvelo' package. "
            "Please install it with: pip install scvelo"
        ) from e

    info: Dict[str, Any] = {}
    had_moments_before = "Mu" in adata_comp.layers and "Ms" in adata_comp.layers

    if use_precomputed and "Mu" in adata_comp.layers and "Ms" in adata_comp.layers:
        logger.info("Using precomputed Mu/Ms layers.")
        used_precomputed = True
        neighbors_source = "precomputed"
        n_neighbors_eff = None
        n_pcs_eff = None

        # Shape check for precomputed layers
        if adata_comp.layers["Mu"].shape != adata_comp.shape:
            raise ValueError("Layer 'Mu' shape does not match adata shape.")
        if adata_comp.layers["Ms"].shape != adata_comp.shape:
            raise ValueError("Layer 'Ms' shape does not match adata shape.")
    else:
        used_precomputed = False
        if "spliced" not in adata_comp.layers or "unspliced" not in adata_comp.layers:
            raise ValueError(
                "Advanced mode requires 'spliced' and 'unspliced' layers "
                "(after any automatic kb_python 'mature'/'nascent' remapping)."
            )

        n_obs = adata_comp.n_obs
        if n_obs < 5:
            raise ValueError(f"Advanced mode requires at least 5 observations (got {n_obs}).")

        if n_obs < 30:
            logger.warning(
                f"Advanced mode has only {n_obs} observations; scVelo moments may be unstable."
            )

        n_neighbors_eff = min(n_neighbors, max(2, n_obs - 1))
        n_pcs_eff = min(n_pcs, max(1, n_obs - 1), max(1, adata_comp.n_vars - 1))

        if n_pcs_eff < 2:
            raise ValueError("Too few observations/genes for reliable PCA in advanced mode.")

        if recompute_neighbors or "neighbors" not in adata_comp.uns:
            sc.pp.neighbors(
                adata_comp,
                n_neighbors=n_neighbors_eff,
                n_pcs=n_pcs_eff,
                random_state=random_state,
            )
            neighbors_source = "computed_by_scatrans"
        else:
            neighbors_source = "preexisting_neighbors"

        scv.pp.moments(adata_comp, n_pcs=n_pcs_eff, n_neighbors=n_neighbors_eff)

        # Shape check after moments
        if adata_comp.layers["Mu"].shape != adata_comp.shape:
            raise ValueError("Layer 'Mu' shape does not match adata shape after moments.")
        if adata_comp.layers["Ms"].shape != adata_comp.shape:
            raise ValueError("Layer 'Ms' shape does not match adata shape after moments.")

        info.update({
            "n_neighbors_effective": n_neighbors_eff,
            "n_pcs_effective": n_pcs_eff,
        })

    delta_velocity, total_us = _compute_velocity_delta(
        adata_comp.layers["Mu"],
        adata_comp.layers["Ms"],
        t_mask,
        r_mask,
        prior_weight,
    )

    info.update({
        "used_precomputed_moments": used_precomputed,
        "neighbors_source": neighbors_source,
        "had_moments_before": had_moments_before,
        "n_neighbors_requested": n_neighbors,
        "n_pcs_requested": n_pcs,
        "recompute_neighbors": recompute_neighbors,
        "random_state": random_state,
        "Mu_sparse": sparse.issparse(adata_comp.layers["Mu"]),
        "Ms_sparse": sparse.issparse(adata_comp.layers["Ms"]),
        "moments_shape": tuple(adata_comp.shape),
    })

    return delta_velocity, total_us, info


# ==================== COMPLETE active_score FUNCTION ====================

def active_score(
    adata_input,
    groupby: str = "condition",
    target_group: str = "GA",
    reference_group: str = "Ctrl",
    subset_col: Optional[str] = None,
    subset_values: Optional[Union[str, List[str], Tuple[str, ...]]] = None,
    weight_fc: float = 1.0,
    weight_unspliced: float = 1.0,
    weight_pval: float = 1.0,
    pval_cutoff: float = 0.05,
    logfc_cutoff: float = 0.5,
    active_fdr_cutoff: float = 0.05,
    de_method: str = "t-test_overestim_var",
    pseudobulk_de_backend: str = "pydeseq2",
    use_permutation: bool = False,
    perm_de_backend: str = "fast",
    n_perm: int = 100,
    n_jobs: int = -1,
    de_preprocess: str = "auto",
    gene_type_filter: Optional[str] = None,
    use_pseudobulk: bool = False,
    sample_col: Optional[str] = None,
    min_cells: int = 10,
    min_counts: int = 1000,
    pb_x_layer: str = "spliced",
    pb_use_total_for_x: bool = True,
    min_total_counts: int = 50,
    random_seed: int = 42,
    show_plot: bool = True,
    auto_adjust_n_perm: bool = True,
    prior_weight: float = 5.0,
    strict_pydeseq2_counts: bool = True,
    # ==================== DUAL-TRACK PARAMETERS ====================
    mode: str = "heuristic",
    advanced_fallback: bool = True,
    advanced_n_neighbors: int = 30,
    advanced_n_pcs: int = 30,
    advanced_use_precomputed: bool = False,
    allow_advanced_pseudobulk: bool = False,
    advanced_recompute_neighbors: bool = True,
    spliced_layer: str = "spliced",
    unspliced_layer: str = "unspliced",
):
    """
    Identify actively transcribed genes.

    Dual-Track Design:
    - mode="heuristic": Original fast group-wise method (default, recommended for most cases)
    - mode="advanced": Uses scVelo moments for neighborhood smoothing on contrast subset,
      then applies the same group-wise ratio delta + Huber bias correction.
      Experimental but more robust to noise.

    Layer name support:
    - Standard Scanpy/velocyto: 'spliced' and 'unspliced' (default)
    - kb_python (--lamanno / velocity mode): 'mature' (spliced) and 'nascent' (unspliced)
      → Automatically detected and remapped if standard names are missing.
    You can also manually specify custom layer names via `spliced_layer` / `unspliced_layer`.
    """
    # ==================== EARLY VALIDATION ====================
    if mode not in {"heuristic", "advanced"}:
        raise ValueError("mode must be either 'heuristic' or 'advanced'")

    if not isinstance(advanced_fallback, bool):
        raise ValueError("advanced_fallback must be boolean.")
    if not isinstance(advanced_n_neighbors, int) or advanced_n_neighbors < 2:
        raise ValueError("advanced_n_neighbors must be integer >= 2.")
    if not isinstance(advanced_n_pcs, int) or advanced_n_pcs < 2:
        raise ValueError("advanced_n_pcs must be integer >= 2.")
    if not isinstance(advanced_use_precomputed, bool):
        raise ValueError("advanced_use_precomputed must be boolean.")
    if not isinstance(allow_advanced_pseudobulk, bool):
        raise ValueError("allow_advanced_pseudobulk must be boolean.")
    if not isinstance(advanced_recompute_neighbors, bool):
        raise ValueError("advanced_recompute_neighbors must be boolean.")

    if mode == "advanced" and use_pseudobulk and not allow_advanced_pseudobulk:
        raise ValueError(
            "mode='advanced' is not supported with use_pseudobulk=True by default. "
            "Set allow_advanced_pseudobulk=True if you really want to try it."
        )

    if mode == "advanced" and use_pseudobulk and allow_advanced_pseudobulk:
        logger.warning(
            "Advanced mode on pseudobulk was explicitly enabled. "
            "This is experimental and may over-smooth sample-level replicates."
        )

    if mode == "advanced":
        logger.info(
            "Advanced mode is experimental and uses scVelo moments for smoothing, "
            "not scVelo's stochastic/dynamical velocity model."
        )

    logger.info(f"scATrans {VERSION} Analysis started. Mode: {mode}")

    # ==================== SUBSET & VALIDATION (original logic) ====================
    if subset_col is not None:
        if subset_col not in adata_input.obs.columns:
            raise ValueError(f"subset_col='{subset_col}' not found in adata.obs.columns")
        if subset_values is None:
            raise ValueError("subset_values must be provided when subset_col is specified")
        if isinstance(subset_values, (str, int, float)):
            subset_values_list = [str(subset_values)]
        else:
            subset_values_list = [str(v) for v in subset_values]
        subset_mask = adata_input.obs[subset_col].astype(str).isin(subset_values_list)
        n_before = adata_input.n_obs
        adata_input = adata_input[subset_mask].copy()
        n_after = adata_input.n_obs
        if n_after == 0:
            raise ValueError(f"No cells remain after subsetting {subset_col}")
        logger.info(f"Subsetted by {subset_col} ({n_after}/{n_before} cells remaining)")

    if not adata_input.var_names.is_unique:
        raise ValueError("adata.var_names must be unique.")

    target_group = str(target_group)
    reference_group = str(reference_group)

    if target_group == reference_group:
        raise ValueError("target_group and reference_group must be different.")

    if groupby not in adata_input.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found.")

    if target_group not in adata_input.obs[groupby].astype(str).unique():
        raise ValueError(f"target_group '{target_group}' not found.")
    if reference_group not in adata_input.obs[groupby].astype(str).unique():
        raise ValueError(f"reference_group '{reference_group}' not found.")

    # ====================== LAYER NAME HANDLING (kb_python support) ======================
    available_layers = list(adata_input.layers.keys())

    # Auto-detect kb_python layers if standard names are missing
    if spliced_layer not in available_layers or unspliced_layer not in available_layers:
        if "mature" in available_layers and "nascent" in available_layers:
            logger.warning(
                "Standard 'spliced'/'unspliced' layers not found in adata.layers. "
                "Auto-detected kb_python layers → using 'mature' as spliced (mature mRNA) "
                "and 'nascent' as unspliced (nascent pre-mRNA). "
                "All internal processing will use standard names after remapping. "
                "You can override with spliced_layer= / unspliced_layer= if needed."
            )
            spliced_layer = "mature"
            unspliced_layer = "nascent"
        else:
            raise ValueError(
                f"Required layers not found. "
                f"Expected '{spliced_layer}' + '{unspliced_layer}' (or kb_python 'mature' + 'nascent'). "
                f"Available layers: {available_layers}"
            )

    # Normalize to standard internal layer names so all existing code works unchanged
    if spliced_layer != "spliced" or unspliced_layer != "unspliced":
        if spliced_layer in adata_input.layers and unspliced_layer in adata_input.layers:
            adata_input.layers["spliced"] = adata_input.layers[spliced_layer].copy()
            adata_input.layers["unspliced"] = adata_input.layers[unspliced_layer].copy()
            logger.info(
                f"Layer remapping applied: '{spliced_layer}' → 'spliced', "
                f"'{unspliced_layer}' → 'unspliced' (internal use only)"
            )

    if "spliced" not in adata_input.layers or "unspliced" not in adata_input.layers:
        raise ValueError("Both 'spliced' and 'unspliced' layers are required after layer handling.")

    keep_mask = adata_input.obs[groupby].astype(str).isin([target_group, reference_group])
    adata = adata_input[keep_mask].copy()

    if gene_type_filter:
        if "gene_type" not in adata.var.columns:
            raise ValueError("'gene_type_filter' provided but 'gene_type' column is missing.")
        adata = adata[:, adata.var["gene_type"] == gene_type_filter].copy()

    if adata.n_vars == 0:
        raise ValueError("No genes remain after filtering.")

    if "gene_length" not in adata.var.columns:
        adata.var["gene_length"] = np.nan
    if "intron_number" not in adata.var.columns:
        adata.var["intron_number"] = np.nan

    gene_length = pd.to_numeric(adata.var["gene_length"], errors="coerce").to_numpy()
    intron_number = pd.to_numeric(adata.var["intron_number"], errors="coerce").to_numpy()

    adata.var["gene_length"] = gene_length
    adata.var["intron_number"] = intron_number

    valid_feat = (
        np.isfinite(gene_length) &
        np.isfinite(intron_number) &
        (gene_length >= 0) &
        (intron_number >= 0)
    )

    is_pseudobulk = False
    if use_pseudobulk:
        if sample_col is None:
            raise ValueError("sample_col must be provided when use_pseudobulk=True")
        logger.info("Performing pseudobulk aggregation...")
        adata = _pseudobulk_with_layers(
            adata, sample_col, groupby, x_layer=pb_x_layer,
            use_total_for_x=pb_use_total_for_x, min_cells=min_cells, min_counts=min_counts
        )
        is_pseudobulk = True
        adata.obs[groupby] = pd.Categorical(adata.obs[groupby].astype(str), categories=[reference_group, target_group])

    if de_preprocess == "normalize_log1p":
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif de_preprocess == "auto" and not (is_pseudobulk and pseudobulk_de_backend == "pydeseq2"):
        if "log1p" not in adata.uns:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
    elif de_preprocess == "none":
        pass

    X_features = np.column_stack([
        np.log1p(gene_length[valid_feat]),
        np.log1p(intron_number[valid_feat])
    ]) if valid_feat.sum() >= 50 else None

    effective_n_jobs = joblib.cpu_count() if n_jobs == -1 else max(1, n_jobs)

    logger.info("Performing differential expression analysis...")
    de_df = _run_de_wrapper(
        adata, groupby, target_group, reference_group,
        de_method=de_method, is_pseudobulk=is_pseudobulk,
        pb_backend=pseudobulk_de_backend, n_jobs=effective_n_jobs,
        strict_pydeseq2_counts=strict_pydeseq2_counts
    )

    adata.var["logFC"] = de_df["logFC"]
    adata.var["p_val"] = de_df["p_val"]
    adata.var["p_adj"] = de_df["p_adj"]

    uns_layer_raw = adata.layers["unspliced"]
    spl_layer_raw = adata.layers["spliced"]

    if is_pseudobulk:
        uns_layer, spl_layer, _, _ = _normalize_velocity_layers_by_size_factor(uns_layer_raw, spl_layer_raw)
    else:
        uns_layer, spl_layer = uns_layer_raw, spl_layer_raw

    obs_labels = adata.obs[groupby].astype(str).values
    t_mask = obs_labels == target_group
    r_mask = obs_labels == reference_group

    # ==================== DUAL-TRACK DELTA COMPUTATION ====================
    moments_info: Dict[str, Any] = {}
    velocity_layer_for_perm_uns = uns_layer
    velocity_layer_for_perm_spl = spl_layer

    if mode == "heuristic":
        delta_velocity, total_us_velocity = _compute_velocity_delta(
            uns_layer, spl_layer, t_mask, r_mask, prior_weight
        )
        velocity_source = "heuristic_global_ratio"

    elif mode == "advanced":
        adata_comp = adata.copy()
        if is_pseudobulk:
            adata_comp.layers["unspliced"] = uns_layer.copy()
            adata_comp.layers["spliced"] = spl_layer.copy()

        try:
            delta_velocity, total_us_velocity, moments_info = _compute_moments_velocity_delta(
                adata_comp,
                t_mask,
                r_mask,
                prior_weight=prior_weight,
                n_neighbors=advanced_n_neighbors,
                n_pcs=advanced_n_pcs,
                use_precomputed=advanced_use_precomputed,
                recompute_neighbors=advanced_recompute_neighbors,
                random_state=random_seed,
            )
            velocity_source = "scvelo_moments_groupwise_ratio"
            velocity_layer_for_perm_uns = adata_comp.layers["Mu"].copy()
            velocity_layer_for_perm_spl = adata_comp.layers["Ms"].copy()
            moments_info["advanced_failed"] = False
        except Exception as e:
            if advanced_fallback:
                logger.warning(f"Advanced mode failed: {e}. Falling back to heuristic.")
                delta_velocity, total_us_velocity = _compute_velocity_delta(
                    uns_layer, spl_layer, t_mask, r_mask, prior_weight
                )
                velocity_source = "heuristic_fallback_from_advanced"
                moments_info = {
                    "advanced_failed": True,
                    "failure_reason": str(e),
                }
            else:
                raise

    # ==================== BIAS CORRECTION ====================
    total_us_raw = (
        np.asarray(uns_layer_raw.sum(axis=0)).ravel() +
        np.asarray(spl_layer_raw.sum(axis=0)).ravel()
    )
    total_us_raw = np.nan_to_num(total_us_raw)

    valid_expr = total_us_raw >= min_total_counts
    residual = np.zeros(adata.n_vars, dtype=float)

    if X_features is not None and (valid_feat & valid_expr).sum() >= 30:
        try:
            fit_mask = valid_feat & valid_expr
            X_fit = np.column_stack([
                np.log1p(gene_length[fit_mask]),
                np.log1p(intron_number[fit_mask])
            ])
            weights = np.clip(
                total_us_raw[fit_mask],
                a_min=None,
                a_max=np.percentile(total_us_raw[fit_mask], 95)
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = HuberRegressor(epsilon=1.35, max_iter=500).fit(
                    X_fit, delta_velocity[fit_mask], sample_weight=weights)
            pred = model.predict(X_features)
            residual[valid_feat] = delta_velocity[valid_feat] - pred
        except Exception as e:
            logger.warning(f"Bias correction failed. Falling back to median. Reason: {e}")

    if np.all(residual == 0) and valid_expr.sum() > 0:
        residual[valid_expr] = delta_velocity[valid_expr] - np.nanmedian(delta_velocity[valid_expr])
    residual[~valid_expr] = 0.0

    adata.var["velocity_delta_raw"] = delta_velocity
    adata.var["velocity_residual"] = residual
    adata.var["total_us_counts"] = total_us_raw
    adata.var["total_us_counts_raw"] = total_us_raw
    adata.var["total_us_counts_velocity_layer"] = total_us_velocity
    adata.var["valid_expr"] = valid_expr
    adata.var["velocity_source"] = velocity_source

    # ==================== SCORING ====================
    lambda_fc = max(_get_exponential_scale_lambda(adata.var["logFC"].values), 0.25)
    lambda_res = max(_get_exponential_scale_lambda(residual), 1e-8)
    lambda_pval = max(
        _get_exponential_scale_lambda(-np.log10(adata.var["p_adj"].values + 1e-300)), 1.0
    )

    s1 = _soft_scale(adata.var["logFC"].values, lambda_fc)
    s2 = _soft_scale(residual, lambda_res)
    s3 = _soft_scale(-np.log10(adata.var["p_adj"].values + 1e-300), lambda_pval)

    total_w = weight_fc + weight_unspliced + weight_pval
    real_score = (weight_fc * s1 + weight_unspliced * s2 + weight_pval * s3) / total_w * 100.0
    adata.var["active_score"] = real_score

    # ==================== PERMUTATION (with consistent layers) ====================
    current_max_perm = None
    use_fdr_for_significance = True
    active_fdr_disabled_reason = None

    if use_permutation:
        if is_pseudobulk:
            n_t, n_r = t_mask.sum(), r_mask.sum()
            current_max_perm = float('inf') if n_t + n_r > 30 else max(1, comb(n_t + n_r, n_t) - 1)

        if perm_de_backend == "fast":
            perm_pb_backend, perm_de_method = "scanpy", "t-test_overestim_var"
        elif perm_de_backend == "same":
            perm_pb_backend, perm_de_method = pseudobulk_de_backend, de_method
        else:
            raise ValueError("perm_de_backend must be 'fast' or 'same'")

        if is_pseudobulk and auto_adjust_n_perm and np.isfinite(current_max_perm) and current_max_perm < n_perm:
            n_perm = int(current_max_perm)

        logger.info(f"Running parallel permutation testing ({n_perm} iterations)...")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            perm_results = Parallel(n_jobs=effective_n_jobs, backend="loky")(
                delayed(_single_permutation_task)(
                    i + random_seed,
                    obs_labels,
                    target_group,
                    reference_group,
                    adata,
                    X_features,
                    valid_feat,
                    velocity_layer_for_perm_uns,   # <-- consistent with advanced
                    velocity_layer_for_perm_spl,   # <-- consistent with advanced
                    total_us_raw,
                    min_total_counts,
                    weight_fc,
                    weight_unspliced,
                    weight_pval,
                    lambda_fc,
                    lambda_res,
                    lambda_pval,
                    is_pseudobulk,
                    perm_pb_backend,
                    perm_de_method,
                    prior_weight,
                    de_preprocess,
                    strict_pydeseq2_counts,
                ) for i in range(n_perm)
            )

        perm_scores_matrix = np.vstack(perm_results)
        exceed_count = np.sum(perm_scores_matrix >= real_score.reshape(1, -1), axis=0)
        pvals = (1.0 + exceed_count) / (n_perm + 1.0)
        adata.var["active_score_pval"] = pvals

        adata.var["active_score_fdr"] = np.ones(adata.n_vars)
        if valid_expr.sum() > 0:
            adata.var.loc[valid_expr, "active_score_fdr"] = multipletests(pvals[valid_expr], method="fdr_bh")[1]

        if current_max_perm is not None and current_max_perm < 100:
            use_fdr_for_significance = False
            active_fdr_disabled_reason = "small_permutation_space"

    # ==================== METADATA ====================
    velocity_delta_layer = (
        "scvelo_Mu_Ms_moments" if velocity_source.startswith("scvelo_moments")
        else ("size_factor_normalized_spliced_unspliced" if is_pseudobulk else "raw_spliced_unspliced")
    )

    adata.uns["scatrans"] = {
        "version": VERSION,
        "groupby": groupby,
        "target_group": target_group,
        "reference_group": reference_group,
        "mode": mode,
        "velocity_source": velocity_source,
        "velocity_delta_layer": velocity_delta_layer,
        "advanced_fallback": advanced_fallback,
        "advanced_use_precomputed": advanced_use_precomputed,
        "allow_advanced_pseudobulk": allow_advanced_pseudobulk,
        "advanced_recompute_neighbors": advanced_recompute_neighbors,
        "advanced_experimental": mode == "advanced",
        "moments_info": moments_info if mode == "advanced" else None,
        "advanced_neighbor_graph_basis": "adata.X_after_de_preprocess" if mode == "advanced" else None,
        "advanced_layer_preprocessing": (
            "existing_spliced_unspliced_layers_no_scv_filter_and_normalize"
            if mode == "advanced" else None
        ),
        "use_permutation": use_permutation,
        "n_perm": int(n_perm) if use_permutation else 0,
        "weight_fc": weight_fc,
        "weight_unspliced": weight_unspliced,
        "weight_pval": weight_pval,
        "pval_cutoff": pval_cutoff,
        "logfc_cutoff": logfc_cutoff,
        "active_fdr_cutoff": active_fdr_cutoff,
        "de_method": de_method,
        "pseudobulk_de_backend": pseudobulk_de_backend,
        "perm_de_backend": perm_de_backend if use_permutation else None,
        "prior_weight": prior_weight,
        "min_total_counts": min_total_counts,
        "random_seed": random_seed,
    }

    # ==================== SIGNIFICANT GENES ====================
    cols = [
        "active_score", "velocity_delta_raw", "velocity_residual",
        "logFC", "p_val", "p_adj",
        "total_us_counts", "total_us_counts_raw", "total_us_counts_velocity_layer",
        "valid_expr", "gene_length", "intron_number"
    ]
    if use_permutation:
        cols.extend(["active_score_pval", "active_score_fdr"])
    cols = [c for c in cols if c in adata.var.columns]

    mask = (
        (adata.var["p_adj"] < pval_cutoff)
        & (adata.var["logFC"] > logfc_cutoff)
        & (adata.var["velocity_residual"] > 0)
        & (adata.var["valid_expr"])
        & (adata.var["active_score"] > 0)
    )

    if use_permutation and use_fdr_for_significance:
        mask = mask & (adata.var["active_score_fdr"] < active_fdr_cutoff)

    significant = adata.var[mask][cols].copy().sort_values("active_score", ascending=False)
    all_results = adata.var[cols].copy().sort_values("active_score", ascending=False)

    logger.info(f"Analysis completed in {mode} mode! Significant active genes: {len(significant)}")

    # ==================== PLOTTING ====================
    if show_plot:
        try:
            from .pl import set_style
            set_style()
        except Exception:
            pass

        if velocity_source.startswith("scvelo_moments"):
            ylabel = "Velocity Residual (scVelo moments Mu/Ms)"
        elif is_pseudobulk:
            ylabel = "Velocity Residual (size-factor normalized U/S)"
        else:
            ylabel = "Velocity Residual (raw U/S)"

        display_mode = mode
        if velocity_source == "heuristic_fallback_from_advanced":
            display_mode = "advanced→heuristic fallback"

        fig, ax = plt.subplots(figsize=(7, 5.5), dpi=150)
        scatter = ax.scatter(
            adata.var["logFC"],
            adata.var["velocity_residual"],
            c=adata.var["active_score"],
            s=6 + adata.var["active_score"] * 0.4,
            cmap="viridis",
            alpha=0.75,
            edgecolors="none"
        )
        sig_bool = adata.var_names.isin(significant.index)
        if sig_bool.sum() > 0:
            ax.scatter(
                adata.var.loc[sig_bool, "logFC"],
                adata.var.loc[sig_bool, "velocity_residual"],
                facecolors="none",
                edgecolors="#d62728",
                s=55,
                linewidths=1.4,
                label="Significant active drivers",
                zorder=5
            )
            ax.legend(frameon=False, fontsize=9)

        ax.axvline(logfc_cutoff, color="#d62728", linestyle="--", linewidth=1.2, alpha=0.75)
        ax.axhline(0, color="#7f7f7f", linestyle=":", linewidth=1.0, alpha=0.8)

        ax.set_xlabel("logFC", fontweight="bold")
        ax.set_ylabel(ylabel, fontweight="bold")
        ax.set_title(f"scATrans Active Landscape ({display_mode})", fontweight="bold", pad=10)

        cbar = plt.colorbar(scatter, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label("Active Score (0-100)", fontsize=9, fontweight="bold", rotation=270, labelpad=12)

        plt.tight_layout()
        plt.show()
        plt.close(fig)

    return adata, significant, all_results

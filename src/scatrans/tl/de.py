"""scatrans.tl.de — internal package module."""

from __future__ import annotations

import logging
import warnings
from typing import Any

import anndata as ad
import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sparse

from .._de import _run_de_wrapper
from .._utils import (
    _apply_de_preprocess,
    _as_contrast_categorical,
    _clear_log_preprocess_metadata,
    _is_integer_counts_like,
    _pseudobulk_with_layers,
    _resolve_aligned_raw_counts,
    _subset_obs_mask,
    _validate_group_contrast,
)
from ._common import (
    VERSION,
    _coerce_memento_de_preprocess,
    _materialize_if_view,
    _require_explicit_groups,
    _select_obs,
    _select_var,
    _validate_de_common_options,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _record_raw_counts_metadata(
    adata: Any, *, save_raw: bool = False, overwrite: bool = False
) -> None:
    """Write scatrans metadata, optional adata.raw, and raw_* velocity layers."""
    if "scatrans" not in adata.uns:
        adata.uns["scatrans"] = {}
    prev_list = adata.uns["scatrans"].get("raw_gene_list")
    n_genes = int(adata.n_vars)
    if prev_list is not None and len(prev_list) != n_genes:
        logger.warning(
            "Updating raw_gene_list (%d → %d genes). If you subsetted after an earlier "
            "store_raw_counts(), enrichment universe will reflect the current gene set only. "
            "For the full pre-HVG universe, keep a separate full-gene AnnData or use save_raw=True "
            "on the original object before subsetting.",
            len(prev_list),
            n_genes,
        )
    # Preserve the original full universe under a sticky key when shrinking.
    if (
        prev_list is not None
        and len(prev_list) > n_genes
        and "raw_gene_list_full" not in adata.uns["scatrans"]
    ):
        adata.uns["scatrans"]["raw_gene_list_full"] = list(prev_list)
        logger.info(
            "Preserved previous full raw_gene_list (%d genes) as "
            "adata.uns['scatrans']['raw_gene_list_full'] for enrichment "
            "background if you re-store after HVG.",
            len(prev_list),
        )
    adata.uns["scatrans"]["raw_gene_list"] = list(adata.var_names)
    adata.uns["scatrans"]["store_raw_counts_n_genes"] = n_genes
    logger.info(
        "Saved the current gene list as the measured universe for enrichment "
        "(in adata.uns['scatrans']['raw_gene_list'])."
    )

    if save_raw:
        if getattr(adata, "raw", None) is not None and not overwrite:
            logger.debug("adata.raw already exists; skipping (pass overwrite=True to replace).")
        else:
            adata.raw = adata.copy()
            logger.info("Set adata.raw to preserve full data.")

    for vel_name in ("spliced", "unspliced", "mature", "nascent"):
        if vel_name in adata.layers:
            raw_vel_name = f"raw_{vel_name}"
            if raw_vel_name not in adata.layers or overwrite:
                adata.layers[raw_vel_name] = adata.layers[vel_name].copy()
                logger.info(f"Saved original {vel_name} to adata.layers['{raw_vel_name}'].")


def ensure_raw_counts(
    adata: Any, layer: str = "counts", save_raw: bool = False, overwrite: bool = False
) -> None:
    """
    Ensure raw integer counts are available in ``adata.layers[layer]``.

    Convenience wrapper around :func:`store_raw_counts` that also tries to recover
    counts from ``adata.raw`` when ``adata.X`` is already normalized or log-transformed
    (common after HVG + ``sc.pp.log1p``).

    Resolution order:
    1. Existing ``layers[layer]`` if it already looks like integer counts
    2. Current ``adata.X`` if integer counts-like
    3. ``adata.raw.X`` when gene names/order match ``adata.var_names``

    Always updates ``adata.uns['scatrans']['raw_gene_list']`` and velocity ``raw_*`` layers
    via the same metadata path as :func:`store_raw_counts`.
    """
    if (
        layer in adata.layers
        and not overwrite
        and _is_integer_counts_like(adata.layers[layer])
        and adata.layers[layer].shape[1] == adata.n_vars
    ):
        _record_raw_counts_metadata(adata, save_raw=save_raw, overwrite=overwrite)
        logger.debug("Layer '%s' already holds aligned integer counts.", layer)
        return

    if _is_integer_counts_like(adata.X):
        store_raw_counts(adata, layer=layer, save_raw=save_raw, overwrite=overwrite)
        return

    raw = getattr(adata, "raw", None)
    if raw is not None and _is_integer_counts_like(raw.X) and raw.shape[1] == adata.n_vars:
        if hasattr(raw, "var_names") and np.array_equal(raw.var_names, adata.var_names):
            adata.layers[layer] = raw.X.copy()
            logger.info(
                "ensure_raw_counts: recovered raw counts from adata.raw into layers['%s'].",
                layer,
            )
            _record_raw_counts_metadata(adata, save_raw=save_raw, overwrite=overwrite)
            return
        logger.warning(
            "adata.raw exists but gene names/order do not match current adata.var_names. "
            "Cannot recover counts automatically."
        )

    logger.warning(
        "ensure_raw_counts: adata.X does not look like raw counts and adata.raw could not be used. "
        "Falling back to store_raw_counts (may warn again)."
    )
    store_raw_counts(adata, layer=layer, save_raw=save_raw, overwrite=overwrite)


def store_raw_counts(
    adata: Any, layer: str = "counts", save_raw: bool = False, overwrite: bool = False
) -> None:
    """
    Store raw counts and the original spliced/unspliced (or mature/nascent) layers
    early in the analysis, right after loading and basic QC, but BEFORE HVG selection,
    normalization, or log1p.

    This is critical for scATrans because:
    - Memento and PyDESeq2 need raw counts for proper modeling.
    - Velocity / active transcription analysis (active_score) needs the original
      spliced/unspliced matrices on as many genes as possible.

    By default we only save to the given layer (save_raw defaults to False so we do
    not automatically touch adata.raw unless you explicitly ask for it).

    We automatically save any existing velocity layers under "raw_spliced",
    "raw_unspliced" (or "raw_mature", "raw_nascent"). These raw_* layers are
    subject to the normal AnnData behavior: if you later gene-subset the object
    (e.g. to HVGs), the layers are subsetted as well. They do **not** magically
    retain the original full-gene matrices after subsetting.

    If you need the full-gene raw velocity data after HVG-based visualization,
    either:
      - call store_raw_counts() on the full object and keep the full object for
        DE / active_score / enrichment while using a copy for visualization, or
      - use save_raw=True (which sets adata.raw).

    Recommended early call:
        scat.store_raw_counts(adata, layer="counts", save_raw=False)
    """
    if layer in adata.layers and not overwrite:
        if _is_integer_counts_like(adata.layers[layer]):
            logger.debug(f"Layer '{layer}' already exists with integer counts; not overwriting.")
        else:
            logger.warning(
                f"Existing layer '{layer}' does not look like raw counts; "
                "pass overwrite=True to replace it."
            )
    else:
        if not _is_integer_counts_like(adata.X):
            logger.warning(
                "Current adata.X does not look like raw integer counts. "
                "store_raw_counts should be called early (after basic QC, before normalize/log1p/HVG)."
            )
        mat = adata.X.copy()
        if mat.shape[1] != adata.n_vars:
            raise ValueError(
                f"Cannot store raw counts: matrix has {mat.shape[1]} columns "
                f"but adata has {adata.n_vars} genes."
            )
        adata.layers[layer] = mat
        logger.info(f"Saved raw counts to adata.layers['{layer}'].")

    if layer in adata.layers and adata.layers[layer].shape[1] != adata.n_vars:
        raise ValueError(
            f"Layer '{layer}' has {adata.layers[layer].shape[1]} columns but adata has "
            f"{adata.n_vars} genes. Pass overwrite=True after fixing alignment."
        )

    _record_raw_counts_metadata(adata, save_raw=save_raw, overwrite=overwrite)


def restore_raw_counts(adata: Any, layer: str = "counts", inplace: bool = False) -> Any | None:
    """
    Restore raw counts from the stored layer (preferred) or adata.raw back into .X.

    This is useful when you have done HVG + log1p on .X for visualization,
    but want to work with (or pass to other tools) the raw counts for the
    genes currently in the adata (or the preserved set).

    It only uses explicitly stored raw data (from store_raw_counts), never
    attempts to recover from log-transformed data.

    Parameters
    ----------
    adata : AnnData
        The AnnData object.
    layer : str
        The layer name where raw counts were stored (default "counts").
    inplace : bool
        If True, modify adata in place and return None.
        If False (default), return a new AnnData with .X set to raw counts.

    Returns
    -------
    AnnData or None
        If not inplace, a copy of adata with raw counts in .X.
    """
    if layer in adata.layers:
        raw = adata.layers[layer].copy()
        source = f"layers['{layer}']"
    elif getattr(adata, "raw", None) is not None:
        raw = adata.raw.X.copy()
        source = "adata.raw"
    else:
        raise ValueError(
            f"No raw counts found in layer '{layer}' or adata.raw. "
            "Call scat.store_raw_counts(adata) early to preserve them."
        )

    if raw.shape[1] != adata.n_vars:
        raise ValueError(
            f"Stored raw counts have {raw.shape[1]} genes, but current adata has {adata.n_vars} genes. "
            "Cannot restore into .X without explicit gene reindexing. "
            "Use the object before gene subsetting, or call store_raw_counts() again on the current object."
        )

    # Guard against same-dimension but permuted gene order (layers and adata.raw).
    if source == "adata.raw" and hasattr(adata.raw, "var_names"):
        if not np.array_equal(adata.raw.var_names, adata.var_names):
            raise ValueError(
                "adata.raw has the same number of genes as current adata, but gene names/order differ. "
                "Cannot restore into .X without explicit gene reindexing."
            )
    elif source.startswith("layers"):
        raw_gene_list = adata.uns.get("scatrans", {}).get("raw_gene_list")
        if (
            raw_gene_list is not None
            and len(raw_gene_list) == adata.n_vars
            and not np.array_equal(np.asarray(raw_gene_list), adata.var_names.to_numpy())
        ):
            raise ValueError(
                f"Stored counts in {source} match n_vars but stored raw_gene_list order differs "
                "from current adata.var_names. Cannot restore into .X without explicit gene "
                "reindexing. Re-run store_raw_counts() on the current object."
            )

    target = adata if inplace else adata.copy()
    target.X = raw
    _clear_log_preprocess_metadata(target)
    if inplace:
        logger.info(f"Restored raw counts from {source} into adata.X (inplace).")
        return None
    logger.info(f"Created copy with raw counts from {source} in .X.")
    return target


def differential_expression(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str | None = None,
    reference_group: str | None = None,
    subset_col: str | None = None,
    subset_values: str | list[str] | tuple[str, ...] | None = None,
    de_method: str = "t-test_overestim_var",
    pseudobulk_de_backend: str = "pydeseq2",
    pydeseq2_min_counts: int = 10,
    use_pseudobulk: bool = False,
    sample_col: str | None = None,
    min_cells: int = 10,
    min_counts: int = 1000,
    pb_x_layer: str = "X",  # for pseudobulk, what to aggregate (usually the count matrix)
    # Default False for DE-only: do not silently sum U+S (total RNA) when velocity
    # layers exist. Prefer .X / counts. active_score keeps pb_use_total_for_x=True.
    pb_use_total_for_x: bool = False,
    de_preprocess: str = "auto",
    min_total_counts: int = 50,  # reserved for API compatibility / future use; currently not enforced in DE path
    strict_pydeseq2_counts: bool = True,
    use_mixed_model: bool = False,
    use_delta_variance_pval: bool = False,
    delta_var_pval_cutoff: float = 0.05,
    mixed_model_pval: str = "wald",
    paired_replicates: bool = False,
    # Memento support (first-class, integrated backend)
    use_memento_de: bool = False,
    memento_capture_rate: float = 0.07,
    memento_num_boot: int = 5000,
    memento_n_cpus: int = -1,
    n_jobs: int = -1,
    gene_type_filter: str | None = None,
    # Allow providing raw counts separately when adata.X is already HVG+log (very common)
    counts: str | np.ndarray | sparse.spmatrix | pd.DataFrame | ad.AnnData | None = None,
    copy_input: bool = True,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """
    Standalone differential expression (DE) using the same flexible backends
    as scATrans (scanpy methods, PyDESeq2 pseudobulk, mixed linear models,
    and Memento -- the Cell 2024 method-of-moments framework).

    This function does **not** require spliced/unspliced (velocity) layers.
    It is intended for users who want high-quality DE (especially via Memento),
    followed by scATrans' downstream tools:

        candidates = scat.filter_active_genes(de_results, pval_cutoff=0.05, logfc_cutoff=0.3)  # upregulated
        # down or both directions:
        # down_cands = scat.filter_active_genes(de_results, pval_cutoff=0.05, logfc_cutoff=0.3, logfc_direction="down")
        # For enrichment, pass adata= (if store_raw_counts was used) so it uses
        # the preserved full measured gene set as universe, not just current HVGs.
        enrich = scat.run_enrichment(candidates.index.tolist(), ..., adata=adata)
        scat.pl.volcano_plot(de_results, ...)
        scat.pl.enrich_dotplot(enrich, ...)

    All DE-related options from `active_score` are supported here
    (pseudobulk, mixed models, Memento, etc.), except permutation-based FDR
    (use ``active_score(..., use_permutation=True)`` when velocity layers are available).
    For a minimal-parameter entry point see ``active_score_simple`` or ``run_default_pipeline``.

    copy_input : bool, default True
        Same semantics as :func:`active_score`: one combined obs-filter copy when
        True; zero ``AnnData.copy()`` calls when False and no obs filtering is needed.

    Returns
    -------
    (adata_with_results, results_df)
        - results_df is sorted by p_adj (ascending; most significant first) and
          contains at minimum: logFC, p_val, p_adj, and (when use_memento_de) the
          native memento_de_* / memento_dv_* columns.
        - adata.var is updated with the same columns for convenience.
        - Metadata is stored under adata.uns["scatrans"].
    """
    _require_explicit_groups(target_group, reference_group, func_name="differential_expression")

    # --- minimal shared validation (subset + group checks) ---
    obs_filter = pd.Series(True, index=adata_input.obs_names)
    if subset_col is not None:
        if subset_col not in adata_input.obs.columns:
            raise ValueError(f"subset_col='{subset_col}' not found in adata.obs.columns")
        if subset_values is None:
            raise ValueError("subset_values must be provided when subset_col is specified")
        subset_mask = _subset_obs_mask(adata_input.obs[subset_col], subset_values)
        if int(subset_mask.sum()) == 0:
            raise ValueError("No cells remain after subsetting.")
        obs_filter &= subset_mask

    if not adata_input.var_names.is_unique:
        raise ValueError("adata.var_names must be unique.")

    if groupby not in adata_input.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found.")

    target_group, reference_group, norm_groups = _validate_group_contrast(
        adata_input.obs[groupby],
        groupby=groupby,
        target_group=str(target_group),
        reference_group=str(reference_group),
    )

    obs_filter &= norm_groups.isin([target_group, reference_group])
    _caller_adata = adata_input
    adata_input = _select_obs(adata_input, obs_filter, copy_input=copy_input)
    if adata_input.n_obs == 0:
        raise ValueError(
            "No cells match target/reference groups after filtering. "
            f"Check target_group='{target_group}' and reference_group='{reference_group}' "
            f"against adata.obs['{groupby}'] (missing labels are excluded)."
        )
    adata_input = _materialize_if_view(adata_input)
    # Never mutate the caller's AnnData (labels / preprocess / .var).
    if adata_input is _caller_adata:
        adata_input = adata_input.copy()
        if not copy_input:
            logger.info(
                "copy_input=False: isolated a working copy before mutation so the "
                "caller's AnnData is left unchanged."
            )
    adata_input.obs[groupby] = norm_groups.loc[obs_filter].values

    if gene_type_filter:
        if "gene_type" not in adata_input.var.columns:
            raise ValueError("'gene_type_filter' provided but 'gene_type' column is missing.")
        adata_input = _select_var(
            adata_input, adata_input.var["gene_type"] == gene_type_filter, copy_input=copy_input
        )
        adata_input = _materialize_if_view(adata_input)

    if adata_input.n_vars == 0:
        raise ValueError("No genes remain after filtering.")

    if use_mixed_model and sample_col is None:
        raise ValueError("sample_col must be provided when use_mixed_model=True")

    if use_pseudobulk and sample_col is None:
        raise ValueError("sample_col must be provided when use_pseudobulk=True")

    # Memento-specific guard (same as in active_score)
    if use_memento_de and use_pseudobulk:
        raise ValueError(
            "use_memento_de=True is not supported with use_pseudobulk=True "
            "(Memento is a cell-level method-of-moments estimator)."
        )
    if use_memento_de and use_mixed_model:
        raise ValueError(
            "use_mixed_model=True and use_memento_de=True are incompatible. "
            "Choose one cell-level DE backend."
        )

    # Memento requires count data; force no log-norm preprocess for the DE leg
    de_preprocess = _coerce_memento_de_preprocess(use_memento_de, de_preprocess)

    # Shared DE option validation (deduplicated via helper)
    _validate_de_common_options(
        de_preprocess=de_preprocess,
        pseudobulk_de_backend=pseudobulk_de_backend,
        n_jobs=n_jobs,
        use_permutation=False,
        n_perm=0,
        use_mixed_model=use_mixed_model,
        mixed_model_pval=mixed_model_pval,
        paired_replicates=paired_replicates,
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        min_cells=min_cells,
        min_counts=min_counts,
    )

    if min_total_counts != 50:
        logger.warning(
            "differential_expression: min_total_counts=%s is not enforced in the DE-only path "
            "(reserved for API compatibility). It affects gene filtering in active_score() only. "
            "Use min_cells / min_counts for pseudobulk filtering instead.",
            min_total_counts,
        )

    if use_delta_variance_pval:
        logger.warning(
            "differential_expression: use_delta_variance_pval=True is not enforced in the DE-only "
            "path (this function returns the full ranked results table, not a significant-gene "
            "subset). Use active_score() for delta-variance filtering, or filter manually via "
            "results['delta_var_pval'] < delta_var_pval_cutoff (currently %.4g).",
            delta_var_pval_cutoff,
        )

    # --- prepare data (pseudobulk if requested) ---
    adata = adata_input

    if use_pseudobulk:
        # sample_col is required above when use_pseudobulk=True
        assert sample_col is not None
        logger.info("Performing pseudobulk aggregation for DE...")
        available_layers = [
            layer for layer in ("spliced", "unspliced", "counts") if layer in adata.layers
        ]
        x_layer_eff = pb_x_layer if pb_x_layer != "X" else None
        if (
            x_layer_eff is None
            and not (
                pb_use_total_for_x and "spliced" in adata.layers and "unspliced" in adata.layers
            )
            and pseudobulk_de_backend == "pydeseq2"
            and "counts" in adata.layers
        ):
            x_layer_eff = "counts"
            logger.info(
                "Pseudobulk: aggregating layers['counts'] into .X for PyDESeq2 "
                "(adata.X may be log-normalized)."
            )
        adata = _pseudobulk_with_layers(
            adata,
            sample_col,
            groupby,
            layers=available_layers,
            x_layer=x_layer_eff,
            use_total_for_x=pb_use_total_for_x
            and ("spliced" in adata.layers and "unspliced" in adata.layers),
            min_cells=min_cells,
            min_counts=min_counts,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            adata.obs[groupby] = _as_contrast_categorical(
                adata.obs[groupby], reference_group, target_group
            )

    # Resolve counts AFTER any pseudobulk so the matrix matches current obs.
    if counts is None and (
        use_memento_de or (use_pseudobulk and pseudobulk_de_backend == "pydeseq2")
    ):
        counts = _resolve_aligned_raw_counts(adata, layer="counts", require_integer=True)
        if counts is None and ("counts" in adata.layers or getattr(adata, "raw", None)):
            logger.warning(
                "Count-based DE was requested but no safely aligned raw counts were found. "
                "Call scat.store_raw_counts(adata) before HVG/normalize, or pass counts= explicitly."
            )
    elif counts is not None and use_pseudobulk:
        # Cell-level counts= matrix is invalid after aggregation; prefer pb layer / .X.
        if hasattr(counts, "shape") and getattr(counts, "shape", (0, 0))[0] != adata.n_obs:
            logger.info(
                "differential_expression: counts= matrix is cell-level but data is now "
                "pseudobulk; using aggregated layers['counts']/.X instead."
            )
            counts = _resolve_aligned_raw_counts(adata, layer="counts", require_integer=True)

    # DE preprocess
    # (Memento coercion to 'none' already performed early via _coerce_memento_de_preprocess)

    _apply_de_preprocess(
        adata,
        de_preprocess,
        skip_auto=use_pseudobulk and pseudobulk_de_backend == "pydeseq2",
    )

    effective_n_jobs = joblib.cpu_count() if n_jobs == -1 else max(1, n_jobs)

    # --- run DE via the shared engine (Memento, scanpy, DESeq2, mixedlm all supported) ---
    logger.info("Performing differential expression analysis (differential_expression mode)...")
    de_df = _run_de_wrapper(
        adata,
        groupby,
        target_group,
        reference_group,
        de_method=de_method,
        is_pseudobulk=use_pseudobulk,
        pb_backend=pseudobulk_de_backend,
        n_jobs=effective_n_jobs,
        strict_pydeseq2_counts=strict_pydeseq2_counts,
        use_mixed_model=use_mixed_model,
        sample_col=sample_col if use_mixed_model else None,
        mixed_model_pval=mixed_model_pval,
        paired_replicates=paired_replicates,
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        memento_n_cpus=memento_n_cpus,
        counts=counts,
        min_counts_per_gene=pydeseq2_min_counts,
    )

    # Store results
    adata.var["logFC"] = de_df["logFC"]
    adata.var["p_val"] = de_df["p_val"]
    adata.var["p_adj"] = de_df["p_adj"]

    for extra in [
        "mixedlm_coef",
        "delta_variance",
        "delta_var_pval",
        "memento_de_se",
        "memento_dv_coef",
        "memento_dv_se",
        "memento_dv_pval",
        "memento_p_adj_native",
    ]:
        if extra in de_df.columns:
            adata.var[extra] = de_df[extra]

    n_mixed_failed_de = (
        int(de_df.attrs.get("n_genes_failed_fit", 0))
        if (use_mixed_model and hasattr(de_df, "attrs"))
        else 0
    )
    mixed_failed_rate_de = (
        float(de_df.attrs.get("failed_fit_rate", 0.0))
        if (use_mixed_model and hasattr(de_df, "attrs"))
        else 0.0
    )
    mixedlm_logfc_method_de = (
        de_df.attrs.get("logFC_method") if (use_mixed_model and hasattr(de_df, "attrs")) else None
    )
    n_mixed_sign_discordant_de = (
        int(de_df.attrs.get("n_genes_logFC_mixedlm_sign_discordant", 0) or 0)
        if (use_mixed_model and hasattr(de_df, "attrs"))
        else 0
    )
    pydeseq2_diag = {}
    if use_pseudobulk and pseudobulk_de_backend == "pydeseq2" and hasattr(de_df, "attrs"):
        pydeseq2_diag = {
            "n_genes_filtered_low_count": int(de_df.attrs.get("n_genes_filtered_low_count", 0)),
            "n_genes_nan_from_deseq2": int(de_df.attrs.get("n_genes_nan_from_deseq2", 0)),
            "neutral_fill": bool(de_df.attrs.get("pydeseq2_neutral_fill", True)),
            "note": (
                "Genes filtered by min_counts or marked NaN by DESeq2 independent filtering "
                "appear as logFC=0, p_adj=1 and are not 'tested and non-significant'."
            ),
        }
    n_memento_not_returned = (
        int(de_df.attrs.get("n_genes_not_returned_by_memento", 0))
        if (use_memento_de and hasattr(de_df, "attrs"))
        else 0
    )

    # Build clean results table (no velocity columns)
    cols = ["logFC", "p_val", "p_adj"]
    for c in [
        "mixedlm_coef",
        "delta_variance",
        "delta_var_pval",
        "memento_de_se",
        "memento_dv_coef",
        "memento_dv_se",
        "memento_dv_pval",
        "memento_p_adj_native",
    ]:
        if c in adata.var.columns:
            cols.append(c)

    # Add a simple base expression measure when possible
    if "total_us_counts" in adata.var.columns:
        cols.append("total_us_counts")
    else:
        # fallback: mean of current X (after any preprocess the user chose)
        try:
            means = np.asarray(adata.X.mean(axis=0)).ravel()
            adata.var["baseMean"] = means
            cols.append("baseMean")
        except Exception:
            pass

    results = adata.var[cols].copy()
    results = results.sort_values("p_adj", ascending=True)
    # Preserve backend diagnostics attrs (pandas often drops them on reindex/copy)
    if hasattr(de_df, "attrs") and de_df.attrs:
        results.attrs.update(dict(de_df.attrs))

    # Metadata — merge to preserve raw_gene_list etc. from store_raw_counts()
    existing = dict(adata.uns.get("scatrans", {}))
    history = existing.get("history", [])
    if "analysis" in existing:
        prev = {
            k: existing.get(k)
            for k in ("analysis", "mode", "target_group", "reference_group")
            if k in existing
        }
        if prev:
            history.append(prev)
            if len(history) > 5:
                history = history[-5:]
    existing["history"] = history

    de_diagnostics: dict[str, Any] = {
        "n_cells": int(adata.n_obs),
        "n_genes_input": int(adata.n_vars),
        "mixed_model": {
            "used": bool(use_mixed_model),
            "sample_col": sample_col if use_mixed_model else None,
            "paired_replicates": paired_replicates if use_mixed_model else None,
            "n_samples": int(adata.obs[sample_col].nunique())
            if (use_mixed_model and sample_col and sample_col in adata.obs.columns)
            else None,
            "mixedlm_grouping": (
                de_df.attrs.get("mixedlm_grouping") if hasattr(de_df, "attrs") else None
            )
            if use_mixed_model
            else None,
            "delta_variance_available": "delta_variance" in adata.var.columns,
            "median_delta_variance": float(np.nanmedian(adata.var["delta_variance"]))
            if "delta_variance" in adata.var.columns
            else np.nan,
            "n_genes_failed_fit": n_mixed_failed_de if use_mixed_model else 0,
            "failed_fit_rate": mixed_failed_rate_de if use_mixed_model else 0.0,
            "logFC_method": mixedlm_logfc_method_de if use_mixed_model else None,
            "n_genes_logFC_mixedlm_sign_discordant": (
                n_mixed_sign_discordant_de if use_mixed_model else 0
            ),
            "note": (
                "Lightweight LMM analogue (log1p + Wald/LRT); not NB-GLMM/voom. "
                "logFC is sample-aware mean-of-means log2FC; p_val tests mixedlm_coef. "
                "Inspect failed_fit_rate and n_genes_logFC_mixedlm_sign_discordant "
                "before publication claims."
                if use_mixed_model
                else None
            ),
        },
        "pydeseq2": pydeseq2_diag or {"used": False},
        "memento": {
            "used": bool(use_memento_de),
            "n_genes_not_returned": n_memento_not_returned,
            "note": (
                "Genes dropped by memento internal filters appear as logFC=0, p_adj=1 "
                "after reindexing and were not tested."
                if use_memento_de
                else None
            ),
        },
    }

    from .._utils import _merge_scatrans_uns

    de_meta = {
        "mode": "differential_expression",
        "analysis": "differential_expression",
        "version": VERSION,
        "groupby": groupby,
        "target_group": target_group,
        "reference_group": reference_group,
        "use_pseudobulk": use_pseudobulk,
        "use_mixed_model": use_mixed_model,
        "use_memento_de": use_memento_de,
        "memento_capture_rate": memento_capture_rate if use_memento_de else None,
        "de_method": de_method,
        "pseudobulk_de_backend": pseudobulk_de_backend,
        "de_preprocess": de_preprocess,
        "strict_pydeseq2_counts": strict_pydeseq2_counts,
        "min_cells": min_cells,
        "min_counts": min_counts,
        "min_total_counts": min_total_counts,
        "sample_col": sample_col if (use_mixed_model or use_pseudobulk) else None,
        "pb_x_layer": pb_x_layer if use_pseudobulk else None,
        "pb_use_total_for_x": pb_use_total_for_x if use_pseudobulk else None,
        "use_delta_variance_pval": use_delta_variance_pval,
        "delta_var_pval_cutoff": delta_var_pval_cutoff,
        "mixed_model_pval": mixed_model_pval if use_mixed_model else None,
        "n_genes_failed_mixed_fit": n_mixed_failed_de if use_mixed_model else None,
        "failed_fit_rate_mixed": mixed_failed_rate_de if use_mixed_model else None,
        "memento_has_native_padj": bool(use_memento_de and "memento_p_adj_native" in de_df.columns)
        if use_memento_de
        else None,
        "memento_num_boot": memento_num_boot if use_memento_de else None,
        "memento_n_cpus": memento_n_cpus if use_memento_de else None,
        "n_genes_not_returned_by_memento": n_memento_not_returned if use_memento_de else None,
        "pydeseq2": pydeseq2_diag or None,
        "n_jobs": n_jobs,
        "gene_type_filter": gene_type_filter,
        "diagnostics": de_diagnostics,
        "history": existing.get("history", []),
    }
    adata.uns["scatrans"] = _merge_scatrans_uns(existing, de_meta)

    logger.info("DE completed. %d genes in results table.", len(results))
    return adata, results

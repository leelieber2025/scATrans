"""scatrans.tl.de — internal package module."""

from __future__ import annotations

import logging
import os
import warnings
from typing import Any

import anndata as ad
import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sparse

from .._de import _run_de_wrapper
from .._utils import (
    _align_snapshot_bundle,
    _apply_de_preprocess,
    _as_contrast_categorical,
    _clear_log_preprocess_metadata,
    _get_raw_snapshot,
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
    """Write scatrans metadata (raw_gene_list) and, optionally, adata.raw."""
    if "scatrans" not in adata.uns:
        adata.uns["scatrans"] = {}
    prev_list = adata.uns["scatrans"].get("raw_gene_list")
    n_genes = int(adata.n_vars)
    if prev_list is not None and len(prev_list) != n_genes:
        logger.warning(
            "Updating raw_gene_list (%d → %d genes). If you subsetted after an earlier "
            "store_raw_counts(), enrichment universe will reflect the current gene set only. "
            "For the full pre-HVG universe, call store_raw_counts(sidecar=True) on the "
            "original object before subsetting and recover it with "
            "restore_raw_counts(full_genes=True).",
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


def _store_raw_snapshot(
    adata: Any,
    mat: Any,
    *,
    overwrite: bool = False,
    ondisk: bool = False,
    snapshot_path: str | None = None,
) -> None:
    """Write a label-indexed raw count snapshot into ``uns['scatrans']``.

    The snapshot stores the full obs x var count matrix together with the complete
    ``obs_names`` and ``var_names``. Because it lives in ``uns`` (not an axis-aligned
    slot) it survives gene subsetting (HVG), cell subsetting, ``copy()`` and
    ``write_h5ad()``. Restoration re-aligns it by label, so the full-gene universe
    and the correct cells can be recovered even after subsetting or reordering.

    Any velocity layers present (``spliced``/``unspliced`` or ``mature``/``nascent``)
    are captured into the snapshot too, so :func:`restore_raw_counts` can recover the
    full-gene velocity matrices for active-transcription analysis after HVG subsetting.

    When ``ondisk`` is True the full matrix is written to ``snapshot_path`` (a
    standalone ``.h5ad``) and only a lightweight pointer (path + names) is kept in
    ``uns`` — useful to avoid doubling the in-memory / in-file count matrix for large
    datasets. The referenced file must remain reachable for later restoration.
    """
    if "scatrans" not in adata.uns:
        adata.uns["scatrans"] = {}
    if adata.uns["scatrans"].get("raw_snapshot") is not None and not overwrite:
        logger.debug("raw_snapshot already exists; skipping (pass overwrite=True to replace).")
        return
    if mat.shape[0] != adata.n_obs or mat.shape[1] != adata.n_vars:
        raise ValueError(
            f"Cannot store raw_snapshot: matrix shape {mat.shape} does not match "
            f"adata shape ({adata.n_obs}, {adata.n_vars})."
        )
    # Preserve the original matrix format (dense stays dense, sparse stays sparse)
    # so restoration returns the same type the caller started with.
    stored = mat.copy()
    obs_names = adata.obs_names.to_numpy().astype(str)
    var_names = adata.var_names.to_numpy().astype(str)
    is_integer = bool(_is_integer_counts_like(stored))

    # Capture velocity layers so they survive HVG/cell subsetting like the counts.
    vel_layers = {
        name: adata.layers[name].copy()
        for name in ("spliced", "unspliced", "mature", "nascent")
        if name in adata.layers
    }

    if ondisk:
        if not snapshot_path:
            raise ValueError(
                "sidecar='ondisk' requires snapshot_path=<file.h5ad> to write the full "
                "count matrix to."
            )
        path = os.path.abspath(snapshot_path)
        if os.path.exists(path) and not overwrite:
            logger.debug("On-disk snapshot %s already exists; reusing (pass overwrite=True).", path)
        else:
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
            snap_ad = ad.AnnData(
                X=stored,
                obs=pd.DataFrame(index=pd.Index(obs_names)),
                var=pd.DataFrame(index=pd.Index(var_names)),
            )
            for name, layer_mat in vel_layers.items():
                snap_ad.layers[name] = layer_mat
            snap_ad.write_h5ad(path)
            logger.info("Wrote on-disk raw snapshot to %s.", path)
        adata.uns["scatrans"]["raw_snapshot"] = {
            "backend": "ondisk",
            "path": path,
            "obs_names": obs_names,
            "var_names": var_names,
            "is_integer": is_integer,
            "velocity_layers": list(vel_layers),
            "version": VERSION,
        }
        logger.info(
            "Recorded on-disk raw snapshot pointer (%d cells x %d genes%s) in "
            "adata.uns['scatrans']['raw_snapshot']; keep %s reachable to restore.",
            adata.n_obs,
            adata.n_vars,
            f", velocity layers: {list(vel_layers)}" if vel_layers else "",
            path,
        )
        return

    adata.uns["scatrans"]["raw_snapshot"] = {
        "backend": "inline",
        "X": stored,
        "obs_names": obs_names,
        "var_names": var_names,
        "is_integer": is_integer,
        "layers": vel_layers,
        "version": VERSION,
    }
    logger.info(
        "Stored label-indexed raw snapshot (%d cells x %d genes%s) in "
        "adata.uns['scatrans']['raw_snapshot']; it survives HVG/cell subsetting.",
        adata.n_obs,
        adata.n_vars,
        f", velocity layers: {list(vel_layers)}" if vel_layers else "",
    )


def _sidecar_mode(sidecar: bool | str) -> tuple[bool, bool]:
    """Normalize the ``sidecar`` argument into ``(enabled, ondisk)``.

    Accepts ``True``/``False`` (in-memory snapshot in ``uns``) or the string
    ``"ondisk"`` (write the full matrix to a standalone ``.h5ad`` and keep only a
    pointer in ``uns``).
    """
    if isinstance(sidecar, str):
        if sidecar == "ondisk":
            return True, True
        raise ValueError(f"Unknown sidecar mode {sidecar!r}; use True, False, or 'ondisk'.")
    return bool(sidecar), False


def ensure_raw_counts(
    adata: Any,
    layer: str = "counts",
    save_raw: bool = False,
    overwrite: bool = False,
    sidecar: bool | str = True,
    snapshot_path: str | None = None,
) -> None:
    """
    Deprecated alias for :func:`store_raw_counts` with ``mode="auto"``.

    ``store_raw_counts(mode="auto")`` performs the same idempotent recovery this
    function used to: it reuses an existing integer counts layer, stores the current
    ``.X`` when it looks like counts, or recovers counts from ``adata.raw`` when
    ``.X`` is already normalized/log-transformed. Prefer calling that directly.
    """
    warnings.warn(
        "ensure_raw_counts() is deprecated; call store_raw_counts(..., mode='auto') instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    store_raw_counts(
        adata,
        layer=layer,
        save_raw=save_raw,
        overwrite=overwrite,
        sidecar=sidecar,
        snapshot_path=snapshot_path,
        mode="auto",
    )


def store_raw_counts(
    adata: Any,
    layer: str = "counts",
    save_raw: bool = False,
    overwrite: bool = False,
    sidecar: bool | str = True,
    snapshot_path: str | None = None,
    mode: str = "force",
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

    Any existing velocity layers (``spliced``/``unspliced`` or ``mature``/``nascent``)
    are captured into the sidecar snapshot (see below), so they survive HVG/cell
    subsetting and can be recovered by :func:`restore_raw_counts`. (Earlier versions
    instead wrote position-aligned ``raw_spliced``/``raw_unspliced`` layers, which were
    trimmed by subsetting and never read back; those are no longer created.)

    Sidecar snapshot (``sidecar=True``, default)
    --------------------------------------------
    In addition to the axis-aligned ``layer``, a label-indexed snapshot of the full
    obs x var count matrix is written to ``adata.uns['scatrans']['raw_snapshot']``.
    Because ``uns`` is not tied to the obs/var axes, this snapshot survives HVG gene
    subsetting, cell subsetting, ``copy()`` and ``write_h5ad()``. Use
    :func:`restore_raw_counts` (optionally with ``full_genes=True``) to recover the
    full-gene universe and the correct cells afterwards, aligned by name. This is the
    recommended way to keep raw counts available for DE / enrichment after HVG or
    after extracting a subset. Pass ``sidecar=False`` to skip it (e.g. to save memory
    when you keep the full object around yourself).

    For large datasets, pass ``sidecar='ondisk'`` together with
    ``snapshot_path='raw_snapshot.h5ad'`` to write the full count matrix to a
    standalone file and keep only a lightweight pointer in ``uns`` (avoids doubling
    the count matrix in memory and inside the main ``.h5ad``). The referenced file
    must stay reachable for :func:`restore_raw_counts` to work.

    If you need the full-gene raw velocity data after HVG-based visualization,
    either:

      - rely on the sidecar snapshot and restore_raw_counts(full_genes=True)
        (recommended), or

      - call store_raw_counts() on the full object and keep the full object for
        DE / active_score / enrichment while using a copy for visualization.

    ``save_raw`` is deprecated: adata.raw is commonly reserved for log-normalized
    data, and the sidecar snapshot already preserves full-gene raw counts.

    mode : {"force", "auto"}, default "force"
        "force" stores the current ``adata.X`` into the counts layer (the classic
        behavior). "auto" is idempotent and recovery-aware: it reuses an existing
        aligned integer counts layer, or recovers counts from ``adata.raw`` when
        ``adata.X`` is already normalized/log-transformed, before falling back to
        storing the current ``.X``. (This is what the deprecated
        :func:`ensure_raw_counts` did.)

    Recommended early call:
        scat.store_raw_counts(adata, layer="counts", save_raw=False)
    """
    if mode not in ("force", "auto"):
        raise ValueError(f"Unknown mode {mode!r}; use 'force' or 'auto'.")

    if save_raw:
        warnings.warn(
            "save_raw=True is deprecated. adata.raw is commonly reserved for "
            "log-normalized data (scanpy convention), so writing raw integer counts "
            "there is ambiguous. The sidecar snapshot (sidecar=True, the default) "
            "already preserves full-gene raw counts across HVG/cell subsetting; use "
            "restore_raw_counts(..., full_genes=True) to recover them.",
            DeprecationWarning,
            stacklevel=2,
        )

    enabled, ondisk = _sidecar_mode(sidecar)

    def _finalize() -> None:
        if enabled and layer in adata.layers:
            _store_raw_snapshot(
                adata,
                adata.layers[layer],
                overwrite=overwrite,
                ondisk=ondisk,
                snapshot_path=snapshot_path,
            )
        _record_raw_counts_metadata(adata, save_raw=save_raw, overwrite=overwrite)

    if mode == "auto":
        # 1. Existing layer already holds aligned integer counts -> just refresh sidecar/metadata.
        if (
            layer in adata.layers
            and not overwrite
            and _is_integer_counts_like(adata.layers[layer])
            and adata.layers[layer].shape[1] == adata.n_vars
        ):
            _finalize()
            logger.debug("Layer '%s' already holds aligned integer counts.", layer)
            return
        # 2. .X is not counts -> try to recover from adata.raw before the force store.
        if not _is_integer_counts_like(adata.X):
            raw = getattr(adata, "raw", None)
            if raw is not None and _is_integer_counts_like(raw.X) and raw.shape[1] == adata.n_vars:
                if hasattr(raw, "var_names") and np.array_equal(raw.var_names, adata.var_names):
                    adata.layers[layer] = raw.X.copy()
                    logger.info(
                        "store_raw_counts(mode='auto'): recovered raw counts from "
                        "adata.raw into layers['%s'].",
                        layer,
                    )
                    _finalize()
                    return
                logger.warning(
                    "adata.raw exists but gene names/order do not match current adata.var_names. "
                    "Cannot recover counts automatically."
                )
            logger.warning(
                "store_raw_counts(mode='auto'): adata.X does not look like raw counts and "
                "adata.raw could not be used. Storing current .X (may warn again)."
            )
        # 3. .X looks like counts -> fall through to the force store below.

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

    _finalize()


def _restore_full_genes_from_snapshot(adata: Any) -> Any:
    """Build a new AnnData with the full pre-HVG gene universe from the uns snapshot.

    Rows are the current (possibly subsetted) cells, aligned by ``obs_names``; columns
    are the complete stored gene set. Obs-aligned annotations (``obs``, ``obsm``) are
    carried over; var-aligned annotations are dropped because the dropped genes have no
    current metadata. ``uns`` is shallow-copied and ``raw_gene_list`` is updated to the
    full gene set so downstream enrichment uses the correct universe.
    """
    aligned = _align_snapshot_bundle(adata, full_genes=True)
    if aligned is None:
        raise ValueError(
            "full_genes=True requires a raw snapshot in adata.uns['scatrans']['raw_snapshot']. "
            "Call scat.store_raw_counts(adata, sidecar=True) before HVG subsetting."
        )
    X_full, var_names, vel_layers = aligned
    var_df = pd.DataFrame(index=pd.Index(var_names, name=adata.var_names.name))
    new = ad.AnnData(X=X_full, obs=adata.obs.copy(), var=var_df)
    for name, layer_mat in vel_layers.items():
        new.layers[name] = layer_mat
    for key in adata.obsm:
        new.obsm[key] = adata.obsm[key].copy()

    new_uns = dict(adata.uns)
    scat = new_uns.get("scatrans")
    if isinstance(scat, dict):
        scat = dict(scat)
        scat["raw_gene_list"] = list(var_names)
        scat["store_raw_counts_n_genes"] = int(len(var_names))
        new_uns["scatrans"] = scat
    new.uns = new_uns

    _clear_log_preprocess_metadata(new)
    logger.info(
        "Restored full-gene raw counts from snapshot into a new AnnData (%d cells x %d genes).",
        new.n_obs,
        new.n_vars,
    )
    return new


def restore_raw_counts(
    adata: Any,
    layer: str = "counts",
    inplace: bool = False,
    full_genes: bool = False,
    prefer_snapshot: bool = True,
) -> Any | None:
    """
    Restore raw counts from the uns snapshot (preferred), the stored layer, or
    adata.raw back into .X.

    This is useful when you have done HVG + log1p on .X for visualization,
    but want to work with (or pass to other tools) the raw counts for the
    genes currently in the adata (or the full pre-HVG universe).

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
        Ignored (forced False) when ``full_genes=True``, since the gene axis changes.
    full_genes : bool
        If True, recover the complete pre-HVG gene set from the uns snapshot and
        return a new AnnData with all genes (for DE / enrichment on the full
        universe). Requires ``store_raw_counts(sidecar=True)`` to have been called
        before subsetting. Cannot be combined with ``inplace=True``.
    prefer_snapshot : bool
        If True (default), use the label-indexed uns snapshot when present. The
        snapshot aligns by cell/gene name, so it also recovers counts correctly after
        cell subsetting or gene reordering. Set to False to force the legacy
        layer/adata.raw path. Note that adata.raw is only used as a fallback when it
        looks like integer counts; a log-normalized adata.raw is never restored as .X.

    Returns
    -------
    AnnData or None
        If not inplace (or full_genes=True), a new AnnData with raw counts in .X.
    """
    if full_genes:
        if inplace:
            raise ValueError(
                "full_genes=True changes the gene axis and cannot be done inplace; "
                "use inplace=False (the default) to get a new full-gene AnnData."
            )
        return _restore_full_genes_from_snapshot(adata)

    if prefer_snapshot and _get_raw_snapshot(adata) is not None:
        aligned = _align_snapshot_bundle(adata, full_genes=False)
        if aligned is not None:
            snap_mat, _, vel_layers = aligned
            target = adata if inplace else adata.copy()
            target.X = snap_mat.copy() if hasattr(snap_mat, "copy") else snap_mat
            for name, layer_mat in vel_layers.items():
                target.layers[name] = layer_mat
            _clear_log_preprocess_metadata(target)
            if inplace:
                logger.info("Restored raw counts from uns raw_snapshot into adata.X (inplace).")
                return None
            logger.info("Created copy with raw counts from uns raw_snapshot in .X.")
            return target

    raw_attr = getattr(adata, "raw", None)
    if layer in adata.layers:
        raw = adata.layers[layer].copy()
        source = f"layers['{layer}']"
    elif raw_attr is not None and _is_integer_counts_like(raw_attr.X):
        # Only trust adata.raw as a counts source when it actually looks like integer
        # counts. By the common scanpy convention adata.raw holds log-normalized data,
        # which must not be restored into .X as if it were raw counts.
        raw = raw_attr.X.copy()
        source = "adata.raw"
    else:
        detail = ""
        if raw_attr is not None:
            detail = (
                " adata.raw is present but does not look like integer counts (it is "
                "assumed to hold log-normalized data) and was not used."
            )
        raise ValueError(
            f"No usable raw counts found in layer '{layer}' or the uns snapshot.{detail} "
            "Call scat.store_raw_counts(adata) early (sidecar=True) to preserve them."
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

"""Permutation-calibrated program placement on the mechanism axis.

``program_mechanism`` and :func:`~scatrans.program_mechanism_induction_matched` are
*relative*: they rank a gene set against a background or against matched non-members.
Neither supplies a trustworthy zero, because ``transcription_support`` is standardized
genome-wide while any particular gene set has its own gene-structure composition
(length, intron count, expression), so its expectation under "no condition effect" is
not exactly zero.

This module supplies that zero empirically. The calibrated program score is

    calibrated = observed program mean − mean program score under shuffled condition labels

where each shuffled replicate is pushed through the **identical** transformation chain as
the observation (same residual estimator, same gene-structure bias correction, same
robust-z, same gene universe). Whatever structural offset the gene set carries is present
in both terms and cancels, so the calibrated value can be read as an absolute displacement
on the mechanism axis rather than merely a rank against another program.

Gene-set *membership* must not vary across replicates, otherwise the null would absorb
DE-selection noise instead of the structural offset. Membership is therefore frozen by
requiring a precomputed DE table (``de=``), which is reused unchanged for every replicate;
only the condition labels move.

Permutation p-values use the (b+1)/(n+1) estimator of Phipson & Smyth (2010), so they are
never zero.

.. note::
   This is a **program-level** procedure. It does not make per-gene labels identifiable;
   see the identifiability limits documented on :func:`annotate_mechanism_class`.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .._utils import _normalize_group_label

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ["program_mechanism_permutation_calibrated"]


def _shuffle_within_blocks(
    labels: np.ndarray,
    blocks: np.ndarray | None,
    rng: np.random.Generator,
    *,
    mask: np.ndarray | None = None,
) -> np.ndarray:
    """Permute ``labels`` among rows in ``mask`` (default: every row).

    When ``blocks`` is given, permute within each block, still restricted to
    ``mask``. Rows outside ``mask`` keep their original labels. Blocking
    preserves condition-by-batch or donor imbalance in the null.
    """
    out = np.array(labels, copy=True)
    if mask is None:
        eligible = np.ones(len(out), dtype=bool)
    else:
        eligible = np.asarray(mask, dtype=bool)
        if eligible.shape != out.shape:
            raise ValueError("mask must have the same length as labels")
    if blocks is None:
        idx = np.flatnonzero(eligible)
        if idx.size > 1:
            out[idx] = rng.permutation(out[idx])
        return out
    # pd.unique keeps NA; ``blocks == na`` is never True, so skip NA keys and
    # leave NA-blocked rows fixed (cannot form an exchangeability class).
    for b in pd.unique(blocks):
        try:
            if b is None or bool(pd.isna(b)):
                continue
        except (TypeError, ValueError):
            pass
        idx = np.flatnonzero((blocks == b) & eligible)
        if idx.size > 1:
            out[idx] = rng.permutation(out[idx])
    return out


def program_mechanism_permutation_calibrated(
    adata: Any,
    gene_sets: Mapping[str, Sequence[str]],
    *,
    de: pd.DataFrame,
    groupby: str = "condition",
    target_group: str = "Disease",
    reference_group: str = "Control",
    n_perm: int = 200,
    random_state: int = 0,
    block_col: str | None = None,
    support_col: str = "transcription_support",
    min_genes: int = 5,
    restrict_to_selected: bool = False,
    **partition_kwargs: Any,
) -> pd.DataFrame:
    """Place gene programs on the mechanism axis against a label-permutation null.

    Parameters
    ----------
    adata
        Object carrying ``spliced``/``unspliced`` layers and ``obs[groupby]``.
    gene_sets
        ``{program_name: [gene, ...]}``. Genes are matched against the gene table index.
    de
        Precomputed DE table. **Required**: it freezes gene-set membership and the gene
        universe across replicates so the null isolates the structural offset rather than
        DE-selection noise. Pass the same table used for the observed partition.
    groupby, target_group, reference_group
        Condition column and the two groups contrasted.
    n_perm
        Number of label shuffles. 200 gives a smallest attainable p of 1/201 ≈ 0.005.
    random_state
        Seed for the permutation stream.
    block_col
        Optional ``obs`` column (e.g. donor, batch) within which labels are shuffled.
        Only cells whose condition is ``target_group`` or ``reference_group`` are
        shuffled; any other ``groupby`` value is left unchanged.
    support_col
        Column holding per-gene transcription support in the partition gene table.
    min_genes
        Skip programs with fewer than this many genes present in the gene table.
    restrict_to_selected
        When True, each program is intersected with the DE-selected genes of the observed
        run, so the calibrated value refers to the same membership the program figures
        show. Default False uses every gene set member with finite support. The
        intersection is taken **once, from the observed run**, and held fixed across
        replicates — recomputing it per replicate would let DE-selection noise leak into
        the null.
    **partition_kwargs
        Forwarded to :func:`~scatrans.partition_de_by_mechanism` (e.g. ``organism``,
        ``logfc_cutoff``). Must be identical between observed and null by construction —
        the same values are used for every replicate.

    Returns
    -------
    pandas.DataFrame
        One row per program with ``n_genes``, ``observed_mean``, ``null_mean``,
        ``null_sd``, ``calibrated`` (observed − null mean), ``z``, ``p_perm`` and
        ``n_perm_effective``. Negative ``calibrated`` = stabilization-weighted,
        positive = transcription-weighted.

    Notes
    -----
    Cost is one residual recomputation per replicate; DE is not re-run, because ``de`` is
    reused. Runtime therefore scales with ``n_perm`` and the number of cells, not with the
    DE backend.
    """
    from .partition import partition_de_by_mechanism

    if de is None or not isinstance(de, pd.DataFrame):
        raise TypeError(
            "de must be a precomputed DE DataFrame so that gene-set membership is "
            "frozen across permutations; pass the table used for the observed partition."
        )
    if groupby not in adata.obs:
        raise KeyError(f"groupby={groupby!r} not in adata.obs")
    # Match active_score / partition label handling (strip, "1.0"↔"1", etc.).
    # Validate on normalized labels so callers can pass target_group=1 against
    # obs values stored as "1.0", but shuffle the *original* labels so the
    # multiset written back into obs is unchanged aside from order.
    target_norm = _normalize_group_label(target_group)
    reference_norm = _normalize_group_label(reference_group)
    if target_norm is None or reference_norm is None:
        raise ValueError("target_group and reference_group must be valid labels")
    if target_norm == reference_norm:
        raise ValueError("target_group and reference_group must be different")
    labels = np.asarray(adata.obs[groupby].values)
    norm_labels = pd.Series(labels).map(_normalize_group_label).to_numpy()
    for g_raw, g_norm in ((target_group, target_norm), (reference_group, reference_norm)):
        if not np.any(norm_labels == g_norm):
            raise ValueError(f"group {g_raw!r} absent from adata.obs[{groupby!r}]")
    # Only the two contrast levels are exchangeable. Other values in groupby
    # (extra conditions, unused cell types) stay fixed.
    contrast_mask = (norm_labels == target_norm) | (norm_labels == reference_norm)
    if block_col is not None and block_col not in adata.obs:
        raise KeyError(f"block_col={block_col!r} not in adata.obs")
    blocks = np.asarray(adata.obs[block_col].values) if block_col else None
    if n_perm < 1:
        raise ValueError("n_perm must be >= 1")

    # Explicit contrast args win over any accidental collisions in partition_kwargs.
    kw = {
        **partition_kwargs,
        "groupby": groupby,
        "target_group": target_group,
        "reference_group": reference_group,
        "de": de,
    }

    def _run(a: Any):
        res = partition_de_by_mechanism(a, **kw)
        gt = res.gene_table
        if support_col not in gt:
            raise KeyError(f"{support_col!r} missing from the partition gene table")
        return res, pd.to_numeric(gt[support_col], errors="coerce")

    observed_res, observed = _run(adata)

    # DE-selected membership is taken once, from the observed run, then frozen: letting it
    # vary per replicate would put DE-selection noise into the null instead of the
    # structural offset the calibration is meant to isolate.
    selected_index = None
    if restrict_to_selected:
        sel = getattr(observed_res, "selected", None)
        if sel is None or not len(sel):
            raise ValueError(
                "restrict_to_selected=True but the observed run selected no genes; "
                "relax the DE gates or set restrict_to_selected=False."
            )
        selected_index = set(sel.index)

    programs: dict[str, list[str]] = {}
    for name, genes in gene_sets.items():
        present = [g for g in dict.fromkeys(genes) if g in observed.index]
        present = [g for g in present if np.isfinite(observed.get(g, np.nan))]
        if selected_index is not None:
            present = [g for g in present if g in selected_index]
        if len(present) >= min_genes:
            programs[name] = present
        else:
            logger.warning(
                "program %s skipped: %d genes with finite support (< min_genes=%d)",
                name,
                len(present),
                min_genes,
            )
    if not programs:
        return pd.DataFrame(
            columns=[
                "n_genes",
                "observed_mean",
                "null_mean",
                "null_sd",
                "calibrated",
                "z",
                "p_perm",
                "n_perm_effective",
            ]
        )

    obs_mean = {k: float(observed.loc[v].mean()) for k, v in programs.items()}

    rng = np.random.default_rng(random_state)
    null: dict[str, list[float]] = {k: [] for k in programs}
    n_failed = 0
    perm_ad = adata.copy()
    for i in range(n_perm):
        # Preserve original label dtype/representation; only order changes.
        perm_ad.obs[groupby] = _shuffle_within_blocks(labels, blocks, rng, mask=contrast_mask)
        try:
            _, sup = _run(perm_ad)
        except Exception as exc:  # a degenerate shuffle must not kill the run
            n_failed += 1
            logger.warning("permutation %d failed (%s); skipped", i, type(exc).__name__)
            continue
        for k, v in programs.items():
            vals = sup.reindex(v).to_numpy(dtype=float)
            # Same gene list as the observation; all-NaN → drop this replicate.
            null[k].append(float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else np.nan)

    if n_failed == n_perm:
        logger.warning("all %d permutations failed; calibrated scores will be NaN", n_perm)

    rows = []
    for k, v in programs.items():
        arr = np.asarray([x for x in null[k] if np.isfinite(x)], dtype=float)
        n_eff = int(arr.size)
        o = obs_mean[k]
        if n_eff == 0:
            rows.append(
                {
                    "program": k,
                    "n_genes": len(v),
                    "observed_mean": o,
                    "null_mean": np.nan,
                    "null_sd": np.nan,
                    "calibrated": np.nan,
                    "z": np.nan,
                    "p_perm": np.nan,
                    "n_perm_effective": 0,
                }
            )
            continue
        m = float(arr.mean())
        s = float(arr.std(ddof=1)) if n_eff > 1 else np.nan
        calibrated = o - m
        # two-sided empirical p, centered on the null mean; Phipson and Smyth (2010)
        b = int(np.sum(np.abs(arr - m) >= abs(calibrated)))
        z = calibrated / s if (np.isfinite(s) and s > 0) else np.nan
        rows.append(
            {
                "program": k,
                "n_genes": len(v),
                "observed_mean": o,
                "null_mean": m,
                "null_sd": s,
                "calibrated": calibrated,
                "z": z,
                "p_perm": (b + 1) / (n_eff + 1),
                "n_perm_effective": n_eff,
            }
        )
    out = pd.DataFrame(rows).set_index("program")
    out.attrs["n_perm_requested"] = n_perm
    out.attrs["n_perm_failed"] = n_failed
    out.attrs["block_col"] = block_col
    return out

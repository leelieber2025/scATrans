"""Pseudobulk nascent-activity score (variance-stabilized) + DE reproducibility.

The default ``unspliced_excess_residual`` residualizes nascent excess on log(gene
length) and log(intron number). Length/intron correction can deflate long or
high-intron genes that are genuinely transcriptionally active. This module
provides an additive **detection** signal that:

* is computed on **pseudobulk** (target vs reference group sums);
* uses a **variance-stabilized Poisson-z** of the nascent increase
  (``nascent_poisson_z``);
* optionally marks spliced-side **DE reproducibility**
  (``de_reproducible`` / ``de_repro_frac``) — annotation only, never membership
  gating. The flag cannot distinguish a DE false positive from a genuine
  stabilization-driven gene (both can show near-zero nascent excess).

Detection is a different question from the transcription-vs-stabilization
**mechanism** partition. The Poisson-z is an absolute, induction-coupled nascent
increase, so it must **not** drive mechanism labels (highly induced
stabilization targets would be mis-read as transcription-driven). Mechanism
annotation stays on the induction-normalized residual
(:func:`~scatrans.tl.annotate_mechanism_class`). Opt in via
``partition_de_by_mechanism(add_nascent_score=True)`` (appends columns only).

On low-capture data the signal is capture-limited; inspect
:func:`~scatrans.qc.regime_diagnosis` first.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd

from .._utils import _validate_group_contrast
from ._common import _resolve_velocity_layer_keys

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

NASCENT_Z_COL = "nascent_poisson_z"
DLOG_U_COL = "dlog_unspliced"
DLOG_S_COL = "dlog_spliced"
REPRO_COL = "de_reproducible"
REPRO_FRAC_COL = "de_repro_frac"

_CPM = 1e6


def _layer(adata: Any, name: str) -> Any:
    if name not in adata.layers:
        raise KeyError(
            f"layer {name!r} not found; nascent_activity_score needs spliced/unspliced "
            f"velocity layers (available: {list(adata.layers)})"
        )
    return adata.layers[name]


def _col_sum(mat: Any, mask: np.ndarray) -> np.ndarray:
    """Per-gene sum over the masked rows (works for sparse or dense)."""
    sub = mat[mask]
    s = sub.sum(axis=0)
    return np.asarray(s).ravel()


def _dlog(u_t: np.ndarray, u_r: np.ndarray) -> np.ndarray:
    """log fold-change of CPM-normalized pseudobulk (target vs reference)."""
    nt = u_t.sum()
    nr = u_r.sum()
    ct = u_t / nt * _CPM if nt > 0 else u_t
    cr = u_r / nr * _CPM if nr > 0 else u_r
    return np.log((ct + 1.0) / (cr + 1.0))


def _fold_labels(
    n: int, groups: np.ndarray | None, n_splits: int, rng: np.random.Generator
) -> np.ndarray:
    """Assign each cell to one of ``n_splits`` folds — by sample when available."""
    if groups is not None:
        uniq = pd.unique(groups)
        if len(uniq) >= n_splits:
            fold_of = {g: i % n_splits for i, g in enumerate(uniq)}
            return np.array([fold_of[g] for g in groups])
    return rng.integers(0, n_splits, size=n)


def nascent_activity_score(
    adata: Any,
    groupby: str = "condition",
    target_group: str | None = None,
    reference_group: str | None = None,
    *,
    unspliced_layer: str | None = None,
    spliced_layer: str | None = None,
    sample_col: str | None = None,
    n_splits: int = 2,
    random_state: int = 0,
) -> pd.DataFrame:
    """Pseudobulk variance-stabilized nascent-activity score + DE reproducibility.

    Parameters
    ----------
    adata
        AnnData with spliced/unspliced (velocity) layers and a ``groupby`` obs column.
    groupby, target_group, reference_group
        Contrast definition (same convention as :func:`~scatrans.tl.active_score`).
        Both groups are required — treatment direction is never guessed (a wrong
        guess would flip the sign of the nascent increase).
    unspliced_layer, spliced_layer
        Velocity layer names. When ``None`` (default) they are auto-resolved
        (``spliced``/``unspliced`` or kb_python ``mature``/``nascent``), matching
        the rest of the package.
    sample_col
        Optional obs column of biological samples/replicates. When given it must
        exist (a wrong name raises, not silently falls back). With at least
        ``n_splits`` samples per group the reproducibility folds are the samples
        themselves (proper cross-replicate check); otherwise cells are split randomly.
    n_splits
        Number of reproducibility folds (default 2).
    random_state
        Seed for the random fold fallback.

    Returns
    -------
    pandas.DataFrame
        Indexed by ``adata.var_names`` with columns:

        - ``nascent_poisson_z`` — variance-stabilized nascent increase as a two-sample
          Poisson / conditional-binomial score ``(U_t - n*f) / sqrt(n*f*(1-f) + 1)``,
          ``n = U_t + U_r``, ``f`` = target depth-exposure fraction
          ``tot_t / (tot_t + tot_r)``. Symmetric and bounded (a gene absent in the
          reference does not get an unbounded z). The active-transcription signal
          (higher = more nascent gain).
        - ``dlog_unspliced`` / ``dlog_spliced`` — CPM log fold-changes (diagnostic).
        - ``de_reproducible`` — bool: the spliced fold-change keeps its sign across all
          reproducibility folds (proxy-independent DE-false-positive flag).
        - ``de_repro_frac`` — fraction of folds agreeing with the overall spliced sign.
    """
    if n_splits < 1:
        raise ValueError(f"n_splits must be >= 1, got {n_splits}")
    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby={groupby!r} not in adata.obs")
    if target_group is None or reference_group is None:
        raise ValueError(
            "target_group and reference_group are required (treatment direction "
            "cannot be inferred safely — a wrong guess flips the nascent sign)"
        )
    if sample_col is not None and sample_col not in adata.obs.columns:
        raise KeyError(
            f"sample_col={sample_col!r} not in adata.obs (pass None for the "
            "within-sample random-fold fallback)"
        )

    # normalized group matching + clear label errors, consistent with active_score/DE
    target_norm, reference_norm, norm_groups = _validate_group_contrast(
        adata.obs[groupby],
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
    )
    labels = norm_groups.to_numpy()
    mt = labels == target_norm
    mr = labels == reference_norm
    n_t, n_r = int(mt.sum()), int(mr.sum())
    if n_t < 2 or n_r < 2:
        raise ValueError(
            f"need >=2 cells per group (target={n_t}, reference={n_r}) for "
            f"groups {target_group!r}/{reference_group!r} in {groupby!r}"
        )

    # resolve velocity layers (spliced/unspliced or kb_python mature/nascent)
    if spliced_layer is None or unspliced_layer is None:
        resolved = _resolve_velocity_layer_keys(adata)
        if resolved is None:
            raise KeyError(
                "no velocity layers found; nascent_activity_score needs "
                "spliced/unspliced (or mature/nascent) layers "
                f"(available: {list(adata.layers)})"
            )
        spliced_layer = spliced_layer or resolved[0]
        unspliced_layer = unspliced_layer or resolved[1]
    U = _layer(adata, unspliced_layer)
    S = _layer(adata, spliced_layer)

    Ut, Ur = _col_sum(U, mt), _col_sum(U, mr)
    St, Sr = _col_sum(S, mt), _col_sum(S, mr)
    tot_t = Ut.sum() + St.sum()
    tot_r = Ur.sum() + Sr.sum()

    # variance-stabilized nascent increase as a two-sample Poisson / conditional-
    # binomial score: given a gene's total unspliced n = U_t + U_r, is the target
    # share higher than its depth-exposure fraction f (null = no nascent change)?
    # This is symmetric and properly bounded — a gene absent in the reference
    # (U_r = 0) does NOT get an unbounded z (the earlier U_t - U_r*d form did),
    # so sparse zero-reference genes cannot be spuriously top-ranked at low capture.
    denom = tot_t + tot_r
    f = (tot_t / denom) if denom > 0 else 0.5
    n_us = Ut + Ur
    E0 = n_us * f
    poisson_z = (Ut - E0) / np.sqrt(n_us * f * (1.0 - f) + 1.0)

    dlog_u = _dlog(Ut, Ur)
    dlog_s = _dlog(St, Sr)

    # DE reproducibility from the SPLICED signal only (proxy-independent). A DE false
    # positive does not keep its spliced fold-change sign across independent folds.
    rng = np.random.default_rng(random_state)
    samp = adata.obs[sample_col].to_numpy() if sample_col is not None else None
    ti = np.where(mt)[0]
    ri = np.where(mr)[0]
    ft = _fold_labels(n_t, samp[mt] if samp is not None else None, n_splits, rng)
    fr = _fold_labels(n_r, samp[mr] if samp is not None else None, n_splits, rng)
    overall_sign = np.sign(dlog_s)
    agree = np.zeros(adata.n_vars, dtype=float)
    used = 0
    for k in range(n_splits):
        tmask = np.zeros(adata.n_obs, dtype=bool)
        rmask = np.zeros(adata.n_obs, dtype=bool)
        tmask[ti[ft == k]] = True
        rmask[ri[fr == k]] = True
        if not tmask.any() or not rmask.any():
            continue
        dls_k = _dlog(_col_sum(S, tmask), _col_sum(S, rmask))
        agree += (np.sign(dls_k) == overall_sign).astype(float)
        used += 1
    repro_frac = agree / used if used else np.full(adata.n_vars, np.nan)
    # a flat gene (overall spliced dlog == 0) has no reproducible direction — do not
    # let sign==0 matches count as "reproducible" (the flag is a directional DE check).
    reproducible = (
        (repro_frac >= 1.0) & (overall_sign != 0) if used else np.zeros(adata.n_vars, dtype=bool)
    )

    out = pd.DataFrame(
        {
            NASCENT_Z_COL: poisson_z,
            DLOG_U_COL: dlog_u,
            DLOG_S_COL: dlog_s,
            REPRO_COL: reproducible,
            REPRO_FRAC_COL: repro_frac,
        },
        index=adata.var_names,
    )
    logger.info(
        "nascent_activity_score: %s vs %s, %d/%d genes DE-reproducible (%d folds)",
        target_group,
        reference_group,
        int(np.nansum(reproducible)),
        adata.n_vars,
        used,
    )
    return out

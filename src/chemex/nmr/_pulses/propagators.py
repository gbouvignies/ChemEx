from __future__ import annotations

from collections.abc import Hashable, Iterable

import numpy as np
from cachetools import cached
from cachetools.keys import hashkey
from scipy.linalg import expm, schur, solve_sylvester
from scipy.linalg.lapack import dtrsen

from chemex.nmr._engine.engine import ISLiouvillianEngine
from chemex.typing import Array

DictArray = dict[str, Array]

# Imaginary-part cutoff separating secular/overdamped modes from oscillatory modes
SMALL_VALUE = 1e-6

# Heuristic 1-norm eigenvector condition-number threshold above which the
# complete-dephasing branch switches from batched eigendecomposition to a
# real-Schur modal projection. Well-conditioned NMR spin-lock Liouvillians
# typically have cond_1(V) < 20; this is a practical numerical heuristic rather
# than a strict accuracy guarantee across arbitrary matrix scales.
_EIGENVECTOR_COND_1_THRESHOLD = 1.0e3


def _cache_key(
    engine: ISLiouvillianEngine, *args: Hashable, **kwargs: Hashable
) -> tuple[Hashable, ...]:
    return tuple(hashkey(engine.basis, *args, **kwargs))


def _parse_delays(delays: float | Iterable[float]) -> tuple[Array, bool]:
    if isinstance(delays, np.ndarray):
        is_scalar = delays.ndim == 0
        delays_array = np.atleast_1d(np.asarray(delays, dtype=np.float64))
    elif isinstance(delays, Iterable) and not isinstance(delays, (str, bytes)):
        is_scalar = False
        delays_array = np.fromiter(delays, dtype=np.float64)
    else:
        is_scalar = True
        delays_array = np.array([delays], dtype=np.float64)

    if delays_array.ndim != 1:
        msg = f"`delays` must be a scalar or a 1D iterable of floats; got shape {delays_array.shape}."
        raise ValueError(msg)

    return delays_array, is_scalar


def _expm_scaled(scaled: Array) -> Array:
    """Evaluate matrix exponentials for a 2D or batched `(..., n, n)` array."""
    if scaled.size == 0:
        return np.empty(scaled.shape, dtype=np.float64)
    if not np.any(scaled):
        size = scaled.shape[-1]
        return np.broadcast_to(np.eye(size, dtype=np.float64), scaled.shape).copy()
    return np.ascontiguousarray(expm(scaled), dtype=np.float64)


def _select_secular_schur_blocks(t_mat: Array) -> Array:
    """Mark diagonal blocks of a real Schur matrix with `|imag| < SMALL_VALUE`."""
    size = t_mat.shape[0]
    frob_t = float(np.linalg.norm(t_mat, ord="fro"))
    real_roundoff_tol = 64.0 * size * np.finfo(np.float64).eps * frob_t
    asymmetry_bound = 1.0 / (_EIGENVECTOR_COND_1_THRESHOLD**2)

    select = np.zeros(size, dtype=np.int32)
    idx = 0
    while idx < size:
        if idx + 1 < size and t_mat[idx + 1, idx] != 0.0:
            # Canonical 2x2 real Schur block [[a, b], [c, a]] with b * c < 0.
            # A real defective block perturbed by O(eps * ||T||_F) roundoff can
            # split by O(sqrt(hi * eps * ||T||_F)); treat the imaginary part as
            # roundoff noise only when the distance `lo` to real eigenvalues is
            # within backward roundoff AND the block is strongly non-normal.
            b = abs(float(t_mat[idx, idx + 1]))
            c = abs(float(t_mat[idx + 1, idx]))
            lo, hi = (b, c) if b <= c else (c, b)
            is_roundoff_split = lo <= real_roundoff_tol and lo <= asymmetry_bound * hi
            imag_mag = 0.0 if is_roundoff_split else float(np.sqrt(b * c))
            if imag_mag < SMALL_VALUE:
                select[idx] = 1
                select[idx + 1] = 1
            idx += 2
        else:
            select[idx] = 1
            idx += 1
    return select


def _dephased_propagators_schur(liouv_array: Array, delays_array: Array) -> Array:
    """Compute modal complete-dephasing propagators via real Schur + Sylvester."""
    batch_shape = liouv_array.shape[:-2]
    size = liouv_array.shape[-1]
    out = np.empty((delays_array.size, *batch_shape, size, size), dtype=np.float64)
    flat_liouv = liouv_array.reshape((-1, size, size))
    flat_out = out.reshape((delays_array.size, -1, size, size))
    eye = np.eye(size, dtype=np.float64)

    zero_mask = delays_array == 0.0
    neg_mask = delays_array < 0.0
    pos_mask = delays_array > 0.0

    for b_idx in range(flat_liouv.shape[0]):
        mat = flat_liouv[b_idx]
        if np.any(zero_mask):
            flat_out[zero_mask, b_idx] = eye
        if np.any(neg_mask):
            flat_out[neg_mask, b_idx] = _expm_scaled(
                delays_array[neg_mask, None, None] * mat[None, :, :]
            )
        if not np.any(pos_mask):
            continue

        t_mat, z_mat = schur(mat, output="real")
        select = _select_secular_schur_blocks(t_mat)
        n_keep = int(np.sum(select))
        pos_delays = delays_array[pos_mask]

        if n_keep == 0:
            flat_out[pos_mask, b_idx] = 0.0
        elif n_keep == size:
            flat_out[pos_mask, b_idx] = _expm_scaled(
                pos_delays[:, None, None] * mat[None, :, :]
            )
        else:
            t_ord, z_ord, _wr, _wi, m_out, _s, _sep, info = dtrsen(
                select, t_mat, z_mat, job="N", wantq=1
            )
            ordered_count = int(m_out)
            if info < 0:
                msg = f"Invalid argument {-int(info)} passed to LAPACK dtrsen."
                raise np.linalg.LinAlgError(msg)
            if info > 0 or ordered_count != n_keep:
                msg = (
                    "Complete-dephasing Schur reordering failed to separate "
                    f"near-degenerate eigenvalues (dtrsen info={int(info)}, "
                    f"selected={n_keep}, ordered={ordered_count})."
                )
                raise np.linalg.LinAlgError(msg)

            t11 = t_ord[:n_keep, :n_keep]
            t12 = t_ord[:n_keep, n_keep:]
            t22 = t_ord[n_keep:, n_keep:]
            x_mat = solve_sylvester(t11, -t22, -t12)
            if not np.isfinite(x_mat).all():
                msg = "Complete-dephasing Sylvester projection produced non-finite values."
                raise np.linalg.LinAlgError(msg)

            z1 = z_ord[:, :n_keep]
            z2 = z_ord[:, n_keep:]
            right = z1.T - x_mat @ z2.T
            exp_t11 = _expm_scaled(pos_delays[:, None, None] * t11[None, :, :])
            flat_out[pos_mask, b_idx] = (z1 @ exp_t11) @ right

    return np.ascontiguousarray(out)


def _dephased_propagators(liouv_array: Array, delays_array: Array) -> Array:
    try:
        eigenvalues, eigenvectors = np.linalg.eig(liouv_array)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            eigenvectors_inv = np.linalg.inv(eigenvectors)
            norm_v = np.max(np.sum(np.abs(eigenvectors), axis=-2), axis=-1)
            norm_vinv = np.max(np.sum(np.abs(eigenvectors_inv), axis=-2), axis=-1)
            cond_1 = norm_v * norm_vinv
        if np.isfinite(eigenvectors_inv).all() and np.all(
            cond_1 <= _EIGENVECTOR_COND_1_THRESHOLD
        ):
            # The complete-dephasing limit removes oscillatory modes for positive
            # durations. At zero duration there is no evolution, so retain identity.
            exponents = np.multiply.outer(delays_array, eigenvalues)
            removed = (delays_array > 0).reshape((-1,) + (1,) * eigenvalues.ndim) & (
                np.abs(eigenvalues.imag) >= SMALL_VALUE
            )
            exp_eigenvalues = np.exp(
                exponents, out=np.zeros_like(exponents), where=~removed
            )
            propagators = (
                eigenvectors[np.newaxis, ...] * exp_eigenvalues[..., np.newaxis, :]
            ) @ eigenvectors_inv[np.newaxis, ...]
            return np.ascontiguousarray(propagators.real)
    except np.linalg.LinAlgError:
        pass

    return _dephased_propagators_schur(liouv_array, delays_array)


def calculate_propagators(
    liouv: Array,
    delays: float | Iterable[float],
    *,
    dephasing: bool = False,
) -> Array:
    """Calculate propagators for one or more delays."""
    liouv_array = np.asarray(liouv, dtype=np.float64)
    if liouv_array.ndim < 2 or liouv_array.shape[-1] != liouv_array.shape[-2]:
        msg = (
            "`liouv` must have square matrix dimensions `(..., n, n)`; "
            f"got shape {liouv_array.shape}."
        )
        raise ValueError(msg)

    delays_array, is_scalar = _parse_delays(delays)

    if not dephasing:
        if is_scalar:
            return _expm_scaled(liouv_array * delays_array[0])
        return _expm_scaled(np.multiply.outer(delays_array, liouv_array))

    if delays_array.size == 0:
        return np.empty((0, *liouv_array.shape), dtype=np.float64)

    propagators = _dephased_propagators(liouv_array, delays_array)
    return propagators[0] if is_scalar else propagators


@cached(cache={}, key=_cache_key)
def make_perfect180(engine: ISLiouvillianEngine, spin: str) -> Array:
    size = engine.size
    identity = np.eye(size).reshape((1, 1, size, size))
    compx, compy, compz = (f"{spin}{axis}" for axis in "xyz")
    perfect180: dict[str, Array] = {comp: identity.copy() for comp in (compx, compy)}
    for comp in engine.basis.components:
        vect = engine.basis.vectors[comp].ravel()
        if compx in comp or compz in comp:
            perfect180[compy] -= 2 * np.diag(vect)
        if compy in comp or compz in comp:
            perfect180[compx] -= 2 * np.diag(vect)
    p180 = [perfect180[comp] for comp in (compx, compy, compx, compy)]
    return np.array(p180)


@cached(cache={}, key=_cache_key)
def make_perfect90(engine: ISLiouvillianEngine, spin: str) -> Array:
    size = engine.size
    zeros = np.zeros((size, size))
    rot = engine.basis.matrices.get(f"b1x_{spin}", zeros)
    return expm(0.25 * rot).reshape(1, 1, size, size)


@cached(cache={}, key=_cache_key)
def get_phases(engine: ISLiouvillianEngine) -> DictArray:
    phases = {}
    size = engine.size
    zeros = np.zeros((size, size))
    for spin in "is":
        l_rotz = engine.basis.matrices.get(f"rotz_{spin}", zeros)
        phases[spin] = np.array([expm(n * 0.5 * np.pi * l_rotz) for n in range(4)])
    return phases

from __future__ import annotations

from functools import cache
from itertools import product
from typing import Literal

import numpy as np
from scipy.linalg import schur, solve_sylvester
from scipy.linalg.lapack import ztrsen

from chemex.nmr._engine.engine import ISLiouvillianEngine
from chemex.nmr._engine.views import reshape_single_liouvillian
from chemex.nmr.basis import Basis
from chemex.typing import Array

# A small value used for numerical stability
SMALL_VALUE = 1e-6

# Inhomogeneous thermal-bath identity coordinates in `*_eq` bases; these form an
# invariant source block with conserved 0 and pure-exchange eigenvalues rather
# than spin-relaxation modes.
_EQUILIBRIUM_COMPONENTS = frozenset(("ie", "se"))

Coherence = Literal["sq_i", "sq_s", "dq_is", "zq_is"]

_INV_SQRT2 = float(np.sqrt(0.5))
# Right-eigenvector basis vectors in ChemEx's Cartesian operator convention,
# chosen so a free-precession frequency Omega = delta * ppm produces eigenvalue
# -R2 + i * Omega (Im(lambda) / ppm == +delta for both gamma > 0 and gamma < 0):
#   - SQ(I):  (I_x - i I_y) / sqrt(2)          -> -R2_I + i Omega_I
#   - SQ(S):  (S_x - i S_y) / sqrt(2)          -> -R2_S + i Omega_S
#   - DQ(IS): (2I_xS_x - i 2I_xS_y - i 2I_yS_x - 2I_yS_y) / 2
#                                              -> -(R2_MQ + mu_IS) + i (Omega_I + Omega_S)
#   - ZQ(IS): (2I_xS_x + i 2I_xS_y - i 2I_yS_x + 2I_yS_y) / 2
#                                              -> -(R2_MQ - mu_IS) + i (Omega_I - Omega_S)
_COHERENCE_OPERATORS: dict[Coherence, tuple[tuple[str, complex], ...]] = {
    "sq_i": (("ix", _INV_SQRT2), ("iy", -1j * _INV_SQRT2)),
    "sq_s": (("sx", _INV_SQRT2), ("sy", -1j * _INV_SQRT2)),
    "dq_is": (
        ("2ixsx", 0.5),
        ("2ixsy", -0.5j),
        ("2iysx", -0.5j),
        ("2iysy", -0.5),
    ),
    "zq_is": (
        ("2ixsx", 0.5),
        ("2ixsy", 0.5j),
        ("2iysx", -0.5j),
        ("2iysy", 0.5),
    ),
}

_EPS = float(np.finfo(np.float64).eps)
_EIG_COND_LIMIT = 1.0e6
_SYLVESTER_NORM_LIMIT = 1.0e6


@cache
def _build_coherence_basis(basis: Basis, coherence: Coherence) -> Array:
    if coherence not in _COHERENCE_OPERATORS:
        msg = f"Unsupported coherence {coherence!r}; expected one of {tuple(_COHERENCE_OPERATORS)}."
        raise ValueError(msg)

    terms = _COHERENCE_OPERATORS[coherence]
    components = basis.components
    missing = [comp for comp, _ in terms if comp not in components]
    if missing:
        msg = (
            f"Basis {basis.type!r} lacks components {missing!r} required for "
            f"coherence {coherence!r}."
        )
        raise ValueError(msg)

    n_states = len(basis.model.states)
    n_comp = len(components)
    u_mat = np.zeros((n_states * n_comp, n_states), dtype=np.complex128)
    comp_offsets = [(components.index(comp), coeff) for comp, coeff in terms]
    for state_idx in range(n_states):
        base = state_idx * n_comp
        for offset, coeff in comp_offsets:
            u_mat[base + offset, state_idx] = coeff
    u_mat.flags.writeable = False
    return u_mat


def _project_coherence_liouvillian(
    liouv: Array,
    basis: Basis,
    coherence: Coherence,
) -> Array:
    u_mat = _build_coherence_basis(basis, coherence)
    liouv_u = liouv @ u_mat
    liouv_sub = u_mat.conj().T @ liouv_u
    leakage = float(np.linalg.norm(liouv_u - u_mat @ liouv_sub, ord="fro"))
    scale = max(1.0, float(np.linalg.norm(liouv, ord="fro")))
    if leakage > 64.0 * liouv.shape[0] * _EPS * scale:
        msg = (
            f"Coherence {coherence!r} is coupled outside its scalar subspace in "
            f"basis {basis.type!r} (leakage={leakage:.3e})."
        )
        raise ValueError(msg)
    return liouv_sub


def _pick_dominant_mode(weights: Array, delta_eigs: Array) -> int:
    abs_weights = np.abs(weights)
    max_abs = float(np.max(abs_weights))
    if max_abs == 0.0:
        return int(np.argmin(np.abs(delta_eigs)))
    candidates = np.flatnonzero(abs_weights >= (1.0 - 16.0 * _EPS) * max_abs)
    if candidates.size == 1:
        return int(candidates[0])
    best_local = int(np.argmin(np.abs(delta_eigs[candidates])))
    return int(candidates[best_local])


def _schur_cluster_projector(
    t_mat: Array,
    z_mat: Array,
    indices: tuple[int, ...],
) -> tuple[Array | None, float, complex]:
    size = t_mat.shape[0]
    count = len(indices)
    if count == size:
        mean_eig = complex(np.trace(t_mat) / size)
        return np.eye(size, dtype=np.complex128), 0.0, mean_eig

    select = np.zeros(size, dtype=np.int32)
    select[list(indices)] = 1
    t_ord, z_ord, _w, m_out, _s, _sep, info = ztrsen(
        select,
        t_mat,
        z_mat,
        job="N",
        wantq=1,
    )
    if info != 0 or int(m_out) != count:
        return None, np.inf, 0j

    t11 = t_ord[:count, :count]
    t12 = t_ord[:count, count:]
    t22 = t_ord[count:, count:]
    mean_eig = complex(np.trace(t11) / count)

    try:
        y_mat = solve_sylvester(t11, -t22, -t12)
    except np.linalg.LinAlgError:
        return None, np.inf, mean_eig

    if not np.isfinite(y_mat).all():
        return None, np.inf, mean_eig

    y_norm = float(np.linalg.norm(y_mat, ord="fro"))
    z1 = z_ord[:, :count]
    z2 = z_ord[:, count:]
    proj = z1 @ (z1.conj().T - y_mat @ z2.conj().T)
    return proj, y_norm, mean_eig


def _try_eig_frequency_offset(
    centered: Array,
    obs_index: int,
    sep_tol: float,
) -> float | None:
    try:
        eigenvalues, eigenvectors = np.linalg.eig(centered)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            inv_v = np.linalg.inv(eigenvectors)
            norm_v = float(np.max(np.sum(np.abs(eigenvectors), axis=0)))
            norm_vinv = float(np.max(np.sum(np.abs(inv_v), axis=0)))
            cond_1 = norm_v * norm_vinv
    except np.linalg.LinAlgError:
        return None

    diffs = np.abs(eigenvalues[:, None] - eigenvalues[None, :])
    np.fill_diagonal(diffs, np.inf)
    if not (
        np.isfinite(inv_v).all()
        and cond_1 <= _EIG_COND_LIMIT
        and float(np.min(diffs)) > sep_tol
    ):
        return None

    weights = eigenvectors.sum(axis=0) * inv_v[:, obs_index]
    best = _pick_dominant_mode(weights, eigenvalues)
    return float(eigenvalues[best].imag)


def _cluster_distance(
    diag: Array,
    c1: tuple[int, ...],
    c2: tuple[int, ...],
) -> float:
    return float(min(abs(diag[i] - diag[j]) for i in c1 for j in c2))


def _initial_schur_clusters(diag: Array, sep_tol: float) -> list[tuple[int, ...]]:
    clusters: list[tuple[int, ...]] = [(idx,) for idx in range(diag.size)]
    while len(clusters) > 1:
        pair = next(
            (
                (i, j)
                for i in range(len(clusters))
                for j in range(i + 1, len(clusters))
                if _cluster_distance(diag, clusters[i], clusters[j]) <= sep_tol
            ),
            None,
        )
        if pair is None:
            break
        lo, hi = pair
        clusters[lo] = tuple(sorted(clusters[lo] + clusters[hi]))
        del clusters[hi]
    return clusters


def _schur_frequency_offset(
    centered: Array,
    obs_index: int,
    sep_tol: float,
) -> float:
    size = centered.shape[0]
    t_mat, z_mat = schur(centered, output="complex")
    diag = np.diag(t_mat)
    clusters = _initial_schur_clusters(diag, sep_tol)

    while len(clusters) > 1:
        bad_idx: int | None = None
        results: list[tuple[Array, complex]] = []
        for idx, cluster in enumerate(clusters):
            proj, y_norm, mean_eig = _schur_cluster_projector(t_mat, z_mat, cluster)
            if proj is None or y_norm > _SYLVESTER_NORM_LIMIT:
                bad_idx = idx
                break
            results.append((proj, mean_eig))

        if bad_idx is None:
            weights = np.array(
                [proj[:, obs_index].sum() for proj, _ in results],
                dtype=np.complex128,
            )
            mean_eigs = np.array(
                [mean_eig for _, mean_eig in results],
                dtype=np.complex128,
            )
            best = _pick_dominant_mode(weights, mean_eigs)
            return float(mean_eigs[best].imag)

        nearest = min(
            (j for j in range(len(clusters)) if j != bad_idx),
            key=lambda j: _cluster_distance(diag, clusters[bad_idx], clusters[j]),
        )
        lo, hi = min(bad_idx, nearest), max(bad_idx, nearest)
        clusters[lo] = tuple(sorted(clusters[lo] + clusters[hi]))
        del clusters[hi]

    return float(complex(np.trace(t_mat) / size).imag)


def _select_observed_frequency_offset(liouv_sub: Array, obs_index: int) -> float:
    """Return Im(lambda_obs - liouv_sub[obs_index, obs_index]) for observed state."""
    size = liouv_sub.shape[0]
    if size == 1:
        return 0.0

    ref_diag = complex(liouv_sub[obs_index, obs_index])
    centered = liouv_sub - ref_diag * np.eye(size, dtype=np.complex128)
    scale = float(np.linalg.norm(centered, ord="fro"))
    if scale == 0.0:
        return 0.0

    sep_tol = 64.0 * size * _EPS * scale
    eig_offset = _try_eig_frequency_offset(centered, obs_index, sep_tol)
    if eig_offset is not None:
        return eig_offset
    return _schur_frequency_offset(centered, obs_index, sep_tol)


class SpectrometerAnalysis:
    """Small analysis surface built on top of a spectrometer engine."""

    def __init__(self, engine: ISLiouvillianEngine) -> None:
        self._engine = engine

    def _homogeneous_spin_submatrix(self, liouv: Array) -> Array:
        basis = self._engine.basis
        if not any(comp in _EQUILIBRIUM_COMPONENTS for comp in basis.components):
            return liouv
        spin_indices = [
            idx
            for idx, (_state, comp) in enumerate(
                product(basis.model.states, basis.components)
            )
            if comp not in _EQUILIBRIUM_COMPONENTS
        ]
        return liouv[np.ix_(spin_indices, spin_indices)]

    def _coherence_liouvillian(self, coherence: Coherence) -> Array:
        liouv = reshape_single_liouvillian(
            self._engine.l_free,
            self._engine.size,
            purpose="Shift eigenvalue calculation",
        )
        return _project_coherence_liouvillian(liouv, self._engine.basis, coherence)

    def calculate_shifts(self, coherence: Coherence = "sq_i") -> Array:
        liouv_sub = self._coherence_liouvillian(coherence)
        return np.linalg.eigvals(liouv_sub).imag

    def calculate_shift(
        self,
        observed_state: str = "a",
        *,
        coherence: Coherence = "sq_i",
        relative: bool = False,
    ) -> float:
        states = self._engine.model.states
        if observed_state not in states:
            msg = (
                f"Unknown observed_state {observed_state!r} for model "
                f"states {states!r}."
            )
            raise ValueError(msg)
        obs_index = states.index(observed_state)
        liouv_sub = self._coherence_liouvillian(coherence)
        delta_omega = _select_observed_frequency_offset(liouv_sub, obs_index)
        if relative:
            return delta_omega
        return float(liouv_sub[obs_index, obs_index].imag) + delta_omega

    def calculate_r1rho(self) -> float:
        liouv = reshape_single_liouvillian(
            self._engine.l_free + self._engine.l_b1x_i,
            self._engine.size,
            purpose="R1rho eigenvalue calculation",
        )
        liouv = self._homogeneous_spin_submatrix(liouv)
        eigenvalues = np.linalg.eigvals(liouv)
        real_eigenvalues = eigenvalues[
            np.isclose(eigenvalues.imag, 0.0, atol=SMALL_VALUE)
        ].real
        if real_eigenvalues.size == 0:
            smallest_imag = float(np.min(np.abs(eigenvalues.imag)))
            msg = (
                "R1rho eigenvalue calculation did not find a nearly real eigenvalue "
                f"within atol={SMALL_VALUE}; smallest |imag| was {smallest_imag:.3e}."
            )
            raise ValueError(msg)
        return -float(np.max(real_eigenvalues))

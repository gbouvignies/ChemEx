from __future__ import annotations

from collections.abc import Hashable, Iterable

import numpy as np
from cachetools import cached
from cachetools.keys import hashkey
from scipy.linalg import expm

from chemex.nmr.is_liouvillian_engine import ISLiouvillianEngine
from chemex.typing import Array

DictArray = dict[str, Array]

# A small value used for numerical stability
SMALL_VALUE = 1e-6


def _cache_key(
    engine: ISLiouvillianEngine, *args: Hashable, **kwargs: Hashable
) -> tuple[Hashable, ...]:
    return tuple(hashkey(engine.basis, *args, **kwargs))


def _dephase_eigenvalues(eigenvalues: Array) -> Array:
    return np.where(
        np.abs(eigenvalues.imag) < SMALL_VALUE,
        eigenvalues,
        eigenvalues * 1e9,
    )


def _as_diagonal_matrices(values: Array) -> Array:
    diagonal_shape = (*values.shape, values.shape[-1])
    diagonal_matrices = np.zeros(diagonal_shape, dtype=np.complex128)
    indices = np.arange(values.shape[-1])
    diagonal_matrices[..., indices, indices] = values
    return diagonal_matrices


def calculate_propagators(
    liouv: Array,
    delays: float | Iterable[float],
    *,
    dephasing: bool = False,
) -> Array:
    """Calculate propagators for one or more delays."""
    liouv_array = np.asarray(liouv, dtype=np.float64)
    delays_array = np.atleast_1d(np.asarray(delays, dtype=np.float64))

    if liouv_array.ndim == 2 and delays_array.size == 1 and not dephasing:
        return expm(liouv_array * delays_array[0])

    eigenvalues, eigenvectors = np.linalg.eig(liouv_array)
    if dephasing:
        eigenvalues = _dephase_eigenvalues(eigenvalues)

    exp_eigenvalues = np.exp(np.multiply.outer(delays_array, eigenvalues))
    diagonal_matrices = _as_diagonal_matrices(exp_eigenvalues)

    vectors = np.expand_dims(eigenvectors, axis=0)
    transformed = vectors @ diagonal_matrices
    propagators_t = np.linalg.solve(
        np.swapaxes(vectors, -1, -2),
        np.swapaxes(transformed, -1, -2),
    )
    propagators = np.swapaxes(propagators_t, -1, -2)

    if propagators.shape[0] == 1:
        propagators = propagators[0]

    return propagators.real


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

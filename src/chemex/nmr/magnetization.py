from __future__ import annotations

from collections.abc import Iterable, Mapping
from itertools import product

import numpy as np

from chemex.nmr.basis import Basis
from chemex.nmr.constants import XI_RATIO
from chemex.parameters.spin_system.nucleus import Nucleus
from chemex.typing import Array


def _zero_magnetization(basis: Basis) -> Array:
    size = len(basis) * len(basis.model.states)
    return np.zeros((size, 1))


def collapse_magnetization(magnetization: Array, weights: Array) -> Array:
    """Collapse a distributed magnetization into a weighted average."""
    if (ndim := magnetization.ndim) < 3:
        return magnetization
    magnetization_weighted = weights * magnetization
    sum_axes = tuple(range(ndim - 2))
    return magnetization_weighted.sum(axis=sum_axes)


def detect_signal(
    detection_vector: Array,
    magnetization: Array,
    weights: Array,
) -> float:
    """Detect a scalar signal from a magnetization vector."""
    collapsed_magnetization = collapse_magnetization(magnetization, weights)
    detected = detection_vector @ collapsed_magnetization
    if np.iscomplexobj(detected):
        detected = np.sign(detected.real) * np.abs(detected)
    return float(detected.item())


def build_equilibrium_magnetization(
    basis: Basis,
    par_values: Mapping[str, float],
) -> Array:
    """Build the thermal equilibrium magnetization for a basis."""
    magnetization = _zero_magnetization(basis)
    for state, (name, nucleus) in product(basis.model.states, basis.nuclei.items()):
        scale = par_values.get(f"p{state}", 0.0) * XI_RATIO.get(nucleus, 1.0)
        magnetization += basis.vectors.get(f"{name}e_{state}", 0.0) * scale
        magnetization += basis.vectors.get(f"{name}z_{state}", 0.0) * scale
    return magnetization


def build_start_magnetization(
    basis: Basis,
    par_values: Mapping[str, float],
    terms: Iterable[str],
    atom: Nucleus = Nucleus.H1,
) -> Array:
    """Build a starting magnetization from selected basis terms."""
    ratio = XI_RATIO.get(atom, 1.0)
    terms_set = set(terms)
    magnetization = _zero_magnetization(basis)

    for component, vector in basis.vectors.items():
        if "_" not in component:
            continue
        state = component[-1]
        for term in terms_set:
            if component.startswith(term.strip("+-")):
                sign = -1.0 if term.startswith("-") else 1.0
                population = par_values.get(f"p{state}", 0.0)
                magnetization += sign * population * ratio * vector

    return magnetization


def keep_components(
    basis: Basis,
    magnetization: Array,
    components: Iterable[str],
) -> Array:
    """Zero all magnetization components except the selected ones."""
    keep_mask = sum(
        (basis.vectors[name] for name in components),
        start=_zero_magnetization(basis),
    )
    keep_mask[keep_mask > 0] = 1.0
    return keep_mask * magnetization

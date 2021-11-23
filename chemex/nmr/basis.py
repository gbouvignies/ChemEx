from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from functools import cache
from itertools import permutations
from itertools import product
from typing import Literal

import numpy as np

from chemex.model import model

_BASES = {
    "ixy": ["ix", "iy"],
    "iz": ["iz"],
    "izsz": ["iz", "2izsz"],
    "iz_eq": ["ie", "iz"],
    "ixyz": ["ix", "iy", "iz"],
    "ixyz_eq": ["ie", "ix", "iy", "iz"],
    "ixysxy": ["2ixsx", "2ixsy", "2iysx", "2iysy"],
    "ixy_ixysxy": ["ix", "iy", "2ixsx", "2ixsy", "2iysx", "2iysy"],
    "ixyzsz": ["ix", "iy", "iz", "2ixsz", "2iysz", "2izsz"],
    "ixyzsz_diff": ["ix", "iy", "iz", "2ixsz", "2iysz", "2izsz"],
    "ixyzsz_eq": ["ie", "ix", "iy", "iz", "2ixsz", "2iysz", "2izsz"],
    "ixyzsxyz": [
        "ix",
        "iy",
        "iz",
        "sx",
        "sy",
        "sz",
        "2ixsz",
        "2iysz",
        "2izsx",
        "2izsy",
        "2ixsx",
        "2ixsy",
        "2iysx",
        "2iysy",
        "2izsz",
    ],
    "ixyzsxyz_eq": [
        "ie",
        "se",
        "ix",
        "iy",
        "iz",
        "sx",
        "sy",
        "sz",
        "2ixsz",
        "2iysz",
        "2izsx",
        "2izsy",
        "2ixsx",
        "2ixsy",
        "2iysx",
        "2iysy",
        "2izsz",
    ],
}
_TRANSITIONS: dict[str, tuple[tuple[str, str, float], ...]] = {
    "r2_i_{state}": (("ix", "ix", -1.0), ("iy", "iy", -1.0)),
    "r2_s_{state}": (("sx", "sx", -1.0), ("sy", "sy", -1.0)),
    "r1_i_{state}": (("iz", "iz", -1.0), ("iz", "ie", +1.0)),
    "r1_s_{state}": (("sz", "sz", -1.0), ("sz", "se", +1.0)),
    "r2a_i_{state}": (("2ixsz", "2ixsz", -1.0), ("2iysz", "2iysz", -1.0)),
    "r2a_s_{state}": (("2izsx", "2izsx", -1.0), ("2izsy", "2izsy", -1.0)),
    "r2mq_is_{state}": (
        ("2ixsx", "2ixsx", -1.0),
        ("2ixsy", "2ixsy", -1.0),
        ("2iysx", "2iysx", -1.0),
        ("2iysy", "2iysy", -1.0),
    ),
    "r1a_is_{state}": (("2izsz", "2izsz", -1.0),),
    "etaxy_i_{state}": (
        ("ix", "2ixsz", -1.0),
        ("iy", "2iysz", -1.0),
        ("2ixsz", "ix", -1.0),
        ("2iysz", "iy", -1.0),
    ),
    "etaxy_s_{state}": (
        ("sx", "2izsx", -1.0),
        ("sy", "2izsy", -1.0),
        ("2izsx", "sx", -1.0),
        ("2izsy", "sy", -1.0),
    ),
    "etaz_i_{state}": (
        ("iz", "2izsz", -1.0),
        ("2izsz", "iz", -1.0),
        ("2izsz", "ie", +1.0),
    ),
    "etaz_s_{state}": (
        ("sz", "2izsz", -1.0),
        ("2izsz", "sz", -1.0),
        ("2izsz", "se", +1.0),
    ),
    "sigma_is_{state}": (
        ("iz", "sz", -1.0),
        ("sz", "iz", -1.0),
        ("sz", "ie", +1.0),
        ("iz", "se", +1.0),
    ),
    "mu_is_{state}": (
        ("2ixsx", "2iysy", +1.0),
        ("2ixsy", "2iysx", -1.0),
        ("2iysx", "2ixsy", -1.0),
        ("2iysy", "2ixsx", +1.0),
    ),
    "rotz_i": (
        ("ix", "iy", -1.0),
        ("iy", "ix", +1.0),
        ("2ixsx", "2iysx", -1.0),
        ("2iysx", "2ixsx", +1.0),
        ("2ixsy", "2iysy", -1.0),
        ("2iysy", "2ixsy", +1.0),
        ("2ixsz", "2iysz", -1.0),
        ("2iysz", "2ixsz", +1.0),
    ),
    "rotz_s": (
        ("sx", "sy", -1.0),
        ("sy", "sx", +1.0),
        ("2ixsx", "2ixsy", -1.0),
        ("2ixsy", "2ixsx", +1.0),
        ("2iysx", "2iysy", -1.0),
        ("2iysy", "2iysx", +1.0),
        ("2izsx", "2izsy", -1.0),
        ("2izsy", "2izsx", +1.0),
    ),
    "cs_i_{state}": (
        ("ix", "iy", -1.0),
        ("iy", "ix", +1.0),
        ("2ixsx", "2iysx", -1.0),
        ("2iysx", "2ixsx", +1.0),
        ("2ixsy", "2iysy", -1.0),
        ("2iysy", "2ixsy", +1.0),
        ("2ixsz", "2iysz", -1.0),
        ("2iysz", "2ixsz", +1.0),
    ),
    "cs_s_{state}": (
        ("sx", "sy", -1.0),
        ("sy", "sx", +1.0),
        ("2ixsx", "2ixsy", -1.0),
        ("2ixsy", "2ixsx", +1.0),
        ("2iysx", "2iysy", -1.0),
        ("2iysy", "2iysx", +1.0),
        ("2izsx", "2izsy", -1.0),
        ("2izsy", "2izsx", +1.0),
    ),
    "carrier_i": (
        ("ix", "iy", +1.0),
        ("iy", "ix", -1.0),
        ("2ixsx", "2iysx", +1.0),
        ("2iysx", "2ixsx", -1.0),
        ("2ixsy", "2iysy", +1.0),
        ("2iysy", "2ixsy", -1.0),
        ("2ixsz", "2iysz", +1.0),
        ("2iysz", "2ixsz", -1.0),
    ),
    "carrier_s": (
        ("sx", "sy", +1.0),
        ("sy", "sx", -1.0),
        ("2ixsx", "2ixsy", +1.0),
        ("2ixsy", "2ixsx", -1.0),
        ("2iysx", "2iysy", +1.0),
        ("2iysy", "2iysx", -1.0),
        ("2izsx", "2izsy", +1.0),
        ("2izsy", "2izsx", -1.0),
    ),
    "offset_i": (
        ("ix", "iy", +2.0 * np.pi),
        ("iy", "ix", -2.0 * np.pi),
        ("2ixsx", "2iysx", +2.0 * np.pi),
        ("2iysx", "2ixsx", -2.0 * np.pi),
        ("2ixsy", "2iysy", +2.0 * np.pi),
        ("2iysy", "2ixsy", -2.0 * np.pi),
        ("2ixsz", "2iysz", +2.0 * np.pi),
        ("2iysz", "2ixsz", -2.0 * np.pi),
    ),
    "offset_s": (
        ("sx", "sy", +2.0 * np.pi),
        ("sy", "sx", -2.0 * np.pi),
        ("2ixsx", "2ixsy", +2.0 * np.pi),
        ("2ixsy", "2ixsx", -2.0 * np.pi),
        ("2iysx", "2iysy", +2.0 * np.pi),
        ("2iysy", "2iysx", -2.0 * np.pi),
        ("2izsx", "2izsy", +2.0 * np.pi),
        ("2izsy", "2izsx", -2.0 * np.pi),
    ),
    "jeff_i": (
        ("ix", "iy", -2.0 * np.pi),
        ("iy", "ix", +2.0 * np.pi),
        ("2ixsx", "2iysx", -2.0 * np.pi),
        ("2iysx", "2ixsx", +2.0 * np.pi),
        ("2ixsy", "2iysy", -2.0 * np.pi),
        ("2iysy", "2ixsy", +2.0 * np.pi),
        ("2ixsz", "2iysz", -2.0 * np.pi),
        ("2iysz", "2ixsz", +2.0 * np.pi),
    ),
    "j_is_{state}": (
        ("ix", "2iysz", -np.pi),
        ("2iysz", "ix", +np.pi),
        ("2ixsz", "iy", -np.pi),
        ("iy", "2ixsz", +np.pi),
        ("sx", "2izsy", -np.pi),
        ("2izsy", "sx", +np.pi),
        ("2izsx", "sy", -np.pi),
        ("sy", "2izsx", +np.pi),
    ),
    "d_{state}": (
        ("ix", "ix", -1.0),
        ("iy", "iy", -1.0),
        ("iz", "iz", -1.0),
        ("sx", "sx", -1.0),
        ("sy", "sy", -1.0),
        ("sz", "sz", -1.0),
        ("2ixsz", "2ixsz", -1.0),
        ("2iysz", "2iysz", -1.0),
        ("2izsx", "2izsx", -1.0),
        ("2izsy", "2izsy", -1.0),
        ("2ixsx", "2ixsx", -1.0),
        ("2ixsy", "2ixsy", -1.0),
        ("2iysx", "2iysx", -1.0),
        ("2iysy", "2iysy", -1.0),
        ("2izsz", "2izsz", -1.0),
    ),
    "b1x_i": (
        ("iy", "iz", -2.0 * np.pi),
        ("iz", "iy", +2.0 * np.pi),
        ("2iysx", "2izsx", -2.0 * np.pi),
        ("2izsx", "2iysx", +2.0 * np.pi),
        ("2iysy", "2izsy", -2.0 * np.pi),
        ("2izsy", "2iysy", +2.0 * np.pi),
        ("2iysz", "2izsz", -2.0 * np.pi),
        ("2izsz", "2iysz", +2.0 * np.pi),
    ),
    "b1y_i": (
        ("iz", "ix", -2.0 * np.pi),
        ("ix", "iz", +2.0 * np.pi),
        ("2izsx", "2ixsx", -2.0 * np.pi),
        ("2ixsx", "2izsx", +2.0 * np.pi),
        ("2izsy", "2ixsy", -2.0 * np.pi),
        ("2ixsy", "2izsy", +2.0 * np.pi),
        ("2izsz", "2ixsz", -2.0 * np.pi),
        ("2ixsz", "2izsz", +2.0 * np.pi),
    ),
    "b1x_s": (
        ("sy", "sz", -2.0 * np.pi),
        ("sz", "sy", +2.0 * np.pi),
        ("2ixsy", "2ixsz", -2.0 * np.pi),
        ("2ixsz", "2ixsy", +2.0 * np.pi),
        ("2iysy", "2iysz", -2.0 * np.pi),
        ("2iysz", "2iysy", +2.0 * np.pi),
        ("2izsy", "2izsz", -2.0 * np.pi),
        ("2izsz", "2izsy", +2.0 * np.pi),
    ),
    "b1y_s": (
        ("sz", "sx", -2.0 * np.pi),
        ("sx", "sz", +2.0 * np.pi),
        ("2ixsz", "2ixsx", -2.0 * np.pi),
        ("2ixsx", "2ixsz", +2.0 * np.pi),
        ("2iysz", "2iysx", -2.0 * np.pi),
        ("2iysx", "2iysz", +2.0 * np.pi),
        ("2izsz", "2izsx", -2.0 * np.pi),
        ("2izsx", "2izsz", +2.0 * np.pi),
    ),
}
_ATOMS = {
    atoms: {"i": atoms[0], "s": atoms[1]} for atoms in ("hn", "hc", "nh", "ch", "cn")
}


@cache
def _build_vectors(basis: Basis) -> dict[str, np.ndarray]:
    size = len(basis) * len(model.states)
    vectors: defaultdict[str, np.ndarray] = defaultdict(lambda: np.zeros((size, 1)))
    for index, (state, name) in enumerate(product(model.states, basis.components)):
        vectors[f"{name}_{state}"][index] = 1.0
        vectors[name][index] = 1.0
    return dict(vectors)


def _get_indices(
    basis: Basis, transition_name: str, state: str
) -> tuple[tuple[list[int], list[int]], list[float]]:
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    offset = model.states.index(state) * len(basis)
    for start, end, value in _TRANSITIONS[transition_name]:
        if {start, end}.issubset(basis.components):
            rows.append(basis.components.index(start) + offset)
            cols.append(basis.components.index(end) + offset)
            vals.append(value)
    return (rows, cols), vals


def _build_spin_matrices(basis: Basis) -> dict[str, np.ndarray]:
    size = len(basis) * len(model.states)
    matrices: dict[str, np.ndarray] = defaultdict(lambda: np.zeros((size, size)))
    for transition_name, state in product(_TRANSITIONS, model.states):
        if not basis.type.endswith("_diff") and transition_name.startswith("d_"):
            continue
        name = transition_name.format(state=state)
        indices, values = _get_indices(basis, transition_name, state)
        if values:
            matrices[name][indices] = values
    return matrices


def _build_exchange_matrices(basis: Basis) -> dict[str, np.ndarray]:
    matrices: dict[str, np.ndarray] = {}
    for (i1, s1), (i2, s2) in permutations(enumerate(model.states), r=2):
        name = f"k{s1}{s2}"
        matrix = np.zeros((len(model.states), len(model.states)))
        matrix[((i1, i2), i1)] = -1.0, 1.0
        matrices[name] = np.kron(matrix, np.eye(len(basis)))
    return matrices


@cache
def _build_matrices(basis: Basis) -> dict[str, np.ndarray]:
    return _build_exchange_matrices(basis) | _build_spin_matrices(basis)


@dataclass(frozen=True)
class Basis:
    type: str
    extension: Literal["", "dq", "tq"] = ""
    spin_system: str = ""

    @property
    def name(self):
        return ".".join([self.type, self.extension, self.spin_system])

    @property
    def components(self):
        return _BASES[self.type]

    @property
    def atoms(self):
        return {
            letter: atom
            for letter, atom in _ATOMS.get(self.spin_system, {}).items()
            if letter in self.type
        }

    @property
    def vectors(self) -> dict[str, np.ndarray]:
        return dict(_build_vectors(self))

    @property
    def matrices(self) -> dict[str, np.ndarray]:
        return _build_matrices(self)

    @property
    def required_names(self) -> set[str]:
        required_names = set(self.matrices)
        required_names |= {f"p{state}" for state in model.states}
        return required_names

    def __len__(self):
        return len(self.components)

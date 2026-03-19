from __future__ import annotations

import numpy as np

from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.constants import XI_RATIO
from chemex.nmr.magnetization import (
    build_equilibrium_magnetization,
    build_start_magnetization,
    collapse_magnetization,
    detect_signal,
    keep_components,
)
from chemex.parameters.spin_system.nucleus import Nucleus


def test_collapse_magnetization_applies_distribution_weights() -> None:
    magnetization = np.array(
        [
            [[1.0], [2.0]],
            [[10.0], [20.0]],
        ]
    )
    weights = np.array(
        [
            [[0.25]],
            [[0.75]],
        ]
    )

    collapsed = collapse_magnetization(magnetization, weights)

    np.testing.assert_allclose(collapsed, [[7.75], [15.5]])


def test_detect_signal_uses_collapsed_magnetization() -> None:
    detection_vector = np.array([[1.0, -1.0]])
    magnetization = np.array(
        [
            [[2.0], [1.0]],
            [[5.0], [2.0]],
        ]
    )
    weights = np.array(
        [
            [[0.25]],
            [[0.75]],
        ]
    )

    detected = detect_signal(detection_vector, magnetization, weights)

    assert detected == 2.5


def test_build_equilibrium_magnetization_uses_populations_and_xi_ratio() -> None:
    basis = Basis(type="iz_eq", spin_system="nh", model=ModelSpec())
    par_values = {"pa": 0.6, "pb": 0.4}

    magnetization = build_equilibrium_magnetization(basis, par_values)
    expected_scale_a = 0.6 * XI_RATIO[Nucleus.N15]
    expected_scale_b = 0.4 * XI_RATIO[Nucleus.N15]
    expected = (
        expected_scale_a * (basis.vectors["ie_a"] + basis.vectors["iz_a"])
        + expected_scale_b * (basis.vectors["ie_b"] + basis.vectors["iz_b"])
    )

    np.testing.assert_allclose(magnetization, expected)


def test_build_start_magnetization_respects_terms_signs_and_populations() -> None:
    basis = Basis(type="iz", spin_system="nh", model=ModelSpec())
    par_values = {"pa": 0.75, "pb": 0.25}

    magnetization = build_start_magnetization(
        basis,
        par_values,
        terms=["iz_a", "-iz_b"],
        atom=Nucleus.H1,
    )
    expected = 0.75 * basis.vectors["iz_a"] - 0.25 * basis.vectors["iz_b"]

    np.testing.assert_allclose(magnetization, expected)


def test_keep_components_masks_unselected_magnetization() -> None:
    basis = Basis(type="izsz", spin_system="nh", model=ModelSpec())
    magnetization = basis.vectors["iz_a"] + 2.0 * basis.vectors["2izsz_b"]

    kept = keep_components(basis, magnetization, components=["iz_a"])

    np.testing.assert_allclose(kept, basis.vectors["iz_a"])

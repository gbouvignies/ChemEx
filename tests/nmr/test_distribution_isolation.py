from __future__ import annotations

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr._engine.engine import ISLiouvillianEngine
from chemex.nmr.basis import Basis
from chemex.nmr.constants import Distribution
from chemex.parameters.spin_system import SpinSystem


def test_distribution_copies_inputs_and_exposes_read_only_arrays() -> None:
    values = np.array([1.0, 2.0])
    weights = np.array([0.25, 0.75])
    distribution = Distribution(values, weights)

    values[0] = 99.0
    weights[0] = 1.0

    np.testing.assert_allclose(distribution.values, [1.0, 2.0])
    np.testing.assert_allclose(distribution.weights, [0.25, 0.75])

    with pytest.raises(ValueError, match="read-only"):
        distribution.values[0] = 3.0


def test_set_b1_distribution_does_not_mutate_input_distribution() -> None:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    liouvillian = ISLiouvillianEngine(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    distribution = Distribution(np.array([10.0, 12.0]), np.array([0.4, 0.6]))

    liouvillian.set_b1_i_distribution(distribution)

    assert distribution.values.shape == (2,)
    assert distribution.weights.shape == (2,)
    assert liouvillian.l_b1x_i.shape[0] == 2


def test_set_jeff_distribution_does_not_mutate_input_distribution() -> None:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    liouvillian = ISLiouvillianEngine(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    distribution = Distribution(np.array([1.0, 2.0]), np.array([0.5, 0.5]))

    liouvillian.jeff_i = distribution

    assert distribution.values.shape == (2,)
    assert distribution.weights.shape == (2,)


def test_distribution_equality_compares_array_contents() -> None:
    left = Distribution(np.array([10.0, 12.0]), np.array([0.4, 0.6]))
    same = Distribution(np.array([10.0, 12.0]), np.array([0.4, 0.6]))
    different = Distribution(np.array([10.0, 12.0]), np.array([0.5, 0.5]))

    assert left == same
    assert left != different

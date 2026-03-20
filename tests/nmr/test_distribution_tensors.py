from __future__ import annotations

import numpy as np
import pytest

from chemex.nmr._engine.tensors import (
    B1DistributionState,
    DistributionAxis,
    JeffDistributionState,
)
from chemex.nmr.constants import Distribution


def test_distribution_axis_rejects_non_positive_ndim() -> None:
    distribution = Distribution(np.array([1.0]), np.array([1.0]))

    with pytest.raises(ValueError, match="must be positive"):
        DistributionAxis.create(distribution, ndim=0)


def test_b1_distribution_state_broadcasts_values_weights_and_rf_terms() -> None:
    distribution = Distribution(np.array([0.8, 1.2]), np.array([0.25, 0.75]))
    matrix_x = np.array([[1.0, 2.0], [3.0, 4.0]])
    matrix_y = np.array([[5.0, 6.0], [7.0, 8.0]])

    state = B1DistributionState.build(
        distribution,
        matrix_x=matrix_x,
        matrix_y=matrix_y,
    )

    assert state.axis.values.shape == (2, 1, 1)
    assert state.axis.weights.shape == (2, 1, 1)
    np.testing.assert_allclose(state.l_x[:, 0, 0], distribution.values)
    np.testing.assert_allclose(state.l_y[:, 0, 0], 5.0 * distribution.values)


def test_jeff_distribution_state_broadcasts_values_weights_and_liouvillian() -> None:
    distribution = Distribution(np.array([1.0, 2.0]), np.array([0.4, 0.6]))
    matrix = np.array([[2.0]])

    state = JeffDistributionState.build(distribution, matrix=matrix)

    assert state.axis.values.shape == (2, 1, 1, 1)
    assert state.axis.weights.shape == (2, 1, 1, 1)
    np.testing.assert_allclose(state.liouvillian[:, 0, 0, 0], 2.0 * distribution.values)

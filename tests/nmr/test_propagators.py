from __future__ import annotations

import numpy as np
from scipy.linalg import expm

from chemex.nmr.propagators import calculate_propagators


def test_calculate_propagators_single_delay_matches_expm() -> None:
    liouv = np.array([[0.0, 1.0], [-2.0, -3.0]])
    delay = 0.25

    propagator = calculate_propagators(liouv, delay)
    expected = expm(liouv * delay)

    np.testing.assert_allclose(propagator, expected)


def test_calculate_propagators_multiple_delays_match_expm() -> None:
    liouv = np.array([[0.0, 1.0], [-2.0, -3.0]])
    delays = np.array([0.1, 0.25, 0.5])

    propagators = calculate_propagators(liouv, delays)
    expected = np.array([expm(liouv * delay) for delay in delays])

    np.testing.assert_allclose(propagators, expected)


def test_calculate_propagators_handles_stacked_liouvillians() -> None:
    liouvillians = np.array(
        [
            [[0.0, 1.0], [-2.0, -3.0]],
            [[-1.0, 0.5], [-0.25, -2.0]],
        ]
    )
    delays = np.array([0.1, 0.25])

    propagators = calculate_propagators(liouvillians, delays)
    expected = np.array(
        [
            [expm(liouv * delay) for liouv in liouvillians]
            for delay in delays
        ]
    )

    np.testing.assert_allclose(propagators, expected)


def test_calculate_propagators_handles_stacked_liouvillians_for_scalar_delay() -> None:
    liouvillians = np.array(
        [
            [[0.0, 1.0], [-2.0, -3.0]],
            [[-1.0, 0.5], [-0.25, -2.0]],
        ]
    )
    delay = 0.25

    propagators = calculate_propagators(liouvillians, delay)
    expected = np.array([expm(liouv * delay) for liouv in liouvillians])

    np.testing.assert_allclose(propagators, expected)


def test_calculate_propagators_dephasing_returns_finite_propagators() -> None:
    liouvillians = np.array(
        [
            [[0.0, 1.0], [-1.0, 0.0]],
            [[0.0, 2.0], [-2.0, 0.0]],
        ]
    )
    delay = 0.25

    propagators = calculate_propagators(liouvillians, delay, dephasing=True)

    assert propagators.shape == (2, 2, 2)
    assert propagators.dtype == np.float64
    assert np.isfinite(propagators).all()
    assert not np.allclose(propagators, calculate_propagators(liouvillians, delay))

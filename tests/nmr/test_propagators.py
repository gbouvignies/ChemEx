from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import expm

from chemex.nmr._pulses.propagators import calculate_propagators


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
        [[expm(liouv * delay) for liouv in liouvillians] for delay in delays]
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


def test_complete_dephasing_removes_undamped_oscillatory_modes() -> None:
    liouvillians = np.array(
        [
            [[0.0, 1.0], [-1.0, 0.0]],
            [[0.0, 2.0], [-2.0, 0.0]],
        ]
    )
    delays = [0.25, 0.5]

    propagators = calculate_propagators(liouvillians, delays, dephasing=True)

    assert propagators.shape == (2, 2, 2, 2)
    assert propagators.dtype == np.float64
    np.testing.assert_array_equal(propagators, np.zeros_like(propagators))


@pytest.mark.parametrize("damping", [0.0, 1.0e-6, 3.0])
def test_complete_dephasing_agrees_with_broad_uniform_b1_limit(
    damping: float,
) -> None:
    delay = 0.25
    liouv = np.array([[-damping, 1.0], [-1.0, -damping]])
    # The finite uniform ensemble tends to zero as its B1 width increases.
    frequencies = np.linspace(-1000.0, 1000.0, 20001)
    phases = frequencies * delay
    ensemble = np.exp(-damping * delay) * np.array(
        [
            [np.cos(phases).mean(), np.sin(phases).mean()],
            [-np.sin(phases).mean(), np.cos(phases).mean()],
        ]
    )

    np.testing.assert_array_equal(
        calculate_propagators(liouv, delay, dephasing=True), np.zeros((2, 2))
    )
    assert np.linalg.norm(ensemble) < 0.01


def test_complete_dephasing_matches_damped_legacy_limit() -> None:
    liouv = np.array([[-3.0, 1.0, 0.0], [-1.0, -3.0, 0.0], [0.0, 0.0, -2.0]])
    expected = np.diag([0.0, 0.0, np.exp(-0.5)])

    np.testing.assert_allclose(
        calculate_propagators(liouv, 0.25, dephasing=True),
        expected,
        rtol=0.0,
        atol=2.0e-15,
    )


@pytest.mark.parametrize(("frequency", "retained"), [(0.5e-6, True), (2.0e-6, False)])
def test_complete_dephasing_uses_historical_strict_imaginary_cutoff(
    frequency: float, retained: bool
) -> None:
    liouv = np.array([[0.0, frequency], [-frequency, 0.0]])
    expected = expm(liouv * 0.25) if retained else np.zeros((2, 2))

    np.testing.assert_allclose(
        calculate_propagators(liouv, 0.25, dephasing=True),
        expected,
        rtol=0.0,
        atol=2.0e-15,
    )


def test_complete_dephasing_at_zero_delay_is_identity() -> None:
    liouv = np.array([[0.0, 1.0], [-1.0, 0.0]])

    np.testing.assert_allclose(
        calculate_propagators(liouv, [0.0, 0.25], dephasing=True),
        [np.eye(2), np.zeros((2, 2))],
        rtol=0.0,
        atol=2.0e-15,
    )


def test_complete_dephasing_handles_repeated_modes_and_non_normal_basis() -> None:
    base = np.array(
        [
            [0.0, 2.0, 0.0, 0.0, 0.0],
            [-2.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 2.0, 0.0],
            [0.0, 0.0, -2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, -1.0],
        ]
    )
    similarity = np.eye(5)
    similarity[0, 4] = 3.0
    liouv = similarity @ base @ np.linalg.inv(similarity)
    expected = (
        similarity
        @ np.diag([0.0, 0.0, 0.0, 0.0, np.exp(-0.25)])
        @ np.linalg.inv(similarity)
    )

    np.testing.assert_allclose(
        calculate_propagators(liouv, 0.25, dephasing=True),
        expected,
        rtol=0.0,
        atol=1.0e-14,
    )


def test_complete_dephasing_handles_nearly_defective_oscillatory_block() -> None:
    frequency = 2.0e-6
    liouv = np.array([[0.0, 1.0], [-(frequency**2), 0.0]])
    assert np.linalg.cond(np.linalg.eig(liouv).eigenvectors) > 1.0e5

    np.testing.assert_allclose(
        calculate_propagators(liouv, 0.25, dephasing=True),
        np.zeros((2, 2)),
        rtol=0.0,
        atol=1.0e-12,
    )

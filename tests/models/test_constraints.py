"""Scientific tests for kinetic population constraints."""

from __future__ import annotations

import numpy as np
import pytest

from chemex.models.constraints import pop_3st

Rates3State = tuple[float, float, float, float, float, float]


def _assert_stationary(rates: Rates3State, populations: dict[str, float]) -> None:
    kab, kba, kac, kca, kbc, kcb = rates
    scale = max(rates)
    kab, kba, kac, kca, kbc, kcb = (
        rate / scale for rate in (kab, kba, kac, kca, kbc, kcb)
    )
    pa, pb, pc = (populations[name] for name in ("pa", "pb", "pc"))
    residual = np.array(
        (
            -(kab + kac) * pa + kba * pb + kca * pc,
            kab * pa - (kba + kbc) * pb + kcb * pc,
            kac * pa + kbc * pb - (kca + kcb) * pc,
        )
    )

    assert sum(populations.values()) == pytest.approx(1.0, abs=1e-15)
    assert all(population >= 0.0 for population in populations.values())
    np.testing.assert_allclose(residual, 0.0, atol=1e-15)


def test_positive_connected_edge_is_not_classified_as_absent() -> None:
    rates = (25.0, 25.0, 0.0, 0.0, 4.0e-49, 4.0e-49)

    populations = pop_3st(*rates)

    assert populations == pytest.approx({"pa": 1 / 3, "pb": 1 / 3, "pc": 1 / 3})
    _assert_stationary(rates, populations)


def test_literal_zero_represents_a_structurally_absent_edge() -> None:
    rates = (4.0, 2.0, 0.0, 0.0, 3.0, 6.0)

    populations = pop_3st(*rates)

    assert populations == pytest.approx({"pa": 0.25, "pb": 0.5, "pc": 0.25})
    _assert_stationary(rates, populations)


def test_literal_zeros_preserve_disconnected_state_behavior() -> None:
    assert pop_3st(4.0, 2.0, 0.0, 0.0, 0.0, 0.0) == pytest.approx(
        {"pa": 1 / 3, "pb": 2 / 3, "pc": 0.0}
    )


@pytest.mark.parametrize("scale", [1.0e-250, 1.0, 1.0e250])
def test_stationary_populations_are_invariant_to_rate_scale(scale: float) -> None:
    base_rates: Rates3State = (2.5, 8.0, 1.25, 4.0, 3.0, 7.5)
    rates: Rates3State = (
        base_rates[0] * scale,
        base_rates[1] * scale,
        base_rates[2] * scale,
        base_rates[3] * scale,
        base_rates[4] * scale,
        base_rates[5] * scale,
    )

    populations = pop_3st(*rates)
    reference = pop_3st(*base_rates)

    assert populations == pytest.approx(reference, rel=1e-14)
    _assert_stationary(rates, populations)


def test_ordinary_three_state_populations_match_linear_solve() -> None:
    rates: Rates3State = (2.5, 8.0, 1.25, 4.0, 3.0, 7.5)
    kab, kba, kac, kca, kbc, kcb = rates
    matrix = np.array(
        [
            [-kab - kac, kba, kca],
            [kab, -kba - kbc, kcb],
            [1.0, 1.0, 1.0],
        ]
    )
    expected = np.linalg.solve(matrix, np.array([0.0, 0.0, 1.0]))

    populations = pop_3st(*rates)

    np.testing.assert_allclose(tuple(populations.values()), expected, rtol=1e-14)
    _assert_stationary(rates, populations)

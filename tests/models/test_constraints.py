"""Scientific tests for kinetic population constraints."""

from __future__ import annotations

import sys

import numpy as np
import pytest

from chemex.models.constraints import pop_2st, pop_3st

Rates3State = tuple[float, float, float, float, float, float]


@pytest.mark.parametrize(
    ("kab", "kba"),
    [(1.0, 1.0), (1.0, 3.0), (1.0, 1000.0), (1000.0, 1.0)],
)
def test_two_state_positive_rates_have_stationary_ratio(kab: float, kba: float) -> None:
    populations = pop_2st(kab, kba)
    swapped = pop_2st(kba, kab)

    assert sum(populations.values()) == pytest.approx(1.0, abs=1e-15)
    assert all(population >= 0.0 for population in populations.values())
    assert populations == pytest.approx(
        {"pa": kba / (kab + kba), "pb": kab / (kab + kba)}
    )
    assert swapped == pytest.approx({"pa": populations["pb"], "pb": populations["pa"]})
    assert kab * populations["pa"] == pytest.approx(kba * populations["pb"])


@pytest.mark.parametrize("scale", [1.0, 1.0e-8, 1.0e-20, 1.0e-250])
def test_two_state_populations_are_invariant_to_rate_scale(scale: float) -> None:
    populations = pop_2st(3.0 * scale, 7.0 * scale)

    assert populations == pytest.approx(pop_2st(3.0, 7.0), rel=1e-15)


def test_two_state_equal_rates_are_continuous_across_historical_cutoff() -> None:
    cutoff = 1.0e-8
    rates = (
        np.nextafter(cutoff, 0.0),
        cutoff,
        np.nextafter(cutoff, np.inf),
    )

    for rate in rates:
        assert pop_2st(rate, rate) == pytest.approx({"pa": 0.5, "pb": 0.5})


def test_two_state_asymmetric_rates_avoid_historical_cutoff_plateau() -> None:
    cutoff = 1.0e-8
    scales = (
        np.nextafter(cutoff / 3.0, 0.0),
        cutoff / 3.0,
        np.nextafter(cutoff / 3.0, np.inf),
        np.nextafter(cutoff, 0.0),
        cutoff,
        np.nextafter(cutoff, np.inf),
    )

    for scale in scales:
        assert pop_2st(3.0 * scale, scale) == pytest.approx({"pa": 0.25, "pb": 0.75})


def test_two_state_subnormal_positive_rates_remain_active() -> None:
    smallest_subnormal = np.nextafter(0.0, 1.0)

    assert pop_2st(smallest_subnormal, smallest_subnormal) == pytest.approx(
        {"pa": 0.5, "pb": 0.5}
    )
    assert pop_2st(smallest_subnormal, 2.0 * smallest_subnormal) == pytest.approx(
        {"pa": 2 / 3, "pb": 1 / 3}
    )


def test_two_state_maximum_finite_rates_do_not_overflow() -> None:
    maximum = sys.float_info.max

    assert pop_2st(maximum, maximum) == pytest.approx({"pa": 0.5, "pb": 0.5})
    assert pop_2st(maximum / 2.0, maximum) == pytest.approx({"pa": 2 / 3, "pb": 1 / 3})


@pytest.mark.parametrize("positive", [1.0, np.nextafter(0.0, 1.0)])
def test_two_state_exact_one_way_endpoints(positive: float) -> None:
    assert pop_2st(0.0, positive) == {"pa": 1.0, "pb": 0.0}
    assert pop_2st(positive, 0.0) == {"pa": 0.0, "pb": 1.0}


def test_two_state_both_zero_preserves_compatibility_fallback() -> None:
    # This historical fallback is not a unique stationary distribution.
    assert pop_2st(0.0, 0.0) == {"pa": 1.0, "pb": 0.0}


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


def test_slow_reduced_two_state_edge_preserves_rate_ratio() -> None:
    assert pop_3st(4.0e-9, 2.0e-9, 0.0, 0.0, 0.0, 0.0) == pytest.approx(
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

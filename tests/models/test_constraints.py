"""Scientific tests for kinetic population constraints."""

from __future__ import annotations

import math
import sys

import numpy as np
import pytest

from chemex.models.constraints import pop_2st, pop_3st, pop_4st

Rates3State = tuple[float, float, float, float, float, float]
Rates4State = tuple[
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
]


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


def _assert_stationary_2st(
    rates: tuple[float, float],
    populations: dict[str, float],
) -> None:
    kab, kba = rates
    scale = max(rates)
    if scale > 0.0:
        kab /= scale
        kba /= scale
    pa, pb = populations["pa"], populations["pb"]
    residual = np.array((-kab * pa + kba * pb, kab * pa - kba * pb))

    assert sum(populations.values()) == pytest.approx(1.0, abs=1e-15)
    assert all(population >= 0.0 for population in populations.values())
    np.testing.assert_allclose(residual, 0.0, atol=2e-15)


def _assert_stationary(rates: Rates3State, populations: dict[str, float]) -> None:
    kab, kba, kac, kca, kbc, kcb = rates
    scale = max(rates)
    if scale > 0.0:
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
    np.testing.assert_allclose(residual, 0.0, atol=2e-15)


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


def _assert_stationary_4st(
    rates: Rates4State,
    populations: dict[str, float],
) -> None:
    kab, kba, kac, kca, kad, kda, kbc, kcb, kbd, kdb, kcd, kdc = rates
    scale = max(rates)
    if scale > 0.0:
        kab, kba, kac, kca, kad, kda, kbc, kcb, kbd, kdb, kcd, kdc = (
            rate / scale
            for rate in (kab, kba, kac, kca, kad, kda, kbc, kcb, kbd, kdb, kcd, kdc)
        )
    pa, pb, pc, pd = (populations[name] for name in ("pa", "pb", "pc", "pd"))
    residual = np.array(
        (
            -(kab + kac + kad) * pa + kba * pb + kca * pc + kda * pd,
            kab * pa - (kba + kbc + kbd) * pb + kcb * pc + kdb * pd,
            kac * pa + kbc * pb - (kca + kcb + kcd) * pc + kdc * pd,
            kad * pa + kbd * pb + kcd * pc - (kda + kdb + kdc) * pd,
        )
    )

    assert sum(populations.values()) == pytest.approx(1.0, abs=1e-15)
    assert all(population >= 0.0 for population in populations.values())
    np.testing.assert_allclose(residual, 0.0, atol=2e-15)


def test_four_state_populations_match_independent_linear_solve() -> None:
    rates: Rates4State = (
        2.5,
        8.0,
        1.25,
        4.0,
        0.75,
        3.5,
        3.0,
        7.5,
        1.5,
        5.0,
        2.0,
        6.0,
    )
    kab, kba, kac, kca, kad, kda, kbc, kcb, kbd, kdb, kcd, kdc = rates
    matrix = np.array(
        [
            [-kab - kac - kad, kba, kca, kda],
            [kab, -kba - kbc - kbd, kcb, kdb],
            [kac, kbc, -kca - kcb - kcd, kdc],
            [1.0, 1.0, 1.0, 1.0],
        ]
    )
    expected = np.linalg.solve(matrix, np.array([0.0, 0.0, 0.0, 1.0]))

    populations = pop_4st(*rates)

    np.testing.assert_allclose(tuple(populations.values()), expected, rtol=1e-14)
    _assert_stationary_4st(rates, populations)


@pytest.mark.parametrize(
    "scale",
    [
        np.nextafter(0.0, 1.0),
        sys.float_info.min,
        1.0e-250,
        1.0e-8,
        5.0e-4,
        1.0,
        1.0e250,
        sys.float_info.max / 32.0,
    ],
)
def test_four_state_populations_are_invariant_to_rate_scale(scale: float) -> None:
    # Exact small integers avoid binary64 subnormal input quantization when scale = nextafter(0, 1).
    base_rates: Rates4State = (
        10.0,
        32.0,
        5.0,
        16.0,
        3.0,
        14.0,
        12.0,
        30.0,
        6.0,
        20.0,
        8.0,
        24.0,
    )
    rates = tuple(rate * scale for rate in base_rates)

    populations = pop_4st(*rates)
    reference = pop_4st(*base_rates)

    assert populations == pytest.approx(reference, rel=1e-14)
    _assert_stationary_4st(rates, populations)


@pytest.mark.parametrize(
    "rate",
    [np.nextafter(0.0, 1.0), 1.0e-8, 5.0e-4, 1.0, sys.float_info.max],
)
def test_four_state_equal_rates_avoid_determinant_cutoff(rate: float) -> None:
    rates: Rates4State = (rate,) * 12

    populations = pop_4st(*rates)

    assert populations == pytest.approx(
        {"pa": 0.25, "pb": 0.25, "pc": 0.25, "pd": 0.25}
    )
    _assert_stationary_4st(rates, populations)


@pytest.mark.parametrize(
    ("rates", "expected"),
    [
        (
            (1.0e-9, 1.0e-9, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0),
            {"pa": 0.25, "pb": 0.25, "pc": 0.25, "pd": 0.25},
        ),
        (
            (4.0e-49, 4.0e-49, 0.0, 0.0, 0.0, 0.0, 25.0, 25.0, 0.0, 0.0, 25.0, 25.0),
            {"pa": 0.25, "pb": 0.25, "pc": 0.25, "pd": 0.25},
        ),
        (
            (
                2.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                3.0e-200,
                1.0e-200,
                0.0,
                0.0,
                4.0e-200,
                1.0e-200,
            ),
            {"pa": 1 / 33, "pb": 2 / 33, "pc": 6 / 33, "pd": 24 / 33},
        ),
    ],
)
def test_four_state_slow_positive_connected_edges_are_not_classified_as_absent(
    rates: Rates4State,
    expected: dict[str, float],
) -> None:
    populations = pop_4st(*rates)

    assert populations == pytest.approx(expected, rel=1e-13)
    _assert_stationary_4st(rates, populations)


def test_four_state_literal_zero_isolated_states_reduce_cleanly() -> None:
    # State A isolated, B-C-D active with equal rates.
    assert pop_4st(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0) == (
        pytest.approx({"pa": 0.0, "pb": 1 / 3, "pc": 1 / 3, "pd": 1 / 3})
    )
    # State B isolated, A-C-D active.
    assert pop_4st(0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0) == (
        pytest.approx({"pa": 1 / 3, "pb": 0.0, "pc": 1 / 3, "pd": 1 / 3})
    )
    # State C isolated, A-B-D active.
    assert pop_4st(2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0) == (
        pytest.approx({"pa": 1 / 3, "pb": 1 / 3, "pc": 0.0, "pd": 1 / 3})
    )
    # State D isolated, A-B-C active.
    assert pop_4st(2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0) == (
        pytest.approx({"pa": 1 / 3, "pb": 1 / 3, "pc": 1 / 3, "pd": 0.0})
    )
    # States C and D isolated, slow active A-B edge preserves ratio.
    assert pop_4st(
        4.0e-9, 2.0e-9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    ) == pytest.approx({"pa": 1 / 3, "pb": 2 / 3, "pc": 0.0, "pd": 0.0})


def test_all_zero_rates_preserve_stationary_pa_compatibility_convention() -> None:
    pop2 = pop_2st()
    pop3 = pop_3st()
    pop4 = pop_4st()

    assert pop2 == {"pa": 1.0, "pb": 0.0}
    assert pop3 == {"pa": 1.0, "pb": 0.0, "pc": 0.0}
    assert pop4 == {"pa": 1.0, "pb": 0.0, "pc": 0.0, "pd": 0.0}
    _assert_stationary_2st((0.0, 0.0), pop2)
    _assert_stationary((0.0,) * 6, pop3)
    _assert_stationary_4st((0.0,) * 12, pop4)


@pytest.mark.parametrize(
    "rates",
    [
        # B <- A -> C: two disjoint absorbing sinks {B} and {C}
        (1.0, 0.0, 1.0, 0.0, 0.0, 0.0),
        # A <- B -> C: two disjoint absorbing sinks {A} and {C}
        (0.0, 1.0, 0.0, 0.0, 1.0, 0.0),
        # A <- C -> B: two disjoint absorbing sinks {A} and {B}
        (0.0, 0.0, 0.0, 1.0, 0.0, 1.0),
    ],
)
def test_three_state_multi_sink_networks_raise_domain_error(
    rates: Rates3State,
) -> None:
    with pytest.raises(ValueError, match="unique stationary distribution"):
        pop_3st(*rates)


@pytest.mark.parametrize(
    "rates",
    [
        # Disconnected {A <-> B} and {C <-> D}: two closed communicating classes
        (2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 1.0),
        # Disconnected {A <-> C} and {B <-> D}: two closed communicating classes
        (0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 3.0, 1.0, 0.0, 0.0),
        # Disconnected {A <-> D} and {B <-> C}: two closed communicating classes
        (0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0),
        # Disconnected one-way A -> B and C -> D: two absorbing sinks {B} and {D}
        (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        # A drains into both {B} and {C <-> D}: two closed communicating classes
        (1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 1.0),
        # State D isolated while active {A, B, C} has two sinks B <- A -> C
        (1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    ],
)
def test_four_state_disconnected_or_multi_sink_networks_raise_domain_error(
    rates: Rates4State,
) -> None:
    with pytest.raises(ValueError, match="unique stationary distribution"):
        pop_4st(*rates)


def test_returned_equilibrium_populations_always_satisfy_stationarity() -> None:
    two_state_cases: tuple[tuple[float, float], ...] = (
        (0.0, 0.0),
        (0.0, 5.0),
        (5.0, 0.0),
        (2.0, 7.0),
        (1.0e-200, 3.0e-200),
    )
    for rates_2st in two_state_cases:
        _assert_stationary_2st(rates_2st, pop_2st(*rates_2st))

    three_state_cases: tuple[Rates3State, ...] = (
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        # Single absorbing sink A (B -> A <- C)
        (0.0, 2.0, 0.0, 3.0, 0.0, 0.0),
        # Single absorbing sink C (A -> B -> C)
        (2.0, 0.0, 0.0, 0.0, 5.0, 0.0),
        # Irreversible cycle A -> B -> C -> A
        (2.0, 0.0, 0.0, 4.0, 3.0, 0.0),
        # Isolated state C with active A <-> B
        (4.0, 2.0, 0.0, 0.0, 0.0, 0.0),
        # Isolated state A with single-sink B -> C
        (0.0, 0.0, 0.0, 0.0, 5.0, 0.0),
    )
    for rates_3st in three_state_cases:
        _assert_stationary(rates_3st, pop_3st(*rates_3st))

    four_state_cases: tuple[Rates4State, ...] = (
        (0.0,) * 12,
        # Single absorbing sink A (B -> A, C -> A, D -> A)
        (0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        # Single absorbing sink D along chain A -> B -> C -> D
        (2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0, 0.0),
        # Irreversible 4-cycle A -> B -> C -> D -> A
        (2.0, 0.0, 0.0, 0.0, 0.0, 5.0, 3.0, 0.0, 0.0, 0.0, 4.0, 0.0),
        # Transient A feeding closed reversible {C <-> D}, with B -> A
        (0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 1.0),
        # State C isolated while B <-> A -> D drains into unique sink D
        (2.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        # States C and D isolated with active A <-> B
        (4.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    )
    for rates_4st in four_state_cases:
        _assert_stationary_4st(rates_4st, pop_4st(*rates_4st))


@pytest.mark.parametrize(
    "invalid_rate",
    [-1.0, -np.nextafter(0.0, 1.0), math.nan, math.inf, -math.inf],
)
def test_population_solvers_reject_negative_or_non_finite_rates(
    invalid_rate: float,
) -> None:
    with pytest.raises(ValueError, match="finite and non-negative"):
        pop_2st(invalid_rate, 1.0)
    with pytest.raises(ValueError, match="finite and non-negative"):
        pop_3st(1.0, invalid_rate, 1.0, 1.0, 1.0, 1.0)
    with pytest.raises(ValueError, match="finite and non-negative"):
        pop_4st(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, invalid_rate)

"""Fast scientific-domain tests for the generic N-state model authority."""

from __future__ import annotations

import math
import struct
from decimal import Decimal, localcontext

import pytest

from chemex.models.kinetic.settings_nst import (
    calculate_pair_rates,
    calculate_population_complement,
    register,
)
from chemex.parameters.userfunctions import (
    AnalyticFunctionLinearization,
    function_linearization_registry,
)


def _positive_ulp_distance(left: float, right: float) -> int:
    """Return the binary64 ULP distance between finite nonnegative values."""
    assert math.isfinite(left) and left >= 0.0
    assert math.isfinite(right) and right >= 0.0
    left_bits = struct.unpack(">Q", struct.pack(">d", left))[0]
    right_bits = struct.unpack(">Q", struct.pack(">d", right))[0]
    return abs(left_bits - right_bits)


def test_population_complement_accepts_closed_simplex_boundaries() -> None:
    assert calculate_population_complement(0.0, 0.0)["pa"] == 1.0
    assert calculate_population_complement(0.4, 0.6)["pa"] == 0.0


@pytest.mark.parametrize(
    "populations",
    (
        (-0.1, 0.2),
        (0.7, 0.7),
        (math.inf, 0.0),
        (math.nan, 0.0),
    ),
)
def test_population_complement_rejects_values_outside_simplex(
    populations: tuple[float, ...],
) -> None:
    with pytest.raises(ValueError, match="population"):
        calculate_population_complement(*populations)


def test_pair_rates_match_independent_ordinary_oracle() -> None:
    rates = calculate_pair_rates(360.0, 0.2, 0.6)

    assert rates == pytest.approx({"forward": 270.0, "reverse": 90.0})
    assert rates["forward"] + rates["reverse"] == 360.0
    assert 0.2 * rates["forward"] == pytest.approx(0.6 * rates["reverse"])


@pytest.mark.parametrize(
    ("p_i", "p_j", "expected"),
    (
        (0.0, 0.4, {"forward": 125.0, "reverse": 0.0}),
        (0.4, 0.0, {"forward": 0.0, "reverse": 125.0}),
    ),
)
def test_pair_rates_preserve_exact_zero_endpoint_semantics(
    p_i: float,
    p_j: float,
    expected: dict[str, float],
) -> None:
    assert calculate_pair_rates(125.0, p_i, p_j) == expected


def test_pair_rates_allow_zero_kex_with_two_zero_endpoints() -> None:
    assert calculate_pair_rates(0.0, 0.0, 0.0) == {
        "forward": 0.0,
        "reverse": 0.0,
    }


def test_pair_rates_reject_positive_kex_with_two_zero_endpoints() -> None:
    with pytest.raises(ValueError, match="both endpoint populations are zero"):
        calculate_pair_rates(1.0, 0.0, 0.0)


@pytest.mark.parametrize(
    ("kex", "p_i", "p_j"),
    (
        (-1.0, 0.2, 0.3),
        (math.inf, 0.2, 0.3),
        (1.0, -0.2, 0.3),
        (1.0, math.nan, 0.3),
    ),
)
def test_pair_rates_reject_invalid_inputs(kex: float, p_i: float, p_j: float) -> None:
    with pytest.raises(ValueError):
        calculate_pair_rates(kex, p_i, p_j)


@pytest.mark.parametrize(
    ("p_i", "p_j", "kex", "expected_forward"),
    (
        (1.0e-120, 2.0e-120, 300.0, 200.0),
        (1.0e-200, 3.0e-200, 400.0, 300.0),
        (math.ulp(0.0), 1.0, 1.0e6, 1.0e6),
    ),
)
def test_pair_rates_have_no_population_denominator_floor(
    p_i: float,
    p_j: float,
    kex: float,
    expected_forward: float,
) -> None:
    rates = calculate_pair_rates(kex, p_i, p_j)

    assert rates["forward"] == pytest.approx(expected_forward, rel=1.0e-15)
    assert rates["forward"] + rates["reverse"] == kex
    assert p_i * rates["forward"] == pytest.approx(
        p_j * rates["reverse"],
        rel=2.0e-15,
        abs=math.ulp(0.0),
    )


def test_minimum_subnormal_pair_rate_matches_high_precision_oracle() -> None:
    p_i = 0.7
    p_j = math.ulp(0.0)
    kex = 1.0e6

    with localcontext() as context:
        context.prec = 200
        expected_forward = float(
            Decimal.from_float(kex)
            * Decimal.from_float(p_j)
            / (Decimal.from_float(p_i) + Decimal.from_float(p_j))
        )
        expected_flux = float(
            Decimal.from_float(p_i)
            * Decimal.from_float(kex)
            * Decimal.from_float(p_j)
            / (Decimal.from_float(p_i) + Decimal.from_float(p_j))
        )

    rates = calculate_pair_rates(kex, p_i, p_j)
    forward_flux = p_i * rates["forward"]
    reverse_flux = p_j * rates["reverse"]

    assert rates["forward"] > 0.0
    assert _positive_ulp_distance(rates["forward"], expected_forward) <= 1
    assert _positive_ulp_distance(forward_flux, expected_flux) <= 1
    assert _positive_ulp_distance(reverse_flux, expected_flux) <= 1
    assert _positive_ulp_distance(forward_flux, reverse_flux) <= 1
    assert rates["forward"] + rates["reverse"] == kex


def test_highly_biased_normal_pair_rate_matches_high_precision_oracle() -> None:
    p_i = 0.7
    p_j = 1.0e-200
    kex = 1.0e6
    with localcontext() as context:
        context.prec = 200
        expected_forward = float(
            Decimal.from_float(kex)
            * Decimal.from_float(p_j)
            / (Decimal.from_float(p_i) + Decimal.from_float(p_j))
        )

    rates = calculate_pair_rates(kex, p_i, p_j)

    assert _positive_ulp_distance(rates["forward"], expected_forward) <= 1
    assert rates["forward"] + rates["reverse"] == kex
    assert p_i * rates["forward"] == pytest.approx(
        p_j * rates["reverse"],
        rel=2.0e-15,
    )


def test_pair_rates_reject_unrepresentable_positive_direction() -> None:
    with pytest.raises(ValueError, match="cannot be represented"):
        calculate_pair_rates(0.5, math.ulp(0.0), 1.0)


def test_registered_pair_rate_partials_match_independent_oracle() -> None:
    register()
    capabilities = {
        (item.function_id, item.component): item
        for item in function_linearization_registry.get("4st")
        if isinstance(item, AnalyticFunctionLinearization)
    }
    forward = capabilities[("pair_rates", "forward")]
    reverse = capabilities[("pair_rates", "reverse")]
    arguments = (360.0, 0.2, 0.6)

    assert tuple(partial(*arguments) for partial in forward.partials) == pytest.approx(
        (0.75, -337.5, 112.5)
    )
    assert tuple(partial(*arguments) for partial in reverse.partials) == pytest.approx(
        (0.25, 337.5, -112.5)
    )


def test_pair_rate_partials_are_finite_at_one_zero_endpoint() -> None:
    register()
    capabilities = {
        (item.function_id, item.component): item
        for item in function_linearization_registry.get("4st")
        if isinstance(item, AnalyticFunctionLinearization)
    }

    forward = capabilities[("pair_rates", "forward")]
    reverse = capabilities[("pair_rates", "reverse")]
    arguments = (125.0, 0.0, 0.4)
    assert tuple(partial(*arguments) for partial in forward.partials) == pytest.approx(
        (1.0, -312.5, 0.0)
    )
    assert tuple(partial(*arguments) for partial in reverse.partials) == pytest.approx(
        (0.0, 312.5, 0.0)
    )


def test_pair_rate_partials_fail_closed_at_a_zero_pair() -> None:
    register()
    capabilities = {
        (item.function_id, item.component): item
        for item in function_linearization_registry.get("4st")
        if isinstance(item, AnalyticFunctionLinearization)
    }

    for component in ("forward", "reverse"):
        for partial in capabilities[("pair_rates", component)].partials:
            with pytest.raises(ValueError, match="undefined at a zero pair"):
                partial(0.0, 0.0, 0.0)


def test_registered_population_complement_partials_preserve_coupled_pa() -> None:
    register()
    capability = next(
        item
        for item in function_linearization_registry.get("6st")
        if isinstance(item, AnalyticFunctionLinearization)
        and (item.function_id, item.component) == ("population_complement", "pa")
    )

    assert (
        tuple(partial(0.1, 0.2, 0.3, 0.1, 0.1) for partial in capability.partials)
        == (-1.0,) * 5
    )

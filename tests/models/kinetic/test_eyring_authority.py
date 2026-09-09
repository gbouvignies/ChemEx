"""Independent numerical tests for the shared Eyring scientific authority."""

from __future__ import annotations

import math
import sys
from decimal import Decimal, localcontext

import pytest
from scipy import constants

from chemex.models.kinetic._eyring import (
    ThermodynamicCoordinate,
    calculate_directional_rate,
    directional_log_rate,
    rate_from_log,
    temperature_to_kelvin,
    thermodynamic_population_partials,
    thermodynamic_populations,
)

ZERO = ThermodynamicCoordinate(enthalpy=0.0, entropy=0.0)


def _decimal_eyring_rate(
    kelvin: str,
    activation_enthalpy: str,
    activation_entropy: str = "0",
) -> Decimal:
    """Evaluate the Eyring equation independently with exact SI literals."""
    with localcontext() as context:
        context.prec = 100
        temperature = Decimal(kelvin)
        gas_constant = Decimal("8.31446261815324")
        frequency_factor = Decimal("1.380649e-23") / Decimal("6.62607015e-34")
        exponent = Decimal(activation_entropy) / gas_constant - Decimal(
            activation_enthalpy
        ) / (gas_constant * temperature)
        return frequency_factor * temperature * exponent.exp()


def _decimal_two_state_population(
    kelvin: Decimal,
    state_enthalpy: float,
    state_entropy: float = 0.0,
) -> Decimal:
    """Return P(B) from an independent exact-SI Boltzmann calculation."""
    with localcontext() as context:
        context.prec = 100
        gas_constant = Decimal("8.31446261815324")
        log_weight = Decimal(str(state_entropy)) / gas_constant - Decimal(
            str(state_enthalpy)
        ) / (gas_constant * kelvin)
        weight = log_weight.exp()
        return weight / (Decimal(1) + weight)


def _enthalpy_for_population(target: Decimal, kelvin: Decimal) -> float:
    """Choose a binary64 H coordinate whose Decimal population is near target."""
    with localcontext() as context:
        context.prec = 100
        gas_constant = Decimal("8.31446261815324")
        log_weight = (target / (Decimal(1) - target)).ln()
        return float(-gas_constant * kelvin * log_weight)


@pytest.mark.parametrize(
    ("celsius", "kelvin"),
    ((25.0, 298.15), (50.0, 323.15), (0.0, 273.15)),
)
def test_public_celsius_temperature_is_converted_to_kelvin(
    celsius: float,
    kelvin: float,
) -> None:
    assert temperature_to_kelvin(celsius) == kelvin


def test_ordinary_rates_match_high_precision_oracle_at_two_temperatures() -> None:
    transition_state = ThermodynamicCoordinate(enthalpy=65_000.0, entropy=0.0)

    rate_298 = calculate_directional_rate(ZERO, transition_state, 25.0)
    rate_323 = calculate_directional_rate(ZERO, transition_state, 50.0)

    assert rate_298 == pytest.approx(float(_decimal_eyring_rate("298.15", "65000")))
    assert rate_323 == pytest.approx(float(_decimal_eyring_rate("323.15", "65000")))
    assert rate_298 == pytest.approx(25.45392415, rel=5.0e-10)
    assert rate_323 == pytest.approx(209.74952472, rel=5.0e-10)


def test_log_k_over_temperature_is_linear_in_inverse_temperature() -> None:
    transition_state = ThermodynamicCoordinate(enthalpy=65_000.0, entropy=12.0)
    temperatures = (280.0, 298.15, 323.15, 360.0)
    ordinates = tuple(
        directional_log_rate(ZERO, transition_state, kelvin) - math.log(kelvin)
        for kelvin in temperatures
    )
    slopes = tuple(
        (right_y - left_y) / (1.0 / right_t - 1.0 / left_t)
        for left_t, right_t, left_y, right_y in zip(
            temperatures[:-1],
            temperatures[1:],
            ordinates[:-1],
            ordinates[1:],
            strict=True,
        )
    )

    assert max(slopes) - min(slopes) <= 2.0e-9 * abs(slopes[0])
    assert slopes[0] == pytest.approx(-65_000.0 / constants.R, rel=2.0e-15)


def test_low_positive_kelvin_is_valid_without_an_arbitrary_upper_limit() -> None:
    lowest_public_celsius = math.nextafter(-constants.zero_Celsius, math.inf)
    low_rate = calculate_directional_rate(ZERO, ZERO, lowest_public_celsius)

    assert temperature_to_kelvin(lowest_public_celsius) > 0.0
    assert 0.0 < low_rate < 1.0
    assert temperature_to_kelvin(1.0e300) == 1.0e300


def test_high_finite_temperature_rate_remains_available_when_representable() -> None:
    celsius = 1.0e100
    transition_state = ThermodynamicCoordinate(enthalpy=0.0, entropy=-2_000.0)

    rate = calculate_directional_rate(ZERO, transition_state, celsius)
    expected = _decimal_eyring_rate("1e100", "0", "-2000")

    assert math.isfinite(rate)
    assert rate == pytest.approx(float(expected), rel=3.0e-14)


def test_direct_factorization_underflow_does_not_erase_representable_rate() -> None:
    kelvin = Decimal("32")
    with localcontext() as context:
        context.prec = 100
        barrier = Decimal("750") * Decimal("8.31446261815324") * kelvin
        expected = _decimal_eyring_rate("32", str(barrier))

    rate = calculate_directional_rate(
        ZERO,
        ThermodynamicCoordinate(enthalpy=float(barrier), entropy=0.0),
        -241.15,
    )

    assert expected > Decimal.from_float(math.nextafter(0.0, 1.0))
    assert expected < Decimal.from_float(sys.float_info.min)
    assert rate > 0.0
    expected_float = float(expected)
    # Four ulps accommodates small cross-libm log/exp differences while remaining
    # billions of ulps away from structural zero for this fixture.
    assert abs(rate - expected_float) <= 4 * math.ulp(expected_float)


def test_binary64_rate_classification_includes_minimum_subnormal() -> None:
    minimum = math.nextafter(0.0, 1.0)
    log_minimum = math.log(minimum)

    assert Decimal.from_float(log_minimum).exp() >= Decimal.from_float(minimum)
    assert rate_from_log(log_minimum) == minimum


def test_binary64_rate_classification_rejects_positive_underflow() -> None:
    minimum = math.nextafter(0.0, 1.0)
    log_below = math.nextafter(math.log(minimum), -math.inf)

    assert Decimal.from_float(log_below).exp() < Decimal.from_float(minimum)
    with pytest.raises(ValueError, match="below binary64 representability"):
        rate_from_log(log_below)


def test_binary64_rate_classification_includes_near_maximum() -> None:
    log_maximum = math.log(sys.float_info.max)
    expected = math.exp(log_maximum)

    assert Decimal.from_float(log_maximum).exp() <= Decimal.from_float(
        sys.float_info.max
    )
    assert rate_from_log(log_maximum) == expected
    assert math.isfinite(expected)


def test_binary64_rate_classification_rejects_overflow() -> None:
    log_above = math.nextafter(math.log(sys.float_info.max), math.inf)

    assert Decimal.from_float(log_above).exp() > Decimal.from_float(sys.float_info.max)
    with pytest.raises(ValueError, match="exceeds maximum finite binary64"):
        rate_from_log(log_above)


def test_state_populations_match_high_precision_boltzmann_oracle() -> None:
    states = {
        "a": ZERO,
        "b": ThermodynamicCoordinate(enthalpy=8_000.0, entropy=10.0),
        "c": ThermodynamicCoordinate(enthalpy=12_000.0, entropy=-5.0),
        "d": ThermodynamicCoordinate(enthalpy=15_000.0, entropy=15.0),
    }
    populations = thermodynamic_populations(states, 25.0)
    with localcontext() as context:
        context.prec = 100
        temperature = Decimal("298.15")
        gas_constant = Decimal("8.31446261815324")
        weights = {
            state: (
                Decimal(str(coordinate.entropy)) / gas_constant
                - Decimal(str(coordinate.enthalpy)) / (gas_constant * temperature)
            ).exp()
            for state, coordinate in states.items()
        }
        total = sum(weights.values())
        expected = {state: float(weight / total) for state, weight in weights.items()}

    assert populations == pytest.approx(expected, rel=3.0e-15)
    assert math.fsum(populations.values()) == pytest.approx(1.0, rel=0.0, abs=2e-16)


def test_representable_subnormal_population_remains_positive() -> None:
    kelvin = Decimal("0.01")
    minimum = Decimal.from_float(math.nextafter(0.0, 1.0))
    enthalpy = _enthalpy_for_population(Decimal(8) * minimum, kelvin)
    expected = _decimal_two_state_population(kelvin, enthalpy)

    populations = thermodynamic_populations(
        {"a": ZERO, "b": ThermodynamicCoordinate(enthalpy, 0.0)},
        float(kelvin) - constants.zero_Celsius,
    )

    assert minimum < expected < Decimal.from_float(sys.float_info.min)
    assert 0.0 < populations["b"] < sys.float_info.min
    assert abs(populations["b"] - float(expected)) <= 4 * math.ulp(float(expected))


def test_minimum_positive_population_boundary_is_representable() -> None:
    kelvin = Decimal("0.01")
    minimum_float = math.nextafter(0.0, 1.0)
    minimum = Decimal.from_float(minimum_float)
    enthalpy = _enthalpy_for_population(Decimal("1.1") * minimum, kelvin)
    expected = _decimal_two_state_population(kelvin, enthalpy)

    populations = thermodynamic_populations(
        {"a": ZERO, "b": ThermodynamicCoordinate(enthalpy, 0.0)},
        float(kelvin) - constants.zero_Celsius,
    )

    assert minimum <= expected < Decimal("1.5") * minimum
    assert populations["b"] == minimum_float


def test_positive_population_just_below_binary64_fails_closed() -> None:
    kelvin = Decimal("0.01")
    minimum = Decimal.from_float(math.nextafter(0.0, 1.0))
    enthalpy = _enthalpy_for_population(Decimal("0.9") * minimum, kelvin)
    expected = _decimal_two_state_population(kelvin, enthalpy)

    assert Decimal(0) < expected < minimum
    with pytest.raises(ValueError, match="population.*below binary64 representability"):
        thermodynamic_populations(
            {"a": ZERO, "b": ThermodynamicCoordinate(enthalpy, 0.0)},
            float(kelvin) - constants.zero_Celsius,
        )


def test_audited_positive_population_underflow_reproducer_fails_closed() -> None:
    kelvin = Decimal("0.01")
    enthalpy = 61.959375430421595
    expected = _decimal_two_state_population(kelvin, enthalpy)

    assert Decimal(0) < expected < Decimal.from_float(math.nextafter(0.0, 1.0))
    assert Decimal("2.3107e-324") < expected < Decimal("2.3108e-324")
    with pytest.raises(ValueError, match="population.*below binary64 representability"):
        thermodynamic_populations(
            {"a": ZERO, "b": ThermodynamicCoordinate(enthalpy, 0.0)},
            float(kelvin) - constants.zero_Celsius,
        )


def test_tiny_representable_population_has_nonzero_normalized_partial() -> None:
    kelvin = Decimal("0.01")
    minimum = Decimal.from_float(math.nextafter(0.0, 1.0))
    enthalpy = _enthalpy_for_population(Decimal(8) * minimum, kelvin)
    expected_population = _decimal_two_state_population(kelvin, enthalpy)
    with localcontext() as context:
        context.prec = 100
        expected_partial = (
            -expected_population
            * (Decimal(1) - expected_population)
            / (Decimal("8.31446261815324") * kelvin)
        )
    partials = thermodynamic_population_partials(("pa", "pb"))

    derivative_a = partials["pa"][0](enthalpy, 0.0, -273.14)
    derivative_b = partials["pb"][0](enthalpy, 0.0, -273.14)

    assert derivative_b < 0.0
    assert derivative_b != -0.0
    assert abs(derivative_b - float(expected_partial)) <= 4 * math.ulp(
        float(expected_partial)
    )
    assert math.fsum((derivative_a, derivative_b)) == 0.0


def test_unrepresentable_population_fails_before_exposing_false_zero_partial() -> None:
    enthalpy = 61.959375430421595
    partial = thermodynamic_population_partials(("pa", "pb"))["pb"][0]

    with pytest.raises(ValueError, match="population.*below binary64 representability"):
        partial(enthalpy, 0.0, -273.14)


def test_population_partial_underflow_does_not_invalidate_population() -> None:
    kelvin = Decimal("298.15")
    minimum = Decimal.from_float(math.nextafter(0.0, 1.0))
    enthalpy = _enthalpy_for_population(Decimal(1000) * minimum, kelvin)
    populations = thermodynamic_populations(
        {"a": ZERO, "b": ThermodynamicCoordinate(enthalpy, 0.0)},
        25.0,
    )
    partial = thermodynamic_population_partials(("pa", "pb"))["pb"][0]

    assert populations["b"] > 0.0
    with pytest.raises(
        ValueError,
        match="population derivative.*below binary64 representability",
    ):
        partial(enthalpy, 0.0, 25.0)


def test_population_partials_preserve_normalized_softmax_coupling() -> None:
    arguments = (8_000.0, 10.0, 12_000.0, -5.0, 25.0)
    partials = thermodynamic_population_partials(("pa", "pb", "pc"))
    jacobian = tuple(
        tuple(partial(*arguments) for partial in partials[component])
        for component in ("pa", "pb", "pc")
    )

    for column in zip(*jacobian, strict=True):
        # Each derivative is O(1e-2); this is below two binary64 ulps of their sum.
        assert math.fsum(column) == pytest.approx(0.0, abs=4.0e-18)

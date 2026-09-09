"""Shared scientific authority for temperature-dependent Eyring models."""

from __future__ import annotations

import math
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

from scipy import constants

MIN_POSITIVE_FLOAT = math.nextafter(0.0, 1.0)
_LOG_MIN_POSITIVE_FLOAT = math.log(MIN_POSITIVE_FLOAT)
_LOG_MAX_FLOAT = math.log(sys.float_info.max)
_LOG_EYRING_FREQUENCY_FACTOR = math.log(constants.k / constants.h)


@dataclass(frozen=True, slots=True)
class ThermodynamicCoordinate:
    """Enthalpy and entropy coordinates relative to reference state A."""

    enthalpy: float
    entropy: float


def _validate_coordinate(
    coordinate: ThermodynamicCoordinate,
    *,
    description: str,
) -> None:
    if not math.isfinite(coordinate.enthalpy) or not math.isfinite(coordinate.entropy):
        msg = f"{description} thermodynamic coordinates must be finite"
        raise ValueError(msg)


def _validate_kelvin(kelvin: float) -> float:
    if not math.isfinite(kelvin) or kelvin <= 0.0:
        msg = "Eyring absolute temperature must be finite and strictly positive"
        raise ValueError(msg)
    return kelvin


def temperature_to_kelvin(temperature: float) -> float:
    """Validate a public Celsius temperature and return absolute temperature."""
    kelvin = float(temperature) + constants.zero_Celsius
    if not math.isfinite(kelvin) or kelvin <= 0.0:
        msg = (
            "Eyring temperature must be finite and above absolute zero "
            f"(-273.15 degrees Celsius), got {temperature!r}"
        )
        raise ValueError(msg)
    return kelvin


def state_log_boltzmann_weight(
    state: ThermodynamicCoordinate,
    kelvin: float,
) -> float:
    """Return ``log(w_i)`` for one state's coordinate relative to state A."""
    _validate_coordinate(state, description="State")
    absolute_temperature = _validate_kelvin(kelvin)
    log_weight = state.entropy / constants.R - state.enthalpy / (
        constants.R * absolute_temperature
    )
    if not math.isfinite(log_weight):
        msg = "Eyring state log Boltzmann weight is outside the finite binary64 domain"
        raise ValueError(msg)
    return log_weight


def directional_log_rate(
    initial_state: ThermodynamicCoordinate,
    transition_state: ThermodynamicCoordinate,
    kelvin: float,
) -> float:
    """Return the Eyring log rate from one state through a shared transition state."""
    _validate_coordinate(initial_state, description="Initial-state")
    _validate_coordinate(transition_state, description="Transition-state")
    absolute_temperature = _validate_kelvin(kelvin)
    activation_enthalpy = transition_state.enthalpy - initial_state.enthalpy
    activation_entropy = transition_state.entropy - initial_state.entropy
    return (
        _LOG_EYRING_FREQUENCY_FACTOR
        + math.log(absolute_temperature)
        + activation_entropy / constants.R
        - activation_enthalpy / (constants.R * absolute_temperature)
    )


def _positive_float_from_log(log_value: float, *, quantity: str) -> float:
    """Convert one positive Eyring quantity from log space or fail closed."""
    if math.isnan(log_value):
        msg = f"Eyring log {quantity} must not be NaN"
        raise ValueError(msg)
    if log_value > _LOG_MAX_FLOAT:
        msg = f"Eyring {quantity} exceeds maximum finite binary64 value"
        raise ValueError(msg)
    if log_value < _LOG_MIN_POSITIVE_FLOAT:
        msg = f"Positive Eyring {quantity} is below binary64 representability"
        raise ValueError(msg)

    value = math.exp(log_value)
    if value == 0.0:
        msg = f"Positive Eyring {quantity} is below binary64 representability"
        raise ValueError(msg)
    if not math.isfinite(value):
        msg = f"Eyring {quantity} exceeds maximum finite binary64 value"
        raise ValueError(msg)
    return value


def rate_from_log(log_rate: float) -> float:
    """Convert a positive mathematical rate to binary64 or fail closed."""
    return _positive_float_from_log(log_rate, quantity="rate")


def calculate_directional_rate(
    initial_state: ThermodynamicCoordinate,
    transition_state: ThermodynamicCoordinate,
    temperature: float,
) -> float:
    """Calculate one directional Eyring rate from public Celsius temperature."""
    kelvin = temperature_to_kelvin(temperature)
    return rate_from_log(directional_log_rate(initial_state, transition_state, kelvin))


def calculate_rate_from_coordinates(
    initial_enthalpy: float,
    initial_entropy: float,
    transition_enthalpy: float,
    transition_entropy: float,
    temperature: float,
) -> float:
    """Scalar constraint-expression adapter for one directional Eyring rate."""
    return calculate_directional_rate(
        ThermodynamicCoordinate(initial_enthalpy, initial_entropy),
        ThermodynamicCoordinate(transition_enthalpy, transition_entropy),
        temperature,
    )


def calculate_rate_component(
    initial_enthalpy: float,
    initial_entropy: float,
    transition_enthalpy: float,
    transition_entropy: float,
    temperature: float,
) -> dict[str, float]:
    """Constraint-expression mapping adapter for one directional Eyring rate."""
    return {
        "rate": calculate_rate_from_coordinates(
            initial_enthalpy,
            initial_entropy,
            transition_enthalpy,
            transition_entropy,
            temperature,
        )
    }


def rate_partial_initial_enthalpy(
    initial_enthalpy: float,
    initial_entropy: float,
    transition_enthalpy: float,
    transition_entropy: float,
    temperature: float,
) -> float:
    rate = calculate_rate_from_coordinates(
        initial_enthalpy,
        initial_entropy,
        transition_enthalpy,
        transition_entropy,
        temperature,
    )
    kelvin = temperature_to_kelvin(temperature)
    return rate / (constants.R * kelvin)


def rate_partial_initial_entropy(
    initial_enthalpy: float,
    initial_entropy: float,
    transition_enthalpy: float,
    transition_entropy: float,
    temperature: float,
) -> float:
    rate = calculate_rate_from_coordinates(
        initial_enthalpy,
        initial_entropy,
        transition_enthalpy,
        transition_entropy,
        temperature,
    )
    return -rate / constants.R


def rate_partial_transition_enthalpy(
    initial_enthalpy: float,
    initial_entropy: float,
    transition_enthalpy: float,
    transition_entropy: float,
    temperature: float,
) -> float:
    return -rate_partial_initial_enthalpy(
        initial_enthalpy,
        initial_entropy,
        transition_enthalpy,
        transition_entropy,
        temperature,
    )


def rate_partial_transition_entropy(
    initial_enthalpy: float,
    initial_entropy: float,
    transition_enthalpy: float,
    transition_entropy: float,
    temperature: float,
) -> float:
    return -rate_partial_initial_entropy(
        initial_enthalpy,
        initial_entropy,
        transition_enthalpy,
        transition_entropy,
        temperature,
    )


def rate_partial_temperature(
    initial_enthalpy: float,
    initial_entropy: float,
    transition_enthalpy: float,
    transition_entropy: float,
    temperature: float,
) -> float:
    rate = calculate_rate_from_coordinates(
        initial_enthalpy,
        initial_entropy,
        transition_enthalpy,
        transition_entropy,
        temperature,
    )
    kelvin = temperature_to_kelvin(temperature)
    activation_enthalpy = transition_enthalpy - initial_enthalpy
    return rate * (1.0 / kelvin + activation_enthalpy / (constants.R * kelvin * kelvin))


EYRING_RATE_PARTIALS = (
    rate_partial_initial_enthalpy,
    rate_partial_initial_entropy,
    rate_partial_transition_enthalpy,
    rate_partial_transition_entropy,
    rate_partial_temperature,
)


def thermodynamic_populations(
    states: Mapping[str, ThermodynamicCoordinate],
    temperature: float,
) -> dict[str, float]:
    """Normalize equilibrium populations from thermodynamic state coordinates."""
    if not states:
        msg = "Eyring population calculation requires at least one state"
        raise ValueError(msg)
    kelvin = temperature_to_kelvin(temperature)
    log_weights = {
        name: state_log_boltzmann_weight(coordinate, kelvin)
        for name, coordinate in states.items()
    }
    reference = max(log_weights.values())
    relative_log_weights = {
        name: log_weight - reference for name, log_weight in log_weights.items()
    }
    log_scaled_total = math.log(
        math.fsum(math.exp(value) for value in relative_log_weights.values())
    )
    # Classify normalized log probabilities before exponentiation so a finite
    # thermodynamic state cannot silently disappear as structural zero.
    return {
        name: _positive_float_from_log(
            relative_log_weight - log_scaled_total,
            quantity="population",
        )
        for name, relative_log_weight in relative_log_weights.items()
    }


def _population_partial(
    state_count: int,
    population_index: int,
    argument_index: int,
) -> Callable[..., float]:
    """Build one exact softmax partial for a flattened Eyring state call."""

    def partial(*arguments: float) -> float:
        expected_arity = 2 * (state_count - 1) + 1
        if len(arguments) != expected_arity:
            msg = (
                "Eyring population derivative received "
                f"{len(arguments)} arguments instead of {expected_arity}"
            )
            raise ValueError(msg)
        temperature = arguments[-1]
        kelvin = temperature_to_kelvin(temperature)
        states = (
            ThermodynamicCoordinate(0.0, 0.0),
            *(
                ThermodynamicCoordinate(arguments[index], arguments[index + 1])
                for index in range(0, expected_arity - 1, 2)
            ),
        )
        populations = tuple(
            thermodynamic_populations(
                {str(index): state for index, state in enumerate(states)},
                temperature,
            ).values()
        )
        population = populations[population_index]

        if argument_index == expected_arity - 1:
            mean_enthalpy = math.fsum(
                probability * state.enthalpy
                for probability, state in zip(populations, states, strict=True)
            )
            centered_enthalpy = states[population_index].enthalpy - mean_enthalpy
            if centered_enthalpy == 0.0:
                return 0.0
            direct = population * centered_enthalpy / (constants.R * kelvin * kelvin)
            if direct != 0.0:
                return direct
            magnitude = _positive_float_from_log(
                math.log(population)
                + math.log(abs(centered_enthalpy))
                - math.log(constants.R)
                - 2.0 * math.log(kelvin),
                quantity="population derivative",
            )
            return math.copysign(magnitude, centered_enthalpy)

        coordinate_index = argument_index // 2 + 1
        coordinate_population = populations[coordinate_index]
        selected = population_index == coordinate_index
        signed_coupling = (
            coordinate_population - float(selected)
            if argument_index % 2 == 0
            else float(selected) - coordinate_population
        )
        denominator = constants.R * kelvin if argument_index % 2 == 0 else constants.R
        direct = population * signed_coupling / denominator
        if direct != 0.0:
            return direct
        coupling = (
            math.fsum(
                probability
                for index, probability in enumerate(populations)
                if index != coordinate_index
            )
            if selected
            else coordinate_population
        )
        magnitude = _positive_float_from_log(
            math.log(population) + math.log(coupling) - math.log(denominator),
            quantity="population derivative",
        )
        if argument_index % 2 == 0:
            return -magnitude if selected else magnitude
        return magnitude if selected else -magnitude

    return partial


def thermodynamic_population_partials(
    components: Sequence[str],
) -> dict[str, tuple[Callable[..., float], ...]]:
    """Return exact coupled derivatives for named normalized populations."""
    names = tuple(components)
    if len(names) < 2 or len(set(names)) != len(names):
        msg = "Eyring populations require at least two distinct component names"
        raise ValueError(msg)
    arity = 2 * (len(names) - 1) + 1
    return {
        name: tuple(
            _population_partial(len(names), population_index, argument_index)
            for argument_index in range(arity)
        )
        for population_index, name in enumerate(names)
    }

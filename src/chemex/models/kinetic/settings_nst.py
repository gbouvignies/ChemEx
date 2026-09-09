"""Generic N-state topology, population, and pairwise exchange authority."""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations, pairwise

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import (
    AnalyticFunctionLinearization,
    FunctionLinearization,
    function_linearization_registry,
    user_function_registry,
)

TPL = ("temperature", "p_total", "l_total")
type Edge = tuple[str, str]


def complete_edges(states: str) -> tuple[Edge, ...]:
    """Return every unordered pair in a complete N-state topology."""
    return tuple(combinations(states, 2))


def linear_edges(states: str) -> tuple[Edge, ...]:
    """Return consecutive pairs in the A-B-C-... path topology."""
    return tuple(pairwise(states))


def fork_edges(states: str) -> tuple[Edge, ...]:
    """Return the A-centered star topology."""
    return tuple((states[0], state) for state in states[1:])


@dataclass(frozen=True, slots=True)
class Topology:
    """One public generic N-state model and its structural exchange edges."""

    model_name: str
    states: str
    edges: tuple[Edge, ...]


# Unsuffixed models are complete. Named variants declare structural absence
# directly; they are not complete graphs with edges deleted afterwards.
TOPOLOGIES = (
    Topology("3st", "abc", complete_edges("abc")),
    Topology("3st_triangle", "abc", complete_edges("abc")),
    Topology("3st_linear", "abc", linear_edges("abc")),
    Topology("3st_fork", "abc", fork_edges("abc")),
    Topology("4st", "abcd", complete_edges("abcd")),
    Topology("4st_linear", "abcd", linear_edges("abcd")),
    Topology("4st_fork", "abcd", fork_edges("abcd")),
    Topology("5st", "abcde", complete_edges("abcde")),
    Topology("5st_linear", "abcde", linear_edges("abcde")),
    Topology("5st_fork", "abcde", fork_edges("abcde")),
    Topology("6st", "abcdef", complete_edges("abcdef")),
    Topology("6st_linear", "abcdef", linear_edges("abcdef")),
    Topology("6st_fork", "abcdef", fork_edges("abcdef")),
)
_TOPOLOGY_BY_NAME = {topology.model_name: topology for topology in TOPOLOGIES}


def calculate_population_complement(
    *independent_populations: float,
) -> dict[str, float]:
    """Validate the closed population simplex and return its state-A complement."""
    if not independent_populations:
        raise ValueError("at least one independent population is required")
    if any(not math.isfinite(value) for value in independent_populations):
        raise ValueError("independent populations must be finite")
    if any(value < 0.0 for value in independent_populations):
        raise ValueError("independent populations must be nonnegative")
    total = math.fsum(independent_populations)
    if total > 1.0:
        raise ValueError("independent population sum must not exceed one")
    pa = 1.0 - total
    if pa < 0.0:
        raise ValueError("derived state-A population must be nonnegative")
    return {"pa": pa}


def _validate_pair_inputs(kex: float, p_i: float, p_j: float) -> None:
    if not math.isfinite(kex):
        raise ValueError("KEX must be finite")
    if kex < 0.0:
        raise ValueError("KEX must be nonnegative")
    if not math.isfinite(p_i) or not math.isfinite(p_j):
        raise ValueError("endpoint populations must be finite")
    if p_i < 0.0 or p_j < 0.0:
        raise ValueError("endpoint populations must be nonnegative")


def _scaled_positive_product_ratio(
    multiplier: float,
    numerator: float,
    denominator: float,
) -> float:
    """Evaluate ``multiplier * numerator / denominator`` with one final scaling.

    Splitting each positive input into a binary mantissa and exponent prevents a
    subnormal population ratio from rounding before the exchange scale is
    applied. ``ldexp`` then performs the only potentially subnormal rounding.
    """
    multiplier_fraction, multiplier_exponent = math.frexp(multiplier)
    numerator_fraction, numerator_exponent = math.frexp(numerator)
    denominator_fraction, denominator_exponent = math.frexp(denominator)
    fraction = multiplier_fraction * numerator_fraction / denominator_fraction
    exponent = multiplier_exponent + numerator_exponent - denominator_exponent
    return math.ldexp(fraction, exponent)


def calculate_pair_rates(kex: float, p_i: float, p_j: float) -> dict[str, float]:
    """Split total pairwise exchange into detailed-balance directions.

    ``forward`` is i -> j and is therefore proportional to ``p_j``.
    ``reverse`` is j -> i and is proportional to ``p_i``.
    """
    _validate_pair_inputs(kex, p_i, p_j)
    if kex == 0.0:
        return {"forward": 0.0, "reverse": 0.0}
    if p_i == 0.0 and p_j == 0.0:
        raise ValueError(
            "directional exchange is undefined when both endpoint populations are zero"
        )
    if p_i == 0.0:
        return {"forward": kex, "reverse": 0.0}
    if p_j == 0.0:
        return {"forward": 0.0, "reverse": kex}

    smaller, larger = sorted((p_i, p_j))
    denominator = math.fsum((smaller, larger))
    smaller_rate = _scaled_positive_product_ratio(kex, smaller, denominator)
    if smaller_rate == 0.0:
        raise ValueError(
            "a mathematically positive directional rate cannot be represented "
            "as positive binary64"
        )
    larger_rate = kex - smaller_rate
    if larger_rate <= 0.0:
        raise ValueError("directional exchange could not be represented in binary64")
    if p_j <= p_i:
        forward, reverse = smaller_rate, larger_rate
    else:
        forward, reverse = larger_rate, smaller_rate
    return {"forward": forward, "reverse": reverse}


def _population_partial(*independent_populations: float) -> float:
    calculate_population_complement(*independent_populations)
    return -1.0


def _pair_derivative_denominator(kex: float, p_i: float, p_j: float) -> float:
    _validate_pair_inputs(kex, p_i, p_j)
    denominator = math.fsum((p_i, p_j))
    if denominator == 0.0:
        raise ValueError("directional-rate derivatives are undefined at a zero pair")
    return denominator


def _forward_d_kex(kex: float, p_i: float, p_j: float) -> float:
    denominator = _pair_derivative_denominator(kex, p_i, p_j)
    return p_j / denominator


def _forward_d_pi(kex: float, p_i: float, p_j: float) -> float:
    denominator = _pair_derivative_denominator(kex, p_i, p_j)
    return -kex * p_j / (denominator * denominator)


def _forward_d_pj(kex: float, p_i: float, p_j: float) -> float:
    denominator = _pair_derivative_denominator(kex, p_i, p_j)
    return kex * p_i / (denominator * denominator)


def _reverse_d_kex(kex: float, p_i: float, p_j: float) -> float:
    denominator = _pair_derivative_denominator(kex, p_i, p_j)
    return p_i / denominator


def _reverse_d_pi(kex: float, p_i: float, p_j: float) -> float:
    denominator = _pair_derivative_denominator(kex, p_i, p_j)
    return kex * p_j / (denominator * denominator)


def _reverse_d_pj(kex: float, p_i: float, p_j: float) -> float:
    denominator = _pair_derivative_denominator(kex, p_i, p_j)
    return -kex * p_i / (denominator * denominator)


def _population_settings(states: str) -> dict[str, ParamLocalSetting]:
    return {
        f"p{state}": ParamLocalSetting(
            name_setting=NameSetting(f"p{state}", "", TPL),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        )
        for state in states[1:]
    }


def _population_complement_setting(states: str) -> ParamLocalSetting:
    arguments = ", ".join(f"{{p{state}}}" for state in states[1:])
    minimum, maximum = (-math.inf, math.inf) if len(states) == 3 else (0.0, 1.0)
    return ParamLocalSetting(
        name_setting=NameSetting("pa", "", TPL),
        min=minimum,
        max=maximum,
        expr=f"population_complement({arguments})['pa']",
    )


def _kex_settings(edges: tuple[Edge, ...]) -> dict[str, ParamLocalSetting]:
    return {
        f"kex_{left}{right}": ParamLocalSetting(
            name_setting=NameSetting(f"kex_{left}{right}", "", TPL),
            value=200.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        )
        for left, right in edges
    }


def _directional_rate_settings(
    edges: tuple[Edge, ...],
) -> dict[str, ParamLocalSetting]:
    settings: dict[str, ParamLocalSetting] = {}
    for left, right in edges:
        arguments = f"{{kex_{left}{right}}}, {{p{left}}}, {{p{right}}}"
        settings[f"k{left}{right}"] = ParamLocalSetting(
            name_setting=NameSetting(f"k{left}{right}", "", TPL),
            expr=f"pair_rates({arguments})['forward']",
        )
        settings[f"k{right}{left}"] = ParamLocalSetting(
            name_setting=NameSetting(f"k{right}{left}", "", TPL),
            expr=f"pair_rates({arguments})['reverse']",
        )
    return settings


def _make_settings(topology: Topology) -> dict[str, ParamLocalSetting]:
    """Build populations, structural KEX values, and directional rate pairs."""
    return {
        **_population_settings(topology.states),
        **_kex_settings(topology.edges),
        "pa": _population_complement_setting(topology.states),
        **_directional_rate_settings(topology.edges),
    }


def _settings_for(model_name: str) -> dict[str, ParamLocalSetting]:
    return _make_settings(_TOPOLOGY_BY_NAME[model_name])


def make_settings_3st(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("3st")


def make_settings_3st_linear(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("3st_linear")


def make_settings_3st_fork(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("3st_fork")


def make_settings_4st(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("4st")


def make_settings_4st_linear(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("4st_linear")


def make_settings_4st_fork(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("4st_fork")


def make_settings_5st(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("5st")


def make_settings_5st_linear(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("5st_linear")


def make_settings_5st_fork(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("5st_fork")


def make_settings_6st(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("6st")


def make_settings_6st_linear(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("6st_linear")


def make_settings_6st_fork(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return _settings_for("6st_fork")


def _linearizations(state_count: int) -> tuple[FunctionLinearization, ...]:
    return (
        AnalyticFunctionLinearization(
            "population_complement",
            "pa",
            "generic-nstate-population-complement-partials-v1",
            (_population_partial,) * (state_count - 1),
        ),
        AnalyticFunctionLinearization(
            "pair_rates",
            "forward",
            "generic-nstate-pair-rate-partials-v1",
            (_forward_d_kex, _forward_d_pi, _forward_d_pj),
        ),
        AnalyticFunctionLinearization(
            "pair_rates",
            "reverse",
            "generic-nstate-pair-rate-partials-v1",
            (_reverse_d_kex, _reverse_d_pi, _reverse_d_pj),
        ),
    )


def register() -> None:
    user_functions = {
        "population_complement": calculate_population_complement,
        "pair_rates": calculate_pair_rates,
    }
    makers = {
        "3st": make_settings_3st,
        "3st_triangle": make_settings_3st,
        "3st_linear": make_settings_3st_linear,
        "3st_fork": make_settings_3st_fork,
        "4st": make_settings_4st,
        "4st_linear": make_settings_4st_linear,
        "4st_fork": make_settings_4st_fork,
        "5st": make_settings_5st,
        "5st_linear": make_settings_5st_linear,
        "5st_fork": make_settings_5st_fork,
        "6st": make_settings_6st,
        "6st_linear": make_settings_6st_linear,
        "6st_fork": make_settings_6st_fork,
    }
    for topology in TOPOLOGIES:
        model_name = topology.model_name
        model_factory.register(name=model_name, setting_maker=makers[model_name])
        user_function_registry.register(model_name, user_functions)
        function_linearization_registry.register(
            model_name,
            _linearizations(len(topology.states)),
        )

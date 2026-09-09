from __future__ import annotations

from collections.abc import Callable
from functools import lru_cache
from typing import Literal

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.models.kinetic._eyring import (
    EYRING_RATE_PARTIALS,
    ThermodynamicCoordinate,
    calculate_directional_rate,
    calculate_rate_component,
    temperature_to_kelvin,
    thermodynamic_population_partials,
    thermodynamic_populations,
)
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import (
    AnalyticFunctionLinearization,
    FunctionLinearization,
    function_linearization_registry,
    user_function_registry,
)

NAME = "3st_eyring"
LINEAR_NAME = "3st_eyring_linear"
FORK_NAME = "3st_eyring_fork"

PL = ("p_total", "l_total")
TPL = ("temperature", "p_total", "l_total")

Edge = Literal["ab", "ac", "bc"]
LINEAR_EDGES: tuple[Edge, ...] = ("ab", "bc")
FORK_EDGES: tuple[Edge, ...] = ("ab", "ac")
REFERENCE_STATE = ThermodynamicCoordinate(enthalpy=0.0, entropy=0.0)


def _calculate_kij_3st_eyring(
    dh_b: float,
    ds_b: float,
    dh_c: float,
    ds_c: float,
    transition_states: dict[Edge, ThermodynamicCoordinate],
    temperature: float,
) -> dict[str, float]:
    states = {
        "a": REFERENCE_STATE,
        "b": ThermodynamicCoordinate(dh_b, ds_b),
        "c": ThermodynamicCoordinate(dh_c, ds_c),
    }
    rates: dict[str, float] = {}
    for edge, transition_state in transition_states.items():
        for initial, final in (edge, edge[::-1]):
            rates[f"k{initial}{final}"] = calculate_directional_rate(
                states[initial],
                transition_state,
                temperature,
            )
    return rates


@lru_cache(maxsize=100)
def calculate_kij_3st_eyring_linear(
    dh_b: float,
    ds_b: float,
    dh_c: float,
    ds_c: float,
    dh_ab: float,
    ds_ab: float,
    dh_bc: float,
    ds_bc: float,
    temperature: float,
) -> dict[str, float]:
    return _calculate_kij_3st_eyring(
        dh_b=dh_b,
        ds_b=ds_b,
        dh_c=dh_c,
        ds_c=ds_c,
        transition_states={
            "ab": ThermodynamicCoordinate(dh_ab, ds_ab),
            "bc": ThermodynamicCoordinate(dh_bc, ds_bc),
        },
        temperature=temperature,
    )


# The historical public model name has always represented the linear topology.
calculate_kij_3st_eyring = calculate_kij_3st_eyring_linear


@lru_cache(maxsize=100)
def calculate_kij_3st_eyring_fork(
    dh_b: float,
    ds_b: float,
    dh_c: float,
    ds_c: float,
    dh_ab: float,
    ds_ab: float,
    dh_ac: float,
    ds_ac: float,
    temperature: float,
) -> dict[str, float]:
    return _calculate_kij_3st_eyring(
        dh_b=dh_b,
        ds_b=ds_b,
        dh_c=dh_c,
        ds_c=ds_c,
        transition_states={
            "ab": ThermodynamicCoordinate(dh_ab, ds_ab),
            "ac": ThermodynamicCoordinate(dh_ac, ds_ac),
        },
        temperature=temperature,
    )


@lru_cache(maxsize=100)
def calculate_populations_3st_eyring(
    dh_b: float,
    ds_b: float,
    dh_c: float,
    ds_c: float,
    temperature: float,
) -> dict[str, float]:
    populations = thermodynamic_populations(
        {
            "a": REFERENCE_STATE,
            "b": ThermodynamicCoordinate(dh_b, ds_b),
            "c": ThermodynamicCoordinate(dh_c, ds_c),
        },
        temperature,
    )
    return {f"p{state}": population for state, population in populations.items()}


def _thermodynamic_settings(edges: tuple[Edge, ...]) -> dict[str, ParamLocalSetting]:
    settings = {
        "dh_b": ParamLocalSetting(
            name_setting=NameSetting("dh_b", "", PL),
            value=8e3,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        ),
        "dh_c": ParamLocalSetting(
            name_setting=NameSetting("dh_c", "", PL),
            value=8e3,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        ),
        "ds_b": ParamLocalSetting(
            name_setting=NameSetting("ds_b", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
        ),
        "ds_c": ParamLocalSetting(
            name_setting=NameSetting("ds_c", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
        ),
    }
    for edge in edges:
        settings[f"dh_{edge}"] = ParamLocalSetting(
            name_setting=NameSetting(f"dh_{edge}", "", PL),
            value=6.5e4,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        )
        settings[f"ds_{edge}"] = ParamLocalSetting(
            name_setting=NameSetting(f"ds_{edge}", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
        )
    return settings


def _rate_settings(
    edges: tuple[Edge, ...],
    temperature: float,
) -> dict[str, ParamLocalSetting]:
    state_arguments = {
        "a": ("0.0", "0.0"),
        "b": ("{dh_b}", "{ds_b}"),
        "c": ("{dh_c}", "{ds_c}"),
    }
    settings: dict[str, ParamLocalSetting] = {}
    for edge in edges:
        for initial, final in (edge, edge[::-1]):
            state_enthalpy, state_entropy = state_arguments[initial]
            settings[f"k{initial}{final}"] = ParamLocalSetting(
                name_setting=NameSetting(f"k{initial}{final}", "", TPL),
                min=0.0,
                expr=(
                    f"eyring_rate({state_enthalpy},{state_entropy},"
                    f"{{dh_{edge}}},{{ds_{edge}}},{temperature})['rate']"
                ),
            )
    return settings


def _population_settings(temperature: float) -> dict[str, ParamLocalSetting]:
    call = f"pop_3st_eyring({{dh_b}},{{ds_b}},{{dh_c}},{{ds_c}},{temperature})"
    return {
        f"p{state}": ParamLocalSetting(
            name_setting=NameSetting(f"p{state}", "", TPL),
            min=0.0,
            max=1.0,
            expr=f"{call}['p{state}']",
        )
        for state in "abc"
    }


def _make_settings_3st_eyring(
    conditions: Conditions,
    edges: tuple[Edge, ...],
) -> dict[str, ParamLocalSetting]:
    celsius = conditions.temperature
    if celsius is None:
        msg = "The 'temperature' is None"
        raise ValueError(msg)
    temperature_to_kelvin(celsius)
    return {
        **_thermodynamic_settings(edges),
        **_rate_settings(edges, celsius),
        **_population_settings(celsius),
    }


def make_settings_3st_eyring_linear(
    conditions: Conditions,
) -> dict[str, ParamLocalSetting]:
    return _make_settings_3st_eyring(conditions, LINEAR_EDGES)


# Compatibility name: 3st_eyring retains the historical linear topology.
make_settings_3st_eyring = make_settings_3st_eyring_linear


def make_settings_3st_eyring_fork(
    conditions: Conditions,
) -> dict[str, ParamLocalSetting]:
    return _make_settings_3st_eyring(
        conditions,
        FORK_EDGES,
    )


def _user_functions(
    rate_function: Callable[..., dict[str, float]],
    function_name: str,
) -> dict[str, object]:
    return {
        "eyring_rate": calculate_rate_component,
        function_name: rate_function,
        "pop_3st": pop_3st,
        "pop_3st_eyring": calculate_populations_3st_eyring,
    }


def _linearizations() -> tuple[FunctionLinearization, ...]:
    components = ("pa", "pb", "pc")
    population_partials = thermodynamic_population_partials(components)
    return (
        AnalyticFunctionLinearization(
            "eyring_rate",
            "rate",
            "eyring-directional-rate-partials-v1",
            EYRING_RATE_PARTIALS,
        ),
        *(
            AnalyticFunctionLinearization(
                "pop_3st_eyring",
                component,
                "eyring-thermodynamic-population-partials-v1",
                population_partials[component],
            )
            for component in components
        ),
    )


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_3st_eyring)
    model_factory.register(
        name=LINEAR_NAME, setting_maker=make_settings_3st_eyring_linear
    )
    model_factory.register(name=FORK_NAME, setting_maker=make_settings_3st_eyring_fork)
    linear_functions = _user_functions(
        calculate_kij_3st_eyring_linear,
        "kij_3st_eyring",
    )
    user_function_registry.register(name=NAME, user_functions=linear_functions)
    user_function_registry.register(name=LINEAR_NAME, user_functions=linear_functions)
    function_linearization_registry.register(NAME, _linearizations())
    function_linearization_registry.register(LINEAR_NAME, _linearizations())
    user_function_registry.register(
        name=FORK_NAME,
        user_functions=_user_functions(
            calculate_kij_3st_eyring_fork,
            "kij_3st_eyring_fork",
        ),
    )
    function_linearization_registry.register(FORK_NAME, _linearizations())

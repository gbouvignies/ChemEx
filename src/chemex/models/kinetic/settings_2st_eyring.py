from __future__ import annotations

from functools import lru_cache

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_2st
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
    function_linearization_registry,
    user_function_registry,
)

NAME = "2st_eyring"

PL = ("p_total", "l_total")
TPL = ("temperature", "p_total", "l_total")

REFERENCE_STATE = ThermodynamicCoordinate(enthalpy=0.0, entropy=0.0)


@lru_cache(maxsize=100)
def calculate_kij_2st_eyring(
    dh_b: float,
    ds_b: float,
    dh_ab: float,
    ds_ab: float,
    temperature: float,
) -> dict[str, float]:
    state_b = ThermodynamicCoordinate(enthalpy=dh_b, entropy=ds_b)
    transition_ab = ThermodynamicCoordinate(enthalpy=dh_ab, entropy=ds_ab)
    return {
        "kab": calculate_directional_rate(REFERENCE_STATE, transition_ab, temperature),
        "kba": calculate_directional_rate(state_b, transition_ab, temperature),
    }


@lru_cache(maxsize=100)
def calculate_populations_2st_eyring(
    dh_b: float,
    ds_b: float,
    temperature: float,
) -> dict[str, float]:
    populations = thermodynamic_populations(
        {
            "a": REFERENCE_STATE,
            "b": ThermodynamicCoordinate(enthalpy=dh_b, entropy=ds_b),
        },
        temperature,
    )
    return {f"p{state}": population for state, population in populations.items()}


def make_settings_2st_eyring(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    celsius = conditions.temperature
    if celsius is None:
        msg = "The 'temperature' is None"
        raise ValueError(msg)
    temperature_to_kelvin(celsius)
    return {
        "dh_b": ParamLocalSetting(
            name_setting=NameSetting("dh_b", "", PL),
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
            vary=True,
        ),
        "dh_ab": ParamLocalSetting(
            name_setting=NameSetting("dh_ab", "", PL),
            value=6.5e4,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        ),
        "ds_ab": ParamLocalSetting(
            name_setting=NameSetting("ds_ab", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=True,
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            min=0.0,
            expr=f"eyring_rate(0.0,0.0,{{dh_ab}},{{ds_ab}},{celsius})['rate']",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            min=0.0,
            expr=(
                f"eyring_rate({{dh_b}},{{ds_b}},{{dh_ab}},{{ds_ab}},{celsius})['rate']"
            ),
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            min=0.0,
            max=1.0,
            expr=f"pop_2st_eyring({{dh_b}},{{ds_b}},{celsius})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            min=0.0,
            max=1.0,
            expr=f"pop_2st_eyring({{dh_b}},{{ds_b}},{celsius})['pb']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_2st_eyring)
    user_functions = {
        "eyring_rate": calculate_rate_component,
        "kij_2st_eyring": calculate_kij_2st_eyring,
        "pop_2st": pop_2st,
        "pop_2st_eyring": calculate_populations_2st_eyring,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)
    population_partials = thermodynamic_population_partials(("pa", "pb"))
    function_linearization_registry.register(
        NAME,
        (
            AnalyticFunctionLinearization(
                "eyring_rate",
                "rate",
                "eyring-directional-rate-partials-v1",
                EYRING_RATE_PARTIALS,
            ),
            *(
                AnalyticFunctionLinearization(
                    "pop_2st_eyring",
                    component,
                    "eyring-thermodynamic-population-partials-v1",
                    population_partials[component],
                )
                for component in ("pa", "pb")
            ),
        ),
    )

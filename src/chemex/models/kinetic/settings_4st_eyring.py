"""4-state Eyring model for chemical exchange kinetics.

This module implements a 4-state kinetic model using Eyring transition state
theory to calculate temperature-dependent exchange rate constants from
thermodynamic parameters.
"""

from __future__ import annotations

from functools import lru_cache

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_4st
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

NAME = "4st_eyring"

PL = ("p_total", "l_total")
TPL = ("temperature", "p_total", "l_total")

REFERENCE_STATE = ThermodynamicCoordinate(enthalpy=0.0, entropy=0.0)
EDGES = ("ab", "ac", "ad", "bc", "bd", "cd")


@lru_cache(maxsize=100)
def calculate_kij_4st_eyring(
    dh_b: float,
    ds_b: float,
    dh_c: float,
    ds_c: float,
    dh_d: float,
    ds_d: float,
    dh_ab: float,
    ds_ab: float,
    dh_ac: float,
    ds_ac: float,
    dh_ad: float,
    ds_ad: float,
    dh_bc: float,
    ds_bc: float,
    dh_bd: float,
    ds_bd: float,
    dh_cd: float,
    ds_cd: float,
    temperature: float,
) -> dict[str, float]:
    """Calculate exchange rate constants using Eyring transition state theory.

    This function computes temperature-dependent rate constants for a 4-state
    exchange system using thermodynamic parameters (enthalpy and entropy changes).
    State A is used as the reference state with ΔH_A = ΔS_A = 0.

    Parameters
    ----------
    dh_b, dh_c, dh_d : float
        Enthalpy differences (J/mol) of states B, C, D relative to state A.
        Positive values indicate states higher in energy than A.
    ds_b, ds_c, ds_d : float
        Entropy differences (J/mol/K) of states B, C, D relative to state A.
        Positive values indicate states with higher entropy than A.
    dh_ab, dh_ac, dh_ad, dh_bc, dh_bd, dh_cd : float
        Activation enthalpies (J/mol) for transitions between states.
        These represent the enthalpy of the transition state relative to state A.
    ds_ab, ds_ac, ds_ad, ds_bc, ds_bd, ds_cd : float
        Activation entropies (J/mol/K) for transitions between states.
        These represent the entropy of the transition state relative to state A.
    temperature : float
        Temperature in Celsius.

    Returns
    -------
    dict[str, float]
        Dictionary containing rate constants (s⁻¹) for all state transitions.
        Keys are formatted as 'kij' where i and j are states ('a', 'b', 'c', 'd').

    Examples
    --------
    >>> rates = calculate_kij_4st_eyring(
    ...     dh_b=8000, ds_b=0, dh_c=12000, ds_c=0, dh_d=15000, ds_d=0,
    ...     dh_ab=75000, ds_ab=0, dh_ac=80000, ds_ac=0, dh_ad=85000, ds_ad=0,
    ...     dh_bc=70000, ds_bc=0, dh_bd=77000, ds_bd=0, dh_cd=72000, ds_cd=0,
    ...     temperature=25.0
    ... )
    >>> print(f"k_ab = {rates['kab']:.2e} s⁻¹")

    """
    states = {
        "a": REFERENCE_STATE,
        "b": ThermodynamicCoordinate(dh_b, ds_b),
        "c": ThermodynamicCoordinate(dh_c, ds_c),
        "d": ThermodynamicCoordinate(dh_d, ds_d),
    }
    transition_states = {
        "ab": ThermodynamicCoordinate(dh_ab, ds_ab),
        "ac": ThermodynamicCoordinate(dh_ac, ds_ac),
        "ad": ThermodynamicCoordinate(dh_ad, ds_ad),
        "bc": ThermodynamicCoordinate(dh_bc, ds_bc),
        "bd": ThermodynamicCoordinate(dh_bd, ds_bd),
        "cd": ThermodynamicCoordinate(dh_cd, ds_cd),
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
def calculate_populations_4st_eyring(
    dh_b: float,
    ds_b: float,
    dh_c: float,
    ds_c: float,
    dh_d: float,
    ds_d: float,
    temperature: float,
) -> dict[str, float]:
    populations = thermodynamic_populations(
        {
            "a": REFERENCE_STATE,
            "b": ThermodynamicCoordinate(dh_b, ds_b),
            "c": ThermodynamicCoordinate(dh_c, ds_c),
            "d": ThermodynamicCoordinate(dh_d, ds_d),
        },
        temperature,
    )
    return {f"p{state}": population for state, population in populations.items()}


def create_kij_4st_eyring_settings(temperature: float) -> dict[str, ParamLocalSetting]:
    state_arguments = {
        "a": ("0.0", "0.0"),
        "b": ("{dh_b}", "{ds_b}"),
        "c": ("{dh_c}", "{ds_c}"),
        "d": ("{dh_d}", "{ds_d}"),
    }
    settings: dict[str, ParamLocalSetting] = {}
    for edge in EDGES:
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


def create_pop_4st_eyring_settings(
    temperature: float,
) -> dict[str, ParamLocalSetting]:
    call = (
        "pop_4st_eyring("
        f"{{dh_b}},{{ds_b}},{{dh_c}},{{ds_c}},{{dh_d}},{{ds_d}},{temperature})"
    )
    return {
        f"p{state}": ParamLocalSetting(
            name_setting=NameSetting(f"p{state}", "", TPL),
            min=0.0,
            max=1.0,
            expr=f"{call}['p{state}']",
        )
        for state in "abcd"
    }


def make_settings_4st_eyring(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    """Create parameter settings for 4-state Eyring kinetic model.

    Parameters
    ----------
    conditions : Conditions
        Experimental conditions including temperature

    Returns
    -------
    dict[str, ParamLocalSetting]
        Dictionary of parameter settings for the model

    Raises
    ------
    ValueError
        If temperature is None

    """
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
        "dh_c": ParamLocalSetting(
            name_setting=NameSetting("dh_c", "", PL),
            value=8e3,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        ),
        "dh_d": ParamLocalSetting(
            name_setting=NameSetting("dh_d", "", PL),
            value=8e3,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        ),
        "dh_ab": ParamLocalSetting(
            name_setting=NameSetting("dh_ab", "", PL),
            value=7.5e4,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        ),
        "dh_ac": ParamLocalSetting(
            name_setting=NameSetting("dh_ac", "", PL),
            value=7.5e4,
            min=-2.0e5,
            max=2.0e5,
            vary=False,
        ),
        "dh_ad": ParamLocalSetting(
            name_setting=NameSetting("dh_ad", "", PL),
            value=7.5e4,
            min=-2.0e5,
            max=2.0e5,
            vary=False,
        ),
        "dh_bc": ParamLocalSetting(
            name_setting=NameSetting("dh_bc", "", PL),
            value=7.5e4,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        ),
        "dh_bd": ParamLocalSetting(
            name_setting=NameSetting("dh_bd", "", PL),
            value=7.5e4,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        ),
        "dh_cd": ParamLocalSetting(
            name_setting=NameSetting("dh_cd", "", PL),
            value=7.5e4,
            min=-2.0e5,
            max=2.0e5,
            vary=True,
        ),
        "ds_b": ParamLocalSetting(
            name_setting=NameSetting("ds_b", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=False,
        ),
        "ds_c": ParamLocalSetting(
            name_setting=NameSetting("ds_c", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=False,
        ),
        "ds_d": ParamLocalSetting(
            name_setting=NameSetting("ds_d", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=False,
        ),
        "ds_ab": ParamLocalSetting(
            name_setting=NameSetting("ds_ab", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=False,
        ),
        "ds_ac": ParamLocalSetting(
            name_setting=NameSetting("ds_ac", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=False,
        ),
        "ds_ad": ParamLocalSetting(
            name_setting=NameSetting("ds_ad", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=False,
        ),
        "ds_bc": ParamLocalSetting(
            name_setting=NameSetting("ds_bc", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=False,
        ),
        "ds_bd": ParamLocalSetting(
            name_setting=NameSetting("ds_bd", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=False,
        ),
        "ds_cd": ParamLocalSetting(
            name_setting=NameSetting("ds_cd", "", PL),
            value=0.0,
            min=-5.0e2,
            max=5.0e2,
            vary=False,
        ),
        **create_kij_4st_eyring_settings(celsius),
        **create_pop_4st_eyring_settings(celsius),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_4st_eyring)
    user_functions = {
        "eyring_rate": calculate_rate_component,
        "kij_4st_eyring": calculate_kij_4st_eyring,
        "pop_4st": pop_4st,
        "pop_4st_eyring": calculate_populations_4st_eyring,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)
    components = ("pa", "pb", "pc", "pd")
    population_partials = thermodynamic_population_partials(components)
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
                    "pop_4st_eyring",
                    component,
                    "eyring-thermodynamic-population-partials-v1",
                    population_partials[component],
                )
                for component in components
            ),
        ),
    )

"""4-state Eyring model for chemical exchange kinetics.

This module implements a 4-state kinetic model using Eyring transition state
theory to calculate temperature-dependent exchange rate constants from
thermodynamic parameters.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import permutations

import numpy as np
from scipy.constants import R, h, k

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_4st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import ArrayFloat

NAME = "4st_eyring"

PL = ("p_total", "l_total")
TPL = ("temperature", "p_total", "l_total")

# Physical constants
CELSIUS_TO_KELVIN = 273.15  # Temperature conversion from Celsius to Kelvin
MAX_RATE_CONSTANT = 1e16  # Maximum rate constant (s⁻¹) for numerical stability
MIN_TEMPERATURE = -20.0  # Minimum temperature (°C) for validation
MAX_TEMPERATURE = 100.0  # Maximum temperature (°C) for validation


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

    Notes
    -----
    The rate constants are calculated using Eyring equation:
    k_ij = (k_B*T/h) * exp(-ΔG‡_ij / RT)

    Where ΔG‡_ij is the activation free energy for transition from state i to j:
    ΔG‡_ij = ΔH‡_ij - T*ΔS‡_ij

    Rate constants are clipped to [0, 1e16] s⁻¹ for numerical stability.


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
    # Validate temperature range
    if not MIN_TEMPERATURE <= temperature <= MAX_TEMPERATURE:
        msg = (
            f"Temperature {temperature}°C is outside the valid range "
            f"[{MIN_TEMPERATURE}, {MAX_TEMPERATURE}]°C"
        )
        raise ValueError(msg)

    kelvin = temperature + CELSIUS_TO_KELVIN
    kbt_h = k * kelvin / h
    rt = R * kelvin

    # Reference state A has zero enthalpy and entropy
    dh_a = ds_a = 0.0

    # Calculate activation free energies for all transitions
    ddg_ij = np.array(
        (
            dh_ab - dh_a - kelvin * (ds_ab - ds_a),
            dh_ac - dh_a - kelvin * (ds_ac - ds_a),
            dh_ad - dh_a - kelvin * (ds_ad - ds_a),
            dh_ab - dh_b - kelvin * (ds_ab - ds_b),
            dh_bc - dh_b - kelvin * (ds_bc - ds_b),
            dh_bd - dh_b - kelvin * (ds_bd - ds_b),
            dh_ac - dh_c - kelvin * (ds_ac - ds_c),
            dh_bc - dh_c - kelvin * (ds_bc - ds_c),
            dh_cd - dh_c - kelvin * (ds_cd - ds_c),
            dh_ad - dh_d - kelvin * (ds_ad - ds_d),
            dh_bd - dh_d - kelvin * (ds_bd - ds_d),
            dh_cd - dh_d - kelvin * (ds_cd - ds_d),
        ),
    )

    # Apply Eyring equation
    kij_values: ArrayFloat = kbt_h * np.exp(-ddg_ij / rt)

    # Clip values for numerical stability
    kij_values = np.clip(kij_values, 0.0, MAX_RATE_CONSTANT)

    # Generate rate constant names in same order as permutations
    kij_names = (f"k{i}{j}" for i, j in permutations("abcd", 2))

    return dict(zip(kij_names, kij_values, strict=True))


def create_kij_4st_eyring_settings(temperature: float) -> dict[str, ParamLocalSetting]:
    return {
        f"k{i}{j}": ParamLocalSetting(
            name_setting=NameSetting(f"k{i}{j}", "", TPL),
            min=0.0,
            expr=(
                f"kij_4st_eyring("
                f"{{dh_b}}, {{ds_b}}, "
                f"{{dh_c}}, {{ds_c}}, "
                f"{{dh_d}}, {{ds_d}},"
                f"{{dh_ab}}, {{ds_ab}}, "
                f"{{dh_ac}}, {{ds_ac}}, "
                f"{{dh_ad}}, {{ds_ad}}, "
                f"{{dh_bc}}, {{ds_bc}}, "
                f"{{dh_bd}}, {{ds_bd}}, "
                f"{{dh_cd}}, {{ds_cd}},"
                f" {temperature}"
                f")['k{i}{j}']"
            ),
        )
        for i, j in permutations("abcd", 2)
    }


def create_pop_4st_eyring_settings() -> dict[str, ParamLocalSetting]:
    return {
        f"p{state}": ParamLocalSetting(
            name_setting=NameSetting(f"p{state}", "", TPL),
            min=0.0,
            max=1.0,
            expr=f"pop_4st("
            f"{{kab}},{{kba}},"
            f"{{kac}},{{kca}},"
            f"{{kad}},{{kda}},"
            f"{{kbc}},{{kcb}},"
            f"{{kbd}},{{kdb}},"
            f"{{kcd}},{{kdc}})['p{state}']",
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
    # Early validation of temperature range (same policy as in calculator)
    if not MIN_TEMPERATURE <= celsius <= MAX_TEMPERATURE:
        msg = (
            f"Temperature {celsius}°C is outside the valid range "
            f"[{MIN_TEMPERATURE}, {MAX_TEMPERATURE}]°C"
        )
        raise ValueError(msg)
    return {
        "dh_b": ParamLocalSetting(
            name_setting=NameSetting("dh_b", "", PL),
            value=8e3,
            vary=True,
        ),
        "dh_c": ParamLocalSetting(
            name_setting=NameSetting("dh_c", "", PL),
            value=8e3,
            vary=True,
        ),
        "dh_d": ParamLocalSetting(
            name_setting=NameSetting("dh_d", "", PL),
            value=8e3,
            vary=True,
        ),
        "dh_ab": ParamLocalSetting(
            name_setting=NameSetting("dh_ab", "", PL),
            value=7.5e4,
            vary=True,
        ),
        "dh_ac": ParamLocalSetting(
            name_setting=NameSetting("dh_ac", "", PL),
            value=7.5e4,
            vary=False,
        ),
        "dh_ad": ParamLocalSetting(
            name_setting=NameSetting("dh_ad", "", PL),
            value=7.5e4,
            vary=False,
        ),
        "dh_bc": ParamLocalSetting(
            name_setting=NameSetting("dh_bc", "", PL),
            value=7.5e4,
            vary=True,
        ),
        "dh_bd": ParamLocalSetting(
            name_setting=NameSetting("dh_bd", "", PL),
            value=7.5e4,
            vary=True,
        ),
        "dh_cd": ParamLocalSetting(
            name_setting=NameSetting("dh_cd", "", PL),
            value=7.5e4,
            vary=True,
        ),
        "ds_b": ParamLocalSetting(
            name_setting=NameSetting("ds_b", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_c": ParamLocalSetting(
            name_setting=NameSetting("ds_c", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_d": ParamLocalSetting(
            name_setting=NameSetting("ds_d", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_ab": ParamLocalSetting(
            name_setting=NameSetting("ds_ab", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_ac": ParamLocalSetting(
            name_setting=NameSetting("ds_ac", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_ad": ParamLocalSetting(
            name_setting=NameSetting("ds_ad", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_bc": ParamLocalSetting(
            name_setting=NameSetting("ds_bc", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_bd": ParamLocalSetting(
            name_setting=NameSetting("ds_bd", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_cd": ParamLocalSetting(
            name_setting=NameSetting("ds_cd", "", PL),
            value=0.0,
            vary=False,
        ),
        **create_kij_4st_eyring_settings(celsius),
        **create_pop_4st_eyring_settings(),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_4st_eyring)
    user_functions = {
        "kij_4st_eyring": calculate_kij_4st_eyring,
        "pop_4st": pop_4st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

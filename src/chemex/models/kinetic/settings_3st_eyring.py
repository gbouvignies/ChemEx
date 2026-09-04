from __future__ import annotations

from collections.abc import Callable
from functools import lru_cache
from typing import Literal

import numpy as np
from scipy import constants

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry

NAME = "3st_eyring"
LINEAR_NAME = "3st_eyring_linear"
FORK_NAME = "3st_eyring_fork"

PL = ("p_total", "l_total")
TPL = ("temperature", "p_total", "l_total")

MAX_RATE_CONSTANT = 1e16

Edge = Literal["ab", "ac", "bc"]
LINEAR_EDGES: tuple[Edge, ...] = ("ab", "bc")
FORK_EDGES: tuple[Edge, ...] = ("ab", "ac")


def _calculate_kij_3st_eyring(
    dh_b: float,
    ds_b: float,
    dh_c: float,
    ds_c: float,
    transition_terms: tuple[tuple[Edge, float, float], ...],
    temperature: float,
) -> dict[str, float]:
    kelvin = temperature + constants.zero_Celsius
    kbt_h = constants.k * kelvin / constants.h
    rt = constants.R * kelvin
    state_terms = {
        "a": (0.0, 0.0),
        "b": (dh_b, ds_b),
        "c": (dh_c, ds_c),
    }
    rates: dict[str, float] = {}
    for edge, dh_transition, ds_transition in transition_terms:
        for initial, final in (edge, edge[::-1]):
            dh_initial, ds_initial = state_terms[initial]
            activation_free_energy = (
                dh_transition - dh_initial - kelvin * (ds_transition - ds_initial)
            )
            rate = kbt_h * np.exp(-activation_free_energy / rt)
            rates[f"k{initial}{final}"] = float(
                np.clip(rate, 0.0, MAX_RATE_CONSTANT),
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
        dh_b,
        ds_b,
        dh_c,
        ds_c,
        (("ab", dh_ab, ds_ab), ("bc", dh_bc, ds_bc)),
        temperature,
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
        dh_b,
        ds_b,
        dh_c,
        ds_c,
        (("ab", dh_ab, ds_ab), ("ac", dh_ac, ds_ac)),
        temperature,
    )


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
    function_name: str,
    temperature: float,
) -> dict[str, ParamLocalSetting]:
    arguments = ["{dh_b}", "{ds_b}", "{dh_c}", "{ds_c}"]
    for edge in edges:
        arguments.extend((f"{{dh_{edge}}}", f"{{ds_{edge}}}"))
    arguments.append(str(temperature))
    call = f"{function_name}({','.join(arguments)})"
    return {
        f"k{initial}{final}": ParamLocalSetting(
            name_setting=NameSetting(f"k{initial}{final}", "", TPL),
            min=0.0,
            expr=f"{call}['k{initial}{final}']",
        )
        for edge in edges
        for initial, final in (edge, edge[::-1])
    }


def _population_settings(edges: tuple[Edge, ...]) -> dict[str, ParamLocalSetting]:
    active_rates = {
        f"k{initial}{final}" for edge in edges for initial, final in (edge, edge[::-1])
    }
    arguments = [
        f"{{{rate}}}" if rate in active_rates else "0.0"
        for rate in ("kab", "kba", "kac", "kca", "kbc", "kcb")
    ]
    call = f"pop_3st({','.join(arguments)})"
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
    function_name: str,
) -> dict[str, ParamLocalSetting]:
    celsius = conditions.temperature
    if celsius is None:
        msg = "The 'temperature' is None"
        raise ValueError(msg)
    return {
        **_thermodynamic_settings(edges),
        **_rate_settings(edges, function_name, celsius),
        **_population_settings(edges),
    }


def make_settings_3st_eyring_linear(
    conditions: Conditions,
) -> dict[str, ParamLocalSetting]:
    return _make_settings_3st_eyring(conditions, LINEAR_EDGES, "kij_3st_eyring")


# Compatibility name: 3st_eyring retains the historical linear topology.
make_settings_3st_eyring = make_settings_3st_eyring_linear


def make_settings_3st_eyring_fork(
    conditions: Conditions,
) -> dict[str, ParamLocalSetting]:
    return _make_settings_3st_eyring(
        conditions,
        FORK_EDGES,
        "kij_3st_eyring_fork",
    )


def _user_functions(
    rate_function: Callable[..., dict[str, float]],
    function_name: str,
) -> dict[str, object]:
    return {function_name: rate_function, "pop_3st": pop_3st}


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
    user_function_registry.register(
        name=FORK_NAME,
        user_functions=_user_functions(
            calculate_kij_3st_eyring_fork,
            "kij_3st_eyring_fork",
        ),
    )

from __future__ import annotations

from functools import lru_cache
from itertools import permutations

import numpy as np
from scipy.constants import constants

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import ArrayFloat

NAME = "3st_eyring"

PL = ("p_total", "l_total")
TPL = ("temperature", "p_total", "l_total")


@lru_cache(maxsize=100)
def calculate_kij_3st_eyring(
    dh_b: float,
    ds_b: float,
    dh_c: float,
    ds_c: float,
    dh_ab: float,
    ds_ab: float,
    dh_ac: float,
    ds_ac: float,
    dh_bc: float,
    ds_bc: float,
    temperature: float,
) -> dict[str, float]:
    kelvin = temperature + 273.15
    kbt_h = constants.k * kelvin / constants.h
    rt = constants.R * kelvin
    dh_a = ds_a = 0.0
    ddg_ij = np.array(
        (
            dh_ab - dh_a - kelvin * (ds_ab - ds_a),
            dh_ac - dh_a - kelvin * (ds_ac - ds_a),
            dh_ab - dh_b - kelvin * (ds_ab - ds_b),
            dh_bc - dh_b - kelvin * (ds_bc - ds_b),
            dh_ac - dh_c - kelvin * (ds_ac - ds_c),
            dh_bc - dh_c - kelvin * (ds_bc - ds_c),
        ),
    )
    kij_values: ArrayFloat = kbt_h * np.exp(-ddg_ij / rt)
    kij_values = np.clip(kij_values, 0.0, 1e16)
    kij_names = (f"k{i}{j}" for i, j in permutations("abc", 2))
    return dict(zip(kij_names, kij_values, strict=True))


def create_kij_3st_eyring_settings(temperature: float) -> dict[str, ParamLocalSetting]:
    return {
        f"k{i}{j}": ParamLocalSetting(
            name_setting=NameSetting(f"k{i}{j}", "", TPL),
            expr=(
                f"calculate_kij_3st_eyring({{dh_b}}, {{ds_b}}, {{dh_c}}, {{ds_c}}, "
                f"{{dh_ab}}, {{ds_ab}}, {{dh_ac}}, {{ds_ac}}, {{dh_bc}}, {{ds_bc}},"
                f" {temperature})['k{i}{j}']"
            ),
        )
        for i, j in permutations("abc", 2)
    }


def create_pop_3st_eyring_settings() -> dict[str, ParamLocalSetting]:
    return {
        f"p{state}": ParamLocalSetting(
            name_setting=NameSetting(f"p{state}", "", TPL),
            expr=f"pop_3st({{kab}},{{kba}},{{kac}},{{kca}},{{kbc}},{{kcb}})['p{state}']",
        )
        for state in "abc"
    }


def make_settings_3st_eyring(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    celsius = conditions.temperature
    if celsius is None:
        msg = "The 'temperature' is None"
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
        "dh_ab": ParamLocalSetting(
            name_setting=NameSetting("dh_ab", "", PL),
            value=6.5e4,
            vary=True,
        ),
        "dh_bc": ParamLocalSetting(
            name_setting=NameSetting("dh_bc", "", PL),
            value=6.5e4,
            vary=True,
        ),
        "dh_ac": ParamLocalSetting(
            name_setting=NameSetting("dh_ac", "", PL),
            value=1.0e10,
            vary=False,
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
        "ds_ab": ParamLocalSetting(
            name_setting=NameSetting("ds_ab", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_bc": ParamLocalSetting(
            name_setting=NameSetting("ds_bc", "", PL),
            value=0.0,
            vary=False,
        ),
        "ds_ac": ParamLocalSetting(
            name_setting=NameSetting("ds_ac", "", PL),
            value=0.0,
            vary=False,
        ),
        **create_kij_3st_eyring_settings(celsius),
        **create_pop_3st_eyring_settings(),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_3st_eyring)
    user_functions = {
        "calculate_kij_3st_eyring": calculate_kij_3st_eyring,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

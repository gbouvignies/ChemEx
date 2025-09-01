from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.constants import constants

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_2st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry

NAME = "2st_eyring"

PL = ("p_total", "l_total")
TPL = ("temperature", "p_total", "l_total")


@lru_cache(maxsize=100)
def calculate_kij_2st_eyring(
    dh_b: float,
    ds_b: float,
    dh_ab: float,
    ds_ab: float,
    temperature: float,
) -> dict[str, float]:
    kelvin = temperature + 273.15
    kbt_h = constants.k * kelvin / constants.h
    rt = constants.R * kelvin
    dh_a = ds_a = 0.0
    kab = kbt_h * np.exp(-(dh_ab - dh_a - kelvin * (ds_ab - ds_a)) / rt)
    kba = kbt_h * np.exp(-(dh_ab - dh_b - kelvin * (ds_ab - ds_b)) / rt)
    return {"kab": kab, "kba": kba}


def make_settings_2st_eyring(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    celsius = conditions.temperature
    return {
        "dh_b": ParamLocalSetting(
            name_setting=NameSetting("dh_b", "", PL),
            value=8e3,
            vary=True,
        ),
        "ds_b": ParamLocalSetting(
            name_setting=NameSetting("ds_b", "", PL),
            value=0.0,
            vary=True,
        ),
        "dh_ab": ParamLocalSetting(
            name_setting=NameSetting("dh_ab", "", PL),
            value=6.5e4,
            vary=True,
        ),
        "ds_ab": ParamLocalSetting(
            name_setting=NameSetting("ds_ab", "", PL),
            value=0.0,
            vary=True,
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            min=0.0,
            expr=f"kij_2st_eyring({{dh_b}},{{ds_b}},{{dh_ab}},{{ds_ab}},{celsius})['kab']",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            min=0.0,
            expr=f"kij_2st_eyring({{dh_b}},{{ds_b}},{{dh_ab}},{{ds_ab}},{celsius})['kba']",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            min=0.0,
            max=1.0,
            expr="pop_2st({kab},{kba})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            min=0.0,
            max=1.0,
            expr="pop_2st({kab},{kba})['pb']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_2st_eyring)
    user_functions = {"kij_2st_eyring": calculate_kij_2st_eyring, "pop_2st": pop_2st}
    user_function_registry.register(name=NAME, user_functions=user_functions)

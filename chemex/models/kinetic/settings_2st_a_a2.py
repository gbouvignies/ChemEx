from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import root

from chemex.models.constraints import pop_2st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting
from chemex.parameters.setting import ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry

if TYPE_CHECKING:
    from chemex.configuration.conditions import Conditions

NAME = "2st_a_a2"

TPL = ("temperature", "p_total", "l_total")


def calculate_concentrations(
    concentrations: np.ndarray, p_total: float, k1: float
) -> np.ndarray:
    p_monomer, p_dimer = concentrations
    return np.array(
        [
            p_monomer + 2.0 * p_dimer - p_total,
            k1 * p_monomer**2 - p_dimer,
        ]
    )


@lru_cache(maxsize=100)
def calculate_monomer_concentration(p_total: float, k1: float) -> float:
    concentrations_start = (p_total, 0.0)
    results = root(calculate_concentrations, concentrations_start, args=(p_total, k1))
    p_monomer, _p_dimer = results["x"]
    return p_monomer


def make_settings_2st_binding(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    p_total = conditions.p_total
    if p_total is None:
        raise ValueError("'p_total' must be specified to use the '2st_a_a2' model")
    return {
        "k_a2_a": ParamLocalSetting(
            name_setting=NameSetting("k_a2_a", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "k1": ParamLocalSetting(
            name_setting=NameSetting("k1", "", ("temperature",)),
            value=2.0,
            min=0.0,
            vary=True,
        ),
        "k_a_a2": ParamLocalSetting(
            name_setting=NameSetting("k_a_a2", "", ("temperature",)),
            min=0.0,
            expr="{k_a2_a} * {k1}",
        ),
        "p_monomer": ParamLocalSetting(
            name_setting=NameSetting("p_monomer", "", TPL),
            expr=f"calc_p_monomer({p_total}, {{k1}})",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="2.0 * {k_a_a2} * {p_monomer}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{k_a2_a}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr="pop_2st({kab}, {kba})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr="pop_2st({kab}, {kba})['pb']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_2st_binding)
    user_functions = {
        "calc_p_monomer": calculate_monomer_concentration,
        "pop_2st": pop_2st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

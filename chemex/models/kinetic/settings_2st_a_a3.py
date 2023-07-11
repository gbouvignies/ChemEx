from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import root

from chemex.models.constraints import pop_2st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry

if TYPE_CHECKING:
    from chemex.configuration.conditions import Conditions
    from chemex.typing import ArrayFloat

NAME = "2st_a_a3"

TPL = ("temperature", "p_total", "l_total")


def calculate_concentrations(
    concentrations: ArrayFloat, p_total: float, k1: float
) -> ArrayFloat:
    """Calculate the concentrations of monomer and trimer."""
    p_monomer, p_trimer = concentrations
    return np.array(
        [
            p_monomer + 3.0 * p_trimer - p_total,
            k1 * p_monomer**3 - p_trimer,
        ]
    )


@lru_cache(maxsize=100)
def calculate_monomer_concentration(p_total: float, k1: float) -> float:
    concentrations_start = (p_total, 0.0)
    results = root(calculate_concentrations, concentrations_start, args=(p_total, k1))
    p_monomer, _p_trimer = results["x"]
    return p_monomer


def make_settings_2st_binding(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    p_total = conditions.p_total
    if p_total is None:
        msg = "'p_total' must be specified to use the '2st_a_a2' model"
        raise ValueError(msg)
    return {
        "k1": ParamLocalSetting(
            name_setting=NameSetting("k1", "", ("temperature",)),
            value=2.0,
            min=0.0,
            vary=True,
        ),
        "k_a3_a": ParamLocalSetting(
            name_setting=NameSetting("k_a3_a", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "k_a_a3": ParamLocalSetting(
            name_setting=NameSetting("k_a_a3", "", ("temperature",)),
            min=0.0,
            expr="{k_a3_a} * {k1}",
        ),
        "p_monomer": ParamLocalSetting(
            name_setting=NameSetting("p_monomer", "", TPL),
            expr=f"calc_p_monomer({p_total}, {{k1}})",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="3.0 * {k_a_a3} * {p_monomer}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{k_a3_a}",
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

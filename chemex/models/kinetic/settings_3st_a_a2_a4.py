from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import root

from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting
from chemex.parameters.setting import ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry

if TYPE_CHECKING:
    from chemex.configuration.conditions import Conditions

NAME = "3st_a_a2_a4"

TPL = ("temperature", "p_total", "l_total")


def calculate_residuals(
    concentrations: np.ndarray,
    p_total: float,
    k1: float,
    k2: float,
) -> np.ndarray:
    p_monomer, p_dimer, p_tetramer = concentrations
    return np.array(
        [
            p_monomer + 2.0 * p_dimer + 4.0 * p_tetramer - p_total,
            k1 * p_monomer**2 - p_dimer,
            k2 * p_dimer**2 - p_tetramer,
        ]
    )


@lru_cache(maxsize=100)
def calculate_concentrations(p_total: float, k1: float, k2: float) -> dict[str, float]:
    concentrations_start = (p_total, 0.0, 0.0)
    results = root(calculate_residuals, concentrations_start, args=(p_total, k1, k2))
    return dict(zip(("p_monomer", "p_dimer", "p_tetramer", results["x"])))


def make_settings_3st_a_a2_a4(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    p_total = conditions.p_total
    if p_total is None:
        raise ValueError(f"'p_total' must be specified to use the '{NAME}' model")
    return {
        "k1": ParamLocalSetting(
            name_setting=NameSetting("k1", "", ("temperature",)),
            value=2.0,
            min=0.0,
            vary=True,
        ),
        "k2": ParamLocalSetting(
            name_setting=NameSetting("k2", "", ("temperature",)),
            value=2.0,
            min=0.0,
            vary=True,
        ),
        "k_a2_a": ParamLocalSetting(
            name_setting=NameSetting("k_a2_a", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "k_a4_a2": ParamLocalSetting(
            name_setting=NameSetting("k_a4_a2", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "k_a_a2": ParamLocalSetting(
            name_setting=NameSetting("k_a_a4", "", ("temperature",)),
            min=0.0,
            expr="{k_a2_a} * {k1}",
        ),
        "k_a2_a4": ParamLocalSetting(
            name_setting=NameSetting("k_a_a4", "", ("temperature",)),
            min=0.0,
            expr="{k_a4_a2} * {k2}",
        ),
        "p_monomer": ParamLocalSetting(
            name_setting=NameSetting("p_monomer", "", TPL),
            expr=f"calc_conc({p_total}, {{k1}}, {{k2}})['p_monomer']",
        ),
        "p_dimer": ParamLocalSetting(
            name_setting=NameSetting("p_dimer", "", TPL),
            expr=f"calc_conc({p_total}, {{k1}}, {{k2}})['p_dimer']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="2.0 * {k_a_a2} * {p_monomer}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{k_a2_a}",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="2.0 * {k_a2_a4} * {p_dimer}",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{k_a4_a2}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr="pop_3st({kab}, {kba}, 0.0, 0.0, {kbc}, {kcb})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr="pop_3st({kab}, {kba}, 0.0, 0.0, {kbc}, {kcb})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr="pop_3st({kab}, {kba}, 0.0, 0.0, {kbc}, {kcb})['pc']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_3st_a_a2_a4)
    user_functions = {
        "calc_conc": calculate_concentrations,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

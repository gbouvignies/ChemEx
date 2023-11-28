from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import root

from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry

if TYPE_CHECKING:
    from chemex.configuration.conditions import Conditions
    from chemex.typing import ArrayFloat

NAME = "3st_monomer_dimer_tetramer"

TP = ("temperature", "p_total")


def calculate_residuals(
    concentrations: ArrayFloat,
    p_total: float,
    kd1: float,
    kd2: float,
) -> ArrayFloat:
    monomer, dimer, tetramer = concentrations
    return np.array(
        [
            monomer + 2.0 * dimer + 4.0 * tetramer - p_total,
            kd1 * dimer - monomer**2,
            kd2 * tetramer - dimer**2,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    kd1: float,
    kd2: float,
) -> dict[str, float]:
    concentrations_start = (p_total, 0.0, 0.0)
    results = root(calculate_residuals, concentrations_start, args=(p_total, kd1, kd2))
    return {
        "monomer": results["x"][0],
        "dimer": results["x"][1],
        "tetramer": results["x"][2],
    }


def make_settings_3st_monomer_dimer_tetramer(
    conditions: Conditions,
) -> dict[str, ParamLocalSetting]:
    p_total = conditions.p_total
    if p_total is None:
        msg = f"'p_total' must be specified to use the '{NAME}' model"
        raise ValueError(msg)
    return {
        "kd1": ParamLocalSetting(
            name_setting=NameSetting("kd1", "", ("temperature",)),
            value=1e-6,
            min=0.0,
            vary=True,
        ),
        "kd2": ParamLocalSetting(
            name_setting=NameSetting("kd2", "", ("temperature",)),
            value=1e-6,
            min=0.0,
            vary=True,
        ),
        "koff1": ParamLocalSetting(
            name_setting=NameSetting("koff1", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "koff2": ParamLocalSetting(
            name_setting=NameSetting("koff2", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kon1": ParamLocalSetting(
            name_setting=NameSetting("kon1", "", ("temperature",)),
            min=0.0,
            expr="{koff1} / max({kd1}, 1e-32)",
        ),
        "kon2": ParamLocalSetting(
            name_setting=NameSetting("kon2", "", ("temperature",)),
            min=0.0,
            expr="{koff2} / max({kd2}, 1e-32)",
        ),
        "c_monomer": ParamLocalSetting(
            name_setting=NameSetting("monomer", "", TP),
            expr=f"concentrations({p_total}, {{kd1}}, {{kd2}})['monomer']",
        ),
        "c_dimer": ParamLocalSetting(
            name_setting=NameSetting("dimer", "", TP),
            expr=f"concentrations({p_total}, {{kd1}}, {{kd2}})['dimer']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TP),
            expr="2.0 * {kon1} * {monomer}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TP),
            expr="{koff1}",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TP),
            expr="2.0 * {kon2} * {dimer}",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TP),
            expr="{koff2}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TP),
            expr="pop_3st({kab}, {kba}, 0.0, 0.0, {kbc}, {kcb})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TP),
            expr="pop_3st({kab}, {kba}, 0.0, 0.0, {kbc}, {kcb})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TP),
            expr="pop_3st({kab}, {kba}, 0.0, 0.0, {kbc}, {kcb})['pc']",
        ),
    }


def register() -> None:
    model_factory.register(
        name=NAME,
        setting_maker=make_settings_3st_monomer_dimer_tetramer,
    )
    user_functions = {
        "concentrations": calculate_concentrations,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

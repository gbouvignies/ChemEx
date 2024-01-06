from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.optimize import root

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import ArrayFloat

NAME = "3st_monomer_dimer_trimer"

TP = ("temperature", "p_total")


def calculate_residuals(
    concentrations: ArrayFloat,
    p_total: float,
    kd1: float,
    kd2: float,
) -> ArrayFloat:
    monomer, dimer, trimer = concentrations
    return np.array(
        [
            monomer + 2.0 * dimer + 3.0 * trimer - p_total,
            kd1 * dimer - monomer**2,
            kd2 * trimer - monomer * dimer,
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
        "trimer": results["x"][2],
    }


def make_settings_3st_monomer_dimer_trimer(
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
            name_setting=NameSetting("c_monomer", "", TP),
            expr=f"concetrations({p_total}, {{kd1}}, {{kd2}})['monomer']",
        ),
        "c_dimer": ParamLocalSetting(
            name_setting=NameSetting("c_dimer", "", TP),
            expr=f"concetrations({p_total}, {{kd1}}, {{kd2}})['dimer']",
        ),
        "c_trimer": ParamLocalSetting(
            name_setting=NameSetting("c_trimer", "", TP),
            expr=f"concetrations({p_total}, {{kd1}}, {{kd2}})['trimer']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TP),
            expr="2.0 * {kon1} * {c_monomer}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TP),
            expr="{koff1}",
        ),
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "", TP),
            expr="{kon2} * {c_dimer}",
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "", TP),
            expr="{koff2} / 3.0",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", TP),
            expr="{kon2} * {c_monomer}",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", TP),
            expr="2.0 * {koff2} / 3.0",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TP),
            expr="pop_3st({kab}, {kba}, {kac}, {kca}, {kbc}, {kcb})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TP),
            expr="pop_3st({kab}, {kba}, {kac}, {kca}, {kbc}, {kcb})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TP),
            expr="pop_3st({kab}, {kba}, {kac}, {kca}, {kbc}, {kcb})['pc']",
        ),
    }


def register() -> None:
    model_factory.register(
        name=NAME,
        setting_maker=make_settings_3st_monomer_dimer_trimer,
    )
    user_functions = {
        "concetrations": calculate_concentrations,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

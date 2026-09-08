from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.optimize import root

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_2st
from chemex.models.factory import model_factory
from chemex.models.kinetic._oligomerization import solve_oligomerization_fractions
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import Array

NAME = "2st_monomer_trimer"

TP = ("temperature", "p_total")


def calculate_residuals(
    concentrations: Array,
    p_total: float,
    kd: float,
) -> Array:
    monomer, trimer = concentrations
    return np.array(
        [
            monomer + 3.0 * trimer - p_total,
            kd * trimer - monomer**3,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(p_total: float, kd: float) -> dict[str, float]:
    if math.isfinite(p_total) and p_total > 0.0 and math.isfinite(kd) and kd >= 1e-32:
        try:
            coefficient = p_total**2 / kd
        except OverflowError:
            msg = "Oligomerization concentration solver produced an invalid coefficient"
            raise RuntimeError(msg) from None
        monomer_fraction, (trimer_fraction,) = solve_oligomerization_fractions(
            ((3, coefficient),),
        )
        return {
            "monomer": p_total * monomer_fraction,
            "trimer": p_total * trimer_fraction,
        }

    concentrations_start = (p_total, 0.0)
    results = root(calculate_residuals, concentrations_start, args=(p_total, kd))
    return {"monomer": results["x"][0], "trimer": results["x"][1]}


def make_settings_2st_monomer_trimer(
    conditions: Conditions,
) -> dict[str, ParamLocalSetting]:
    p_total = conditions.p_total
    if p_total is None:
        msg = f"'p_total' must be specified to use the '{NAME}' model"
        raise ValueError(msg)
    return {
        "koff": ParamLocalSetting(
            name_setting=NameSetting("koff", "", ("temperature",)),
            value=100.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "kd": ParamLocalSetting(
            name_setting=NameSetting("kd", "", ("temperature",)),
            value=1e-6,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "kon": ParamLocalSetting(
            name_setting=NameSetting("kon", "", ("temperature",)),
            min=0.0,
            expr="{koff} / max({kd}, 1e-32)",
        ),
        "c_monomer": ParamLocalSetting(
            name_setting=NameSetting("c_monomer", "", TP),
            expr=f"concentrations({p_total}, {{kd}})['monomer']",
        ),
        "c_trimer": ParamLocalSetting(
            name_setting=NameSetting("c_trimer", "", TP),
            expr=f"concentrations({p_total}, {{kd}})['trimer']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TP),
            expr="3.0 * {kon} * {c_monomer} * {c_monomer}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TP),
            expr="{koff}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TP),
            expr="pop_2st({kab}, {kba})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TP),
            expr="pop_2st({kab}, {kba})['pb']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_2st_monomer_trimer)
    user_functions = {
        "concentrations": calculate_concentrations,
        "pop_2st": pop_2st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

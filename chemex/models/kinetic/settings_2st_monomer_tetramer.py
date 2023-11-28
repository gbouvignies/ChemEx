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

NAME = "2st_monomer_tetramer"

TP = ("temperature", "p_total")


def calculate_residuals(
    concentrations: ArrayFloat,
    p_total: float,
    kd: float,
) -> ArrayFloat:
    monomer, tetramer = concentrations
    return np.array(
        [
            monomer + 4.0 * tetramer - p_total,
            kd * tetramer - monomer**4,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(p_total: float, kd: float) -> dict[str, float]:
    concentrations_start = (p_total, 0.0)
    results = root(calculate_residuals, concentrations_start, args=(p_total, kd))
    return {"monomer": results["x"][0], "tetramer": results["x"][1]}


def make_settings_2st_monomer_tetramer(
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
            vary=True,
        ),
        "kd": ParamLocalSetting(
            name_setting=NameSetting("kd", "", ("temperature",)),
            value=1e-6,
            min=0.0,
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
        "c_tetramer": ParamLocalSetting(
            name_setting=NameSetting("c_tetramer", "", TP),
            expr=f"concentrations({p_total}, {{kd}})['tetramer']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TP),
            expr="4.0 * {kon} * {c_monomer}",
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
    model_factory.register(name=NAME, setting_maker=make_settings_2st_monomer_tetramer)
    user_functions = {
        "concentrations": calculate_concentrations,
        "pop_2st": pop_2st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

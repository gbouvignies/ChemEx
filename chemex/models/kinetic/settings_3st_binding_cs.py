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

NAME = "3st_binding_cs"

TPL = ("temperature", "p_total", "l_total")


def calculate_residuals(
    concentrations: ArrayFloat,
    p_total: float,
    l_total: float,
    kab: float,
    kba: float,
    kd: float,
) -> ArrayFloat:
    p1, p2, pl, l_ = concentrations
    return np.array(
        [
            p_total - p1 - p2 - pl,
            l_total - l_ - pl,
            kab * p1 - kba * p2,
            kd * pl - p2 * l_,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd: float,
    kab: float,
    kba: float,
) -> dict[str, float]:
    concentrations_start = (p_total, 0.0, 0.0, l_total)
    results = root(
        calculate_residuals,
        concentrations_start,
        args=(p_total, l_total, kab, kba, kd),
    )
    p1, p2, pl, l_ = results["x"]
    return {"p1": p1, "p2": p2, "pl": pl, "l": l_}


def make_settings_3st_binding_cs(
    conditions: Conditions,
) -> dict[str, ParamLocalSetting]:
    p0 = conditions.p_total
    l0 = conditions.l_total
    if p0 is None:
        msg = f"'p_total' must be specified to use the '{NAME}' model"
        raise ValueError(msg)
    if l0 is None:
        msg = f"'l_total' must be specified to use the '{NAME}' model"
        raise ValueError(msg)
    return {
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "koff": ParamLocalSetting(
            name_setting=NameSetting("koff", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kd": ParamLocalSetting(
            name_setting=NameSetting("kd", "", ("temperature",)),
            value=1e-3,
            min=0.0,
            vary=True,
        ),
        "kon": ParamLocalSetting(
            name_setting=NameSetting("kon", "", ("temperature",)),
            expr="{koff} / max({kd}, 1e-100)",
        ),
        "c_l": ParamLocalSetting(
            name_setting=NameSetting("c_l", "", TPL),
            expr=f"concentrations({p0}, {l0}, {{kab}}, {{kba}}, {{kd}})['l']",
        ),
        "c_p1": ParamLocalSetting(
            name_setting=NameSetting("c_p1", "", TPL),
            expr=f"concentrations({p0}, {l0}, {{kab}}, {{kba}}, {{kd}})['p1']",
        ),
        "c_p2": ParamLocalSetting(
            name_setting=NameSetting("c_p2", "", TPL),
            expr=f"concentrations({p0}, {l0}, {{kab}}, {{kba}}, {{kd}})['p2']",
        ),
        "c_pl": ParamLocalSetting(
            name_setting=NameSetting("c_pl", "", TPL),
            expr=f"concentrations({p0}, {l0}, {{kab}}, {{kba}}, {{kd}})['pl']",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", TPL),
            expr="{kon} * {c_l}",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", TPL),
            expr="{koff}",
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
    model_factory.register(name=NAME, setting_maker=make_settings_3st_binding_cs)
    user_functions = {
        "concentrations": calculate_concentrations,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

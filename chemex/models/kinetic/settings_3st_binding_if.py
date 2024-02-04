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

NAME = "3st_binding_if"

TPL = ("temperature", "p_total", "l_total")


def calculate_residuals(
    populations: ArrayFloat,
    p_total: float,
    l_total: float,
    kd_ab: float,
    kbc: float,
    kcb: float,
) -> ArrayFloat:
    pfree, pl1, pl2 = populations
    lfree = l_total - p_total + pfree
    return np.array(
        [
            pfree + pl1 + pl2 - p_total,
            lfree * pfree - kd_ab * pl1,
            kbc * pl1 - kcb * pl2,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd_ab: float,
    kbc: float,
    kcb: float,
) -> dict[str, float]:
    concentrations_start = (p_total, 0.0, 0.0)
    results = root(
        calculate_residuals,
        concentrations_start,
        args=(p_total, l_total, kd_ab, kbc, kcb),
    )
    return dict(zip(("pfree", "pl1", "pl2"), results["x"], strict=True))


def make_settings_3st_induced_fit(
    conditions: Conditions,
) -> dict[str, ParamLocalSetting]:
    p_total = conditions.p_total
    l_total = conditions.l_total
    if p_total is None:
        msg = f"'p_total' must be specified to use the '{NAME}' model"
        raise ValueError(msg)
    if l_total is None:
        msg = f"'l_total' must be specified to use the '{NAME}' model"
        raise ValueError(msg)
    return {
        "kd_app": ParamLocalSetting(
            name_setting=NameSetting("kd_app", "", ("temperature",)),
            value=1e-3,
            min=0.0,
            vary=True,
        ),
        "koff_ab": ParamLocalSetting(
            name_setting=NameSetting("koff_ab", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kd_ab": ParamLocalSetting(
            name_setting=NameSetting("kd_ab", "", ("temperature",)),
            expr="{kd_app} * (1 + {kbc} / {kcb})",
        ),
        "kon": ParamLocalSetting(
            name_setting=NameSetting("kon", "", ("temperature",)),
            expr="{koff_ab} / max({kd_ab}, 1e-100)",
        ),
        "pfree": ParamLocalSetting(
            name_setting=NameSetting("pfree", "", TPL),
            expr=(f"calc_conc({p_total},{l_total},{{kd_ab}},{{kbc}},{{kcb}})['pfree']"),
        ),
        "l_free": ParamLocalSetting(
            name_setting=NameSetting("l_free", "", TPL),
            expr=f"{l_total} - {p_total} + {{pfree}}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="{kon} * {l_free}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff_ab}",
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
    model_factory.register(name=NAME, setting_maker=make_settings_3st_induced_fit)
    user_functions = {
        "calc_conc": calculate_concentrations,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

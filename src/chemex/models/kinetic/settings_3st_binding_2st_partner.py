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

NAME = "3st_binding_partner_2st"

TPL = ("temperature", "p_total", "l_total")


def calculate_residuals(
    concentrations: ArrayFloat,
    p_total: float,
    l_total: float,
    kd1: float,
    kd2: float,
    keq: float,
) -> ArrayFloat:
    p, l1, l2, pl1, pl2 = concentrations
    return np.array(
        [
            l_total - (l1 + l2 + pl1 + pl2),
            p_total - (p + pl1 + pl2),
            kd1 * pl1 - p * l1,
            kd2 * pl2 - p * l2,
            keq * l1 - l2,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd1: float,
    kd2: float,
    keq: float,
) -> dict[str, float]:
    concentrations_start = (p_total, l_total / 2, l_total / 2, 0.0, 0.0)
    results = root(
        calculate_residuals,
        concentrations_start,
        args=(p_total, l_total, kd1, kd2, keq),
    )
    return dict(zip(("p", "l1", "l2", "pl1", "pl2"), results["x"], strict=True))


def make_settings_3st_binding_partner_2st(
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
        "koff_ab": ParamLocalSetting(
            name_setting=NameSetting("koff_ab", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kd_ab": ParamLocalSetting(
            name_setting=NameSetting("kd_ab", "", ("temperature",)),
            value=1e-3,
            min=0.0,
            vary=True,
        ),
        "koff_ac": ParamLocalSetting(
            name_setting=NameSetting("koff_ac", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kd_ac": ParamLocalSetting(
            name_setting=NameSetting("kd_ac", "", ("temperature",)),
            value=1e-3,
            min=0.0,
            vary=True,
        ),
        "keq": ParamLocalSetting(
            name_setting=NameSetting("keq", "", ("temperature",)),
            value=1.0,
            min=0.0,
            vary=True,
        ),
        "kex_bc": ParamLocalSetting(
            name_setting=NameSetting("kex_bc", "", ("temperature",)),
            value=1e3,
            min=0.0,
            vary=True,
        ),
        "kon_ab": ParamLocalSetting(
            name_setting=NameSetting("kon_ab", "", ("temperature",)),
            expr="{koff_ab} / max({kd_ab}, 1e-100)",
        ),
        "kon_ac": ParamLocalSetting(
            name_setting=NameSetting("kon_ac", "", ("temperature",)),
            expr="{koff_ac} / max({kd_ac}, 1e-100)",
        ),
        "l1_free": ParamLocalSetting(
            name_setting=NameSetting("l1_free", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq}})['l1']",
        ),
        "l2_free": ParamLocalSetting(
            name_setting=NameSetting("l1_free", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq}})['l2']",
        ),
        "pl1": ParamLocalSetting(
            name_setting=NameSetting("pl1", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq}})['pl1']",
        ),
        "pl2": ParamLocalSetting(
            name_setting=NameSetting("pl2", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq}})['pl2']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="{kon_ab} * {l1_free}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff_ab}",
        ),
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "", TPL),
            expr="{kon_ac} * {l2_free}",
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "", TPL),
            expr="{koff_ac}",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", TPL),
            expr="{kex_bc} * {pl2} / max({pl1} + {pl2}, 1e-100)",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", TPL),
            expr="{kex_bc} * {pl1} / max({pl1} + {pl2}, 1e-100)",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr="pop_3st({kab},{kba},{kac},{kca},{kbc},{kcb})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr="pop_3st({kab},{kba},{kac},{kca},{kbc},{kcb})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr="pop_3st({kab},{kba},{kac},{kca},{kbc},{kcb})['pc']",
        ),
    }


def register() -> None:
    model_factory.register(
        name=NAME,
        setting_maker=make_settings_3st_binding_partner_2st,
    )
    user_functions = {
        "calc_conc": calculate_concentrations,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.optimize import root

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import ArrayFloat

NAME = "4st_binding_partner_2st"

TPL = ("temperature", "p_total", "l_total")


def calculate_residuals(
    concentrations: ArrayFloat,
    p_total: float,
    l_total: float,
    kd1: float,
    kd2: float,
    keq_l: float,
    keq_pl: float,
) -> ArrayFloat:
    p, l1, l2, pl1, pl2, pl3 = concentrations
    return np.array(
        [
            l_total - (l1 + l2 + pl1 + pl2 + pl3),
            p_total - (p + pl1 + pl2 + pl3),
            kd1 * pl1 - p * l1,
            kd2 * pl2 - p * l2,
            keq_l * l1 - l2,
            keq_pl * pl2 - pl3,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd1: float,
    kd2: float,
    keq_l: float,
    keq_pl: float,
) -> dict[str, float]:
    concentrations_start = (p_total, l_total / 2, l_total / 2, 0.0, 0.0, 0.0)
    results = root(
        calculate_residuals,
        concentrations_start,
        args=(p_total, l_total, kd1, kd2, keq_l, keq_pl),
    )
    return dict(zip(("p", "l1", "l2", "pl1", "pl2", "pl3"), results["x"], strict=True))


def make_settings_4st_binding_partner_2st(
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
        "keq_l": ParamLocalSetting(
            name_setting=NameSetting("keq_l", "", ("temperature",)),
            value=1.0,
            min=0.0,
            vary=True,
        ),
        "keq_pl": ParamLocalSetting(
            name_setting=NameSetting("keq_pl", "", ("temperature",)),
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
        "kex_cd": ParamLocalSetting(
            name_setting=NameSetting("kex_cd", "", ("temperature",)),
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
        "p_free": ParamLocalSetting(
            name_setting=NameSetting("p_free", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq_l}},{{keq_pl}})['p']",
        ),
        "l1_free": ParamLocalSetting(
            name_setting=NameSetting("l1_free", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq_l}},{{keq_pl}})['l1']",
        ),
        "l2_free": ParamLocalSetting(
            name_setting=NameSetting("l1_free", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq_l}},{{keq_pl}})['l2']",
        ),
        "pl1": ParamLocalSetting(
            name_setting=NameSetting("pl1", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq_l}},{{keq_pl}})['pl1']",
        ),
        "pl2": ParamLocalSetting(
            name_setting=NameSetting("pl2", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq_l}},{{keq_pl}})['pl2']",
        ),
        "pl3": ParamLocalSetting(
            name_setting=NameSetting("pl3", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq_l}},{{keq_pl}})['pl3']",
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
        "kcd": ParamLocalSetting(
            name_setting=NameSetting("kcd", "", TPL),
            expr="{kex_cd} * {pl3} / max({pl2} + {pl3}, 1e-100)",
        ),
        "kdc": ParamLocalSetting(
            name_setting=NameSetting("kdc", "", TPL),
            expr="{kex_cd} * {pl2} / max({pl2} + {pl3}, 1e-100)",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr=f"{{p_free}} / {p_total}",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr=f"{{pl1}} / {p_total}",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr=f"{{pl2}} / {p_total}",
        ),
        "pd": ParamLocalSetting(
            name_setting=NameSetting("pd", "", TPL),
            expr=f"{{pl3}} / {p_total}",
        ),
    }


def register() -> None:
    model_factory.register(
        name=NAME,
        setting_maker=make_settings_4st_binding_partner_2st,
    )
    user_functions = {
        "calc_conc": calculate_concentrations,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

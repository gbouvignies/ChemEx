from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.optimize import root

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import ArrayFloat

NAME = "4st_binding_3_bound_states"

TPL = ("temperature", "p_total", "l_total")


def calculate_residuals(
    concentrations: ArrayFloat,
    p_total: float,
    l_total: float,
    kd_ab: float,
    keq_bc: float,
    keq_cd: float,
) -> ArrayFloat:
    p_, l_, pl1, pl2, pl3 = concentrations
    return np.array(
        [
            l_total - (l_ + pl1 + pl2 + pl3),
            p_total - (p_ + pl1 + pl2 + pl3),
            kd_ab * pl1 - p_ * l_,
            keq_bc * pl1 - pl2,
            keq_cd * pl2 - pl3,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd_ab: float,
    keq_bc: float,
    keq_cd: float,
) -> dict[str, float]:
    concentrations_start = (p_total, l_total, 0.0, 0.0, 0.0)
    results = root(
        calculate_residuals,
        concentrations_start,
        args=(p_total, l_total, kd_ab, keq_bc, keq_cd),
    )
    p_, l_, pl1, pl2, pl3 = results["x"]
    return {"p": p_, "l": l_, "pl1": pl1, "pl2": pl2, "pl3": pl3}


def make_settings_4st_binding_3_bound_states(
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
            value=1e-6,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "koff_ab": ParamLocalSetting(
            name_setting=NameSetting("koff_ab", "", ("temperature",)),
            value=100.0,
            min=0.0,
            max=1e6,
            vary=True,
        ),
        "kex_bc": ParamLocalSetting(
            name_setting=NameSetting("kex_bc", "", ("temperature",)),
            value=1000.0,
            min=0.0,
            max=1e6,
            vary=True,
        ),
        "keq_bc": ParamLocalSetting(
            name_setting=NameSetting("keq_bc", "", ("temperature",)),
            value=1.0,
            min=0.0,
            max=100.0,
            vary=True,
        ),
        "kex_cd": ParamLocalSetting(
            name_setting=NameSetting("kex_cd", "", ("temperature",)),
            value=1000.0,
            min=0.0,
            max=1e6,
            vary=True,
        ),
        "keq_cd": ParamLocalSetting(
            name_setting=NameSetting("keq_cd", "", ("temperature",)),
            value=1.0,
            min=0.0,
            max=100.0,
            vary=True,
        ),
        "kd_ab": ParamLocalSetting(
            name_setting=NameSetting("kd_ab", "", ("temperature",)),
            expr="{kd_app} * (1 + {keq_bc} + {keq_bc} * {keq_cd})",
        ),
        "kon_ab": ParamLocalSetting(
            name_setting=NameSetting("kon_ab", "", ("temperature",)),
            expr="{koff_ab} / max({kd_ab}, 1e-32)",
        ),
        "c_p": ParamLocalSetting(
            name_setting=NameSetting("c_p", "", TPL),
            expr=f"concentrations({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}})['p']",
        ),
        "c_l": ParamLocalSetting(
            name_setting=NameSetting("c_l", "", TPL),
            expr=f"concentrations({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}})['l']",
        ),
        "c_pl1": ParamLocalSetting(
            name_setting=NameSetting("c_pl1", "", TPL),
            expr=f"concentrations({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}})['pl1']",
        ),
        "c_pl2": ParamLocalSetting(
            name_setting=NameSetting("c_pl2", "", TPL),
            expr=f"concentrations({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}})['pl2']",
        ),
        "c_pl3": ParamLocalSetting(
            name_setting=NameSetting("c_pl3", "", TPL),
            expr=f"concentrations({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}})['pl3']",
        ),
        "c_pl": ParamLocalSetting(
            name_setting=NameSetting("c_pl", "", TPL),
            expr="{c_pl1}+{c_pl2}+{c_pl3}",
        ),
        "kd_eff": ParamLocalSetting(
            name_setting=NameSetting("kd_eff", "", TPL),
            expr="({c_p}+{c_l}) / max({c_pl}, 1e-32)",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="{kon_ab} * {c_l}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff_ab}",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", TPL),
            expr="{keq_bc}*{kex_bc}/(1.0+{keq_bc})",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", TPL),
            expr="{kex_bc}/(1.0+{keq_bc})",
        ),
        "kcd": ParamLocalSetting(
            name_setting=NameSetting("kcd", "", TPL),
            expr="{keq_cd}*{kex_cd}/(1.0+{keq_cd})",
        ),
        "kdc": ParamLocalSetting(
            name_setting=NameSetting("kdc", "", TPL),
            expr="{kex_cd}/(1.0+{keq_cd})",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr=f"{{c_p}} / {p_total}",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr=f"{{c_pl1}} / {p_total}",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr=f"{{c_pl2}} / {p_total}",
        ),
        "pd": ParamLocalSetting(
            name_setting=NameSetting("pd", "", TPL),
            expr=f"{{c_pl3}} / {p_total}",
        ),
    }


def register() -> None:
    model_factory.register(
        name=NAME,
        setting_maker=make_settings_4st_binding_3_bound_states,
    )
    user_functions = {
        "concentrations": calculate_concentrations,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

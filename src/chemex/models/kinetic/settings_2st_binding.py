from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.optimize import root

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import ArrayFloat

NAME = "2st_binding"

TPL = ("temperature", "p_total", "l_total")


def calculate_residuals(
    concentrations: ArrayFloat,
    p_total: float,
    l_total: float,
    kd: float,
) -> ArrayFloat:
    p_free, l_free, pl = concentrations
    return np.array(
        [
            l_total - l_free - pl,
            p_total - p_free - pl,
            kd * pl - p_free * l_free,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd: float,
) -> dict[str, float]:
    concentrations_start = (p_total, l_total, 0.0)
    results = root(
        calculate_residuals,
        concentrations_start,
        args=(p_total, l_total, kd),
    )
    p_free, l_free, pl = results["x"]
    return {"p_free": p_free, "l_free": l_free, "pl": pl}


def make_settings_2st_binding(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    p_total = conditions.p_total
    l_total = conditions.l_total
    if p_total is None or l_total is None:
        msg = "'p_total' and 'l_total' must be specified to use the '2st_binding' model"
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
            value=1e-3,
            min=0.0,
            vary=True,
        ),
        "kon": ParamLocalSetting(
            name_setting=NameSetting("kon", "", ("temperature",)),
            expr="{koff} / max({kd}, 1e-100)",
        ),
        "p_free": ParamLocalSetting(
            name_setting=NameSetting("p_free", "", TPL),
            expr=f"calc_conc({p_total}, {l_total}, {{kd}})['p_free']",
        ),
        "l_free": ParamLocalSetting(
            name_setting=NameSetting("l_free", "", TPL),
            expr=f"calc_conc({p_total}, {l_total}, {{kd}})['l_free']",
        ),
        "pl": ParamLocalSetting(
            name_setting=NameSetting("pl", "", TPL),
            expr=f"calc_conc({p_total}, {l_total}, {{kd}})['pl']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="{kon} * {l_free}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr=f"{{p_free}} / {p_total}",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr=f"{{pl}} / {p_total}",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_2st_binding)
    user_functions = {"calc_conc": calculate_concentrations}
    user_function_registry.register(name=NAME, user_functions=user_functions)

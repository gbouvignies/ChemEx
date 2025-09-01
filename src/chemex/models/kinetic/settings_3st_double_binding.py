from __future__ import annotations

"""This code imports the necessary modules and functions to define a 3-state
double binding model. The model is named "3st_double_binding" and it is
registered in the model factory. The model requires two parameters,
"p_total" and "l_total", to be specified in the conditions. The code
also defines a function, "calculate_concentrations", which calculates
the concentrations of free protein and protein bound to each ligand
given the total protein and ligand concentrations, as well as the
dissociation constants for each binding site.
"""

from functools import lru_cache

import numpy as np
from scipy.optimize import root

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import ArrayFloat

NAME = "3st_double_binding"

TPL = ("temperature", "p_total", "l_total")


def calculate_residuals(
    populations: ArrayFloat,
    p_total: float,
    l_total: float,
    kd_ab: float,
    kd_ac: float,
) -> ArrayFloat:
    pfree, pl1, pl2 = populations
    lfree = l_total - p_total + pfree
    return np.array(
        [
            pfree + pl1 + pl2 - p_total,
            lfree * pfree - kd_ab * pl1,
            lfree * pfree - kd_ac * pl2,
        ],
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd_ab: float,
    kd_ac: float,
) -> dict[str, float]:
    concentrations_start = (p_total, 0.0, 0.0)
    results = root(
        calculate_residuals,
        concentrations_start,
        args=(p_total, l_total, kd_ab, kd_ac),
    )
    return dict(zip(("pfree", "pl1", "pl2"), results["x"], strict=True))


def make_settings_3st_double_binding(
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
        "kon_ab": ParamLocalSetting(
            name_setting=NameSetting("kon_ab", "", ("temperature",)),
            expr="{koff_ab} / max({kd_ab}, 1e-100)",
        ),
        "kon_ac": ParamLocalSetting(
            name_setting=NameSetting("kon_ac", "", ("temperature",)),
            expr="{koff_ac} / max({kd_ac}, 1e-100)",
        ),
        "pfree": ParamLocalSetting(
            name_setting=NameSetting("pfree", "", TPL),
            expr=f"calc_conc({p_total}, {l_total}, {{kd_ab}}, {{kd_ac}})['pfree']",
        ),
        "l_free": ParamLocalSetting(
            name_setting=NameSetting("l_free", "", TPL),
            expr=f"{l_total} - {p_total} + {{pfree}}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="{kon_ab} * {l_free}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff_ab}",
        ),
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "", TPL),
            expr="{kon_ac} * {l_free}",
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "", TPL),
            expr="{koff_ac}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr="pop_3st({kab}, {kba}, {kac}, {kca}, 0.0, 0.0)['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr="pop_3st({kab}, {kba}, {kac}, {kca}, 0.0, 0.0)['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr="pop_3st({kab}, {kba}, {kac}, {kca}, 0.0, 0.0)['pc']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_3st_double_binding)
    user_functions = {
        "calc_conc": calculate_concentrations,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

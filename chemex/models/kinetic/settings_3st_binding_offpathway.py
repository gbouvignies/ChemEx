from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import root

from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry

if TYPE_CHECKING:
    from chemex.configuration.conditions import Conditions

NAME = "3st_binding_offpathway"

TPL = ("temperature", "p_total", "l_total")


def calculate_residuals(
    populations: np.ndarray,
    p_total: float,
    l_total: float,
    kd_ab: float,
    kac: float,
    kca: float,
) -> np.ndarray:
    pfree_a, pfree_c, pl_b = populations
    lfree = l_total - pl_b
    return np.array(
        [
            pfree_a + pfree_c + pl_b - p_total,
            lfree * pfree_a - kd_ab * pl_b,
            kac * pfree_a - kca * pfree_c,
        ]
    )


@lru_cache(maxsize=100)
def calculate_lfree(
    p_total: float, l_total: float, kd_ab: float, kac: float, kca: float
) -> float:
    concentrations_start = (p_total, 0.0, 0.0)
    results = root(
        calculate_residuals,
        concentrations_start,
        args=(p_total, l_total, kd_ab, kac, kca),
    )
    _pfree_a, _pfree_c, pl_b = results["x"]
    return l_total - pl_b


def make_settings_3st_binding_offpathway(
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
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kon_ab": ParamLocalSetting(
            name_setting=NameSetting("kon_ab", "", ("temperature",)),
            expr="{koff_ab} / max({kd_ab}, 1e-100)",
        ),
        "l_free": ParamLocalSetting(
            name_setting=NameSetting("l_free", "", TPL),
            expr=f"calculate_lfree({p_total}, {l_total}, {{kd_ab}}, {{kac}}, {{kca}})",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="{kon_ab} * {l_free}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff_ab}",
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
    model_factory.register(
        name=NAME, setting_maker=make_settings_3st_binding_offpathway
    )
    user_functions = {
        "calculate_lfree": calculate_lfree,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.optimize import root

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.models.kinetic._oligomerization import (
    MIN_POSITIVE_FLOAT,
    OligomerizationEquilibrium,
    concentrations_from_log_fractions,
    detailed_balance_forward_rate,
    log_equilibrium_coefficient,
    solve_oligomerization_fractions,
    validate_oligomerization_kd,
)
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry
from chemex.typing import Array

NAME = "3st_monomer_dimer_tetramer"

TP = ("temperature", "p_total")


def calculate_residuals(
    concentrations: Array,
    p_total: float,
    kd1: float,
    kd2: float,
) -> Array:
    monomer, dimer, tetramer = concentrations
    return np.array(
        [
            monomer + 2.0 * dimer + 4.0 * tetramer - p_total,
            kd1 * dimer - monomer**2,
            kd2 * tetramer - dimer**2,
        ],
    )


@lru_cache(maxsize=100)
def _calculate_equilibrium(
    p_total: float,
    kd1: float,
    kd2: float,
) -> OligomerizationEquilibrium:
    return solve_oligomerization_fractions(
        (
            (2, log_equilibrium_coefficient(p_total, 1, kd1)),
            (4, log_equilibrium_coefficient(p_total, 3, kd1, kd1, kd2)),
        ),
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    kd1: float,
    kd2: float,
) -> dict[str, float]:
    validate_oligomerization_kd(kd1)
    validate_oligomerization_kd(kd2)
    if math.isfinite(p_total) and p_total > 0.0:
        equilibrium = _calculate_equilibrium(p_total, kd1, kd2)
        monomer, dimer, tetramer = concentrations_from_log_fractions(
            p_total,
            (
                (1, equilibrium.log_monomer_fraction),
                (2, equilibrium.log_oligomer_fractions[0]),
                (4, equilibrium.log_oligomer_fractions[1]),
            ),
        )
        return {
            "monomer": monomer,
            "dimer": dimer,
            "tetramer": tetramer,
        }

    concentrations_start = (p_total, 0.0, 0.0)
    results = root(calculate_residuals, concentrations_start, args=(p_total, kd1, kd2))
    return {
        "monomer": results["x"][0],
        "dimer": results["x"][1],
        "tetramer": results["x"][2],
    }


@lru_cache(maxsize=100)
def calculate_rates(
    p_total: float,
    kd1: float,
    kd2: float,
    koff1: float,
    koff2: float,
) -> dict[str, float]:
    validate_oligomerization_kd(kd1)
    validate_oligomerization_kd(kd2)
    if p_total == 0.0:
        detailed_balance_forward_rate(koff1, 0.0, 0.0)
        detailed_balance_forward_rate(koff2, 0.0, 0.0)
        return {"kab": 0.0, "kbc": 0.0}

    equilibrium = _calculate_equilibrium(p_total, kd1, kd2)
    log_pa, log_pb, log_pc = equilibrium.log_tagged_fractions
    return {
        "kab": detailed_balance_forward_rate(koff1, log_pa, log_pb),
        "kbc": detailed_balance_forward_rate(koff2, log_pb, log_pc),
    }


def make_settings_3st_monomer_dimer_tetramer(
    conditions: Conditions,
) -> dict[str, ParamLocalSetting]:
    p_total = conditions.p_total
    if p_total is None:
        msg = f"'p_total' must be specified to use the '{NAME}' model"
        raise ValueError(msg)
    return {
        "kd1": ParamLocalSetting(
            name_setting=NameSetting("kd1", "", ("temperature",)),
            value=1e-6,
            min=MIN_POSITIVE_FLOAT,
            max=1.0,
            vary=True,
        ),
        "kd2": ParamLocalSetting(
            name_setting=NameSetting("kd2", "", ("temperature",)),
            value=1e-6,
            min=MIN_POSITIVE_FLOAT,
            max=1.0,
            vary=True,
        ),
        "koff1": ParamLocalSetting(
            name_setting=NameSetting("koff1", "", ("temperature",)),
            value=100.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "koff2": ParamLocalSetting(
            name_setting=NameSetting("koff2", "", ("temperature",)),
            value=100.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "c_monomer": ParamLocalSetting(
            name_setting=NameSetting("c_monomer", "", TP),
            expr=f"concentrations({p_total}, {{kd1}}, {{kd2}})['monomer']",
        ),
        "c_dimer": ParamLocalSetting(
            name_setting=NameSetting("c_dimer", "", TP),
            expr=f"concentrations({p_total}, {{kd1}}, {{kd2}})['dimer']",
        ),
        "c_tetramer": ParamLocalSetting(
            name_setting=NameSetting("c_tetramer", "", TP),
            expr=f"concentrations({p_total}, {{kd1}}, {{kd2}})['tetramer']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TP),
            expr=f"rates({p_total}, {{kd1}}, {{kd2}}, {{koff1}}, {{koff2}})['kab']",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TP),
            expr="{koff1}",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", TP),
            expr=f"rates({p_total}, {{kd1}}, {{kd2}}, {{koff1}}, {{koff2}})['kbc']",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", TP),
            expr="{koff2}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TP),
            expr="pop_3st({kab}, {kba}, 0.0, 0.0, {kbc}, {kcb})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TP),
            expr="pop_3st({kab}, {kba}, 0.0, 0.0, {kbc}, {kcb})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TP),
            expr="pop_3st({kab}, {kba}, 0.0, 0.0, {kbc}, {kcb})['pc']",
        ),
    }


def register() -> None:
    model_factory.register(
        name=NAME,
        setting_maker=make_settings_3st_monomer_dimer_tetramer,
    )
    user_functions = {
        "concentrations": calculate_concentrations,
        "rates": calculate_rates,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)

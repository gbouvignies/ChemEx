from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.optimize import root

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_2st
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
from chemex.parameters.userfunctions import (
    function_linearization_registry,
    population_linearizations,
    user_function_registry,
)
from chemex.typing import Array

NAME = "2st_monomer_tetramer"

TP = ("temperature", "p_total")


def calculate_residuals(
    concentrations: Array,
    p_total: float,
    kd: float,
) -> Array:
    monomer, tetramer = concentrations
    return np.array(
        [
            monomer + 4.0 * tetramer - p_total,
            kd * tetramer - monomer**4,
        ],
    )


@lru_cache(maxsize=100)
def _calculate_equilibrium(p_total: float, kd: float) -> OligomerizationEquilibrium:
    return solve_oligomerization_fractions(
        ((4, log_equilibrium_coefficient(p_total, 3, kd)),),
    )


@lru_cache(maxsize=100)
def calculate_concentrations(p_total: float, kd: float) -> dict[str, float]:
    validate_oligomerization_kd(kd)
    if math.isfinite(p_total) and p_total > 0.0:
        equilibrium = _calculate_equilibrium(p_total, kd)
        monomer, tetramer = concentrations_from_log_fractions(
            p_total,
            (
                (1, equilibrium.log_monomer_fraction),
                (4, equilibrium.log_oligomer_fractions[0]),
            ),
        )
        return {
            "monomer": monomer,
            "tetramer": tetramer,
        }

    concentrations_start = (p_total, 0.0)
    results = root(calculate_residuals, concentrations_start, args=(p_total, kd))
    return {"monomer": results["x"][0], "tetramer": results["x"][1]}


@lru_cache(maxsize=100)
def calculate_populations(p_total: float, kd: float) -> dict[str, float]:
    validate_oligomerization_kd(kd)
    if p_total == 0.0:
        return {"pa": 1.0, "pb": 0.0}
    equilibrium = _calculate_equilibrium(p_total, kd)
    return {
        "pa": equilibrium.monomer_fraction,
        "pb": 4.0 * equilibrium.oligomer_fractions[0],
    }


@lru_cache(maxsize=100)
def calculate_rates(p_total: float, kd: float, koff: float) -> dict[str, float]:
    validate_oligomerization_kd(kd)
    if p_total == 0.0:
        detailed_balance_forward_rate(koff, 0.0, 0.0)
        return {"kab": 0.0}
    equilibrium = _calculate_equilibrium(p_total, kd)
    return {
        "kab": detailed_balance_forward_rate(
            koff,
            equilibrium.log_tagged_fractions[0],
            equilibrium.log_tagged_fractions[1],
        ),
    }


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
            max=1.0e6,
            vary=True,
        ),
        "kd": ParamLocalSetting(
            name_setting=NameSetting("kd", "", ("temperature",)),
            value=1e-6,
            min=MIN_POSITIVE_FLOAT,
            max=1.0,
            vary=True,
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
            expr=f"rates({p_total}, {{kd}}, {{koff}})['kab']",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TP),
            expr="{koff}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TP),
            expr=f"populations({p_total}, {{kd}})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TP),
            expr=f"populations({p_total}, {{kd}})['pb']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_2st_monomer_tetramer)
    user_functions = {
        "concentrations": calculate_concentrations,
        "populations": calculate_populations,
        "rates": calculate_rates,
        "pop_2st": pop_2st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)
    function_linearization_registry.register(
        NAME,
        population_linearizations(
            (1.0e-3, 1.0e-3),
            ("nonnegative", "positive"),
            "pa",
            "pb",
        ),
    )

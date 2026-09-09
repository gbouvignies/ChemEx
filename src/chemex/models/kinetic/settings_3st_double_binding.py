"""This code imports the necessary modules and functions to define a 3-state
double binding model. The model is named "3st_double_binding" and it is
registered in the model factory. The model requires two parameters,
"p_total" and "l_total", to be specified in the conditions. The code
also defines a function, "calculate_concentrations", which calculates
the concentrations of free protein and protein bound to each ligand
given the total protein and ligand concentrations, as well as the
dissociation constants for each binding site.
"""

from __future__ import annotations

from functools import lru_cache

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.models.kinetic._binding import (
    MIN_POSITIVE_FLOAT,
    BindingEquilibrium,
    detailed_balance_rate,
    log_binding_weight,
    solve_binding_equilibrium,
)
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import (
    function_linearization_registry,
    population_linearizations,
    user_function_registry,
)

NAME = "3st_double_binding"

TPL = ("temperature", "p_total", "l_total")


@lru_cache(maxsize=100)
def _calculate_equilibrium(
    p_total: float,
    l_total: float,
    kd_ab: float,
    kd_ac: float,
) -> BindingEquilibrium:
    return solve_binding_equilibrium(
        p_total,
        l_total,
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(
            log_binding_weight(kd_ab),
            log_binding_weight(kd_ac),
        ),
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd_ab: float,
    kd_ac: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd_ab, kd_ac)
    return {
        "pfree": equilibrium.protein_free,
        "lfree": equilibrium.ligand_free,
        "pl1": equilibrium.complexes[0],
        "pl2": equilibrium.complexes[1],
    }


@lru_cache(maxsize=100)
def calculate_populations(
    p_total: float,
    l_total: float,
    kd_ab: float,
    kd_ac: float,
) -> dict[str, float]:
    pa, pb, pc = _calculate_equilibrium(
        p_total,
        l_total,
        kd_ab,
        kd_ac,
    ).populations
    return {"pa": pa, "pb": pb, "pc": pc}


@lru_cache(maxsize=100)
def calculate_rates(
    p_total: float,
    l_total: float,
    kd_ab: float,
    kd_ac: float,
    koff_ab: float,
    koff_ac: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd_ab, kd_ac)
    log_pa, log_pb, log_pc = equilibrium.log_populations
    return {
        "kab": detailed_balance_rate(koff_ab, log_pa, log_pb),
        "kac": detailed_balance_rate(koff_ac, log_pa, log_pc),
    }


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
            max=1.0e6,
            vary=True,
        ),
        "kd_ab": ParamLocalSetting(
            name_setting=NameSetting("kd_ab", "", ("temperature",)),
            value=1e-3,
            min=MIN_POSITIVE_FLOAT,
            max=1.0,
            vary=True,
        ),
        "koff_ac": ParamLocalSetting(
            name_setting=NameSetting("koff_ac", "", ("temperature",)),
            value=100.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "kd_ac": ParamLocalSetting(
            name_setting=NameSetting("kd_ac", "", ("temperature",)),
            value=1e-3,
            min=MIN_POSITIVE_FLOAT,
            max=1.0,
            vary=True,
        ),
        "kon_ab": ParamLocalSetting(
            name_setting=NameSetting("kon_ab", "", ("temperature",)),
            expr="{koff_ab} / {kd_ab}",
            report_only=True,
        ),
        "kon_ac": ParamLocalSetting(
            name_setting=NameSetting("kon_ac", "", ("temperature",)),
            expr="{koff_ac} / {kd_ac}",
            report_only=True,
        ),
        "pfree": ParamLocalSetting(
            name_setting=NameSetting("pfree", "", TPL),
            expr=f"calc_conc({p_total}, {l_total}, {{kd_ab}}, {{kd_ac}})['pfree']",
        ),
        "l_free": ParamLocalSetting(
            name_setting=NameSetting("l_free", "", TPL),
            expr=f"calc_conc({p_total}, {l_total}, {{kd_ab}}, {{kd_ac}})['lfree']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr=(
                f"rates({p_total}, {l_total}, {{kd_ab}}, {{kd_ac}}, "
                "{koff_ab}, {koff_ac})['kab']"
            ),
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff_ab}",
        ),
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "", TPL),
            expr=(
                f"rates({p_total}, {l_total}, {{kd_ab}}, {{kd_ac}}, "
                "{koff_ab}, {koff_ac})['kac']"
            ),
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "", TPL),
            expr="{koff_ac}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr=f"populations({p_total}, {l_total}, {{kd_ab}}, {{kd_ac}})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr=f"populations({p_total}, {l_total}, {{kd_ab}}, {{kd_ac}})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr=f"populations({p_total}, {l_total}, {{kd_ab}}, {{kd_ac}})['pc']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_3st_double_binding)
    user_functions = {
        "calc_conc": calculate_concentrations,
        "populations": calculate_populations,
        "rates": calculate_rates,
        "pop_3st": pop_3st,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)
    function_linearization_registry.register(
        NAME,
        population_linearizations(
            (1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3),
            ("nonnegative", "nonnegative", "positive", "positive"),
            "pa",
            "pb",
            "pc",
        ),
    )

from __future__ import annotations

from functools import lru_cache

from chemex.configuration.conditions import Conditions
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

NAME = "2st_binding"

TPL = ("temperature", "p_total", "l_total")


@lru_cache(maxsize=100)
def _calculate_equilibrium(
    p_total: float,
    l_total: float,
    kd: float,
) -> BindingEquilibrium:
    return solve_binding_equilibrium(
        p_total,
        l_total,
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(log_binding_weight(kd),),
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd)
    return {
        "p_free": equilibrium.protein_free,
        "l_free": equilibrium.ligand_free,
        "pl": equilibrium.complexes[0],
    }


@lru_cache(maxsize=100)
def calculate_rates(
    p_total: float,
    l_total: float,
    kd: float,
    koff: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd)
    return {
        "kab": detailed_balance_rate(
            koff,
            equilibrium.log_populations[0],
            equilibrium.log_populations[1],
        ),
    }


@lru_cache(maxsize=100)
def calculate_populations(
    p_total: float,
    l_total: float,
    kd: float,
) -> dict[str, float]:
    pa, pb = _calculate_equilibrium(p_total, l_total, kd).populations
    return {"pa": pa, "pb": pb}


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
            max=1.0e6,
            vary=True,
        ),
        "kd": ParamLocalSetting(
            name_setting=NameSetting("kd", "", ("temperature",)),
            value=1e-3,
            min=MIN_POSITIVE_FLOAT,
            max=1.0,
            vary=True,
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
        "kon": ParamLocalSetting(
            name_setting=NameSetting("kon", "", ("temperature",)),
            expr="{koff} / {kd}",
            report_only=True,
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr=f"rates({p_total}, {l_total}, {{kd}}, {{koff}})['kab']",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr=f"populations({p_total}, {l_total}, {{kd}})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr=f"populations({p_total}, {l_total}, {{kd}})['pb']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_2st_binding)
    user_functions = {
        "calc_conc": calculate_concentrations,
        "populations": calculate_populations,
        "rates": calculate_rates,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)
    function_linearization_registry.register(
        NAME,
        population_linearizations(
            (1.0e-3, 1.0e-3, 1.0e-3),
            ("nonnegative", "nonnegative", "positive"),
            "pa",
            "pb",
        ),
    )

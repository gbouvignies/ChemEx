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
    log_equilibrium_ratio,
    solve_binding_equilibrium,
    split_exchange_rate,
)
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import (
    function_linearization_registry,
    population_linearizations,
    user_function_registry,
)

NAME = "3st_binding_partner_2st"

TPL = ("temperature", "p_total", "l_total")


@lru_cache(maxsize=100)
def _calculate_equilibrium(
    p_total: float,
    l_total: float,
    kd1: float,
    kd2: float,
    keq: float,
) -> BindingEquilibrium:
    log_keq = log_equilibrium_ratio(keq)
    log_pl1 = log_binding_weight(kd1)
    log_pl2 = log_keq + log_binding_weight(kd2)
    return solve_binding_equilibrium(
        p_total,
        l_total,
        free_ligand_log_weights=(0.0, log_keq),
        complex_log_weights=(
            log_pl1,
            log_pl2,
        ),
        edge_log_ratios=(log_pl2 - log_pl1,),
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd1: float,
    kd2: float,
    keq: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd1, kd2, keq)
    l1, l2 = equilibrium.free_ligands
    pl1, pl2 = equilibrium.complexes
    return {
        "p": equilibrium.protein_free,
        "l1": l1,
        "l2": l2,
        "pl1": pl1,
        "pl2": pl2,
    }


@lru_cache(maxsize=100)
def calculate_populations(
    p_total: float,
    l_total: float,
    kd1: float,
    kd2: float,
    keq: float,
) -> dict[str, float]:
    pa, pb, pc = _calculate_equilibrium(
        p_total,
        l_total,
        kd1,
        kd2,
        keq,
    ).populations
    return {"pa": pa, "pb": pb, "pc": pc}


@lru_cache(maxsize=100)
def calculate_rates(
    p_total: float,
    l_total: float,
    kd1: float,
    kd2: float,
    keq: float,
    koff_ab: float,
    koff_ac: float,
    kex_bc: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd1, kd2, keq)
    log_pa, log_pb, log_pc = equilibrium.log_populations
    kbc, kcb = split_exchange_rate(
        kex_bc,
        0.0,
        equilibrium.edge_log_ratios[0],
    )
    return {
        "kab": detailed_balance_rate(koff_ab, log_pa, log_pb),
        "kac": detailed_balance_rate(koff_ac, log_pa, log_pc),
        "kbc": kbc,
        "kcb": kcb,
    }


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
    rates = (
        f"rates({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq}},"
        "{koff_ab},{koff_ac},{kex_bc})"
    )
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
        "keq": ParamLocalSetting(
            name_setting=NameSetting("keq", "", ("temperature",)),
            value=1.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "kex_bc": ParamLocalSetting(
            name_setting=NameSetting("kex_bc", "", ("temperature",)),
            value=1e3,
            min=0.0,
            max=1.0e6,
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
        "l1_free": ParamLocalSetting(
            name_setting=NameSetting("l1_free", "", TPL),
            expr=f"calc_conc({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq}})['l1']",
        ),
        "l2_free": ParamLocalSetting(
            name_setting=NameSetting("l2_free", "", TPL),
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
            expr=f"{rates}['kab']",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff_ab}",
        ),
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "", TPL),
            expr=f"{rates}['kac']",
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "", TPL),
            expr="{koff_ac}",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", TPL),
            expr=f"{rates}['kbc']",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", TPL),
            expr=f"{rates}['kcb']",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq}})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq}})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_ab}},{{kd_ac}},{{keq}})['pc']",
        ),
    }


def register() -> None:
    model_factory.register(
        name=NAME,
        setting_maker=make_settings_3st_binding_partner_2st,
    )
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
            (1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0),
            (
                "nonnegative",
                "nonnegative",
                "positive",
                "positive",
                "nonnegative",
            ),
            "pa",
            "pb",
            "pc",
        ),
    )

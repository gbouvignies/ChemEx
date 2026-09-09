from __future__ import annotations

from functools import lru_cache

from chemex.configuration.conditions import Conditions
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

NAME = "4st_binding_3_bound_states"

TPL = ("temperature", "p_total", "l_total")


@lru_cache(maxsize=100)
def _calculate_equilibrium(
    p_total: float,
    l_total: float,
    kd_ab: float,
    keq_bc: float,
    keq_cd: float,
) -> BindingEquilibrium:
    log_keq_bc = log_equilibrium_ratio(keq_bc)
    log_keq_cd = log_equilibrium_ratio(keq_cd)
    log_pl1 = log_binding_weight(kd_ab)
    return solve_binding_equilibrium(
        p_total,
        l_total,
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(
            log_pl1,
            log_pl1 + log_keq_bc,
            log_pl1 + log_keq_bc + log_keq_cd,
        ),
        edge_log_ratios=(log_keq_bc, log_keq_cd),
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd_ab: float,
    keq_bc: float,
    keq_cd: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(
        p_total,
        l_total,
        kd_ab,
        keq_bc,
        keq_cd,
    )
    pl1, pl2, pl3 = equilibrium.complexes
    return {
        "p": equilibrium.protein_free,
        "l": equilibrium.ligand_free,
        "pl1": pl1,
        "pl2": pl2,
        "pl3": pl3,
    }


@lru_cache(maxsize=100)
def calculate_rates(
    p_total: float,
    l_total: float,
    kd_ab: float,
    keq_bc: float,
    keq_cd: float,
    koff_ab: float,
    kex_bc: float,
    kex_cd: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(
        p_total,
        l_total,
        kd_ab,
        keq_bc,
        keq_cd,
    )
    kbc, kcb = split_exchange_rate(
        kex_bc,
        0.0,
        equilibrium.edge_log_ratios[0],
    )
    kcd, kdc = split_exchange_rate(
        kex_cd,
        0.0,
        equilibrium.edge_log_ratios[1],
    )
    return {
        "kab": detailed_balance_rate(
            koff_ab,
            equilibrium.log_populations[0],
            equilibrium.log_populations[1],
        ),
        "kbc": kbc,
        "kcb": kcb,
        "kcd": kcd,
        "kdc": kdc,
    }


@lru_cache(maxsize=100)
def calculate_populations(
    p_total: float,
    l_total: float,
    kd_ab: float,
    keq_bc: float,
    keq_cd: float,
) -> dict[str, float]:
    pa, pb, pc, pd = _calculate_equilibrium(
        p_total,
        l_total,
        kd_ab,
        keq_bc,
        keq_cd,
    ).populations
    return {"pa": pa, "pb": pb, "pc": pc, "pd": pd}


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
    rates = (
        f"rates({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}},"
        "{koff_ab},{kex_bc},{kex_cd})"
    )
    return {
        "kd_app": ParamLocalSetting(
            name_setting=NameSetting("kd_app", "", ("temperature",)),
            value=1e-6,
            min=MIN_POSITIVE_FLOAT,
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
            expr="{koff_ab} / {kd_ab}",
            report_only=True,
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
            expr="{kd_app}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr=f"{rates}['kab']",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{koff_ab}",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", TPL),
            expr=f"{rates}['kbc']",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", TPL),
            expr=f"{rates}['kcb']",
        ),
        "kcd": ParamLocalSetting(
            name_setting=NameSetting("kcd", "", TPL),
            expr=f"{rates}['kcd']",
        ),
        "kdc": ParamLocalSetting(
            name_setting=NameSetting("kdc", "", TPL),
            expr=f"{rates}['kdc']",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}})['pc']",
        ),
        "pd": ParamLocalSetting(
            name_setting=NameSetting("pd", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_ab}},{{keq_bc}},{{keq_cd}})['pd']",
        ),
    }


def register() -> None:
    model_factory.register(
        name=NAME,
        setting_maker=make_settings_4st_binding_3_bound_states,
    )
    user_functions = {
        "concentrations": calculate_concentrations,
        "populations": calculate_populations,
        "rates": calculate_rates,
    }
    user_function_registry.register(name=NAME, user_functions=user_functions)
    function_linearization_registry.register(
        NAME,
        population_linearizations(
            (1.0e-3, 1.0e-3, 1.0e-3, 1.0, 1.0),
            (
                "nonnegative",
                "nonnegative",
                "positive",
                "nonnegative",
                "nonnegative",
            ),
            "pa",
            "pb",
            "pc",
            "pd",
        ),
    )

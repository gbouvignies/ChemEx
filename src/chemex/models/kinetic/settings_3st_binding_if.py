"""Three-state induced-fit binding equilibrium and kinetics."""

from __future__ import annotations

import math
from functools import lru_cache

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.models.kinetic._binding import (
    MIN_POSITIVE_FLOAT,
    BindingEquilibrium,
    detailed_balance_rate,
    log_binding_weight,
    log_equilibrium_ratio,
    report_only_positive_log_value,
    solve_binding_equilibrium,
    split_exchange_rate,
)
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import (
    NumericalFunctionLinearization,
    function_linearization_registry,
    population_linearizations,
    user_function_registry,
)

NAME = "3st_binding_if"
TPL = ("temperature", "p_total", "l_total")


@lru_cache(maxsize=100)
def _calculate_equilibrium(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_bc: float,
) -> BindingEquilibrium:
    log_keq_bc = log_equilibrium_ratio(keq_bc)
    log_bound_partition = math.log1p(keq_bc)
    log_first_complex = log_binding_weight(kd_app) - log_bound_partition
    return solve_binding_equilibrium(
        p_total,
        l_total,
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(
            log_first_complex,
            -math.inf if log_keq_bc == -math.inf else log_first_complex + log_keq_bc,
        ),
        edge_log_ratios=(log_keq_bc,),
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_bc: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd_app, keq_bc)
    return {
        "a": equilibrium.free_proteins[0],
        "b": equilibrium.complexes[0],
        "c": equilibrium.complexes[1],
        "l": equilibrium.ligand_free,
    }


@lru_cache(maxsize=100)
def calculate_populations(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_bc: float,
) -> dict[str, float]:
    pa, pb, pc = _calculate_equilibrium(p_total, l_total, kd_app, keq_bc).populations
    return {"pa": pa, "pb": pb, "pc": pc}


@lru_cache(maxsize=100)
def calculate_conformational_rates(
    kex_bc: float,
    keq_bc: float,
) -> dict[str, float]:
    kbc, kcb = split_exchange_rate(kex_bc, 0.0, log_equilibrium_ratio(keq_bc))
    return {"kbc": kbc, "kcb": kcb}


@lru_cache(maxsize=100)
def calculate_binding_rates(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_bc: float,
    koff_ab: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd_app, keq_bc)
    return {
        "kab": detailed_balance_rate(
            koff_ab,
            equilibrium.log_populations[0],
            equilibrium.log_populations[1],
        ),
        "kba": koff_ab,
    }


def calculate_rates(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_bc: float,
    kex_bc: float,
    koff_ab: float,
) -> dict[str, float]:
    return calculate_binding_rates(
        p_total, l_total, kd_app, keq_bc, koff_ab
    ) | calculate_conformational_rates(kex_bc, keq_bc)


def calculate_intrinsic_kd(kd_app: float, keq_bc: float) -> float:
    return report_only_positive_log_value(math.log(kd_app) + math.log1p(keq_bc))


def calculate_kon(koff_ab: float, kd_app: float, keq_bc: float) -> float:
    if koff_ab == 0.0:
        return 0.0
    log_kd_ab = math.log(kd_app) + math.log1p(keq_bc)
    return report_only_positive_log_value(math.log(koff_ab) - log_kd_ab)


def calculate_intrinsic_values(kd_app: float, keq_bc: float) -> dict[str, float]:
    return {"kd": calculate_intrinsic_kd(kd_app, keq_bc)}


def calculate_kon_values(
    koff_ab: float,
    kd_app: float,
    keq_bc: float,
) -> dict[str, float]:
    return {"kon": calculate_kon(koff_ab, kd_app, keq_bc)}


def make_settings_3st_induced_fit(
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
    equilibrium = f"equilibrium({p_total},{l_total},{{kd_app}},{{keq_bc}})"
    binding_rates = (
        f"binding_rates({p_total},{l_total},{{kd_app}},{{keq_bc}},{{koff_ab}})"
    )
    conformational_rates = "conformational_rates({kex_bc},{keq_bc})"
    return {
        "kd_app": ParamLocalSetting(
            name_setting=NameSetting("kd_app", "", ("temperature",)),
            value=1.0e-3,
            min=MIN_POSITIVE_FLOAT,
            max=1.0,
            vary=True,
        ),
        "koff_ab": ParamLocalSetting(
            name_setting=NameSetting("koff_ab", "", ("temperature",)),
            value=100.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "keq_bc": ParamLocalSetting(
            name_setting=NameSetting("keq_bc", "", ("temperature",)),
            value=1.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "kex_bc": ParamLocalSetting(
            name_setting=NameSetting("kex_bc", "", ("temperature",)),
            value=200.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", ("temperature",)),
            expr=f"{conformational_rates}['kbc']",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", ("temperature",)),
            expr=f"{conformational_rates}['kcb']",
        ),
        "kd_ab": ParamLocalSetting(
            name_setting=NameSetting("kd_ab", "", ("temperature",)),
            expr="intrinsic_values({kd_app},{keq_bc})['kd']",
            report_only=True,
        ),
        "kon_ab": ParamLocalSetting(
            name_setting=NameSetting("kon_ab", "", ("temperature",)),
            expr="kon_values({koff_ab},{kd_app},{keq_bc})['kon']",
            report_only=True,
        ),
        "c_l": ParamLocalSetting(
            name_setting=NameSetting("c_l", "", TPL),
            expr=f"{equilibrium}['l']",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr=f"{binding_rates}['kab']",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr=f"{binding_rates}['kba']",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_app}},{{keq_bc}})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_app}},{{keq_bc}})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_app}},{{keq_bc}})['pc']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_3st_induced_fit)
    user_function_registry.register(
        name=NAME,
        user_functions={
            "equilibrium": calculate_concentrations,
            "populations": calculate_populations,
            "conformational_rates": calculate_conformational_rates,
            "binding_rates": calculate_binding_rates,
            "intrinsic_values": calculate_intrinsic_values,
            "kon_values": calculate_kon_values,
        },
    )
    function_linearization_registry.register(
        NAME,
        (
            *population_linearizations(
                (1.0e-3, 1.0e-3, 1.0e-3, 1.0),
                ("nonnegative", "nonnegative", "positive", "nonnegative"),
                "pa",
                "pb",
                "pc",
            ),
            NumericalFunctionLinearization(
                "conformational_rates",
                "kbc",
                (200.0, 1.0),
                ("nonnegative", "nonnegative"),
                output_scale=MIN_POSITIVE_FLOAT,
            ),
            NumericalFunctionLinearization(
                "conformational_rates",
                "kcb",
                (200.0, 1.0),
                ("nonnegative", "nonnegative"),
                output_scale=MIN_POSITIVE_FLOAT,
            ),
            NumericalFunctionLinearization(
                "binding_rates",
                "kab",
                (1.0e-3, 1.0e-3, 1.0e-3, 1.0, 100.0),
                (
                    "nonnegative",
                    "nonnegative",
                    "positive",
                    "nonnegative",
                    "nonnegative",
                ),
                output_scale=MIN_POSITIVE_FLOAT,
            ),
            NumericalFunctionLinearization(
                "binding_rates",
                "kba",
                (1.0e-3, 1.0e-3, 1.0e-3, 1.0, 100.0),
                (
                    "nonnegative",
                    "nonnegative",
                    "positive",
                    "nonnegative",
                    "nonnegative",
                ),
                output_scale=MIN_POSITIVE_FLOAT,
            ),
            NumericalFunctionLinearization(
                "intrinsic_values",
                "kd",
                (1.0e-3, 1.0),
                ("positive", "nonnegative"),
            ),
            NumericalFunctionLinearization(
                "kon_values",
                "kon",
                (100.0, 1.0e-3, 1.0),
                ("nonnegative", "positive", "nonnegative"),
            ),
        ),
    )

"""Three-state conformational-selection binding equilibrium and kinetics."""

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

NAME = "3st_binding_cs"
TPL = ("temperature", "p_total", "l_total")


def _log_keq_ab(keq_ab: float) -> float:
    if not math.isfinite(keq_ab) or keq_ab <= 0.0:
        msg = "3st_binding_cs KEQ_AB must be finite and strictly positive"
        raise ValueError(msg)
    return math.log(keq_ab)


@lru_cache(maxsize=100)
def _calculate_equilibrium(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_ab: float,
) -> BindingEquilibrium:
    log_keq_ab = _log_keq_ab(keq_ab)
    log_apo_weight = math.log1p(keq_ab)
    return solve_binding_equilibrium(
        p_total,
        l_total,
        free_protein_log_weights=(0.0, log_keq_ab),
        free_ligand_log_weights=(0.0,),
        complex_log_weights=(log_binding_weight(kd_app) + log_apo_weight,),
        edge_log_ratios=(log_keq_ab,),
    )


@lru_cache(maxsize=100)
def calculate_concentrations(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_ab: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd_app, keq_ab)
    a, b = equilibrium.free_proteins
    return {
        "a": a,
        "b": b,
        "c": equilibrium.complexes[0],
        "l": equilibrium.ligand_free,
    }


@lru_cache(maxsize=100)
def calculate_populations(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_ab: float,
) -> dict[str, float]:
    pa, pb, pc = _calculate_equilibrium(p_total, l_total, kd_app, keq_ab).populations
    return {"pa": pa, "pb": pb, "pc": pc}


@lru_cache(maxsize=100)
def calculate_conformational_rates(
    kex_ab: float,
    keq_ab: float,
) -> dict[str, float]:
    kab, kba = split_exchange_rate(kex_ab, 0.0, _log_keq_ab(keq_ab))
    return {"kab": kab, "kba": kba}


@lru_cache(maxsize=100)
def calculate_binding_rates(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_ab: float,
    koff_bc: float,
) -> dict[str, float]:
    equilibrium = _calculate_equilibrium(p_total, l_total, kd_app, keq_ab)
    return {
        "kbc": detailed_balance_rate(
            koff_bc,
            equilibrium.log_populations[1],
            equilibrium.log_populations[2],
        ),
        "kcb": koff_bc,
    }


def calculate_rates(
    p_total: float,
    l_total: float,
    kd_app: float,
    keq_ab: float,
    kex_ab: float,
    koff_bc: float,
) -> dict[str, float]:
    return calculate_conformational_rates(kex_ab, keq_ab) | calculate_binding_rates(
        p_total, l_total, kd_app, keq_ab, koff_bc
    )


def calculate_intrinsic_kd(kd_app: float, keq_ab: float) -> float:
    """Return KD_BC without silently flooring a positive unrepresentable value."""
    log_value = math.log(kd_app) + _log_keq_ab(keq_ab) - math.log1p(keq_ab)
    return report_only_positive_log_value(log_value)


def calculate_kon(koff_bc: float, kd_app: float, keq_ab: float) -> float:
    if koff_bc == 0.0:
        return 0.0
    log_kd_bc = math.log(kd_app) + _log_keq_ab(keq_ab) - math.log1p(keq_ab)
    return report_only_positive_log_value(math.log(koff_bc) - log_kd_bc)


def calculate_intrinsic_values(kd_app: float, keq_ab: float) -> dict[str, float]:
    return {"kd": calculate_intrinsic_kd(kd_app, keq_ab)}


def calculate_kon_values(
    koff_bc: float,
    kd_app: float,
    keq_ab: float,
) -> dict[str, float]:
    return {"kon": calculate_kon(koff_bc, kd_app, keq_ab)}


def make_settings_3st_binding_cs(
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
    equilibrium = f"equilibrium({p_total},{l_total},{{kd_app}},{{keq_ab}})"
    conformational_rates = "conformational_rates({kex_ab},{keq_ab})"
    binding_rates = (
        f"binding_rates({p_total},{l_total},{{kd_app}},{{keq_ab}},{{koff_bc}})"
    )
    return {
        "kd_app": ParamLocalSetting(
            name_setting=NameSetting("kd_app", "", ("temperature",)),
            value=1.0e-6,
            min=MIN_POSITIVE_FLOAT,
            max=1.0,
            vary=True,
        ),
        "koff_bc": ParamLocalSetting(
            name_setting=NameSetting("koff_bc", "", ("temperature",)),
            value=100.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "keq_ab": ParamLocalSetting(
            name_setting=NameSetting("keq_ab", "", ("temperature",)),
            value=1.0,
            min=MIN_POSITIVE_FLOAT,
            max=1.0e6,
            vary=True,
        ),
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting("kex_ab", "", ("temperature",)),
            value=200.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", ("temperature",)),
            expr=f"{conformational_rates}['kab']",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", ("temperature",)),
            expr=f"{conformational_rates}['kba']",
        ),
        "kd_bc": ParamLocalSetting(
            name_setting=NameSetting("kd_bc", "", ("temperature",)),
            expr="intrinsic_values({kd_app},{keq_ab})['kd']",
            report_only=True,
        ),
        "kon_bc": ParamLocalSetting(
            name_setting=NameSetting("kon_bc", "", ("temperature",)),
            expr="kon_values({koff_bc},{kd_app},{keq_ab})['kon']",
            report_only=True,
        ),
        "c_l": ParamLocalSetting(
            name_setting=NameSetting("c_l", "", TPL),
            expr=f"{equilibrium}['l']",
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", TPL),
            expr=f"{binding_rates}['kbc']",
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", TPL),
            expr=f"{binding_rates}['kcb']",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_app}},{{keq_ab}})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_app}},{{keq_ab}})['pb']",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            expr=f"populations({p_total},{l_total},{{kd_app}},{{keq_ab}})['pc']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_3st_binding_cs)
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
                ("nonnegative", "nonnegative", "positive", "positive"),
                "pa",
                "pb",
                "pc",
            ),
            NumericalFunctionLinearization(
                "conformational_rates",
                "kab",
                (200.0, 1.0),
                ("nonnegative", "positive"),
                output_scale=MIN_POSITIVE_FLOAT,
            ),
            NumericalFunctionLinearization(
                "conformational_rates",
                "kba",
                (200.0, 1.0),
                ("nonnegative", "positive"),
                output_scale=MIN_POSITIVE_FLOAT,
            ),
            NumericalFunctionLinearization(
                "binding_rates",
                "kbc",
                (1.0e-3, 1.0e-3, 1.0e-3, 1.0, 100.0),
                ("nonnegative", "nonnegative", "positive", "positive", "nonnegative"),
                output_scale=MIN_POSITIVE_FLOAT,
            ),
            NumericalFunctionLinearization(
                "binding_rates",
                "kcb",
                (1.0e-3, 1.0e-3, 1.0e-3, 1.0, 100.0),
                ("nonnegative", "nonnegative", "positive", "positive", "nonnegative"),
                output_scale=MIN_POSITIVE_FLOAT,
            ),
            NumericalFunctionLinearization(
                "intrinsic_values",
                "kd",
                (1.0e-3, 1.0),
                ("positive", "positive"),
            ),
            NumericalFunctionLinearization(
                "kon_values",
                "kon",
                (100.0, 1.0e-3, 1.0),
                ("nonnegative", "positive", "positive"),
            ),
        ),
    )

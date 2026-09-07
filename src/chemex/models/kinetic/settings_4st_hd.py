from __future__ import annotations

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting

NAME = "4st_hd"

TPL = ("temperature", "p_total", "l_total")


def create_pop_4st_hd_settings() -> dict[str, ParamLocalSetting]:
    population_conditions = (*TPL, "d2o")
    denominator_a = "(1.0 + {d2o} * ({phi_a} - 1.0))"
    denominator_b = "(1.0 + {d2o} * ({phi_b} - 1.0))"
    return {
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "g", population_conditions),
            min=0.0,
            max=1.0,
            expr=f"(1.0 - {{pop_b}}) * (1.0 - {{d2o}}) / {denominator_a}",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "g", population_conditions),
            min=0.0,
            max=1.0,
            expr=f"{{pop_b}} * (1.0 - {{d2o}}) / {denominator_b}",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "g", population_conditions),
            min=0.0,
            max=1.0,
            expr=f"(1.0 - {{pop_b}}) * {{d2o}} * {{phi_a}} / {denominator_a}",
        ),
        "pd": ParamLocalSetting(
            name_setting=NameSetting("pd", "g", population_conditions),
            min=0.0,
            max=1.0,
            expr=f"{{pop_b}} * {{d2o}} * {{phi_b}} / {denominator_b}",
        ),
    }


def make_settings_4st_hd(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    d2o: float = conditions.d2o if conditions.d2o is not None else 0.1
    return {
        "d2o": ParamLocalSetting(
            name_setting=NameSetting(
                "d2o",
                "",
                ("d2o",),
                allow_residue_specific=False,
            ),
            value=d2o,
            min=0.0,
            max=1.0,
        ),
        "pop_b": ParamLocalSetting(
            name_setting=NameSetting("pop_b", "", TPL),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting("kex_ab", "", TPL),
            min=0.0,
            max=1.0e6,
            value=0.0,
            vary=True,
        ),
        "kdh_a": ParamLocalSetting(
            name_setting=NameSetting("kdh_a", "g", ("temperature",)),
            value=1.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "kdh_b": ParamLocalSetting(
            name_setting=NameSetting("kdh_b", "g", ("temperature",)),
            value=1.0,
            min=0.0,
            max=1.0e6,
            vary=True,
        ),
        "phi_a": ParamLocalSetting(
            name_setting=NameSetting("phi_a", "g", ("temperature",)),
            value=1.1,
            min=0.75,
            max=1.50,
        ),
        "phi_b": ParamLocalSetting(
            name_setting=NameSetting("phi_b", "g", ("temperature",)),
            value=1.1,
            min=0.75,
            max=1.50,
            expr="{phi_a}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            min=0.0,
            expr="{pop_b} * {kex_ab}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            min=0.0,
            expr="(1.0 - {pop_b}) * {kex_ab}",
        ),
        "kcd": ParamLocalSetting(
            name_setting=NameSetting("kcd", "", TPL),
            min=0.0,
            expr="{pop_b} * {kex_ab}",
        ),
        "kdc": ParamLocalSetting(
            name_setting=NameSetting("kdc", "", TPL),
            min=0.0,
            expr="(1.0 - {pop_b}) * {kex_ab}",
        ),
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "g", ("temperature", "d2o")),
            min=0.0,
            expr="{d2o} * {kdh_a} * {phi_a}",
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "g", ("temperature", "d2o")),
            min=0.0,
            expr="(1.0 - {d2o}) * {kdh_a}",
        ),
        "kbd": ParamLocalSetting(
            name_setting=NameSetting("kbd", "g", ("temperature", "d2o")),
            min=0.0,
            expr="{d2o} * {kdh_b} * {phi_b}",
        ),
        "kdb": ParamLocalSetting(
            name_setting=NameSetting("kdb", "g", ("temperature", "d2o")),
            min=0.0,
            expr="(1.0 - {d2o}) * {kdh_b}",
        ),
        **create_pop_4st_hd_settings(),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_4st_hd)

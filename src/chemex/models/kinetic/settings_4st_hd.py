from __future__ import annotations

from itertools import combinations

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_4st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry

NAME = "4st_hd"

TPL = ("temperature", "p_total", "l_total")


def create_pop_4st_eyring_settings() -> dict[str, ParamLocalSetting]:
    arguments = ", ".join(
        f"{{k{i}{j}}}, {{k{j}{i}}}" for i, j in combinations("abcd", 2)
    )
    return {
        f"p{state}": ParamLocalSetting(
            name_setting=NameSetting(f"p{state}", "", TPL),
            expr=(f"pop_4st({arguments})['p{state}']"),
        )
        for state in "abcd"
    }


def make_settings_4st_hd(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    d2o: float = conditions.d2o if conditions.d2o is not None else 0.1
    return {
        "d2o": ParamLocalSetting(
            name_setting=NameSetting("d2o", "", ("d2o",)),
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
            value=0.0,
            vary=True,
        ),
        "kdh_a": ParamLocalSetting(
            name_setting=NameSetting("kdh_a", "g", ("temperature",)),
            value=1.0,
            min=0.0,
            vary=True,
        ),
        "kdh_b": ParamLocalSetting(
            name_setting=NameSetting("kdh_b", "g", ("temperature",)),
            value=1.0,
            min=0.0,
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
            value=0.02,
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
        **create_pop_4st_eyring_settings(),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_4st_hd)
    user_function_registry.register(name=NAME, user_functions={"pop_4st": pop_4st})

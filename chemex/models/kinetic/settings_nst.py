from __future__ import annotations

from itertools import permutations

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting

TPL = ("temperature", "p_total", "l_total")


def kex_p_to_k(states: str) -> str:
    kex = f"{{kex_{min(states)}{max(states)}}}"
    p1 = f"{{p{states[0]}}}"
    p2 = f"{{p{states[1]}}}"
    return f"{kex} * {p2} / max({p1} + {p2}, 1e-100)"


def create_kij_settings(states: str = "abc") -> dict[str, ParamLocalSetting]:
    if len(states) <= 1:
        msg = "'state_nb' should be larger than 2"
        raise ValueError(msg)
    return {
        f"k{i}{j}": ParamLocalSetting(
            name_setting=NameSetting(f"k{i}{j}", "", TPL),
            expr=kex_p_to_k(f"{i}{j}"),
        )
        for i, j in permutations(states, 2)
    }


def make_settings_3st_linear(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    kij_settings = create_kij_settings("abc")
    del kij_settings["kac"]
    del kij_settings["kca"]
    return {
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting("kex_ab", "", TPL),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "kex_bc": ParamLocalSetting(
            name_setting=NameSetting("kex_bc", "", TPL),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr="1.0 - {pb} - {pc}",
        ),
        **kij_settings,
    }


def make_settings_3st_fork(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    kij_settings = create_kij_settings("abc")
    del kij_settings["kbc"]
    del kij_settings["kcb"]
    return {
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", TPL),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting("kex_ab", "", TPL),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "kex_ac": ParamLocalSetting(
            name_setting=NameSetting("kex_ac", "", TPL),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr="1.0 - {pb} - {pc}",
        ),
        **kij_settings,
    }


def make_settings_3st(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {
        **make_settings_3st_fork(conditions),
        **make_settings_3st_linear(conditions),
    }


def make_settings_4st(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {
        **make_settings_3st(conditions),
        "pd": ParamLocalSetting(
            name_setting=NameSetting("pd", "", TPL),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "kex_ad": ParamLocalSetting(
            name_setting=NameSetting("kex_ad", "", TPL),
            min=0.0,
            value=0.0,
        ),
        "kex_bd": ParamLocalSetting(
            name_setting=NameSetting("kex_bd", "", TPL),
            min=0.0,
            value=0.0,
        ),
        "kex_cd": ParamLocalSetting(
            name_setting=NameSetting("kex_cd", "", TPL),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            min=0.0,
            max=1.0,
            expr="1.0 - {pb} - {pc} - {pd}",
        ),
        **create_kij_settings("abcd"),
    }


def register() -> None:
    model_factory.register(name="3st", setting_maker=make_settings_3st)
    model_factory.register(name="3st_triangle", setting_maker=make_settings_3st)
    model_factory.register(name="3st_linear", setting_maker=make_settings_3st_linear)
    model_factory.register(name="3st_fork", setting_maker=make_settings_3st_fork)
    model_factory.register(name="4st", setting_maker=make_settings_4st)

from __future__ import annotations

from functools import cache
from typing import TYPE_CHECKING

from scipy import constants

from chemex.parameters.setting import NameSetting
from chemex.parameters.setting import ParamLocalSetting

if TYPE_CHECKING:
    from chemex.configuration.conditions import Conditions


def make_settings_1st(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa"), value=1.0, min=0.99, max=1.0
        ),
    }


def make_settings_2st(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_ab", "", ("temperature", "p_total", "l_total")
            ),
            value=200.0,
            min=0.0,
            vary=True,
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", ("temperature", "p_total", "l_total")),
            value=0.05,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr="1.0 - {pb}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr="{kex_ab} * {pb}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr="{kex_ab} * {pa}",
        ),
    }


def make_settings_2st_rs(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_ab", "g", ("temperature", "p_total", "l_total")
            ),
            value=200.0,
            min=0.0,
            vary=True,
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "g", ("temperature", "p_total", "l_total")),
            value=0.05,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "g", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr="1.0 - {pb}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "g", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr="{kex_ab} * {pb}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "g", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr="{kex_ab} * {pa}",
        ),
    }


def make_settings_2st_hd(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    d2o = conditions.d2o if conditions.d2o is not None else 0.1
    return {
        "d2o": ParamLocalSetting(
            name_setting=NameSetting("d2o", "", ("d2o",)),
            value=d2o,
            min=0.0,
            max=1.0,
        ),
        "kdh": ParamLocalSetting(
            name_setting=NameSetting("kdh", "g", ("temperature",)),
            value=1.0,
            min=0.0,
            vary=True,
        ),
        "phi": ParamLocalSetting(
            name_setting=NameSetting("phi", "g", ("temperature",)),
            value=1.1,
            min=0.75,
            max=1.50,
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "g", ("temperature", "d2o")),
            min=0.0,
            expr="{d2o} * {kdh} * {phi}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "g", ("temperature", "d2o")),
            min=0.0,
            expr="(1.0 - {d2o}) * {kdh}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting(
                "pa", "g", ("temperature", "p_total", "l_total", "d2o")
            ),
            min=0.0,
            max=1.0,
            expr="(1.0 - {d2o}) / (1.0 + {d2o} * ({phi} - 1.0))",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting(
                "pb", "g", ("temperature", "p_total", "l_total", "d2o")
            ),
            min=0.0,
            max=1.0,
            expr="{d2o} * {phi} / (1.0 + {d2o} * ({phi} - 1.0))",
        ),
    }


def hs_to_k(states, celsius):
    kelvin = celsius + 273.15
    kbt_h = constants.k * kelvin / constants.h
    rt = constants.R * kelvin
    expr_dh = f"{{dh_{min(states)}{max(states)}}}"
    expr_ds = f"{{ds_{min(states)}{max(states)}}}"
    if states[0] != "a":
        expr_dh += f" - {{dh_{states[0]}}}"
        expr_ds += f" - {{ds_{states[0]}}}"
    return f"{kbt_h} * exp(-(({expr_dh}) - {kelvin} * ({expr_ds})) / {rt})"


def hs_to_p(state, states, celsius):
    kelvin = celsius + 273.15
    rt = constants.R * kelvin
    dg = {
        a_state: f"exp(" f"-({{dh_{a_state}}} - {kelvin} * {{ds_{a_state}}}) / {rt})"
        if a_state != "a"
        else "1.0"
        for a_state in states
    }

    return f"{dg[state]} / ({' + '.join(dg.values())})"


def make_settings_2st_eyring(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    celsius = conditions.temperature
    return {
        "dh_b": ParamLocalSetting(
            name_setting=NameSetting("dh_b", "", ("p_total", "l_total")),
            value=8e3,
            vary=True,
        ),
        "ds_b": ParamLocalSetting(
            name_setting=NameSetting("ds_b", "", ("p_total", "l_total")),
            value=0.0,
            vary=False,
        ),
        "dh_ab": ParamLocalSetting(
            name_setting=NameSetting("dh_ab", "", ("p_total", "l_total")),
            value=6.5e4,
            vary=True,
        ),
        "ds_ab": ParamLocalSetting(
            name_setting=NameSetting("ds_ab", "", ("p_total", "l_total")),
            value=0.0,
            vary=False,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr=hs_to_p("a", "ab", celsius),
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr=hs_to_p("b", "ab", celsius),
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=hs_to_k("ab", celsius),
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=hs_to_k("ba", celsius),
        ),
    }


def make_settings_2st_binding(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    p_total = conditions.p_total
    l_total = conditions.l_total
    if p_total is None or l_total is None:
        raise ValueError(
            "'p_total' and 'l_total' must be specified to use the '2st_binding' model"
        )
    delta = l_total - p_total
    expr_kab = (
        f"{{kon}} * 0.5 * ({delta} - {{kd}} "
        f"+ sqrt(({delta} - {{kd}}) ** 2 + 4.0 * {{kd}} * {l_total}))"
    )
    return {
        "koff": ParamLocalSetting(
            name_setting=NameSetting("koff", "", ("temperature",)),
            value=100.0,
            min=0.0,
            vary=True,
        ),
        "kd": ParamLocalSetting(
            name_setting=NameSetting("kd", "", ("temperature",)),
            min=0.0,
            vary=True,
        ),
        "kon": ParamLocalSetting(
            name_setting=NameSetting("kon", "", ("temperature",)),
            value=1.0e7,
            min=0.0,
            expr="{koff} / max({kd}, 1e-100)",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=expr_kab,
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr="{koff}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr="{kba} / max({kba} + {kab}, 1e-100)",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr="{kab} / max({kba} + {kab}, 1e-100)",
        ),
    }


def kex_p_to_k(states):
    kex = f"{{kex_{min(states)}{max(states)}}}"
    p1 = f"{{p{states[0]}}}"
    p2 = f"{{p{states[1]}}}"
    return f"{kex} * {p2} / max({p1} + {p2}, 1e-100)"


def make_settings_3st(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", ("temperature", "p_total", "l_total")),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", ("temperature", "p_total", "l_total")),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_ab", "", ("temperature", "p_total", "l_total")
            ),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "kex_ac": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_ac", "", ("temperature", "p_total", "l_total")
            ),
            min=0.0,
            value=0.0,
        ),
        "kex_bc": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_bc", "", ("temperature", "p_total", "l_total")
            ),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr="1.0 - {pb} - {pc}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("ab"),
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("ba"),
        ),
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("ac"),
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("ca"),
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("bc"),
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("cb"),
        ),
    }


def make_settings_3st_eyring(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    celsius = conditions.temperature
    return {
        "dh_b": ParamLocalSetting(
            name_setting=NameSetting("dh_b", "", ("p_total", "l_total")),
            value=8e3,
            vary=True,
        ),
        "dh_c": ParamLocalSetting(
            name_setting=NameSetting("dh_c", "", ("p_total", "l_total")),
            value=8e3,
            vary=True,
        ),
        "dh_ab": ParamLocalSetting(
            name_setting=NameSetting("dh_ab", "", ("p_total", "l_total")),
            value=6.5e4,
            vary=True,
        ),
        "dh_bc": ParamLocalSetting(
            name_setting=NameSetting("dh_bc", "", ("p_total", "l_total")),
            value=6.5e4,
            vary=True,
        ),
        "dh_ac": ParamLocalSetting(
            name_setting=NameSetting("dh_ac", "", ("p_total", "l_total")),
            value=1.0e10,
            vary=False,
        ),
        "ds_b": ParamLocalSetting(
            name_setting=NameSetting("ds_b", "", ("p_total", "l_total")),
            value=0.0,
            vary=False,
        ),
        "ds_c": ParamLocalSetting(
            name_setting=NameSetting("ds_c", "", ("p_total", "l_total")),
            value=0.0,
            vary=False,
        ),
        "ds_ab": ParamLocalSetting(
            name_setting=NameSetting("ds_ab", "", ("p_total", "l_total")),
            value=0.0,
            vary=False,
        ),
        "ds_bc": ParamLocalSetting(
            name_setting=NameSetting("ds_bc", "", ("p_total", "l_total")),
            value=0.0,
            vary=False,
        ),
        "ds_ac": ParamLocalSetting(
            name_setting=NameSetting("ds_ac", "", ("p_total", "l_total")),
            value=0.0,
            vary=False,
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=hs_to_k("ab", celsius),
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=hs_to_k("ba", celsius),
        ),
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=hs_to_k("ac", celsius),
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=hs_to_k("ca", celsius),
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=hs_to_k("bc", celsius),
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=hs_to_k("cb", celsius),
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr=hs_to_p("a", "abc", celsius),
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr=hs_to_p("b", "abc", celsius),
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr=hs_to_p("c", "abc", celsius),
        ),
    }


def make_settings_4st(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", ("temperature", "p_total", "l_total")),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting("pc", "", ("temperature", "p_total", "l_total")),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "pd": ParamLocalSetting(
            name_setting=NameSetting("pd", "", ("temperature", "p_total", "l_total")),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_ab", "", ("temperature", "p_total", "l_total")
            ),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "kex_ac": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_ac", "", ("temperature", "p_total", "l_total")
            ),
            min=0.0,
            value=0.0,
        ),
        "kex_ad": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_ad", "", ("temperature", "p_total", "l_total")
            ),
            min=0.0,
            value=0.0,
        ),
        "kex_bc": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_bc", "", ("temperature", "p_total", "l_total")
            ),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "kex_bd": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_bd", "", ("temperature", "p_total", "l_total")
            ),
            min=0.0,
            value=0.0,
        ),
        "kex_cd": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_cd", "", ("temperature", "p_total", "l_total")
            ),
            min=0.0,
            value=200.0,
            vary=True,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            max=1.0,
            expr="1.0 - {pb} - {pc} - {pd}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("ab"),
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("ba"),
        ),
        "kac": ParamLocalSetting(
            name_setting=NameSetting("kac", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("ac"),
        ),
        "kca": ParamLocalSetting(
            name_setting=NameSetting("kca", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("ca"),
        ),
        "kad": ParamLocalSetting(
            name_setting=NameSetting("kad", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("ad"),
        ),
        "kda": ParamLocalSetting(
            name_setting=NameSetting("kda", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("da"),
        ),
        "kbc": ParamLocalSetting(
            name_setting=NameSetting("kbc", "", ("temperature", "p_total", "l_total")),
            expr=kex_p_to_k("bc"),
        ),
        "kcb": ParamLocalSetting(
            name_setting=NameSetting("kcb", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("cb"),
        ),
        "kbd": ParamLocalSetting(
            name_setting=NameSetting("kbd", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("bd"),
        ),
        "kdb": ParamLocalSetting(
            name_setting=NameSetting("kdb", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("db"),
        ),
        "kcd": ParamLocalSetting(
            name_setting=NameSetting("kcd", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("cd"),
        ),
        "kdc": ParamLocalSetting(
            name_setting=NameSetting("kdc", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr=kex_p_to_k("dc"),
        ),
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
            name_setting=NameSetting(
                "pop_b", "", ("temperature", "p_total", "l_total")
            ),
            value=0.02,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting(
                "kex_ab", "", ("temperature", "p_total", "l_total")
            ),
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
            name_setting=NameSetting("kab", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr="{pop_b} * {kex_ab}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr="(1.0 - {pop_b}) * {kex_ab}",
        ),
        "kcd": ParamLocalSetting(
            name_setting=NameSetting("kcd", "", ("temperature", "p_total", "l_total")),
            min=0.0,
            expr="{pop_b} * {kex_ab}",
        ),
        "kdc": ParamLocalSetting(
            name_setting=NameSetting("kdc", "", ("temperature", "p_total", "l_total")),
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
        "pa": ParamLocalSetting(
            name_setting=NameSetting(
                "pa", "g", ("temperature", "p_total", "l_total", "d2o")
            ),
            min=0.0,
            max=1.0,
            expr="(1.0 - {pop_b}) * (1.0 - {d2o}) / (1.0 + {d2o} * ({phi_a} - 1.0))",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting(
                "pb", "g", ("temperature", "p_total", "l_total", "d2o")
            ),
            min=0.0,
            max=1.0,
            expr="{pop_b} * (1.0 - {d2o}) / (1.0 + {d2o} * ({phi_b} - 1.0))",
        ),
        "pc": ParamLocalSetting(
            name_setting=NameSetting(
                "pc", "g", ("temperature", "p_total", "l_total", "d2o")
            ),
            min=0.0,
            max=1.0,
            expr="(1.0 - {pop_b}) * {d2o} * {phi_a} / (1.0 + {d2o} * ({phi_a} - 1.0))",
        ),
        "pd": ParamLocalSetting(
            name_setting=NameSetting(
                "pd", "g", ("temperature", "p_total", "l_total", "d2o")
            ),
            min=0.0,
            max=1.0,
            expr="{pop_b} * {d2o} * {phi_a} / (1.0 + {d2o} * ({phi_b} - 1.0))",
        ),
    }


kinetic_setting_makers = {
    "1st": make_settings_1st,
    "2st": make_settings_2st,
    "3st": make_settings_3st,
    "4st": make_settings_4st,
    "2st_rs": make_settings_2st_rs,
    "2st_hd": make_settings_2st_hd,
    "2st_eyring": make_settings_2st_eyring,
    "3st_eyring": make_settings_3st_eyring,
    "2st_binding": make_settings_2st_binding,
    "4st_hd": make_settings_4st_hd,
}


@cache
def build_kinetic_param_settings(
    model_name: str, conditions: Conditions
) -> dict[str, ParamLocalSetting]:
    return kinetic_setting_makers[model_name](conditions)

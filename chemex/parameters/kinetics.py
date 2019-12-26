import collections
import functools as ft
import string

import numpy as np
import scipy.constants as cst


Model = collections.namedtuple("Model", ["name", "states", "model_free"])


def parse_model(name):
    name_, *ext = name.split(".")
    name_ = _check_model_name(name_)
    state_nb = int(name_[0])
    states = _get_state_names(state_nb)
    make_settings_free = ext[0] == "mf" if ext else False
    return Model(name_, states, make_settings_free)


def _get_state_names(state_nb):
    return string.ascii_lowercase[:state_nb]


def _check_model_name(name):
    if name not in make_settings:
        print("Warning: The 'model' option should either be:")
        for make_settings_name in make_settings:
            print(f"    - '{make_settings_name}'")
        print("Set to the default value: '2st'.")
        return "2st"
    return name


@ft.lru_cache()
def make_settings_2st(conditions, spin_system):
    return {
        "kex_ab": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 200.0,
            "min": 0.0,
            "max": np.inf,
            "vary": True,
        },
        "pb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 0.05,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "pa": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "1.0 - {pb}",
        },
        "kab": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": np.inf,
            "expr": "{kex_ab} * {pb}",
        },
        "kba": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": np.inf,
            "expr": "{kex_ab} * {pa}",
        },
    }


@ft.lru_cache()
def make_settings_2st_rs(conditions, spin_system):
    return {
        "kex_ab": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "value": 200.0,
            "min": 0.0,
            "max": np.inf,
            "vary": True,
        },
        "pb": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "value": 0.05,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "pa": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "1.0 - {pb}",
        },
        "kab": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": np.inf,
            "expr": "{kex_ab} * {pb}",
        },
        "kba": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": np.inf,
            "expr": "{kex_ab} * {pa}",
        },
    }


@ft.lru_cache()
def make_settings_2st_hd(conditions, spin_system):
    return {
        "kex_ab": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "value": 200.0,
            "min": 0.0,
            "vary": True,
        },
        "pb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 0.05,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "pa": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "1.0 - {pb}",
        },
        "kab": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "expr": "{kex_ab} * {pb}",
        },
        "kba": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "expr": "{kex_ab} * {pa}",
        },
    }


@ft.lru_cache()
def make_settings_2st_eyring(conditions, spin_system):
    celsius = conditions.temperature
    return {
        "dh_b": {"attributes": ("p_total", "l_total"), "value": 8e3, "vary": True},
        "ds_b": {"attributes": ("p_total", "l_total"), "value": 0e2, "vary": False},
        "dh_ab": {"attributes": ("p_total", "l_total"), "value": 6.5e4, "vary": True},
        "ds_ab": {"attributes": ("p_total", "l_total"), "value": 0e4, "vary": False},
        "pa": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": hs_to_p("a", "ab", celsius),
        },
        "pb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": hs_to_p("b", "ab", celsius),
        },
        "kab": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": hs_to_k("ab", celsius),
        },
        "kba": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": hs_to_k("ba", celsius),
        },
    }


@ft.lru_cache()
def make_settings_2st_binding(conditions, spin_system):
    p_total = conditions.p_total
    l_total = conditions.l_total
    delta = l_total - p_total
    expr_kab = (
        f"{{kon}} * 0.5 * ({delta} - {{kd}} "
        f"+ sqrt(({delta} - {{kd}}) ** 2 + 4.0 * {{kd}} * {l_total}))"
    )
    return {
        "kon": {
            "attributes": ("temperature",),
            "value": 1.0e7,
            "min": 0.0,
            "vary": True,
        },
        "koff": {
            "attributes": ("temperature",),
            "value": 10.0,
            "min": 0.0,
            "vary": True,
        },
        "kd": {
            "attributes": ("temperature",),
            "min": 0.0,
            "expr": "{koff} / max({kon}, 1e-100)",
        },
        "pa": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "{kba} / max({kba} + {kab}, 1e-100)",
        },
        "pb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "{kab} / max({kba} + {kab}, 1e-100)",
        },
        "kab": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": expr_kab,
        },
        "kba": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": "{koff}",
        },
    }


@ft.lru_cache()
def make_settings_3st(conditions, spin_system):
    return {
        "pb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 0.02,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "pc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 0.02,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "kex_ab": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 200.0,
            "vary": True,
        },
        "kex_ac": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 0.0,
            "vary": False,
        },
        "kex_bc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 200.0,
            "vary": True,
        },
        "pa": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "1.0 - {pb} - {pc}",
        },
        "kab": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ab"),
        },
        "kba": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ba"),
        },
        "kac": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ac"),
        },
        "kca": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ca"),
        },
        "kbc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("bc"),
        },
        "kcb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("cb"),
        },
    }


@ft.lru_cache()
def make_settings_3st_eyring(conditions, spin_system):
    celcius = conditions.temperature
    return {
        "dh_b": {"attributes": ("p_total", "l_total"), "value": 8e3, "vary": True},
        "dh_c": {"attributes": ("p_total", "l_total"), "value": 8e3, "vary": True},
        "dh_ab": {"attributes": ("p_total", "l_total"), "value": 6.5e4, "vary": True},
        "dh_bc": {"attributes": ("p_total", "l_total"), "value": 6.5e4, "vary": True},
        "dh_ac": {"attributes": ("p_total", "l_total"), "value": 1.0e10, "vary": False},
        "ds_b": {"attributes": ("p_total", "l_total"), "value": 0.0, "vary": False},
        "ds_c": {"attributes": ("p_total", "l_total"), "value": 0.0, "vary": False},
        "ds_ab": {"attributes": ("p_total", "l_total"), "value": 0.0, "vary": False},
        "ds_bc": {"attributes": ("p_total", "l_total"), "value": 0.0, "vary": False},
        "ds_ac": {"attributes": ("p_total", "l_total"), "value": 0.0, "vary": False},
        "kab": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": hs_to_k("ab", celcius),
        },
        "kba": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": hs_to_k("ba", celcius),
        },
        "kac": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": hs_to_k("ac", celcius),
        },
        "kca": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": hs_to_k("ca", celcius),
        },
        "kbc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": hs_to_k("bc", celcius),
        },
        "kcb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": hs_to_k("cb", celcius),
        },
        "pa": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": hs_to_p("a", "abc", celcius),
        },
        "pb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": hs_to_p("b", "abc", celcius),
        },
        "pc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": hs_to_p("c", "abc", celcius),
        },
    }


@ft.lru_cache()
def make_settings_4st(conditions, spin_system):
    return {
        "pb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 0.02,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "pc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 0.02,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "pd": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 0.02,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "kex_ab": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 200.0,
            "vary": True,
        },
        "kex_ac": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 0.0,
        },
        "kex_ad": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 0.0,
        },
        "kex_bc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 200.0,
            "vary": True,
        },
        "kex_bd": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 0.0,
        },
        "kex_cd": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 200.0,
            "vary": True,
        },
        "pa": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "1.0 - {pb} - {pc} - {pd}",
        },
        "kab": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ab"),
        },
        "kba": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ba"),
        },
        "kac": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ac"),
        },
        "kca": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ca"),
        },
        "kad": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ad"),
        },
        "kda": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("da"),
        },
        "kbc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "expr": kex_p_to_k("bc"),
        },
        "kcb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("cb"),
        },
        "kbd": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("bd"),
        },
        "kdb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("db"),
        },
        "kcd": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("cd"),
        },
        "kdc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("dc"),
        },
    }


@ft.lru_cache()
def make_settings_4st_hd(conditions, spin_system):
    return {
        "frac": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 0.1,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "pc": {
            "attributes": ("temperature", "p_total", "l_total"),
            "value": 0.02,
            "min": 0.0,
            "max": 1.0,
            "vary": True,
        },
        "pa": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "1.0 / (1.0 + {frac}) - {pc}",
        },
        "pb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "{frac} / (1.0 + {frac}) - {frac} * {pc}",
        },
        "pd": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "max": 1.0,
            "expr": "{frac} * {pc}",
        },
        "kex_ab": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 200.0,
            "vary": True,
        },
        "kex_ac": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 0.0,
            "vary": True,
        },
        "kex_bd": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": "{kex_ac}",
        },
        "kex_cd": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "min": 0.0,
            "value": 200.0,
            "vary": True,
        },
        "kab": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ab"),
        },
        "kba": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ba"),
        },
        "kac": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ac"),
        },
        "kca": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("ca"),
        },
        "kbd": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("bd"),
        },
        "kdb": {
            "attributes": ("temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("db"),
        },
        "kcd": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("cd"),
        },
        "kdc": {
            "attributes": ("spin_system", "temperature", "p_total", "l_total"),
            "min": 0.0,
            "expr": kex_p_to_k("dc"),
        },
    }


make_settings = {
    "2st": make_settings_2st,
    "3st": make_settings_2st,
    "4st": make_settings_2st,
    "2st_rs": make_settings_2st_rs,
    "2st_hd": make_settings_2st_hd,
    "2st_eyring": make_settings_2st_eyring,
    "3st_eyring": make_settings_3st_eyring,
    "2st_binding": make_settings_2st_binding,
    "4st_hd": make_settings_4st_hd,
}


def kex_p_to_k(states):
    kex = f"{{kex_{min(states)}{max(states)}}}"
    p1 = f"{{p{states[0]}}}"
    p2 = f"{{p{states[1]}}}"
    return f"{kex} * {p2} / max({p1} + {p2}, 1e-100)"


def hs_to_k(states, celcius):
    kelvin = celcius + 273.15
    kbt_h = cst.k * kelvin / cst.h
    rt = cst.R * kelvin
    expr_dh = f"{{dh_{min(states)}{max(states)}}}"
    expr_ds = f"{{ds_{min(states)}{max(states)}}}"
    if states[0] != "a":
        expr_dh += f" - {{dh_{states[0]}}}"
        expr_ds += f" - {{ds_{states[0]}}}"
    return f"{kbt_h} * exp(-(({expr_dh}) - {kelvin} * ({expr_ds})) / {rt})"


def hs_to_p(state, states, celcius):
    kelvin = celcius + 273.15
    rt = cst.R * kelvin
    dg = {}
    for a_state in states:
        if a_state != "a":
            dg[a_state] = (
                f"exp(" f"-({{dh_{a_state}}} - {kelvin} * {{ds_{a_state}}}) / {rt})"
            )
        else:
            dg[a_state] = "1.0"
    return f"{dg[state]} / ({' + '.join(dg.values())})"

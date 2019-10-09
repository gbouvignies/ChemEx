import collections
import re
import sys

import jsonschema as js
import numpy as np
import scipy.constants as cst

from chemex.parameters import helper as cph


Model = collections.namedtuple("Model", ["name", "state_nb", "kind"])


def parse_model(name):
    name_ = _check_model_name(name)
    match = re.match(r"(\d)st\.(\w+)", name_, re.IGNORECASE)
    if not match:
        raise NameError(f"Impossible to parse the model name '{name_}'.")
    state_nb = int(match.group(1))
    kind = match.group(2)
    return Model(name_, state_nb, kind)


def _check_model_name(name):
    if name not in models:
        print("Warning: The 'model' option should either be:")
        for model_name in models:
            print(f"    - '{model_name}'")
        print("Set to the default value: '2st.pb_kex'.")
        return "2st.pb_kex"
    return name


def create_params_k(model, conditions, spin_system=None):
    """Update the experimental and fitting parameters depending on the model."""
    fnames, params = models[model.name](conditions, spin_system)
    return fnames, params


def model_2st_pb_kex(conditions, spin_system):
    settings = {
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
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def model_2st_pb_kex_rs(conditions, spin_system):
    settings = {
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
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def model_2st_hd_exch(conditions, spin_system):
    settings = {
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
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def model_2st_eyring(conditions, spin_system):
    celsius = conditions.temperature
    settings = {
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
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def model_2st_binding(conditions, spin_system):
    p_total = conditions.p_total
    l_total = conditions.l_total
    delta = l_total - p_total
    expr_kab = (
        f"{{kon}} * 0.5 * ({delta} - {{kd}} "
        f"+ sqrt(({delta} - {{kd}}) ** 2 + 4.0 * {{kd}} * {l_total}))"
    )
    settings = {
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
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def model_3st_pb_kex(conditions, spin_system):
    settings = {
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
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def model_3st_eyring(conditions, spin_system):
    celcius = conditions.temperature
    settings = {
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
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def model_4st_pb_kex(conditions, spin_system):
    settings = {
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
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def model_4st_hd_exch(conditions, spin_system):
    settings = {
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
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


models = {
    "2st.pb_kex": model_2st_pb_kex,
    "2st.pb_kex_rs": model_2st_pb_kex_rs,
    "2st.hd_exch": model_2st_hd_exch,
    "3st.pb_kex": model_3st_pb_kex,
    "4st.pb_kex": model_4st_pb_kex,
    "2st.eyring": model_2st_eyring,
    "3st.eyring": model_3st_eyring,
    "2st.binding": model_2st_binding,
    "4st.hd_exch": model_4st_hd_exch,
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
                f"exp(" f"-({{dh_{a_state}}} - {kelvin} * {{ds_{a_state}}}) " f"/ {rt})"
            )
        else:
            dg[a_state] = "1.0"
    return f"{dg[state]} / ({' + '.join(dg.values())})"


def validates_conditions(config):
    model = config["model"]
    _schema = {
        "type": "object",
        "properties": {
            "conditions": {
                "type": "object",
                "properties": {
                    "h_larmor_frq": {"type": "number"},
                    "temperatures": {"type": "number"},
                    "p_total": {"type": "number"},
                    "l_total": {"type": "number"},
                },
                "dependencies": {"p_total": ["l_total"], "l_total": ["p_total"]},
                "required": ["h_larmor_frq"],
            }
        },
        "required": ["conditions"],
    }
    if "binding" in model:
        _schema["properties"]["conditions"]["required"].extend(["p_total", "l_total"])

    if "eyring" in model:
        _schema["properties"]["conditions"]["required"].append("temperature")

    try:
        js.validate(config, _schema)
    except js.ValidationError as e:
        filename = config["filename"]
        if len(e.path) == 1:
            sys.exit(
                f"\nerror: The experiment file '{filename}' has no section "
                f"'[conditions]'."
            )
        else:
            sys.exit(
                f"\nerror: The model '{model}' requires the condition '{e.instance}' "
                f"to be defined in the experiment file. '{e.instance}' is missing from"
                f"'{filename}'."
            )

import itertools as it

import numpy as np

import chemex.parameters.helper as cph


def create_params_l(basis, model, conditions, spin_system=None):
    settings = _make_settings(basis, model)
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def _make_settings(basis, model):
    shortnames = _get_shortnames(basis)
    settings = {}
    for shortname, state in it.product(shortnames, model.states):
        name = f"{shortname}_{state}"
        settings[name] = SETTINGS[shortname].copy()
        if state != "a":
            settings[name]["expr"] = f"{{{shortname}_a}}"
        if (
            basis.extension
            and basis.extension in {"dq", "tq"}
            and shortname in _TO_BE_EXTENDED
        ):
            settings[name]["ext"] = basis.extension
    if basis.spin_system in set_constraints:
        settings = set_constraints[basis.spin_system](settings, model.states)
    settings = _add_dw(settings, model.states)
    return settings


def _get_shortnames(basis):
    shortnames = BASIS_TO_SHORTNAMES[basis.type]
    if basis.spin_system in ("hn", "cn", "hc") and "r1_s" not in shortnames:
        shortnames.append("r1_s")
    return shortnames


def _add_dw(settings, states):
    settings_ = settings.copy()
    for spin, state in it.product("is", states[1:]):
        cs_name = f"cs_{spin}_{state}"
        if cs_name in settings:
            name = f"dw_{spin}_a{state}"
            settings_[name] = SETTINGS[f"cs_{spin}"].copy()
            settings_[cs_name]["expr"] = f"{{cs_{spin}_a}} + {{dw_{spin}_a{state}}}"
    return settings_


def _set_constraints_hn(settings, states):
    settings_ = settings.copy()
    for state in states:
        r1_s = f"r1_s_{state}"
        r1_i = f"r1_i_{state}"
        r2_i = f"r2_i_{state}"
        r1a_is = f"r1a_is_{state}"
        r2a_i = f"r2a_i_{state}"
        if {r1_i, r1a_is, r1_s} <= set(settings_):
            settings_[r1_i]["expr"] = f"{{{r1a_is}}} - {{{r1_s}}}"
        if {r2a_i, r2_i, r1_s} <= set(settings_):
            settings_[r2a_i]["expr"] = f"{{{r2_i}}} - {{{r1_s}}}"
    return settings_


def _set_constraints_nh(settings, states):
    settings_ = settings.copy()
    for state in states:
        r1_i = f"r1_i_{state}"
        r2_i = f"r2_i_{state}"
        r2_s = f"r2_s_{state}"
        r1a_is = f"r1a_is_{state}"
        r2a_i = f"r2a_i_{state}"
        r2a_s = f"r2a_s_{state}"
        if {r2a_i, r2_i, r1_i} <= set(settings_):
            settings_[r2a_i]["expr"] = f"{{{r2_i}}} + {{{r1a_is}}} - {{{r1_i}}}"
        if {r2a_s, r2_s, r1_i} <= set(settings_):
            settings_[r2a_s]["expr"] = f"{{{r2_s}}} - {{{r1_i}}}"
    return settings


def _set_constraints_hc(settings, states):
    settings_ = settings.copy()
    for state in states:
        r1_s = f"r1_s_{state}"
        r1_i = f"r1_i_{state}"
        r2_i = f"r2_i_{state}"
        r1a_is = f"r1a_is_{state}"
        r2a_i = f"r2a_i_{state}"
        if {r1a_is, r1_i, r1_s} <= set(settings_):
            settings_[r1a_is]["expr"] = f"{{{r1_i}}} + {{{r1_s}}}"
        if {r2a_i, r2_i, r1_s} <= set(settings_):
            settings_[r2a_i]["expr"] = f"{{{r2_i}}} - {{{r1_s}}}"
    return settings_


set_constraints = {
    "hn": _set_constraints_hn,
    "cn": _set_constraints_hn,
    "nh": _set_constraints_nh,
    "ch": _set_constraints_nh,
    "hc": _set_constraints_hn,
}


SETTINGS = {
    "r2_i": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 10.0,
        "min": 0.0,
        "max": np.inf,
        "vary": False,
    },
    "r2_s": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 10.0,
        "min": 0.0,
        "max": np.inf,
        "vary": False,
    },
    "r1_i": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 1.5,
        "min": 0.0,
        "max": np.inf,
        "vary": False,
    },
    "r1_s": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 1.5,
        "min": 0.0,
        "max": np.inf,
        "vary": False,
    },
    "r2a_i": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 12.0,
        "min": 0.0,
        "max": np.inf,
        "vary": False,
    },
    "r2a_s": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 12.0,
        "min": 0.0,
        "max": np.inf,
        "vary": False,
    },
    "r2mq_is": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 15.0,
        "min": 0.0,
        "max": np.inf,
        "vary": False,
    },
    "r1a_is": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 4.0,
        "min": 0.0,
        "max": np.inf,
        "vary": False,
    },
    "etaxy_i": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 0.0,
        "min": -np.inf,
        "max": np.inf,
        "vary": False,
    },
    "etaxy_s": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 0.0,
        "min": -np.inf,
        "max": np.inf,
        "vary": False,
    },
    "etaz_i": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 0.0,
        "min": -np.inf,
        "max": np.inf,
        "vary": False,
    },
    "etaz_s": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 0.0,
        "min": -np.inf,
        "max": np.inf,
        "vary": False,
    },
    "sigma_is": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 0.0,
        "min": -np.inf,
        "max": np.inf,
        "vary": False,
    },
    "mu_is": {
        "attributes": ("spin_system", "temperature", "h_larmor_frq"),
        "value": 0.0,
        "min": -np.inf,
        "max": np.inf,
        "vary": False,
    },
    "cs_i": {
        "attributes": ("spin_system", "temperature"),
        "value": 0.0,
        "min": -np.inf,
        "max": np.inf,
        "vary": False,
    },
    "cs_s": {
        "attributes": ("spin_system", "temperature"),
        "value": 0.0,
        "min": -np.inf,
        "max": np.inf,
        "vary": False,
    },
    "j_is": {
        "attributes": ("spin_system",),
        "value": 0.0,
        "min": -np.inf,
        "max": np.inf,
        "vary": False,
    },
    "d": {
        "attributes": ("temperature",),
        "value": 0.0,
        "min": 0.0,
        "max": 1.0,
        "vary": False,
    },
}


BASIS_TO_SHORTNAMES = {
    "ixy": ["r2_i", "cs_i"],
    "iz": ["r1_i"],
    "izsz": ["r1_i", "r1a_is"],
    "iz_eq": ["r1_i"],
    "ixyz": ["r2_i", "r1_i", "cs_i"],
    "ixyz_eq": ["r2_i", "r1_i", "cs_i"],
    "ixysxy": ["r2mq_is", "mu_is", "cs_i", "cs_s"],
    "ixy_ixysxy": ["r2_i", "r2mq_is", "mu_is", "cs_i", "cs_s"],
    "ixyzsz": ["r2_i", "r1_i", "r2a_i", "etaxy_i", "etaz_i", "r1a_is", "cs_i", "j_is"],
    "ixyzsz_dif": [
        "r2_i",
        "r1_i",
        "r2a_i",
        "etaxy_i",
        "etaz_i",
        "r1a_is",
        "cs_i",
        "j_is",
        "d",
    ],
    "ixyzsz_eq": [
        "r2_i",
        "r1_i",
        "r2a_i",
        "etaxy_i",
        "etaz_i",
        "r1a_is",
        "cs_i",
        "j_is",
    ],
    "ixyzsxyz": [
        "r2_i",
        "r2_s",
        "r1_i",
        "r1_s",
        "r2a_i",
        "r2a_s",
        "r2mq_is",
        "r1a_is",
        "etaxy_i",
        "etaxy_s",
        "etaz_i",
        "etaz_s",
        "sigma_is",
        "mu_is",
        "cs_i",
        "cs_s",
        "j_is",
    ],
    "ixyzsxyz_eq": [
        "r2_i",
        "r2_s",
        "r1_i",
        "r1_s",
        "r2a_i",
        "r2a_s",
        "r2mq_is",
        "r1a_is",
        "etaxy_i",
        "etaxy_s",
        "etaz_i",
        "etaz_s",
        "sigma_is",
        "mu_is",
        "cs_i",
        "cs_s",
        "j_is",
    ],
}

_TO_BE_EXTENDED = ["r2_i", "r1_i", "r2a_i", "etaxy_i", "etaz_i", "r1a_is"]

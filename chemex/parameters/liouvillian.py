import itertools as it

import numpy as np

import chemex.nmr.spin_system as cns
import chemex.parameters.helper as cph


def create_params_l(basis, model, conditions, spin_system=None, constraints=None):
    if basis not in BASIS_TO_SHORTNAMES:
        message = "The 'basis' option should either be:\n"
        for basis_name in BASIS_TO_SHORTNAMES:
            message += f"    - '{basis_name}'\n"
        raise NameError(message)
    states = cns.get_state_names(model.state_nb)
    shortnames = _get_shortnames(basis, constraints)
    settings = _make_settings(shortnames, states, constraints)
    fnames, params = cph.make_params(settings, conditions, spin_system)
    return fnames, params


def _get_shortnames(basis, constraints=None):
    if constraints is None:
        constraints = []
    shortnames = BASIS_TO_SHORTNAMES[basis]
    for constraint in constraints:
        if constraint in ("hn",) and "r1_s" not in shortnames:
            shortnames.append("r1_s")
    return shortnames


def _make_settings(shortnames, states, constraints=None):
    if constraints is None:
        constraints = []
    settings = {}
    for shortname, state in it.product(shortnames, states):
        name = f"{shortname}_{state}"
        settings[name] = SETTINGS[shortname].copy()
        if state != "a":
            settings[name]["expr"] = f"{{{shortname}_a}}"
    for constraint in constraints:
        if constraint in set_constraints:
            settings = set_constraints[constraint](settings, states)
    settings = _add_dw(settings, states)
    return settings


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


def _set_constraints_ch3_1htq(settings, states):
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
    "nh": _set_constraints_nh,
    "ch3_1htq": _set_constraints_ch3_1htq,
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
    "iz_eq": ["r1_i"],
    "ixyz": ["r2_i", "r1_i", "cs_i"],
    "ixyz_eq": ["r2_i", "r1_i", "cs_i"],
    "ixysxy": ["r2mq_is", "mu_is", "cs_i"],
    "ixyzsz": ["r2_i", "r1_i", "r2a_i", "etaxy_i", "etaz_i", "r1a_is", "cs_i", "j_is"],
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
    "ch3_1htq": ["r2_i", "r1_i", "r1_s", "r2a_i", "r1a_is", "cs_i", "j_is"],
    "ch3_1htq_grad": ["r2_i", "r1_i", "r1_s", "r2a_i", "r1a_is", "cs_i", "j_is", "d"],
}

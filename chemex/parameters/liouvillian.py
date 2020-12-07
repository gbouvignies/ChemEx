import copy
import functools as ft
import itertools as it

import chemex.nmr.constants as cnc
import chemex.nmr.rates as cnr


@ft.lru_cache()
def make_settings(basis, model, conditions):
    settings = {}
    extension = basis.extension if basis.extension in {"dq", "tq"} else ""
    spin_system = basis.spin_system
    for state in model.states:
        settings.update(
            {
                f"tauc_{state}": {
                    "attributes": ("spin_system", "temperature"),
                    "value": 4.0,
                    "min": 0.0,
                },
                f"s2_{state}": {
                    "attributes": ("spin_system", "temperature"),
                    "value": 0.9,
                    "min": 0.0,
                },
                f"khh_{state}": {
                    "attributes": ("spin_system", "temperature"),
                    "value": 0.0,
                },
                f"r2_i_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 10.0,
                    "min": 0.0,
                    "ext": extension,
                },
                f"r2_s_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 10.0,
                    "min": 0.0,
                },
                f"r1_i_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 1.5,
                    "min": 0.0,
                    "ext": extension,
                },
                f"r1_s_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 1.5,
                    "min": 0.0,
                },
                f"r2a_i_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 12.0,
                    "min": 0.0,
                    "ext": extension,
                    "expr": f"{{r2_i_{state}}} - {{r1_s_{state}}}",
                },
                f"r2a_s_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 12.0,
                    "min": 0.0,
                    "expr": f"{{r2_s_{state}}} - {{r1_i_{state}}}",
                },
                f"r2mq_is_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 15.0,
                    "min": 0.0,
                },
                f"r1a_is_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 4.0,
                    "min": 0.0,
                    "ext": extension,
                    "expr": f"{{r1_i_{state}}} + {{r1_s_{state}}}",
                },
                f"etaxy_i_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 0.0,
                    "ext": extension,
                },
                f"etaxy_s_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 0.0,
                },
                f"etaz_i_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 0.0,
                    "ext": extension,
                },
                f"etaz_s_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 0.0,
                },
                f"sigma_is_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 0.0,
                },
                f"mu_is_{state}": {
                    "attributes": ("spin_system", "temperature", "h_larmor_frq"),
                    "value": 0.0,
                },
                f"cs_i_{state}": {
                    "attributes": ("spin_system",),
                    "value": 0.0,
                },
                f"cs_s_{state}": {
                    "attributes": ("spin_system",),
                    "value": 0.0,
                },
                f"j_is_{state}": {
                    "attributes": ("spin_system",),
                    "value": cnc.J_COUPLINGS.get(spin_system, 0.0),
                },
                f"d_{state}": {"attributes": ("temperature",), "value": 0.0},
            }
        )
        # As 1HN sites exchange with water, the "expr" should be modified accordingly
        if spin_system == "hn":
            r2a_s = f"{{r2_s_{state}}} + {{r1a_is_{state}}} - {{r1_s_{state}}}"
            settings[f"r2a_s_{state}"]["expr"] = r2a_s
        elif spin_system == "nh":
            r2a_i = f"{{r2_i_{state}}} + {{r1a_is_{state}}} - {{r1_i_{state}}}"
            settings[f"r2a_i_{state}"]["expr"] = r2a_i
            settings[f"r1a_is_{state}"]["expr"] = ""
    if model.name.startswith("4st_hd"):
        settings = _add_dw_hd(settings, model)
    else:
        settings = _add_dw(settings, model)
    settings = _set_equal_to_a(settings)
    settings_mf = _set_model_free(settings, basis, conditions)
    return settings, settings_mf


def _set_equal_to_a(settings: dict) -> dict:
    result = copy.deepcopy(settings)
    for k, v in result.items():
        if k[-1] != "a" and not v.get("expr"):
            expr = k[:-1] + "a"
            if expr in result:
                v["expr"] = f"{{{expr}}}"
    return result


def _add_dw(settings: dict, model) -> dict:
    result = copy.deepcopy(settings)
    for spin, state in it.product("is", model.states[1:]):
        cs_name = f"cs_{spin}_{state}"
        dw_name = f"dw_{spin}_a{state}"
        result[dw_name] = copy.deepcopy(settings[cs_name])
        result[dw_name]["vary"] = True
        result[cs_name]["expr"] = f"{{cs_{spin}_a}} + {{{dw_name}}}"
    return result


def _add_dw_hd(settings: dict, model) -> dict:
    result = copy.deepcopy(settings)
    result["dw_i_ab"] = copy.deepcopy(settings["cs_i_a"])
    result["dwhd_i_a"] = copy.deepcopy(settings["cs_i_a"])
    result["dwhd_i_b"] = copy.deepcopy(settings["cs_i_a"])
    result["dw_i_ab"]["vary"] = True
    result["dwhd_i_a"]["vary"] = True
    result["dwhd_i_b"]["vary"] = True
    result["cs_i_b"]["expr"] = "{cs_i_a} + {dw_i_ab}"
    result["cs_i_c"]["expr"] = "{cs_i_a} + {dwhd_i_a}"
    result["cs_i_d"]["expr"] = "{cs_i_a} + {dw_i_ab} + {dwhd_i_b}"
    return result


def _set_model_free(settings: dict, basis, conditions) -> dict:
    result = copy.deepcopy(settings)
    spin_system = basis.spin_system
    if "2H" in conditions.label:
        spin_system += "_d"
    rate_constraints = cnr.RATE_CONSTRAINTS.get(spin_system)
    if rate_constraints is None:
        return result
    for rate in set(result) & set(rate_constraints):
        result[rate]["expr"] = rate_constraints[rate]
    return result

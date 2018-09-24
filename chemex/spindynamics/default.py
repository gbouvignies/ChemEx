import itertools
import re
import string

import lmfit as lf
from scipy import constants as cst

from chemex import parameters

OTHER_SPIN = {"i": "s", "s": "i"}

STATES = "abcd"

RE_VARNAME = re.compile(
    "^(?P<name>.+?)(_(?P<spin>[is]{1,2}))?(_(?P<state>[abcd]{1,2}))?$"
)


def create_params(
    basis,
    model,
    nuclei=None,
    conditions=None,
    hn_ap_constraints=False,
    nh_constraints=False,
):
    """TODO"""

    fnames_k, params_k = create_params_k(model=model, conditions=conditions)
    fnames_lr, params_lr = create_params_lr(basis, conditions, nuclei)

    params_lr = set_thermal_factors(fnames_lr, params_lr, nuclei)

    if hn_ap_constraints:
        fnames_lr, params_lr = set_hn_ap_constraints(
            fnames_lr, params_lr, nuclei=nuclei, conditions=conditions
        )

    if nh_constraints:
        fnames_lr, params_lr = set_nh_constraints(fnames_lr, params_lr)

    fnames = dict(**fnames_k, **fnames_lr)

    # hack to add parameters depending on parameters not added yet
    params = lf.Parameters()
    params.add_many(*(lf.Parameter(fname) for fname in fnames.values()))
    params.update(params_lr)
    params.update(params_k)

    return fnames, params


def create_params_lr(basis, conditions, nuclei):
    names = get_default_names(basis)
    fnames = set_fnames(names, nuclei, conditions)
    params = buid_params(fnames)
    for name, fname in fnames.items():
        if name.endswith(("_b", "_c", "_d")):
            name_a = "".join((name[:-1], "a"))
            fname_a = fnames[name_a]
            params[fname].set(expr=fname_a)
    for state in "bcd":
        short_name = "cs_{}".format(state)
        short_name_expr = "cs_a + dw_a{}".format(state)
        parameters.set_param_expr(params, short_name, short_name_expr)
    return fnames, params


def get_default_names(basis):
    """TODO"""

    names = set()
    for name in basis._matrices:

        name_parsed = RE_VARNAME.match(name).groupdict()

        if not name.startswith(("k", "p", "w1", "carrier", "j_eff", "w_")):
            names.add(name)

            if name_parsed["name"] == "cs" and name_parsed["state"] != "a":
                name_dw = "dw_{spin}_a{state}".format(**name_parsed)
                names.add(name_dw)

    return names


def set_fnames(names=None, nuclei=None, conditions=None):
    """TODO"""

    if names is None:
        names = {}

    fnames = {}

    for name in names:
        name_parsed = RE_VARNAME.match(name).groupdict()

        short_name = "_".join(
            name_parsed[key]
            for key in ("name", "state")
            if name_parsed[key] is not None
        )

        definition = {"temperature": conditions.get("temperature", 0.0)}

        spin = name_parsed["spin"]

        if spin is None:
            spin = "is"

        definition["nuclei"] = nuclei[spin]

        if not name.startswith(("dw", "cs", "j")):
            definition["h_larmor_frq"] = conditions.get("h_larmor_frq", 1.0)

        fnames[name] = parameters.ParameterName(short_name, **definition).to_full_name()

    return fnames


def buid_params(fnames):
    """TODO"""

    params = lf.Parameters()

    for name, full_name in fnames.items():
        param = lf.Parameter(full_name, value=0.0, vary=False)

        if name.startswith("p"):
            param.set(min=0.0, max=1.0)
        elif name.startswith("k"):
            param.set(min=0.0)
        elif name.startswith("r2"):
            param.set(value=10.0, min=0.0)
        elif name.startswith("r1"):
            param.set(value=1.0, min=0.0)
        elif name == "j":
            param.set(value=-93.0)

        params[full_name] = param

    return params


def set_thermal_factors(fnames, params, nuclei):
    """TODO"""

    for spin, state in itertools.product("is", STATES):

        theta = "theta_{}_{}".format(spin, state)

        if theta in fnames:

            r1 = "r1_{}_{}".format(spin, state)
            sigma = "sigma_{}".format(state)
            pop = "p{}".format(state)

            terms = []

            terms.append(
                "{g_ratio} * {r1}".format(
                    g_ratio=nuclei.atoms[spin]["g_ratio"], r1=fnames[r1]
                )
            )

            if sigma in fnames:
                terms.append(
                    "{g_ratio} * {sigma}".format(
                        g_ratio=nuclei.atoms[OTHER_SPIN[spin]]["g_ratio"],
                        sigma=fnames[sigma],
                    )
                )

            expr = "{pop} * ({sum})".format(pop=fnames[pop], sum=" + ".join(terms))

            params[fnames[theta]].set(expr=expr)

    for state in STATES:

        theta = "theta_is_{}".format(state)

        if theta in fnames:

            pop = "p{}".format(state)

            terms = []

            for spin in "is":

                etaz = "etaz_{}_{}".format(spin, state)

                if etaz in fnames:
                    terms.append(
                        "{gamma} * {etaz}".format(
                            gamma=nuclei.atoms[spin]["g_ratio"], etaz=fnames[etaz]
                        )
                    )

            expr = "{pop} * ({sum})".format(pop=fnames[pop], sum=" + ".join(terms))

            params[fnames[theta]].set(expr=expr)

    return params


def set_hn_ap_constraints(fnames, params, nuclei=None, conditions=None):

    name_s = nuclei.get("s")

    if name_s is None:
        return fnames, params

    attributes = {"nuclei": name_s}
    attributes.update(conditions)

    names = (
        "r1_s_{}".format(state) for state in STATES if "r1_i_{}".format(state) in fnames
    )
    fnames, params = _add_params(fnames, params, names, attributes)

    expr_r1_i = "{{r1a_{state}}} - {{r1_s_{state}}}"
    expr_r2a_i = "{{r2_i_{state}}} - {{r1_s_{state}}}"

    settings = {}

    for state in STATES:

        r1_s = "r1_s_{}".format(state)
        r1_i = "r1_i_{}".format(state)
        r2a_i = "r2a_i_{}".format(state)

        if r1_s not in fnames:
            continue

        expr_r1_s = None
        if state != "a":
            expr_r1_s = "{r1a_a}".format(**fnames)

        settings[r1_s] = {"value": 1.0, "min": 0.0, "vary": False, "expr": expr_r1_s}
        settings[r1_i] = {"expr": expr_r1_i.format(state=state).format(**fnames)}
        settings[r2a_i] = {"expr": expr_r2a_i.format(state=state).format(**fnames)}

    fnames, params = _apply_settings(fnames, params, settings)

    return fnames, params


def set_nh_constraints(fnames, params):

    expr_r2a_i = "{{r2_i_{state}}} + {{r1a_{state}}} - {{r1_i_{state}}}"
    expr_r2a_s = "{{r2_s_{state}}} - {{r1_i_{state}}}"

    settings = {}

    settings.update(
        {
            "r2a_i_{}".format(state): {
                "expr": expr_r2a_i.format(state=state).format(**fnames)
            }
            for state in STATES
            if "r2a_i_{}".format(state) in fnames
        }
    )

    settings.update(
        {
            "r2a_s_{}".format(state): {
                "expr": expr_r2a_s.format(state=state).format(**fnames)
            }
            for state in STATES
            if "r2a_s_{}".format(state) in fnames
        }
    )

    fnames, params = _apply_settings(fnames, params, settings)

    return fnames, params


def create_params_k(model=None, conditions=None):
    """Update the experimental and fitting parameters depending on the model.

    """

    fnames = {}
    params = lf.Parameters()

    attributes = {
        key: conditions.get(key) for key in ("temperature", "p_total", "l_total")
    }

    states = string.ascii_lowercase[: model.state_nb]
    names = ["k{}{}".format(*pairs) for pairs in itertools.permutations(states, 2)]
    names.extend(["p{}".format(state) for state in states])

    fnames, params = _add_params(fnames, params, names, attributes)

    kinetic_models = {
        "2st.pb_kex": model_2st_pb_kex,
        "3st.pb_kex": model_3st_pb_kex,
        "4st.pb_kex": model_4st_pb_kex,
        "2st.eyring": model_2st_eyring,
        "3st.eyring": model_3st_eyring,
        "2st.binding": model_2st_binding,
    }

    kinetic_model = kinetic_models.get(model.name)

    if kinetic_model is None:
        print("Warning: The 'model' option should either be:")
        for model in sorted(kinetic_models):
            print("    - '{}'".format(model))
        print("Set to the default value: '2st.pb_kex'.")
        kinetic_model = kinetic_models["2st.pb_kex"]

    fnames, params = kinetic_model(conditions, fnames, params)

    return fnames, params


def model_2st_pb_kex(conditions, fnames, params):

    attributes = {
        key: conditions.get(key) for key in ("temperature", "p_total", "l_total")
    }

    names = ("kex_ab",)

    fnames, params = _add_params(fnames, params, names, attributes)

    settings = {
        "kex_ab": {"value": 200.0, "min": 0.0, "vary": True},
        "pa": {"expr": "1.0 - {pb}".format(**fnames)},
        "pb": {"value": 0.05, "min": 0.0, "max": 1.0, "vary": True},
        "kab": {"expr": "{kex_ab} * {pb}".format(**fnames)},
        "kba": {"expr": "{kex_ab} * {pa}".format(**fnames)},
    }

    fnames, params = _apply_settings(fnames, params, settings)

    return fnames, params


def model_2st_eyring(conditions, fnames, params):

    attributes = {key: conditions.get(key) for key in ("p_total", "l_total")}
    names = ("dh_b", "ds_b", "dh_ab", "ds_ab")
    fnames, params = _add_params(fnames, params, names, attributes)

    thermo = _define_thermo(conditions.get("temperature"))

    settings = {
        "dh_b": {"value": 6.5e+03, "vary": True},
        "ds_b": {"value": 1.0e+10},
        "dh_ab": {"value": 6.5e+04, "vary": True},
        "ds_ab": {"value": 1.0e+10},
        "pa": {"expr": "{kba} / ({kba} + {kab})".format(**fnames)},
        "pb": {"expr": "{kab} / ({kba} + {kab})".format(**fnames)},
        "kab": {"expr": _get_k_from_h_s("ab").format(**fnames, **thermo)},
        "kba": {"expr": _get_k_from_h_s("ba").format(**fnames, **thermo)},
    }

    fnames, params = _apply_settings(fnames, params, settings)

    return fnames, params


def model_2st_binding(conditions, fnames, params):

    attributes = {key: conditions.get(key) for key in ("temperature",)}
    names = ("kd", "kon", "koff")
    fnames, params = _add_params(fnames, params, names, attributes)

    p_total = conditions.get("p_total", 0.0)
    l_total = conditions.get("l_total", 0.0)
    delta = l_total - p_total
    extra = {"delta": delta, "l_total": l_total}

    expr_kab = "{kon} * 0.5 * ({delta} - {kd} + sqrt(({delta} - {kd}) ** 2 + 4.0 * {kd} * {l_total}))"

    settings = {
        "kon": {"value": 10.0, "vary": True},
        "koff": {"value": 1.0e7, "vary": True},
        "kd": {"expr": "{koff} / {kon}".format(**fnames)},
        "pa": {"expr": "{kba} / ({kba} + {kab})".format(**fnames)},
        "pb": {"expr": "{kab} / ({kba} + {kab})".format(**fnames)},
        "kab": {"expr": expr_kab.format(**fnames, **extra)},
        "kba": {"expr": "{koff}".format(**fnames)},
    }

    fnames, params = _apply_settings(fnames, params, settings)

    return fnames, params


def model_3st_pb_kex(conditions, fnames, params):

    attributes = {
        key: conditions.get(key) for key in ("temperature", "p_total", "l_total")
    }
    names = ("kex_ab", "kex_ac", "kex_bc")
    fnames, params = _add_params(fnames, params, names, attributes)

    settings = {
        "pa": {"expr": "1.0 - {pb} - {pc}".format(**fnames)},
        "pb": {"value": 0.02, "min": 0.0, "max": 1.0, "vary": True},
        "pc": {"value": 0.02, "min": 0.0, "max": 1.0, "vary": True},
        "kex_ab": {"min": 0.0, "value": 200.0, "vary": True},
        "kex_ac": {"min": 0.0, "value": 0.0},
        "kex_bc": {"min": 0.0, "value": 200.0, "vary": True},
        "kab": {"expr": _get_k_from_kex_p("ab").format(**fnames)},
        "kba": {"expr": _get_k_from_kex_p("ba").format(**fnames)},
        "kac": {"expr": _get_k_from_kex_p("ac").format(**fnames)},
        "kca": {"expr": _get_k_from_kex_p("ca").format(**fnames)},
        "kbc": {"expr": _get_k_from_kex_p("bc").format(**fnames)},
        "kcb": {"expr": _get_k_from_kex_p("cb").format(**fnames)},
    }

    fnames, params = _apply_settings(fnames, params, settings)

    return fnames, params


def model_3st_eyring(conditions, fnames, params):

    attributes = {key: conditions.get(key) for key in ("p_total", "l_total")}
    names = ["dh_b", "dh_c", "dh_ab", "dh_ac", "dh_bc"]
    names.extend(["ds_b", "ds_c", "ds_ab", "ds_ac", "ds_bc"])
    fnames, params = _add_params(fnames, params, names, attributes)

    thermo = _define_thermo(conditions.get("temperature"))

    expr_pa = "{kba} * {kca} / ({kba} * {kca} + {kab} * {kca} + {kac} * {kba})"
    expr_pb = "{kab} * {kcb} / ({kab} * {kcb} + {kba} * {kcb} + {kbc} * {kab})"
    expr_pc = "{kbc} * {kac} / ({kbc} * {kac} + {kcb} * {kab} + {kca} * {kbc})"

    settings = {
        "dh_b": {"value": 6.5e+03, "vary": True},
        "dh_c": {"value": 6.5e+03, "vary": True},
        "dh_ab": {"value": 6.5e+04, "vary": True},
        "dh_bc": {"value": 6.5e+04, "vary": True},
        "dh_ac": {"value": 1.0e+10},
        "ds_b": {"value": 1.0e+10},
        "ds_c": {"value": 1.0e+10},
        "ds_ab": {"value": 1.0e+10},
        "ds_bc": {"value": 1.0e+10},
        "ds_ac": {"value": 1.0e+10},
        "kab": {"expr": _get_k_from_h_s("ab").format(**fnames, **thermo)},
        "kba": {"expr": _get_k_from_h_s("bs").format(**fnames, **thermo)},
        "kac": {"expr": _get_k_from_h_s("ac").format(**fnames, **thermo)},
        "kca": {"expr": _get_k_from_h_s("ca").format(**fnames, **thermo)},
        "kbc": {"expr": _get_k_from_h_s("bc").format(**fnames, **thermo)},
        "kcb": {"expr": _get_k_from_h_s("cb").format(**fnames, **thermo)},
        "pa": {"expr": expr_pa.format(**fnames)},
        "pb": {"expr": expr_pb.format(**fnames)},
        "pc": {"expr": expr_pc.format(**fnames)},
    }

    fnames, params = _apply_settings(fnames, params, settings)

    return fnames, params


def model_4st_pb_kex(conditions, fnames, params):

    attributes = {
        key: conditions.get(key) for key in ("temperature", "p_total", "l_total")
    }
    names = ("kex_ab", "kex_ac", "kex_ad", "kex_bc", "kex_bd", "kex_cd")
    fnames, params = _add_params(fnames, params, names, attributes)

    settings = {
        "pa": {"expr": "1.0 - {pb} - {pc} - {pd}".format(**fnames)},
        "pb": {"value": 0.02, "min": 0.0, "max": 1.0, "vary": True},
        "pc": {"value": 0.02, "min": 0.0, "max": 1.0, "vary": True},
        "pd": {"value": 0.02, "min": 0.0, "max": 1.0, "vary": True},
        "kex_ab": {"min": 0.0, "value": 200.0, "vary": True},
        "kex_ac": {"min": 0.0, "value": 0.0},
        "kex_ad": {"min": 0.0, "value": 0.0},
        "kex_bc": {"min": 0.0, "value": 200.0, "vary": True},
        "kex_bd": {"min": 0.0, "value": 0.0},
        "kex_cd": {"min": 0.0, "value": 200.0, "vary": True},
        "kab": {"expr": _get_k_from_kex_p("ab").format(**fnames)},
        "kba": {"expr": _get_k_from_kex_p("ba").format(**fnames)},
        "kac": {"expr": _get_k_from_kex_p("ac").format(**fnames)},
        "kca": {"expr": _get_k_from_kex_p("ca").format(**fnames)},
        "kad": {"expr": _get_k_from_kex_p("ad").format(**fnames)},
        "kda": {"expr": _get_k_from_kex_p("da").format(**fnames)},
        "kbc": {"expr": _get_k_from_kex_p("bc").format(**fnames)},
        "kcb": {"expr": _get_k_from_kex_p("cb").format(**fnames)},
        "kbd": {"expr": _get_k_from_kex_p("bd").format(**fnames)},
        "kdb": {"expr": _get_k_from_kex_p("db").format(**fnames)},
        "kcd": {"expr": _get_k_from_kex_p("cd").format(**fnames)},
        "kdc": {"expr": _get_k_from_kex_p("dc").format(**fnames)},
    }

    fnames, params = _apply_settings(fnames, params, settings)

    return fnames, params


def _add_params(fnames, params, names, attributes):
    for name in names:
        fname = _get_fullname(name, attributes=attributes)
        fnames[name] = fname
        params.add(name=fname, vary=False)

    return fnames, params


def _get_fullname(shortname, attributes=None):
    if attributes is None:
        attributes = {}
    return parameters.ParameterName(shortname, **attributes).to_full_name()


def _apply_settings(fnames, params, settings):

    settings_filtered = {key: val for key, val in settings.items() if key in fnames}

    for name, setting in settings_filtered.items():
        fname = fnames[name]
        params[fname]._delay_asteval = False
        params[fname].set(**setting)

    for name in settings_filtered:
        fname = fnames[name]
        params[fname]._delay_asteval = True

    return fnames, params


def _get_k_from_kex_p(states):
    kwargs = {
        "kex": "{{kex_{}{}}}".format(*sorted(states)),
        "p1": "{{p{}}}".format(states[0]),
        "psum": "{{p{}}} + {{p{}}}".format(*sorted(states)),
    }
    return "{kex} * {p1} / {psum} if {psum} else 0.0".format(**kwargs)


def _get_k_from_h_s(states):
    dh = "{{dh_{}{}}}".format(*sorted(states))
    ds = "{{ds_{}{}}}".format(*sorted(states))
    if states[0] != "a":
        dh += " - {{dh_{}}}".format(states[0])
        ds += " - {{ds_{}}}".format(states[0])
    return "{{kbt_h}} * exp(-(({dh}) - {{tk}} * ({ds})) / {{rt}})".format(dh=dh, ds=ds)


def _define_thermo(celcius):
    tk = celcius + 273.15
    kbt_h = cst.k * tk / cst.h
    rt = cst.R * tk
    thermo = {"tk": tk, "kbt_h": kbt_h, "rt": rt}
    return thermo

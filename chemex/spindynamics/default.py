import itertools
import re

import lmfit as lf
from scipy import constants as cst

from chemex import parameters
from chemex.spindynamics import util


OTHER_SPIN = {"i": "s", "s": "i"}

RE_VARNAME = re.compile(
    "^(?P<name>.+?)(_(?P<spin>[is]{1,2}))?(_(?P<state>[abcd]{1,2}))?$"
)

DEFAULT_PARAMS = lf.Parameters()

for state1, state2 in itertools.combinations("abcd", 2):
    DEFAULT_PARAMS.add(f"k{state1}{state2}", value=0.0, min=0.0)

for state in "abcd":
    DEFAULT_PARAMS.add(f"p{state}", value=0.0, min=0.0, max=1.0)
    DEFAULT_PARAMS.add(f"j_{state}", value=0.0)
    DEFAULT_PARAMS.add(f"r1a_{state}", value=2.0, min=0.0)
    DEFAULT_PARAMS.add(f"sigma_{state}", value=0.0)
    DEFAULT_PARAMS.add(f"r2_mq_{state}", value=10.0, min=0.0)
    DEFAULT_PARAMS.add(f"mu_mq_{state}", value=0.0)

for spin, state in itertools.product("is", "abcd"):
    DEFAULT_PARAMS.add(f"cs_{spin}_{state}", value=0.0)
    DEFAULT_PARAMS.add(f"r1_{spin}_{state}", value=1.0, min=0.0)
    DEFAULT_PARAMS.add(f"r2_{spin}_{state}", value=10.0, min=0.0)
    DEFAULT_PARAMS.add(f"r2a_{spin}_{state}", value=12.0, min=0.0)
    DEFAULT_PARAMS.add(f"etaxy_{spin}_{state}", value=0.0)
    DEFAULT_PARAMS.add(f"etaz_{spin}_{state}", value=0.0)


def create_params(
    basis,
    model,
    nuclei=None,
    conditions=None,
    hn_ap_constraints=False,
    nh_constraints=False,
):
    """TODO"""

    fnames_k, params_k = create_params_k(model_name=model.name, conditions=conditions)
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
    map_names = set_map_names(names, nuclei, conditions)
    params = buid_params(map_names)
    for name, full_name in map_names.items():
        if name.endswith(("_b", "_c", "_d")):
            name_a = "".join((name[:-1], "a"))
            full_name_a = map_names[name_a]
            params[full_name].set(expr=full_name_a)
    for state in ("b", "c", "d"):
        short_name = "cs_{}".format(state)
        short_name_expr = "cs_a + dw_a{}".format(state)
        parameters.set_param_expr(params, short_name, short_name_expr)
    return map_names, params


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


def set_map_names(names=None, nuclei=None, conditions=None):
    """TODO"""

    if names is None:
        names = {}

    map_names = {}

    for name in names:
        name_parsed = RE_VARNAME.match(name).groupdict()

        short_name = "_".join(
            name_parsed[key]
            for key in ("name", "state")
            if name_parsed[key] is not None
        )

        definition = {"temperature": conditions.get("temperature", 0.0)}

        if name.startswith(("k", "p")):
            definition["p_total"] = conditions.get("p_total", 1.0)
            definition["l_total"] = conditions.get("l_total", 1.0)
        else:
            spin = name_parsed["spin"]
            if spin is None:
                spin = "is"

            definition["nuclei"] = nuclei[spin]

        if not name.startswith(("dw", "cs", "j")):
            definition["h_larmor_frq"] = conditions.get("h_larmor_frq", 1.0)

        map_names[name] = parameters.ParameterName(
            short_name, **definition
        ).to_full_name()

    return map_names


def buid_params(map_names):
    """TODO"""

    params = lf.Parameters()

    for name, full_name in map_names.items():
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

    for spin, state in itertools.product("is", "abcd"):

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

    for state in "abcd":

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


def set_hn_ap_constraints(map_names, params, nuclei=None, conditions=None):

    try:
        name_s = nuclei["s"]
    except IndexError:
        raise

    definition = {"nuclei": name_s}
    definition.update(conditions)

    for state in "abcd":

        r1_i = "r1_i_{}".format(state)
        r2_i = "r2_i_{}".format(state)
        r1_s = "r1_s_{}".format(state)
        r2a_i = "r2a_i_{}".format(state)
        r1a = "r1a_{}".format(state)
        r1 = "r1_{}".format(state)

        if r1_i in map_names:

            if r1_s not in map_names:

                map_names[r1_s] = parameters.ParameterName(
                    r1, **definition
                ).to_full_name()

                if state != "a":
                    expr = map_names["r1_s_a"]
                else:
                    expr = None

                params[map_names[r1_s]] = lf.Parameter(
                    map_names[r1_s], value=1.0, min=0.0, vary=False, expr=expr
                )

            params[map_names[r1_i]].set(
                expr="{} - {}".format(map_names[r1a], map_names[r1_s])
            )

            params[map_names[r2a_i]].set(
                expr="{} - {}".format(map_names[r2_i], map_names[r1_s])
            )

    return map_names, params


def set_nh_constraints(map_names, params):

    for state in "abcd":
        r1_i = "r1_i_{}".format(state)
        r2_i = "r2_i_{}".format(state)
        r2a_i = "r2a_i_{}".format(state)
        r1a = "r1a_{}".format(state)

        if r2a_i in map_names:
            params[map_names[r2a_i]].set(
                expr="{} + {} - {}".format(
                    map_names[r2_i], map_names[r1a], map_names[r1_i]
                )
            )

    return map_names, params


def build_global_name(shortname, attributes=None):
    if attributes is None:
        attributes = {}
    return parameters.ParameterName(shortname, **attributes).to_full_name()


def create_params_k(model_name=None, conditions=None):
    """Update the experimental and fitting parameters depending on the model.

    """

    fnames = {}
    params = lf.Parameters()

    available_model = (
        "2st.pb_kex",
        "2st.eyring",
        "2st.binding",
        "3st.pb_kex",
        "4st.pb_kex",
    )

    if model_name not in available_model:
        print("Warning: The 'model' option should either be:")
        for model_name in available_model:
            print("    - '{}'".format(model_name))
        print("Set to the default value: '2st.pb_kex'.")
        model_name = "2st.pb_kex"

    model = util.parse_model(model_name)
    attributes = {
        key: conditions.get(key, None) for key in {"temperature", "p_total", "l_total"}
    }

    states_all = "abcd"
    states = states_all[: model.state_nb]

    for pairs in itertools.permutations(states, 2):
        name = "k{}{}".format(*pairs)
        name_global = build_global_name(name, attributes)
        fnames[name] = name_global
        params.add(name=name_global, min=0.0, vary=False)

    for state in states:
        name = "p{}".format(state)
        name_global = build_global_name(name, attributes)
        fnames[name] = name_global
        params.add(name=name_global, min=0.0, max=1.0, vary=False)

    value_dict = {}
    vary_list = []
    extra = {}
    expr_dict = {}

    if model_name == "2st.pb_kex":

        fnames["kex_ab"] = build_global_name("kex_ab", attributes)
        params.add(name=fnames["kex_ab"], min=0.0, vary=False)

        value_dict = {"pb": 0.05, "kex_ab": 200.0}

        vary_list = ("pb", "kex_ab")

        expr_dict = {
            "pa": "1.0 - {pb}",
            "kab": "{kex_ab} * {pb}",
            "kba": "{kex_ab} * {pa}",
        }

    elif model_name == "2st.eyring":

        attributes = {key: conditions.get(key, None) for key in {"p_total", "l_total"}}

        for name in ("dh_b", "ds_b", "dh_ab", "ds_ab"):
            fname = build_global_name(name, attributes=attributes)
            fnames[name] = fname
            params.add(name=fname, vary=False)

        value_dict = {
            "dh_b": 6.5e+03,
            "dh_ab": 6.5e+04,
            "ds_b": 1.0e+10,
            "ds_ab": 1.0e+10,
        }

        vary_list = ("dh_b", "dh_ab")

        tk = conditions.get("temperature", 0.0) + 273.15
        kbt_h = cst.k * tk / cst.h
        rt = cst.R * tk

        extra = {"tk": tk, "kbt_h": kbt_h, "rt": rt}

        expr_dict = {
            "kab": "{kbt_h} * exp(-({dh_ab} - {tk} * {ds_ab}) / {rt})",
            "kba": "{kbt_h} * exp(-(({dh_ab} - {dh_b}) - {tk} * ({ds_ab} - {ds_b})) / {rt})",
            "pa": "{kba} / ({kba} + {kab})",
            "pb": "{kab} / ({kba} + {kab})",
        }

    elif model_name == "2st.binding":

        attributes = {"temperature": conditions.get("temperature", None)}

        for name in ("kd", "kon", "koff"):
            fname = build_global_name(name, attributes=attributes)
            fnames[name] = fname
            params.add(name=fname, min=0.0, vary=False)

        value_dict = {"kd": 1e-6, "kon": 1e7, "koff": 10.0}

        vary_list = ("kon", "koff")

        p_total = conditions.get("p_total", 0.0)
        l_total = conditions.get("l_total", 0.0)
        delta = l_total - p_total

        extra = {"delta": delta, "l_total": l_total}

        expr_dict = {
            "kd": "{koff} / {kon}",
            "kab": "{kon} * 0.5 * ({delta} - {kd} + sqrt(({delta} - {kd}) ** 2 + 4.0 * {kd} * {l_total}))",
            "kba": "{koff}",
            "pa": "{kba} / ({kba} + {kab})",
            "pb": "{kab} / ({kba} + {kab})",
        }

    elif model_name == "3st.pb_kex":

        for name in ("kex_ab", "kex_ac", "kex_bc"):
            fname = build_global_name(name, attributes)
            fnames[name] = fname
            params.add(name=fname, min=0.0)

        value_dict = {
            "pb": 0.02,
            "pc": 0.02,
            "kex_ab": 200.0,
            "kex_ac": 0.0,
            "kex_bc": 200.0,
        }

        vary_list = ("pb", "pc", "kex_ab", "kex_bc")

        expr_dict = {
            "pa": "1.0 - {pb} - {pc}",
            "kab": "{kex_ab} * {pb} / ({pa} + {pb}) if {pa} + {pb} else 0.0",
            "kba": "{kex_ab} * {pa} / ({pa} + {pb}) if {pa} + {pb} else 0.0",
            "kac": "{kex_ac} * {pc} / ({pa} + {pc}) if {pa} + {pc} else 0.0",
            "kca": "{kex_ac} * {pa} / ({pa} + {pc}) if {pa} + {pc} else 0.0",
            "kbc": "{kex_bc} * {pc} / ({pb} + {pc}) if {pb} + {pc} else 0.0",
            "kcb": "{kex_bc} * {pb} / ({pb} + {pc}) if {pb} + {pc} else 0.0",
        }

    elif model_name == "3st.eyring":

        attributes = {key: conditions.get(key, None) for key in {"p_total", "l_total"}}

        for name in (
            "dh_b",
            "ds_b",
            "dh_c",
            "ds_c",
            "dh_ab",
            "ds_ab",
            "dh_ac",
            "ds_ac",
            "dh_bc",
            "ds_bc",
        ):
            fname = build_global_name(name, attributes=attributes)
            fnames[name] = fname
            params.add(name=fname, vary=False)

        value_dict = {
            "dh_b": 6.5e+03,
            "dh_c": 6.5e+03,
            "dh_ab": 6.5e+04,
            "dh_ac": 1.0e+10,
            "dh_bc": 6.5e+04,
            "ds_b": 1.0e+10,
            "ds_c": 1.0e+10,
            "ds_ab": 1.0e+10,
            "ds_ac": 1.0e+10,
            "ds_bc": 1.0e+10,
        }

        vary_list = ("dh_b", "dh_c", "dh_ab", "dh_bc")

        tk = conditions.get("temperature", 0.0) + 273.15
        kbt_h = cst.k * tk / cst.h
        rt = cst.R * tk

        extra = {"tk": tk, "kbt_h": kbt_h, "rt": rt}

        expr_dict = {
            "kab": "{kbt_h} * exp(-({dh_ab} - {tk} * {ds_ab}) / {rt})",
            "kba": "{kbt_h} * exp(-(({dh_ab} - {dh_b}) - {tk} * ({ds_ab} - {ds_b})) / {rt})",
            "kac": "{kbt_h} * exp(-(({dh_ac} - {dh_a}) - {tk} * ({ds_ac} - {ds_a})) / {rt})",
            "kca": "{kbt_h} * exp(-(({dh_ac} - {dh_c}) - {tk} * ({ds_ac} - {ds_c})) / {rt})",
            "kbc": "{kbt_h} * exp(-(({dh_bc} - {dh_b}) - {tk} * ({ds_bc} - {ds_b})) / {rt})",
            "kcb": "{kbt_h} * exp(-(({dh_bc} - {dh_c}) - {tk} * ({ds_bc} - {ds_c})) / {rt})",
            "pa": "{kba} * {kca} / ({kba} * {kca} + {kab} * {kca} + {kac} * {kba})",
            "pb": "{kab} * {kcb} / ({kab} * {kcb} + {kba} * {kcb} + {kbc} * {kab})",
            "pc": "{kbc} * {kac} / ({kbc} * {kac} + {kcb} * {kab} + {kca} * {kbc})",
        }

    elif model_name == "4st.pb_kex":

        for name in ("kex_ab", "kex_ac", "kex_ad", "kex_bc", "kex_bd", "kex_cd"):
            fname = build_global_name(name, attributes)
            fnames[name] = fname
            params.add(name=fname, min=0.0)

        value_dict = {
            "pb": 0.02,
            "pc": 0.02,
            "pd": 0.02,
            "kex_ab": 200.0,
            "kex_ac": 0.0,
            "kex_ad": 0.0,
            "kex_bc": 200.0,
            "kex_bd": 0.0,
            "kex_cd": 200.0,
        }

        vary_list = ("pb", "pc", "pd", "kex_ab", "kex_bc", "kex_cd")

        expr_dict = {
            "pa": "1.0 - {pb} - {pc} - {pd}",
            "kab": "{kex_ab} * {pb} / ({pa} + {pb}) if {pa} + {pb} else 0.0",
            "kba": "{kex_ab} * {pa} / ({pa} + {pb}) if {pa} + {pb} else 0.0",
            "kac": "{kex_ac} * {pc} / ({pa} + {pc}) if {pa} + {pc} else 0.0",
            "kca": "{kex_ac} * {pa} / ({pa} + {pc}) if {pa} + {pc} else 0.0",
            "kad": "{kex_ad} * {pd} / ({pa} + {pd}) if {pa} + {pd} else 0.0",
            "kda": "{kex_ad} * {pa} / ({pa} + {pd}) if {pa} + {pd} else 0.0",
            "kbc": "{kex_bc} * {pc} / ({pb} + {pc}) if {pb} + {pc} else 0.0",
            "kcb": "{kex_bc} * {pb} / ({pb} + {pc}) if {pb} + {pc} else 0.0",
            "kbd": "{kex_bd} * {pd} / ({pb} + {pd}) if {pb} + {pd} else 0.0",
            "kdb": "{kex_bd} * {pb} / ({pb} + {pd}) if {pb} + {pd} else 0.0",
            "kcd": "{kex_cd} * {pd} / ({pc} + {pd}) if {pc} + {pd} else 0.0",
            "kdc": "{kex_cd} * {pc} / ({pc} + {pd}) if {pc} + {pd} else 0.0",
        }

    for name, value in value_dict.items():
        fname = fnames[name]
        params[fname].value = value

    for name in vary_list:
        fname = fnames[name]
        params[fname].vary = True

    for name, expr in expr_dict.items():
        fname = fnames[name]
        params[fname].expr = expr.format(**fnames, **extra)

    return fnames, params

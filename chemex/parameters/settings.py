"""The parameters module contains the code for handling of the experimental and
fitting parameters."""
import ast
import collections as cl
import configparser
import difflib as dl
import re

import asteval.astutils as aa

import chemex.helper as ch
import chemex.nmr.helper as cnh
import chemex.nmr.rates as cnr
import chemex.parameters.name as cpn


RE_GROUPNAME = re.compile(r"^[A-Za-z0-9_-]+$")

_SCHEMA_CONFIG_PARAM = {
    "type": "object",
    "properties": {
        "model_free": {
            "type": "object",
            "properties": {
                "tauc": {"type": "number"},
                "taui": {"type": "number"},
                "s2": {"type": "number"},
                "deuterated": {"type": "boolean"},
            },
            "additionalProperties": False,
            "default": {},
        }
    },
    "additionalProperties": {
        "type": "object",
        "additionalProperties": {
            "oneOf": [
                {
                    "type": "array",
                    "items": {"type": "number"},
                    "minItems": 3,
                    "maxItems": 4,
                },
                {"type": "number"},
            ]
        },
    },
}


def set_params_from_files(params, experiments, filenames):
    """Read the parameter file and set initial values and optional bounds and brute
    step size."""
    config = _read_config(filenames)
    matches = cl.Counter()
    cfg_model_free = config.pop("model_free")
    _set_param_mf(params, experiments, cfg_model_free)
    for section, settings in config.items():
        prefix = f"{section}, NUC->" if section != "global" else ""
        for key, values in settings.items():
            name = cpn.ParamName.from_section(f"{prefix}{key}")
            if isinstance(values, float):
                values = [values]
            values_ = dict(zip(("value", "min", "max", "brute_step"), values))
            matched = set_params(params, name, **values_)
            matches.update({name.to_section_name(): len(matched)})
    _print_matches(matches)


def _set_param_mf(params, experiments, cfg_model_free):
    model_free = cnr.ModelFree(**cfg_model_free)
    experiments.set_params(params, model_free)


def _read_config(filenames):
    print(f"\nSetting parameters starting values...")
    config = ch.read_toml_multi(filenames, _SCHEMA_CONFIG_PARAM)
    config_ = {"global": config.pop("global", {})}
    config_.update(config)
    return config_


def set_param_status(params, settings):
    """Set whether or not to vary a fitting parameter or to use a mathemetical
    expression."""
    vary = {"fix": False, "fit": True}
    if settings:
        print("\nSettings parameter status and constraints...")
    matches = cl.Counter()
    for key, status in settings.items():
        name = cpn.ParamName.from_section(key)
        if status in vary:
            matched = set_params(params, name, vary=vary[status], expr="")
        else:
            matched = set_param_expr(params, name, expr=status)
        matches.update({name.to_section_name(): len(matched)})
    _print_matches(matches)


def _print_matches(matches):
    for key, value in matches.items():
        if value > 0:
            print(f"  - [{key}] -> {value} parameter(s) updated")
        else:
            print(f"  - [{key}] -> Nothing to update")


def set_param_expr(params, name, expr=None):
    """Set an optional parameter expression, used to constrain its value during
    the fit."""
    if expr is None:
        expr = ""
    if not isinstance(name, cpn.ParamName):
        name = cpn.ParamName.from_section(name)
    fnames = [fname for fname in params if name.match(fname)]
    names_expr = aa.get_ast_names(ast.parse(expr))
    fnames_expr = {
        name: [
            fname_expr
            for fname_expr in params
            if cpn.ParamName.from_section(name).match(fname_expr)
        ]
        for name in names_expr
    }
    matches = set()
    for fname in fnames:
        expr_ = str(expr)
        for name_expr in names_expr:
            fname_expr = dl.get_close_matches(fname, fnames_expr[name_expr], n=1)[0]
            expr_ = expr_.replace(name_expr, fname_expr)
        param = params[fname]
        repr_ = repr(param)
        param.expr = expr_
        if repr(param) != repr_:
            matches.add(fname)
    return matches


def set_params(
    params,
    name_short,
    value=None,
    vary=None,
    min=None,
    max=None,
    expr=None,
    brute_step=None,
):
    """Set the initial value and (optional) bounds and brute step size for
    parameters."""
    matches = set()
    for name, param in params.items():
        if name_short.match(name):
            repr_ = repr(param)
            if expr is None and param.expr and vary is None:
                param.value = value
            else:
                param.set(value, vary, min, max, expr, brute_step)
            if repr(param) != repr_:
                matches.add(name)
    return matches


def write_par(params, path):
    """Write the fitting parameters and their uncertainties to a file."""
    filename = path / "parameters.toml"
    par_dict = {}
    for name, param in params.items():
        par_name = cpn.ParamName.from_full_name(name)
        if not par_name.spin_system:  # global parameter
            name_print = par_name
            section = "GLOBAL"
        else:  # residue-specific parameter
            name_print = cnh.SpinSystem(par_name.spin_system)
            section = par_name.to_section_name()
        par_dict.setdefault(section, []).append((name_print, _param_to_string(param)))
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    section_global = par_dict.pop("GLOBAL", None)
    if section_global is not None:
        cfg.add_section("GLOBAL")
        for name, val in sorted(section_global):
            cfg.set("GLOBAL", f"{str(name):10s}", val)
    for section, name_vals in sorted(par_dict.items()):
        if not RE_GROUPNAME.match(section):
            section = f'"{section}"'
        cfg.add_section(section)
        for peak, val in sorted(name_vals):
            cfg.set(section, f"{peak.name.upper():10s}", val)
    with open(filename, "w") as f:
        cfg.write(f)


def _param_to_string(param):
    comments = []
    error = ""
    comments_fmtd = ""
    if param.vary or param.expr:
        if param.stderr:
            error += f" ±{param.stderr:.5e}"
        elif not param.expr:
            comments.append("error not calculated")
    if param.expr:
        comments.append("constrained")
    elif not param.vary:
        comments.append("fixed")
    if comments:
        comments_fmtd = f" ({'; '.join(comments)})"
    string_ = f"{param.value: .5e}"
    if error or comments_fmtd:
        string_ += f" #{error}{comments_fmtd}"
    return string_


def write_constraints(params, path):
    """Write the (optional) parameter expression constraints to a file."""
    filename = path / "constraints.fit"
    param_dict = dict()
    for name, param in params.items():
        par_name = cpn.ParamName.from_full_name(name)
        if param.expr:
            name_formatted = "[{}]".format(
                par_name.to_section_name(show_spin_system=True)
            )
            expr_formatted = param.expr
            for name_dep in aa.get_ast_names(ast.parse(param.expr)):
                par_name_dep = cpn.ParamName.from_full_name(name_dep)
                if str(par_name_dep):
                    expr_formatted = expr_formatted.replace(
                        name_dep,
                        "[{}]".format(
                            par_name_dep.to_section_name(show_spin_system=True)
                        ),
                    )
            param_dict[par_name] = f"{name_formatted} = {expr_formatted}\n"
    with open(filename, "w") as f:
        for name, constraint in sorted(param_dict.items()):
            f.write(constraint)

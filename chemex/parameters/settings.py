"""The parameters module contains the code for handling of the experimental and
fitting parameters."""
import ast
import collections as cl
import difflib as dl
import re

import asteval.astutils as aa
import numpy as np

import chemex.helper as ch
import chemex.parameters.name as cpn


RE_GROUPNAME = re.compile(r"^[A-Za-z0-9_-]+$")

_SCHEMA_CONFIG_PARAM = {
    "type": "object",
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


def read_defaults(filenames):
    print("\nReading parameters default values...")
    config = ch.read_toml_multi(filenames, _SCHEMA_CONFIG_PARAM)
    defaults = {"global": config.pop("global", {}), **config}
    return defaults_to_re_list(defaults)


def defaults_to_re_list(defaults):
    defaults_re = []
    for section, settings in defaults.items():
        prefix = f"{section}, NUC->" if section != "global" else ""
        for key, values in settings.items():
            name = cpn.ParamName.from_section(f"{prefix}{key}")
            if not isinstance(values, cl.abc.Iterable):
                values = [values]
            values_ = dict(zip(("value", "min", "max", "brute_step"), values))
            defaults_re.append((name, values_))
    return defaults_re


def set_values(params, defaults):
    """Read the parameter file and set initial values and optional bounds and brute
    step size."""

    fnames = {param.name for param in params.values() if not param.expr}

    for name, values in reversed(defaults):
        matches = set()
        for fname in fnames:
            if name.match(fname):
                params[fname].set(**values)
                matches.add(fname)
        fnames -= matches
    params.update_constraints()
    _check_params(params)
    return params


def _check_params(params):
    messages = []
    for param in params.values():
        pname = param.user_data["pname"]
        sname = pname.name.upper()
        atoms = set(pname.spin_system.atoms.values())
        if sname.startswith("J_") and atoms == {"N", "H"} and param.value > 0.0:
            messages.append(
                "Warning: Some 1J(NH) scalar couplings are set with positive values. \n"
                "This can cause the TROSY and anti-TROSY components to be switched in \n"
                "some experiments."
            )
        if sname.startswith("J_") and atoms == {"C", "H"} and param.value < 0.0:
            messages.append(
                "Warning: Some 1J(CH) scalar couplings are set with negative values. \n"
                "This can cause the TROSY and anti-TROSY components to be switched in \n"
                "some experiments."
            )
    if messages:
        print("")
        print("\n\n".join(np.unique(messages)))


def set_status(params, settings=None, verbose=True):
    """Set whether or not to vary a fitting parameter or to use a mathemetical
    expression."""
    if settings is None:
        settings = {}
    if verbose and settings:
        print("\nSettings parameter status and constraints...")
    vary = {"fix": False, "fit": True}
    matches = cl.Counter()
    for key, status in settings.items():
        name = cpn.ParamName.from_section(key)
        if status in vary:
            matched = _set_params(params, name, {"vary": vary[status], "expr": ""})
        else:
            matched = _set_param_expr(params, name, expr=status)
        matches.update({name.section: len(matched)})
    if verbose:
        _print_matches(matches)


def _print_matches(matches):
    for key, value in matches.items():
        if value > 0:
            print(f"  - [{key}] -> {value} parameter(s) updated")
        else:
            print(f"  - [{key}] -> Nothing to update")


def _set_param_expr(params, name, expr):
    """Set an optional parameter expression, used to constrain its value during
    the fit."""
    if not isinstance(name, cpn.ParamName):
        name = cpn.ParamName.from_section(name)
    pnames = [fname for fname in params if name.match(fname)]
    names_expr = aa.get_ast_names(ast.parse(expr))
    pnames_expr = {
        name: [
            fname_expr
            for fname_expr in params
            if cpn.ParamName.from_section(name).match(fname_expr)
        ]
        for name in names_expr
    }
    matches = set()
    for pname in pnames:
        expr_ = str(expr)
        for name_expr in names_expr:
            pname_expr = dl.get_close_matches(pname, pnames_expr[name_expr], n=1)[0]
            expr_ = expr_.replace(name_expr, pname_expr)
        param = params[pname]
        repr_ = repr(param)
        param.expr = expr_
        print_expr = str(expr_)
        for fname in aa.get_ast_names(ast.parse(expr_)):
            print_expr = print_expr.replace(
                fname, str(params[fname].user_data["pname"])
            )
        param.user_data["print_expr"] = print_expr
        if repr(param) != repr_:
            matches.add(pname)
    return matches


def _set_params(params, name_short, values):
    """Set the initial value and (optional) bounds and brute step size for
    parameters."""
    matches = set()
    value = values.pop("value", None)
    for name, param in params.items():
        if name_short.match(name):
            repr_ = repr(param)
            if value is not None:
                param.value = value
            param.set(**values)
            if repr(param) != repr_:
                matches.add(name)
    return matches


def write_par(params, path):
    """Write the fitting parameters and their uncertainties to a file."""
    fitted_, fixed_, constrained_ = {}, {}, {}
    for name, param in params.items():
        if param.vary:
            fitted_[name] = param
        elif param.expr:
            constrained_[name] = param
        else:
            fixed_[name] = param
    fitted = _params_to_text(fitted_)
    constrained = _params_to_text(constrained_)
    fixed = _params_to_text(fixed_)
    path_ = path / "Parameters"
    path_.mkdir(parents=True, exist_ok=True)
    if fitted:
        (path_ / "fitted.toml").write_text(fitted)
    if constrained:
        (path_ / "constrained.toml").write_text(constrained)
    if fixed:
        (path_ / "fixed.toml").write_text(fixed)


def _params_to_text(params):
    par_string, col_width = _params_to_string(params)
    result = []
    for section, key_values in par_string.items():
        result.append(f"[{_quote(section)}]")
        width = col_width[section]
        for key, value in sorted(key_values):
            result.append(f"{_quote(key):<{width}} = {value}")
        result.append("")
    return "\n".join(result)


def _quote(text):
    text_ = str(text)
    if RE_GROUPNAME.match(text_):
        return text_
    return f'"{text_}"'


def _params_to_string(params):
    col_width = {}
    par_string = {"GLOBAL": []}
    sorted_params = sorted(params.values(), key=lambda x: x.user_data["pname"])
    for param in sorted_params:
        pname = param.user_data["pname"]
        if pname.spin_system:  # residue-specific parameter
            name_print = pname.spin_system
            section = pname.section
        else:  # global parameter
            name_print = pname.section_res
            section = "GLOBAL"
        value_print = _param_to_string(param)
        par_string.setdefault(section, []).append((name_print, value_print))
        col_width[section] = max(10, len(str(name_print)), col_width.get(section, 0))
    if not par_string["GLOBAL"]:
        del par_string["GLOBAL"]
    return par_string, col_width


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
        comments.append(param.user_data["print_expr"])
    elif not param.vary:
        comments.append("fixed")
    if comments:
        comments_fmtd = f" ({'; '.join(comments)})"
    value = param.value if param.value != 0.0 else 0.0
    string_ = f"{value: .5e}"
    if error or comments_fmtd:
        string_ += f" #{error}{comments_fmtd}"
    return string_

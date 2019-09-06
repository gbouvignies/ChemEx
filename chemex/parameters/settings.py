"""The parameters module contains the code for handling of the experimental and
fitting parameters."""
import ast
import collections as cl
import configparser
import difflib as dl

import asteval.astutils as aa
import jsonschema as js

import chemex.helper as ch
import chemex.nmr.helper as cnh
import chemex.parameters.name as cpn


_SCHEMA_CONFIG_PARAM = {
    "type": "object",
    "patternProperties": {
        r"^[A-Za-z_][A-Za-z0-9_,->]*$": {
            "type": "object",
            "patternProperties": {
                r"^[A-Za-z_][A-Za-z0-9_,->]*$": {
                    "oneOf": [
                        {
                            "type": "array",
                            "items": {"type": "number"},
                            "minItems": 3,
                            "maxItems": 4,
                        },
                        {"type": "number"},
                    ]
                }
            },
        }
    },
}


def set_params_from_config_file(params, config_filename):
    """Read the parameter file and set initial values and optional bounds and brute
    step size."""
    print(f"\nReading '{config_filename}'...")
    config = ch.read_toml(config_filename)
    try:
        js.validate(config, _SCHEMA_CONFIG_PARAM)
    except js.ValidationError as err:
        print("Validation error: {0}".format(err))
        raise
    matches = cl.Counter()
    for section in config:
        is_not_global = section.lower() != "global"
        prefix = f"{section}, NUC->" if is_not_global else ""
        for key, values in config[section].items():
            name = cpn.ParamName.from_section(f"{prefix}{key}")
            if isinstance(values, float):
                values = [values]
            values_ = dict(zip(("value", "min", "max", "brute_step"), values))
            matched = set_params(params, name, **values_)
            matches.update({name.to_section_name(): len(matched)})
    for key, value in matches.items():
        print(f"  - [{key}] -> {value} parameter(s) updated")


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
    for key, value in matches.items():
        print(f"  - [{key}] -> {value} parameter(s) updated")


def set_param_expr(params, name, expr=None):
    """Set an optional parameter expression, used to constrain its value during
    the fit."""
    if expr is None:
        expr = ""
    if not isinstance(name, cpn.ParamName):
        name = cpn.ParamName.from_section(name)
    names_full = [name_full for name_full in params if name.match(name_full)]
    names_expr = aa.get_ast_names(ast.parse(expr))
    names_full_expr = {
        name: [
            name_full_expr
            for name_full_expr in params
            if cpn.ParamName.from_section(name).match(name_full_expr)
        ]
        for name in names_expr
    }
    matches = set()
    for name_full in names_full:
        expr_ = expr
        for name_expr in names_expr:
            name_full_expr = dl.get_close_matches(
                name_full, names_full_expr[name_expr], n=1
            )[0]
            expr_ = expr_.replace(name_expr, name_full_expr)
        params[name_full].expr = expr_
        matches.add(name_full)
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
            if expr is None and param.expr and vary is None:
                param.value = value
            else:
                param.set(value, vary, min, max, expr, brute_step)
            matches.add(name)
    return matches


def write_par(params, path):
    """Write the fitting parameters and their uncertainties to a file."""
    filename = path / "parameters.toml"
    print(f"  - {filename}")
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
        if not ch.RE_GROUPNAME.match(section):
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
    print(f"  - {filename}")
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


def remove_comments(line, sep):
    """Remove (optional) comments."""
    for a_sep in sep:
        line = line.split(a_sep)[0]
    return line.strip()

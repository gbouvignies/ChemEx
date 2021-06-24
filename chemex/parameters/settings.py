"""The parameters module contains the code for handling of the experimental and
fitting parameters."""
import collections as cl
import difflib as dl
import re

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
    # This is to put the "GLOBAL" section first
    defaults = {"global": config.pop("global", {}), **config}
    return _defaults_to_list(defaults)


def _defaults_to_list(defaults):
    """Take the TOML dictionary and make it a list of pairs of ParamName
    and keyword arguments for the method lmfit.Parameter.set"""
    defaults_list = []
    for section, settings in defaults.items():
        prefix = f"{section}, NUC->" if section != "global" else ""
        for key, values in settings.items():
            pname = cpn.ParamName.from_section(f"{prefix}{key}")
            values = values if isinstance(values, cl.abc.Iterable) else [values]
            values_ = dict(zip(("value", "min", "max", "brute_step"), values))
            defaults_list.append((pname, values_))
    return defaults_list


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
    """Check wheither the J couplings have the right sign"""
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

    matches_con = _set_expr(params, settings.get("constraints"))
    if verbose:
        _print_matches(matches_con, params, "constrained")

    matches_fix = _set_vary(params, settings.get("fix"), vary=False)
    if verbose:
        _print_matches(matches_fix, params, "fixed")

    matches_fit = _set_vary(params, settings.get("fit"), vary=True)
    if verbose:
        _print_matches(matches_fit, params, "varied")
        print("")


def _set_vary(params, snames, vary=True):
    """Set wheither the parameters corresponding to snames vary or not"""
    matches = set()
    if not snames:
        return matches
    fnames = {param.name for param in params.values() if param.vary != vary}
    for sname in reversed(snames):
        pname = cpn.ParamName.from_section(sname.strip("[] "))
        for fname in fnames:
            if pname.match(fname):
                params[fname].set(vary=vary)
                matches.add(fname)
        fnames -= matches
    return matches


def _set_expr(params, expr_list):
    """Set the constraints given in expr_list"""
    matches = set()
    if not expr_list:
        return matches
    fnames = set(params)
    for expr in reversed(expr_list):
        left, right, *somethingelse = expr.split("=")
        if somethingelse:
            print(f'\nError reading constraints:\n  -> "{expr}"\n\nProgram aborted\n')
            exit()
        fnames_left = _get_fnames_left(left, params)
        fnames_right = _get_fnames_right(right, params)
        for fname in fnames_left:
            expr_new = print_expr = right.strip()
            for sname, fname_set in fnames_right.items():
                fname_replace = dl.get_close_matches(fname, fname_set, n=1).pop()
                sname_replace = str(params[fname_replace].user_data["pname"])
                expr_new = expr_new.replace(sname, fname_replace)
                print_expr = print_expr.replace(sname, sname_replace)
            if expr_new != params[fname].expr:
                params[fname].expr = expr_new
                params[fname].user_data["print_expr"] = print_expr
                matches.add(fname)
        fnames -= matches
    return matches


def _get_fnames_left(left, params):
    pname_left = cpn.ParamName.from_section(left.strip("[] "))
    return {fname for fname in params if pname_left.match(fname)}


def _get_fnames_right(right, params):
    fnames_right = {}
    for match in re.finditer(r"\[(.+?)\]", right.strip()):
        pname = cpn.ParamName.from_section(match.group(1))
        fnames_right[match.group(0)] = {fname for fname in params if pname.match(fname)}
    return fnames_right


def read_grid(grid, params):
    if grid is None:
        return None, params
    re_ = re.compile(
        fr"(lin[(]{cpn._RE_FLOAT},{cpn._RE_FLOAT},\d+[)]$)|"
        fr"(log[(]{cpn._RE_FLOAT},{cpn._RE_FLOAT},\d+[)]$)|"
        fr"([(](({cpn._RE_FLOAT})(,|[)]$))+)"
    )
    fnames_all = set(params)
    grid_values = {}
    for entry in reversed(grid):
        left, right = entry.replace(" ", "").split("=")
        if not re_.match(right):
            print(
                f'\nError reading grid settings:\n  -> "{entry}"\n\nProgram aborted\n'
            )
            exit()
        fnames_left = _get_fnames_left(left, fnames_all)
        expr = right.replace("lin", "np.linspace").replace("log", "np.geomspace")
        values = eval(expr)
        grid_values.update({fname: values for fname in fnames_left})
        fnames_all -= fnames_left
    for fname in grid_values:
        params[fname].vary = False
    return grid_values, params


def _print_matches(matches, params, status):
    count_dict = {}
    for fname in matches:
        sname = params[fname].user_data["pname"].section
        count_dict.setdefault(sname, []).append(fname)
    for key, fnames in count_dict.items():
        number = len(fnames)
        if number > 0:
            print(f"  - [{key}] -> {number} parameter(s) set to be {status}")


def write_par(params, path):
    """Write the model parameter values and their uncertainties to a file"""
    path_ = path / "Parameters"
    path_.mkdir(parents=True, exist_ok=True)
    for status in ("fitted", "constrained", "fixed"):
        selection = _select[status](params)
        if not selection:
            continue
        par_strings = _params_to_strings(selection, status)
        formatted_strings = _format_strings(par_strings)
        filename = path_ / f"{status}.toml"
        filename.write_text(formatted_strings)


def _select_fitted(params):
    return {fname: param for fname, param in params.items() if param.vary}


def _select_constrained(params):
    return {fname: param for fname, param in params.items() if param.expr}


def _select_fixed(params):
    return {
        fname: param
        for fname, param in params.items()
        if not (param.vary or param.expr)
    }


_select = {
    "fitted": _select_fitted,
    "constrained": _select_constrained,
    "fixed": _select_fixed,
}


def _params_to_strings(params, status):
    par_strings = {"GLOBAL": {}}
    params_dict = {param.user_data["pname"]: param for param in params.values()}
    for pname, param in sorted(params_dict.items()):
        if pname.spin_system:  # residue-specific parameter
            section = pname.section
            name_print = str(pname.spin_system)
        else:  # global parameter
            section = "GLOBAL"
            name_print = pname.section_res
        value_print = _format_param[status](param)
        par_strings.setdefault(section, {})[name_print] = value_print
    if not par_strings["GLOBAL"]:
        del par_strings["GLOBAL"]
    return par_strings


def _format_fitted(param):
    error = f"±{param.stderr:.5e}" if param.stderr else "(error not calculated)"
    return f"{param.value: .5e} # {error}"


def _format_constrained(param):
    error = f"±{param.stderr:.5e} " if param.stderr else ""
    constraint = param.user_data["print_expr"]
    return f"{param.value: .5e} # {error}({constraint})"


def _format_fixed(param):
    return f"{param.value: .5e} # (fixed)"


_format_param = {
    "fitted": _format_fitted,
    "constrained": _format_constrained,
    "fixed": _format_fixed,
}


def _format_strings(par_strings):
    result = []
    for section, key_values in par_strings.items():
        result.append(f"[{_quote(section)}]")
        width = len(max(key_values, key=len))
        for key, value in key_values.items():
            result.append(f"{_quote(key):<{width}} = {value}")
        result.append("")
    return "\n".join(result)


def _quote(text):
    text_ = str(text)
    if RE_GROUPNAME.match(text_):
        return text_
    return f'"{text_}"'

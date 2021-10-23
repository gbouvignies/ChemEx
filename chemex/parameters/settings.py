"""The parameters module contains the code for handling of the experimental and
fitting parameters."""
from __future__ import annotations

from collections.abc import Iterable
from difflib import get_close_matches
from re import compile
from re import finditer

import numpy as np
from lmfit import Parameters
from Levenshtein import distance

from chemex.containers.conditions import Conditions
from chemex.helper import read_toml_multi
from chemex.nmr.spin_system import Atom
from chemex.nmr.spin_system import Group
from chemex.nmr.spin_system import Nucleus
from chemex.parameters.name import _RE_FLOAT
from chemex.parameters.name import ParamName


RE_GROUPNAME = compile(r"^[A-Za-z0-9_-]+$")

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


class ParamIndex:
    index: dict[str | Group | Atom | Nucleus | Conditions, set[str]]

    def __init__(self, params: Parameters) -> None:
        self.index = self._index_params(params)

    @staticmethod
    def _index_params(params: Parameters) -> dict:
        index: dict = {}
        for fname, param in params.items():
            search_keys: set = param.user_data["pname"].search_keys
            for key in search_keys:
                index.setdefault(key, set()).add(fname)
        return index

    def match(self, pname: ParamName) -> set[str]:
        search_keys = pname.search_keys
        return set.intersection(*(self.index.get(key, set()) for key in search_keys))


def read_defaults(filenames):
    print("\nReading parameters default values...")
    config = read_toml_multi(filenames, _SCHEMA_CONFIG_PARAM)
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
            pname = ParamName.from_section(f"{prefix}{key}")
            values = values if isinstance(values, Iterable) else [values]
            values_ = dict(zip(("value", "min", "max", "brute_step"), values))
            defaults_list.append((pname, values_))
    return defaults_list


def set_values(params, defaults):
    """Read the parameter file and set initial values and optional bounds and brute
    step size."""

    fnames = {param.name for param in params.values() if not param.expr}

    index = ParamIndex(params)

    for pname, values in reversed(defaults):
        matches = index.match(pname)
        for fname in matches & fnames:
            params[fname].set(**values)
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
        atoms = {atom.name for atom in pname.spin_system.atoms.values()}
        if sname.startswith("J_") and atoms == {"N", "H"} and param.value > 0.0:
            messages.append(
                "Warning: Some 1J(NH) scalar couplings are set with positive values.\n"
                "This can cause the TROSY and anti-TROSY components to be switched in\n"
                "some experiments."
            )
        if sname.startswith("J_") and atoms == {"C", "H"} and param.value < 0.0:
            messages.append(
                "Warning: Some 1J(CH) scalar couplings are set with negative values.\n"
                "This can cause the TROSY and anti-TROSY components to be switched in\n"
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
    all_matches = set()

    if not snames:
        return all_matches

    fnames = {param.name for param in params.values() if param.vary != vary}

    index = ParamIndex(params)

    for sname in reversed(snames):
        pname = ParamName.from_section(sname.strip("[] "))
        matches = index.match(pname)
        for fname in matches & fnames:
            params[fname].set(vary=vary)
        fnames -= matches
        all_matches.update(matches)
    return all_matches


def _set_expr(params, expr_list):
    """Set the constraints given in expr_list"""
    all_matches = set()

    if not expr_list:
        return all_matches

    fnames = set(params)

    index = ParamIndex(params)

    for expr in reversed(expr_list):

        left, right, *somethingelse = expr.split("=")

        if somethingelse:
            print(f'\nError reading constraints:\n  -> "{expr}"\n\nProgram aborted\n')
            exit()

        fnames_left = _get_fnames_left(left, index)
        fnames_right = _get_fnames_right(right, index)

        for fname in fnames_left:
            expr_new = print_expr = right.strip()
            for sname, fname_set in fnames_right.items():
                # String matching was nitially implemented using `difflib` library,
                # but using the `distance` function from the `Levenshtein` package is
                # orders of magnitude faster.
                #
                # Other implementations:
                # fname_replace = difflib.get_close_matches(fname, fname_set, n=1).pop()
                # fname_replace, _ = fuzzywuzzy.process.extractOne(fname, fname_set)
                # fname_replace, _ = thefuzz.process.extractOne(fname, fname_set)
                #
                fname_replace = min(
                    fname_set, key=lambda a_fname: distance(a_fname, fname)
                )
                sname_replace = str(params[fname_replace].user_data["pname"])
                expr_new = expr_new.replace(sname, fname_replace)
                print_expr = print_expr.replace(sname, sname_replace)

            if expr_new != params[fname].expr:
                params[fname].expr = expr_new
                params[fname].user_data["print_expr"] = print_expr
                all_matches.add(fname)

        fnames -= all_matches

    return all_matches


def _get_fnames_left(left, index):
    pname = ParamName.from_section(left.strip("[] "))
    return index.match(pname)


def _get_fnames_right(right, index):
    fnames_right = {}
    for match in finditer(r"\[(.+?)\]", right.strip()):
        pname = ParamName.from_section(match.group(1))
        fnames_right[match.group(0)] = index.match(pname)
    return fnames_right


def read_grid(grid, params):
    if grid is None:
        return None, params
    re_ = compile(
        fr"(lin[(]{_RE_FLOAT},{_RE_FLOAT},\d+[)]$)|"
        fr"(log[(]{_RE_FLOAT},{_RE_FLOAT},\d+[)]$)|"
        fr"([(](({_RE_FLOAT})(,|[)]$))+)"
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
        params_remaining = {fname: params[fname] for fname in fnames_all}
        fnames_left = _get_fnames_left(left, params_remaining)
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

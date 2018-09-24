"""The parameters module contains the code for handling of the experimental and
fitting parameters."""

import ast
import configparser
import copy
import os
import os.path
import re
from difflib import get_close_matches

import lmfit
import numpy as np
from asteval import astutils

from chemex import peaks, util

NAME_MARKERS = {
    "name": "__n_{}_n__",
    "nuclei": "__r_{}_r__",
    "temperature": "__t_{}_t__",
    "h_larmor_frq": "__b_{}_b__",
    "p_total": "__p_{}_p__",
    "l_total": "__l_{}_l__",
}

FRIENDLY_MARKERS = {
    "name": "{}",
    "nuclei": "NUC->{}",
    "temperature": "T->{:.1f}C",
    "h_larmor_frq": "B0->{:.1f}MHz",
    "p_total": "[P]->{:e}M",
    "l_total": "[L]->{:e}M",
}

RE_QUALIFIERS = re.compile(
    """
        (^\s*(?P<name>\w+)) |
        (NUC\s*->\s*(?P<nuclei>(\w|-)+)) |
        (T\s*->s*(?P<temperature>{0})) |
        (B0\s*->\s*(?P<h_larmor_frq>{0})) |
        (\[P\]\s*->\s*(?P<p_total>{0})) |
        (\[L\]\s*->\s*(?P<l_total>{0})) |
    """.format(
        "[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?"
    ),
    re.IGNORECASE | re.VERBOSE,
)

# Regular expression to pick values of the form: intial value [min, max]
RE_VALUE_MIN_MAX = re.compile(
    """
        ^\s*
        (?P<value>{0})?\s*
        (\[\s*(?P<min>({0}|{1}))\s*,\s*(?P<max>({0}|{1}))\s*\]\s*)?
        .*$
    """.format(
        "[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?", "[-+]?inf"
    ),
    re.IGNORECASE | re.VERBOSE,
)

# Regular expression to pick values of the form: intial value [min, max, brute_stepsize]
RE_VALUE_MIN_MAX_BRUTE = re.compile(
    """
        ^\s*
        (?P<value>{0})?\s*
        (\[\s*(?P<min>({0}|{1}))\s*,\s*(?P<max>({0}|{1}))\s*,\s*(?P<brute_step>({0}|None))\s*\]\s*)?
        .*$
    """.format(
        "[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?", "[-+]?inf"
    ),
    re.IGNORECASE | re.VERBOSE,
)

RE_PARNAME = re.compile(
    """
        (__n_(?P<name>\w+)_n__)?
        (__r_(?P<nuclei>(\w|-)+)_r__)?
        (__t_(?P<temperature>{0})_t__)?
        (__b_(?P<h_larmor_frq>{0})_b__)?
        (__p_(?P<p_total>{0})_p__)?
        (__l_(?P<l_total>{0})_l__)?
    """.format(
        "[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?"
    ),
    re.IGNORECASE | re.VERBOSE,
)


class MakeTranslate(object):
    """MakeTranslate class for translating parameter names."""

    def __init__(self, *args, **kwds):
        self.dictionary = dict(*args, **kwds)
        self.rx = self.make_rx()

    def make_rx(self):
        """TODO: method docstring."""
        return re.compile("|".join(map(re.escape, self.dictionary)), re.IGNORECASE)

    def one_xlat(self, match):
        """TODO: method docstring."""
        return self.dictionary[match.group(0)]

    def __call__(self, text):
        return self.rx.sub(self.one_xlat, text)


expand = MakeTranslate({"-": "__minus__", "+": "__plus__", ".": "__point__"})
compress = MakeTranslate({"__minus__": "-", "__plus__": "+", "__point__": "."})


class ParameterName(object):
    """ParameterName class."""

    def __init__(
        self,
        name=None,
        nuclei=None,
        temperature=None,
        h_larmor_frq=None,
        p_total=None,
        l_total=None,
    ):
        self.name = self.nuclei = self.temperature = self.h_larmor_frq = None
        self.p_total = self.l_total = None

        if name is not None:
            self.name = name.lower()

        if nuclei is not None:
            self.nuclei = peaks.Peak(nuclei).assignment

        if temperature is not None:
            self.temperature = round(float(temperature), 1)

        if h_larmor_frq is not None:
            self.h_larmor_frq = round(float(h_larmor_frq), 1)

        if p_total is not None:
            self.p_total = float(p_total)

        if l_total is not None:
            self.l_total = float(l_total)

    @classmethod
    def from_full_name(cls, full_name=None):
        """TODO: method docstring."""
        if full_name is None:
            full_name = ""

        full_name = compress(full_name)

        match = re.match(RE_PARNAME, full_name)
        qualifiers = {}
        if match is not None:
            qualifiers.update(match.groupdict())

        return cls(**qualifiers)

    @classmethod
    def from_section(cls, section=None):
        """TODO: method docstring."""
        if section is None:
            section = ""
        qualifiers = re_to_dict(RE_QUALIFIERS, section)

        return cls(**qualifiers)

    def update_nuclei(self, nuclei=None):
        """TODO: method docstring."""
        if nuclei is not None:
            self.nuclei = peaks.Peak(nuclei).assignment

        return self

    def to_full_name(self):
        """TODO: method docstring."""
        name_components = []

        for attribute, value in vars(self).items():
            if value is not None:
                name_components.append(NAME_MARKERS[attribute].format(value))

        full_name = expand("".join(name_components))

        return full_name

    def to_section_name(self, nuclei=False):
        """TODO: method docstring."""

        name_components = []

        for attribute, value in vars(self).items():
            if (attribute != "nuclei" or nuclei) and value is not None:
                name_components.append(FRIENDLY_MARKERS[attribute].format(value))

        section_name = ", ".join(name_components).upper()

        return section_name

    def to_re(self):
        """TODO: method docstring."""
        re_components = [NAME_MARKERS["name"].format(expand(self.name))]

        if self.nuclei is not None:
            group_name = peaks.Peak(self.nuclei)._resonances["i"]["group"]
            if not group_name:
                all_res = "\D?[0-9]+[abd-gi-mopr-z]*"
            else:
                all_res = ""
            re_components.append(
                NAME_MARKERS["nuclei"].format("".join([all_res, expand(self.nuclei)]))
            )
        else:
            re_components.append(".*")

        if self.temperature is not None:
            re_components.append(
                NAME_MARKERS["temperature"].format(expand(str(self.temperature)))
            )
        elif re_components[-1] != ".*":
            re_components.append(".*")

        if self.h_larmor_frq is not None:
            re_components.append(
                NAME_MARKERS["h_larmor_frq"].format(expand(str(self.h_larmor_frq)))
            )
        elif re_components[-1] != ".*":
            re_components.append(".*")

        if self.p_total is not None:
            re_components.append(
                NAME_MARKERS["p_total"].format(expand(str(self.p_total)))
            )
        elif re_components[-1] != ".*":
            re_components.append(".*")

        if self.l_total is not None:
            re_components.append(
                NAME_MARKERS["l_total"].format(expand(str(self.l_total)))
            )
        elif re_components[-1] != ".*":
            re_components.append(".*")

        re_to_match = re.compile("".join(re_components), re.IGNORECASE)

        return re_to_match

    def __repr__(self):
        return self.to_section_name(nuclei=True)

    def __lt__(self, other):
        tuple_self = (
            self.name,
            self.temperature,
            self.h_larmor_frq,
            self.p_total,
            self.l_total,
            peaks.Peak(self.nuclei),
        )

        tuple_other = (
            other.name,
            other.temperature,
            other.h_larmor_frq,
            other.p_total,
            other.l_total,
            peaks.Peak(other.nuclei),
        )

        return tuple_self < tuple_other

    def intersection(self, other):
        """TODO: method docstring."""
        name = temperature = h_larmor_frq = p_total = l_total = None

        if self.name == other.name:
            name = self.name

        nuclei = (
            peaks.Peak(self.nuclei).intersection(peaks.Peak(other.nuclei)).assignment
        )

        if not nuclei:
            nuclei = None

        if self.temperature == other.temperature:
            temperature = self.temperature

        if self.h_larmor_frq == other.h_larmor_frq:
            h_larmor_frq = self.h_larmor_frq

        if self.p_total == other.p_total:
            p_total = self.p_total

        if self.l_total == other.l_total:
            l_total = self.l_total

        return ParameterName(
            name=name,
            nuclei=nuclei,
            temperature=temperature,
            h_larmor_frq=h_larmor_frq,
            p_total=p_total,
            l_total=l_total,
        )


def create_params(data):
    """Create the array of fitting parameters."""
    params = lmfit.Parameters()

    for profile in data:
        for name, param in profile.params.items():
            if name in params:
                vary = params[name].vary
            else:
                vary = False
            param._delay_asteval = True
            params[name] = param
            params[name].vary = vary | param.vary

    for p in params.values():
        p._delay_asteval = False

    return params


def set_params_from_config_file(params, config_filename):
    """Read the parameter file and set initial values and optional bounds and brute step size."""

    print("File Name: {}".format(config_filename), end="\n\n")

    config = util.read_cfg_file(config_filename)

    print("{:<45s} {:<30s}".format("Section", "Matches"))
    print("{:<45s} {:<30s}".format("-------", "-------"))

    for section in config.sections():
        if section.lower() in ("global", "default"):
            print("{:<45s}".format("[{}]".format(section)))

            for key, value in config.items(section):
                name = ParameterName.from_section(key)
                if value.count(",") == 2:
                    default = re_to_dict(RE_VALUE_MIN_MAX_BRUTE, value)
                else:
                    default = re_to_dict(RE_VALUE_MIN_MAX, value)
                default = {key: np.float64(val) for key, val in default.items()}
                matches = set_params(params, name, **default)

                print("  {:<43s} {:<30d}".format("({})".format(key), len(matches)))

        else:
            name = ParameterName.from_section(section)

            pairs = []

            for key, value in config.items(section):
                if "file" in key:
                    for filename in value.split():
                        filename_ = util.normalize_path(
                            os.path.dirname(config_filename), filename
                        )
                        pairs.extend(get_pairs_from_file(filename_, name))

                elif peaks.re_peak_name.match(key):
                    name.update_nuclei(key)
                    pairs.append((copy.deepcopy(name), value))

            total_matches = set()

            for name, value in pairs:
                if value.count(",") == 2:
                    default = re_to_dict(RE_VALUE_MIN_MAX_BRUTE, value)
                else:
                    default = re_to_dict(RE_VALUE_MIN_MAX, value)
                default = {key: np.float64(val) for key, val in default.items()}
                matches = set_params(params, name, **default)
                total_matches.update(matches)

            print("{:<45s} {:<30d}".format("[{}]".format(section), len(total_matches)))


def get_pairs_from_file(filename, name):
    """Read residue specific values for fitting parameters from a file.

    The file should be formatted like a Sparky peak list.
    Examples:
      * To set G23N to 105.0 and G23H to 8.0:
          G23N-H  105.0  8.0
      * To set a parameter depending on multiple nuclei (e.g., G23N and G23H):
          G23N-H  -93.0

    """
    pairs = []

    with open(filename) as f:
        for line in f:
            if "Assignment" in line:
                continue

            line = remove_comments(line, "#;")
            line = re.sub("\s*\[\s*", "[", line)
            line = re.sub("\s*\]\s*", "]", line)
            elements = line.split()

            if len(elements) > 1:
                peak = peaks.Peak(elements[0])
                n_resonances = len(peak)
                n_cols = len(elements[1:])

                if n_cols == n_resonances:
                    for nuc_name, value in zip(peak.names.values(), elements[1:]):
                        name.update_nuclei(nuc_name)
                        pairs.append((copy.deepcopy(name), value))

                else:
                    name.update_nuclei(peak.assignment)
                    pairs.append((copy.deepcopy(name), elements[1]))

    return pairs


def set_param_status(params, items):
    """Set whether or not to vary a fitting parameter or to use a mathemetical expression."""

    vary = {"fix": False, "fit": True}

    for key, status in items:
        name = ParameterName.from_section(key)

        if status in vary:
            set_params(params, name, vary=vary[status])
        else:
            set_param_expr(params, name, expr=status)


def set_param_expr(params, name, expr=None):
    """Set an optional parameter expression, used to constrain its value during
    the fit."""

    if expr is None:
        expr = ""

    if not isinstance(name, ParameterName):
        name = ParameterName.from_section(name)

    names_full = [name_full for name_full in params if name.to_re().match(name_full)]
    names_expr = astutils.get_ast_names(ast.parse(expr))
    names_full_expr = {
        name: [
            name_full_expr
            for name_full_expr in params
            if ParameterName.from_section(name).to_re().match(name_full_expr)
        ]
        for name in names_expr
    }

    matches = set()

    for name_full in names_full:

        expr_ = expr

        for name_expr in names_expr:
            name_full_expr = get_close_matches(
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
    """Set the initial value and (optional) bounds and brute step size for parameters."""
    matches = set()
    name_short_re = name_short.to_re()

    for name, param in params.items():
        if name_short_re.match(name):
            if expr is None and param.expr and vary is None:
                param.value = value
            else:
                param.set(value, vary, min, max, expr, brute_step)
            matches.add(name)

    return matches


def write_par(params, output_dir="./"):
    """Write the fitting parameters and their uncertainties to a file."""
    filename = os.path.join(output_dir, "parameters.fit")

    print("  * {}".format(filename))

    par_dict = {}

    for name, param in params.items():

        par_name = ParameterName.from_full_name(name)

        if par_name.nuclei is None:  # global parameter
            name_print = par_name
            section = "GLOBAL"

        else:  # residue-specific parameter
            name_print = peaks.Peak(par_name.nuclei)
            section = par_name.to_section_name()

        if not param.vary and param.expr is None:
            val_print = "{: .5e} ; fixed".format(param.value)
        elif param.stderr is None:
            val_print = "{: .5e}  ; error not calculated".format(param.value)
        elif param.expr:
            val_print = "{: .5e} +/- {:.5e} ; constrained".format(
                param.value, param.stderr
            )
        else:
            val_print = "{: .5e} +/- {:.5e}".format(param.value, param.stderr)

        par_dict.setdefault(section, []).append((name_print, val_print))

    cfg = configparser.ConfigParser()
    cfg.optionxform = str

    section_global = par_dict.pop("GLOBAL", None)

    if section_global is not None:
        cfg.add_section("GLOBAL")
        for name, val in sorted(section_global):
            cfg.set("GLOBAL", str(name), val)

    for section, name_vals in sorted(par_dict.items()):
        cfg.add_section(section)
        for peak, val in sorted(name_vals):
            cfg.set(section, peak.assignment.upper(), val)

    with open(filename, "w") as f:
        cfg.write(f)


def write_constraints(params, output_dir="./"):
    """Write the (optional) parameter expression constraints to a file."""
    filename = os.path.join(output_dir, "constraints.fit")

    print("  * {}".format(filename))

    param_dict = dict()

    for name, param in params.items():

        par_name = ParameterName.from_full_name(name)

        if param.expr:
            name_formatted = "[{}]".format(par_name.to_section_name(nuclei=True))
            expr_formatted = param.expr
            for name_dep in astutils.get_ast_names(ast.parse(param.expr)):
                par_name_dep = ParameterName.from_full_name(name_dep)
                if str(par_name_dep):
                    expr_formatted = expr_formatted.replace(
                        name_dep,
                        "[{}]".format(par_name_dep.to_section_name(nuclei=True)),
                    )

            param_dict[par_name] = "{} = {}\n".format(name_formatted, expr_formatted)

    with open(filename, "w") as f:
        for name, constraint in sorted(param_dict.items()):
            f.write(constraint)


def remove_comments(line, sep):
    """Remove (optional) comments."""
    for s in sep:
        line = line.split(s)[0]

    return line.strip()


def re_to_dict(re_to_match, text):
    """TODO: function docstring."""
    return {
        match_key: match_value
        for match in re_to_match.finditer(text)
        for match_key, match_value in match.groupdict().items()
        if match_value is not None
    }

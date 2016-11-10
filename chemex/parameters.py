import ast
import configparser
import copy
import os
import os.path
import re

import lmfit
import numpy as np
from lmfit import astutils

from chemex import peaks, util

name_markers = {
    'name'        : "__n_{}_n__",
    'nuclei'      : "__r_{}_r__",
    'temperature' : "__t_{}_t__",
    'h_larmor_frq': "__b_{}_b__",
    'p_total'     : "__p_{}_p__",
    'l_total'     : "__l_{}_l__",
}

friendly_markers = {
    'name'        : "{}",
    'nuclei'      : "NUC->{}",
    'temperature' : "T->{:.1f}C",
    'h_larmor_frq': "B0->{:.1f}MHz",
    'p_total'     : "[P]->{:e}M",
    'l_total'     : "[L]->{:e}M",
}

re_qualifiers = re.compile(
    '''
        (^\s*(?P<name>\w+)) |
        (NUC\s*->\s*(?P<nuclei>(\w|-)+)) |
        (T\s*->s*(?P<temperature>{0})) |
        (B0\s*->\s*(?P<h_larmor_frq>{0})) |
        (\[P\]\s*->\s*(?P<p_total>{0})) |
        (\[L\]\s*->\s*(?P<l_total>{0})) |
    '''.format('[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?'),
    re.IGNORECASE | re.VERBOSE
)

# Regular expression to pick values of the form: 8.3 [2.0, 10.0]
re_value_min_max = re.compile(
    '''
        ^\s*
        (?P<value>{0})?\s*
        (\[\s*(?P<min>({0}|{1}))\s*,\s*(?P<max>({0}|{1}))\s*\]\s*)?
        .*$
    '''.format('[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?', '[-+]?inf'),
    re.IGNORECASE | re.VERBOSE
)

re_par_name = re.compile(
    '''
        (__n_(?P<name>\w+)_n__)?
        (__r_(?P<nuclei>(\w|-)+)_r__)?
        (__t_(?P<temperature>{0})_t__)?
        (__b_(?P<h_larmor_frq>{0})_b__)?
        (__p_(?P<p_total>{0})_p__)?
        (__l_(?P<l_total>{0})_l__)?
    '''.format('[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?'),
    re.IGNORECASE | re.VERBOSE
)


class MakeTranslate(object):
    def __init__(self, *args, **kwds):
        self.dictionary = dict(*args, **kwds)
        self.rx = self.make_rx()

    def make_rx(self):
        return re.compile('|'.join(map(re.escape, self.dictionary)), re.IGNORECASE)

    def one_xlat(self, match):
        return self.dictionary[match.group(0)]

    def __call__(self, text):
        return self.rx.sub(self.one_xlat, text)


expand = MakeTranslate({
    '-': '__minus__',
    '+': '__plus__',
    '.': '__point__',
})

compress = MakeTranslate({
    '__minus__': '-',
    '__plus__' : '+',
    '__point__': '.',
})


class ParameterName(object):
    def __init__(self, name=None, nuclei=None, temperature=None, h_larmor_frq=None, p_total=None, l_total=None):

        self.name = self.nuclei = self.temperature = self.h_larmor_frq = self.p_total = self.l_total = None

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
        if full_name is None:
            full_name = ''

        full_name = compress(full_name)

        match = re.match(re_par_name, full_name)
        qualifiers = {}
        if match is not None:
            qualifiers.update(match.groupdict())

        return cls(**qualifiers)

    @classmethod
    def from_section(cls, section=None):
        if section is None:
            section = ''
        qualifiers = re_to_dict(re_qualifiers, section)
        return cls(**qualifiers)

    def update(self, other):
        if isinstance(other, ParameterName):
            if other.name is not None:
                self.name = other.name
            if other.nuclei is not None:
                self.nuclei = other.nuclei
            if other.temperature is not None:
                self.temperature = other.temperature
            if other.h_larmor_frq is not None:
                self.h_larmor_frq = other.h_larmor_frq
            if other.p_total is not None:
                self.p_total = other.p_total
            if other.l_total is not None:
                self.l_total = other.l_total
        return self

    def update_nuclei(self, nuclei=None):
        if nuclei is not None:
            self.nuclei = peaks.Peak(nuclei).assignment
        return self

    def to_full_name(self):

        name_components = []

        for name in ('name', 'nuclei', 'temperature', 'h_larmor_frq', 'p_total', 'l_total'):
            attr = getattr(self, name)
            if attr is not None:
                name_components.append(name_markers[name].format(attr))

        full_name = expand(''.join(name_components))

        return full_name

    def to_section_name(self, nuclei=False):

        name_components = []

        if nuclei:
            attributes = ('name', 'nuclei', 'temperature', 'h_larmor_frq', 'p_total', 'l_total')
        else:
            attributes = ('name', 'temperature', 'h_larmor_frq', 'p_total', 'l_total')

        for name in attributes:
            attr = getattr(self, name)
            if attr is not None:
                name_components.append(friendly_markers[name].format(attr))

        section_name = ', '.join(name_components).upper()

        return section_name

    def to_re(self):

        re_components = [name_markers['name'].format(expand(self.name))]

        if self.nuclei is not None:
            group_name = peaks.Peak(self.nuclei).resonances[0]['group']
            if not group_name:
                all_res = '\D?[0-9]+[abd-gi-mopr-z]*'
            else:
                all_res = ''
            re_components.append(name_markers['nuclei'].format(''.join([all_res, expand(self.nuclei)])))
        else:
            re_components.append('.*')

        if self.temperature is not None:
            re_components.append(name_markers['temperature'].format(expand(str(self.temperature))))
        elif re_components[-1] != '.*':
            re_components.append('.*')

        if self.h_larmor_frq is not None:
            re_components.append(name_markers['h_larmor_frq'].format(expand(str(self.h_larmor_frq))))
        elif re_components[-1] != '.*':
            re_components.append('.*')

        if self.p_total is not None:
            re_components.append(name_markers['p_total'].format(expand(str(self.p_total))))
        elif re_components[-1] != '.*':
            re_components.append('.*')

        if self.l_total is not None:
            re_components.append(name_markers['l_total'].format(expand(str(self.l_total))))
        elif re_components[-1] != '.*':
            re_components.append('.*')

        re_to_match = re.compile(''.join(re_components), re.IGNORECASE)

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
            peaks.Peak(self.nuclei)
        )

        tuple_other = (
            other.name,
            other.temperature,
            other.h_larmor_frq,
            other.p_total,
            other.l_total,
            peaks.Peak(other.nuclei)
        )

        return tuple_self < tuple_other

    def intersection(self, other):

        name = temperature = h_larmor_frq = p_total = l_total = None

        if self.name == other.name:
            name = self.name

        nuclei = peaks.Peak(self.nuclei).intersection(peaks.Peak(other.nuclei)).assignment

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

        return ParameterName(name=name, nuclei=nuclei, temperature=temperature, h_larmor_frq=h_larmor_frq,
                             p_total=p_total, l_total=l_total)


def create_params(data):
    """
    Creates the array of parameters that will be used for the fitting
    along with the dictionary that associate the name and the index of each
    parameter in the array.
    """

    params = lmfit.Parameters()

    for profile in data:
        for name, param in profile.default_params.items():
            params[name] = lmfit.Parameter(name)

    for profile in data:
        for name, param in profile.default_params.items():
            params[name] = param

    return params


def set_params_from_config_file(params, config_filename):
    """
    Read the file containing the initial guess for the fitting parameters.
    """

    print("File Name: {:s}".format(config_filename), end='\n\n')

    config = util.read_cfg_file(config_filename)

    print("{:<45s} {:<30s}".format("Section", "Matches"))
    print("{:<45s} {:<30s}".format("-------", "-------"))

    for section in config.sections():

        if section.lower() in ('global', 'default'):

            print("{:<45s}".format("[{}]".format(section)))

            for key, value in config.items(section):
                name = ParameterName.from_section(key)
                default = re_to_dict(re_value_min_max, value)
                default = {key: np.float64(val) for key, val in default.items()}
                matches = set_params(params, name, **default)

                print("  {:<43s} {:<30d}".format("({})".format(key), len(matches)))

        else:

            name = ParameterName.from_section(section)

            pairs = []

            for key, value in config.items(section):

                if 'file' in key:
                    for filename in value.split():
                        filename_ = util.normalize_path(os.path.dirname(config_filename), filename)
                        pairs.extend(get_pairs_from_file(filename_, name))

                elif peaks.re_peak_name.match(key):
                    name.update_nuclei(key)
                    pairs.append((copy.deepcopy(name), value))

            total_matches = set()

            for name, value in pairs:
                default = re_to_dict(re_value_min_max, value)
                default = {key: np.float64(val) for key, val in default.items()}
                matches = set_params(params, name, **default)
                total_matches.update(matches)

            print("{:<45s} {:<30d}".format("[{}]".format(section), len(total_matches)))


def get_pairs_from_file(filename, name):
    """Reads a file containing values associated with a nucleus name.
    The file should be formatted like sparky peak lists.

    For example:
      * To set G23N to 105.0 and G23H to 8.0:
          G23N-H  105.0  8.0
      * To set a parameter depending on multiple nuclei, eg. G23N and G23H, to 1.0:
          G23N-H  -93.0
    """

    pairs = []

    with open(filename) as f:

        for line in f:

            if 'Assignment' in line:
                continue

            line = remove_comments(line, '#;')
            line = re.sub('\s*\[\s*', '[', line)
            line = re.sub('\s*\]\s*', ']', line)
            elements = line.split()

            if len(elements) > 1:

                peak = peaks.Peak(elements[0])
                n_resonances = len(peak.resonances)
                n_cols = len(elements[1:])

                if n_cols == n_resonances:
                    for resonance, value in zip(peak.resonances, elements[1:]):
                        name.update_nuclei(resonance['name'])
                        pairs.append((copy.deepcopy(name), value))

                else:
                    name.update_nuclei(peak.assignment)
                    pairs.append((copy.deepcopy(name), elements[1]))
    return pairs


def set_param_status(params, items):
    """Fix (or not) fit variables according to what set in the protocol file"""

    vary = {'fix': False, 'fit': True}

    for key, status in items:
        name = ParameterName.from_section(key)
        if status in vary:
            set_params(params, name, vary=vary[status])
        else:
            name_short_expr = ParameterName.from_section(status)
            set_param_expr(params, name, name_short_expr=name_short_expr)


def set_param_expr(params, name_short, name_short_expr=None):
    matches = set()
    name_short_re = name_short.to_re()

    for name, param in params.items():

        if name_short_re.match(name):

            expr = None

            if name_short_expr == '':
                expr = name_short_expr

            elif name_short_expr is not None:
                name_expr = ParameterName.from_full_name(name).update(name_short_expr).to_full_name()
                if name_expr != name and name_expr in params:
                    min = param.min
                    max = param.max
                    param.set(expr=name_expr, min=min, max=max)

            if expr is not None:
                matches.add(name)

    return matches


def set_params(params, name_short, value=None, vary=None, min=None, max=None):
    matches = set()
    name_short_re = name_short.to_re()

    for name, param in params.items():

        if name_short_re.match(name):

            if not param.expr or vary is not None:
                min = param.min if min is None else min
                max = param.max if max is None else max
                param.set(value=value, vary=vary, min=min, max=max)

            matches.add(name)

    return matches


def write_par(params, output_dir='./'):
    """Write fitted parameters int a file"""

    filename = os.path.join(output_dir, 'parameters.fit')

    print("  * {}".format(filename))

    par_dict = {}

    for name, param in params.items():

        if not param.vary and param.expr is None:
            val_print = '{: .5e} ; fixed'.format(param.value)
        elif param.stderr is None:
            val_print = '{: .5e}  ; error not calculated'.format(param.value)
        elif param.expr:
            val_print = '{: .5e} +/- {:.5e} ; constrained'.format(param.value, param.stderr)
        else:
            val_print = '{: .5e} +/- {:.5e}'.format(param.value, param.stderr)

        par_name = ParameterName.from_full_name(name)

        if par_name.nuclei is None:

            # This is a non-residue-specific parameter

            name_print = par_name
            section = 'GLOBAL'

        else:

            # This is a residue-specific parameter

            name_print = peaks.Peak(par_name.nuclei)
            section = par_name.to_section_name()

        par_dict.setdefault(section, []).append((name_print, val_print))

    cfg = configparser.ConfigParser()
    cfg.optionxform = str

    section_global = par_dict.pop('GLOBAL', None)

    if section_global is not None:
        cfg.add_section('GLOBAL')
        for name, val in sorted(section_global):
            cfg.set('GLOBAL', str(name), val)

    for section, name_vals in sorted(par_dict.items()):
        cfg.add_section(section)
        for peak, val in sorted(name_vals):
            cfg.set(section, peak.assignment.upper(), val)

    with open(filename, 'w') as f:
        cfg.write(f)


def write_constraints(params, output_dir='./'):
    """Write fitted parameters int a file"""

    filename = os.path.join(output_dir, 'constraints.fit')

    print("  * {}".format(filename))

    param_dict = dict()

    for name, param in params.items():
        if param.expr:
            name_formatted = '[{}]'.format(ParameterName.from_full_name(name).to_section_name(nuclei=True))
            expr_formatted = param.expr
            for name_dep in astutils.get_ast_names(ast.parse(param.expr)):
                expr_formatted = expr_formatted.replace(
                    name_dep,
                    '[{}]'.format(ParameterName.from_full_name(name_dep).to_section_name(nuclei=True))
                )
            param_dict[ParameterName.from_full_name(name)] = '{} = {}\n'.format(name_formatted, expr_formatted)

    with open(filename, 'w') as f:
        for name, constraint in sorted(param_dict.items()):
            f.write(constraint)


def remove_comments(line, sep):
    for s in sep:
        line = line.split(s)[0]
    return line.strip()


def re_to_dict(re_to_match, text):
    return {
        match_key: match_value
        for match in re_to_match.finditer(text)
        for match_key, match_value in match.groupdict().items()
        if match_value is not None
        }


def main():
    pass


if __name__ == '__main__':
    main()

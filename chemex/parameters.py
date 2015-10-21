from __future__ import print_function
from __future__ import absolute_import

import copy
import os
import os.path
import re

import numpy as np
import lmfit

from six.moves import map
from six.moves import zip
from chemex import util, peaks

name_markers = {
    'name': "__n_{}_n__",
    'nuclei': "__r_{}_r__",
    'temperature': "__t_{}_t__",
    'h_larmor_frq': "__b_{}_b__",
}

friendly_markers = {
    'name': "{}",
    'nuclei': "{}",
    'temperature': "{} C",
    'h_larmor_frq': "{} MHz",
}

re_qualifiers = re.compile(
    '''
        (^\s*(?P<name>\w+)) |
        (?P<nuclei>(\D?\d+[abd-gi-mopr-z]*)?[hncq][a-z0-9]*) |
        ((?P<temperature>{0})\s*C) |
        ((?P<h_larmor_frq>{0})\s*MHz)
    '''.format('[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?'),
    re.IGNORECASE | re.VERBOSE
)

# Regular expression to pick values of the form: 8.3 [2.0, 10.0]
re_value_min_max = re.compile(
    '''
        ^\s*
        (?P<value>{0})?\s*
        (\[\s*(?P<min>({0}|{1}))\s*,\s*(?P<max>({0}|{1}))\s*\]\s*)?
        $
    '''.format('[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?', '[-+]?inf'),
    re.IGNORECASE | re.VERBOSE
)

re_par_name = re.compile(
    '''
        (__n_(?P<name>.+)_n__)?
        (__r_(?P<nuclei>.+)_r__)?
        (__t_(?P<temperature>.+)_t__)?
        (__b_(?P<h_larmor_frq>.+)_b__)?
    ''',
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
    '__plus__': '+',
    '__point__': '.',
})


class ParameterName(object):
    def __init__(self, name=None, nuclei=None, temperature=None, h_larmor_frq=None):
        self.name = name
        self.temperature = temperature
        self.h_larmor_frq = h_larmor_frq

        if self.name is not None:
            self.name = self.name.lower()

        if nuclei is not None:
            self.nuclei = peaks.Peak(nuclei).assignment
        else:
            self.nuclei = None

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
        return self

    def update_nuclei(self, nuclei=None):
        if nuclei is not None:
            self.nuclei = peaks.Peak(nuclei).assignment
        return self

    def to_full_name(self):

        name_components = []

        for name in ('name', 'nuclei', 'temperature', 'h_larmor_frq'):
            attr = getattr(self, name)
            if attr is not None:
                name_components.append(name_markers[name].format(attr))

        full_name = ''.join(name_components)

        full_name = expand(full_name)

        return full_name

    def to_section_name(self):

        name_components = []

        for name in ('name', 'temperature', 'h_larmor_frq'):
            attr = getattr(self, name)
            if attr is not None:
                name_components.append(friendly_markers[name].format(attr.upper()))

        section_name = ', '.join(name_components)

        return section_name

    def to_re(self):

        re_components = [name_markers['name'].format(expand(self.name))]
        if self.nuclei is not None:
            group_name = peaks.Peak(self.nuclei).resonances[0]['group']
            if not group_name:
                all_res = '.+'
            else:
                all_res = ''
            re_components.append(name_markers['nuclei'].format(''.join([all_res, expand(self.nuclei)])))
        else:
            re_components.append('.*')
        if self.temperature is not None:
            re_components.append(name_markers['temperature'].format(expand(self.temperature)))
        elif re_components[-1] != '.*':
            re_components.append('.*')
        if self.h_larmor_frq is not None:
            re_components.append(name_markers['h_larmor_frq'].format(expand(self.h_larmor_frq)))
        elif re_components[-1] != '.*':
            re_components.append('.*')

        re_to_match = re.compile(''.join(re_components))

        return re_to_match

    def __repr__(self):
        name_components = []

        for name in ('name', 'nuclei', 'temperature', 'h_larmor_frq'):
            attr = getattr(self, name)
            if attr is not None:
                name_components.append(friendly_markers[name].format(attr.upper()))

        nice_name = ', '.join(name_components)

        return nice_name

    def __lt__(self, other):

        tuple_self = (
            self.name,
            peaks.Peak(self.nuclei),
            float(self.temperature) if self.temperature else self.temperature,
            float(self.h_larmor_frq) if self.h_larmor_frq else self.h_larmor_frq
        )
        tuple_other = (
            other.name,
            peaks.Peak(other.nuclei),
            float(other.temperature) if other.temperature else other.temperature,
            float(other.h_larmor_frq) if other.h_larmor_frq else other.h_larmor_frq
        )
        return tuple_self < tuple_other

    def intersection(self, other):
        name = nuclei = temperature = h_larmor_frq = None
        if self.name == other.name:
            name = self.name
        if self.nuclei == other.nuclei:
            nuclei = self.nuclei
        if self.temperature == other.temperature:
            temperature = self.temperature
        if self.h_larmor_frq == other.h_larmor_frq:
            h_larmor_frq = self.h_larmor_frq
        return ParameterName(name=name, nuclei=nuclei, temperature=temperature, h_larmor_frq=h_larmor_frq)


def create_params(data):
    """
    Creates the array of parameters that will be used for the fitting
    along with the dictionary that associate the name and the index of each
    parameter in the array.
    """

    params = lmfit.Parameters()

    for profile in data:
        for name, param in profile.create_default_parameters().items():
            if name not in params or not params[name].vary:
                params.update({name: param})

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
            pairs = []
            name = ParameterName.from_section(section)
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
            name_other = ParameterName.from_section(status)
            set_params(params, name, name_other=name_other)


def set_params(params, name, value=None, vary=None, min=None, max=None, name_other=None):
    matches = set()
    re_to_match = name.to_re()
    for name_, param in params.items():
        if re_to_match.match(name_):
            expr = None
            if name_other is not None:
                name_updated = ParameterName.from_full_name(name_).update(name_other).to_full_name()
                if name_updated != name_ and name_updated in params:
                    expr = name_updated
            if (value, vary, min, max, expr) is not (None, None, None, None, None):
                param.set(value=value, vary=vary, min=min, max=max, expr=expr)
                matches.add(name_)

    return matches


def write_par(params, output_dir='./'):
    """Write fitted parameters int a file"""

    from six.moves.configparser import ConfigParser

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

    cfg = ConfigParser()
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

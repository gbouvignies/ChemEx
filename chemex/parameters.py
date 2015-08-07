from __future__ import print_function

import os
import os.path
import re
import copy

import lmfit

from chemex import util, peaks

name_markers = {
    'name': "__n_{}_n__",
    'nuclei': "__r_{}_r__",
    'temperature': "__t_{}_t__",
    'h_larmor_frq': "__b_{}_b__",
}

friendly_markers = {
    'name': "{}",
    'temperature': "{} C",
    'h_larmor_frq': "{} MHz",
}

re_qualifiers = re.compile(
    '''
        (^\s*(?P<name>\w+)) |
        (?P<nuclei>(\D?\d+[abd-gi-mopr-z]*)?[hncq][a-z0-9]*) |
        ((?P<temperature>[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?)\s*C) |
        ((?P<h_larmor_frq>[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?)\s*MHz)
    ''',
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
        self.name = name.lower()
        self.temperature = temperature
        self.h_larmor_frq = h_larmor_frq

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
        qualifiers = {
            match_key: match_value
            for match in re_qualifiers.finditer(section)
            for match_key, match_value in match.groupdict().items()
            if match_value is not None}
        return cls(**qualifiers)

    def update_nuclei(self, nuclei):
        if nuclei is not None:
            self.nuclei = peaks.Peak(nuclei).assignment

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

        if section.lower() == 'global':
            print("{:<45s}".format("[global]"))
            for key, value in config.items(section):
                name = ParameterName.from_section(key)
                value_ = float(value.split()[0])
                matches = set_params(params, name, value=value_)

                print("    {:<41s} {:<30d}".format(key, len(matches)))

        else:
            pairs = []
            name = ParameterName.from_section(section)
            for key, value in config.items(section):
                if 'file' in key:
                    for filename in value.split():
                        filename_ = util.normalize_path(os.path.dirname(config_filename), filename)
                        pairs.extend(get_pairs_from_file(filename_, name))
                else:
                    name.update_nuclei(key)
                    value = float(value.split()[0])
                    pairs.append((name, value))

            total_matches = set()
            for name, value in pairs:
                matches = set_params(params, name, value=value)
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
            elements = line.split()
            if len(elements) > 1:
                peak = peaks.Peak(elements[0])
                n_resonances = len(peak.resonances)
                n_cols = len(elements[1:])
                if n_cols == n_resonances:
                    for resonance, value in zip(peak.resonances, elements[1:]):
                        name.update_nuclei(resonance['name'])
                        pairs.append((copy.deepcopy(name), float(value)))
                else:
                    name.update_nuclei(peak.assignment)
                    pairs.append((copy.deepcopy(name), float(elements[1])))
    return pairs


def set_param_status(params, items):
    """Fix (or not) fit variables according to what set in the protocol file"""

    vary = {'fix': False, 'fit': True}

    for key, status in items:
        name = ParameterName.from_section(key)
        set_params(params, name, vary=vary[status])


def set_params(params, name, value=None, vary=None):
    matches = set()
    re_to_match = name.to_re()
    for name_, param in params.items():
        if re_to_match.match(name_):
            if value is not None:
                param.set(value=value)
            if vary is not None:
                param.set(vary=vary)
            matches.add(name_)

    return matches


def write_par(params, output_dir='./'):
    """Write fitted parameters int a file"""

    from ConfigParser import SafeConfigParser

    filename = os.path.join(output_dir, 'parameters.fit')

    print("  * {}".format(filename))

    par_dict = {}

    for name, param in params.items():

        if not param.vary:
            val_print = '{: .5e} fixed'.format(param.value)
        elif param.stderr is not None:
            val_print = '{: .5e} +/- {:.5e}'.format(param.value, param.stderr)
        else:
            val_print = '{: .5e}  ; Error not calculated'.format(param.value)

        par_name = ParameterName.from_full_name(name)

        if par_name.nuclei is None:

            # This is a non-residue-specific parameter

            name_print = par_name.to_section_name()
            section = 'GLOBAL'

        else:

            # This is a residue-specific parameter

            name_print = peaks.Peak(par_name.nuclei)
            section = par_name.to_section_name()

        par_dict.setdefault(section, {})[name_print] = val_print

    cfg = SafeConfigParser()
    cfg.optionxform = str

    section_global = par_dict.pop('GLOBAL', None)

    if section_global is not None:
        cfg.add_section('GLOBAL')
        for name, val in sorted(section_global.items()):
            cfg.set('GLOBAL', name, val)

    for section, name_vals in sorted(par_dict.items()):
        cfg.add_section(section)
        for peak, val in sorted(name_vals.items()):
            cfg.set(section, peak.assignment.upper(), val)

    with open(filename, 'w') as f:
        cfg.write(f)


def remove_comments(line, sep):
    for s in sep:
        line = line.split(s)[0]
    return line.strip()


def main():
    pass


if __name__ == '__main__':
    main()

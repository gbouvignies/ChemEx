from __future__ import print_function

import os
import os.path

import lmfit as lf
import scipy as sp

from chemex import util
from chemex.experiments import sputil


def create_params(data):
    """Creates the array of parameters that will be used for the fitting
    along with the dictionary that associate the name and the index of each
    parameter in the array.
    """

    params = lf.Parameters()

    for profile in data:
        params.update(profile.make_default_parameters())

    return params


def set_params_from_config_file(params, config_filename):
    """Read the file containing the initial guess for the fitting parameters.
    """

    print("\n[{:s}]".format(config_filename))

    config = util.read_cfg_file(config_filename)

    for section in config.sections():

        if section in ['default', 'global']:
            prefix = None
        else:
            prefix = section

        for key, value in config.items(section):

            if 'file' in key:

                filenames = value.split()

                for filename in filenames:
                    filename = util.normalize_path(
                        os.path.dirname(config_filename),
                        filename
                    )

                    set_param_values_from_file(
                        filename, params, prefix=prefix
                    )

            else:

                # Just take the first value in case the file contain a value and
                # an uncertainty

                value = value.split()[0]

                if prefix is not None:
                    key = ','.join([prefix, key])

                set_params(params, key, value)


def set_param_values_from_file(filename, parameters, prefix=None):
    """Reads a file containing values associated with a nucleus name.
    The file should be formatted like sparky peak lists.

    For example:
      * To set G23N to 105.0 and G23H to 8.0:
          G23N-H  105.0  8.0
      * To set a parameter depending on multiple nuclei, eg. G23N and G23H, to 1.0:
          G23N-H  -93.0
    """

    if prefix is None:
        prefix = ''

    try:
        data = sp.genfromtxt(filename, dtype=None)
    except IOError:
        print('The file \'{}\' is empty or does not exist!'
              .format(filename))

    # Hack to solve the problem of 0d-array when 'filename' is a single line
    # file
    if data.ndim < 1:
        data = data.reshape(1, -1)

    for line in data:

        # Don't take into account sparky peak file header
        if line[0] in ['Assignment']:
            continue

        peak = sputil.Peak(line[0])

        if len(peak.resonances) == len(line) - 1:

            for index, resonance in enumerate(peak.resonances, 1):
                subname = resonance.name
                full_name = ','.join([prefix, subname])
                set_params(parameters, full_name, value=line[index])

        elif len(line) == 2:
            subname = sputil.format_assignment(peak.resonances)
            full_name = ','.join([prefix, subname])
            set_params(parameters, full_name, value=line[1])

        else:

            error_message = '\n'.join(
                ["Problem reading the parameter file '{}'.".format(filename),
                 "The number of columns does not match to the number of nuclei"
                 " contained in the assignment:"])
            exit(error_message)

    return parameters


def set_param_status(params, items):
    """Fix (or not) fit variables according to what set in the protocol file"""

    vary = {'fix': False, 'fit': True}

    for key, status in items:
        set_params(params, key, vary=vary[status])


def set_params(parameters, key, value=None, vary=None):
    key_split = set([_.lower() for _ in key.split(',')])

    for name in parameters:

        name_split = set([_.lower().replace('_p_', '.') for _ in name.split('__')])

        if key_split.issubset(name_split):

            if value is not None:
                parameters[name].set(value=float(value))

            if vary is not None:
                parameters[name].set(vary=vary)


def write_par(params, output_dir='./'):
    """Write fitted parameters int a file"""

    from ConfigParser import SafeConfigParser, DuplicateSectionError

    filename = os.path.join(output_dir, 'parameters.fit')

    print("  * {}".format(filename))

    par_name_global = set(['KEX', 'KEX_AB', 'KEX_BC', 'KEX_AC', 'PB', 'PC'])

    par_dict = {}

    for name in params:

        val = params[name].value
        err = params[name].stderr

        if not params[name].vary:
            par_dict[name] = '{: .5e} fixed'.format(val)
        elif err is not None:
            par_dict[name] = '{: .5e} +/- {:.5e}'.format(val, err)
        else:
            par_dict[name] = '{: .5e}   ; Error not calculated'.format(val)

    cfg = SafeConfigParser()
    cfg.optionxform = str

    for name, val in sorted(par_dict.items()):

        name_list = name.replace('_p_', '.').split('__')

        if name_list[0].upper() in par_name_global:
            name_str = ', '.join([str(_).upper() for _ in name_list])
            section = 'global'

        else:
            name_str = str(name_list.pop(1)).upper()
            section = ', '.join([str(_).upper() for _ in name_list])

        try:
            cfg.add_section(section)

        except DuplicateSectionError:
            pass

        cfg.set(section, name_str, val)

    with open(filename, 'w') as f:
        cfg.write(f)


def main():
    pass


if __name__ == '__main__':
    main()

import os.path

import lmfit as lf
import scipy as sp

from chemex import parsing
from chemex import utils


def create_params(data):
    """
    Creates the array of parameters that will be used for the fitting
    along with the dictionary that associate the name and the index of each
    parameter in the array.
    """

    params = lf.Parameters()

    for profile in data:
        params.update(profile.make_default_parameters())

    return params


def set_params_from_config_file(params, config_filename):
    """
    Read the file containing the initial guess for the fitting parameters.
    """

    print("\n[{:s}]".format(config_filename))

    # Parse the config file
    config = utils.read_cfg_file(config_filename)

    for section in config.sections():

        if section in ['default', 'global']:
            prefix = ''
        else:
            prefix = section

        for key, value in config.items(section):

            if 'file' in key:

                filenames = value.split()

                for filename in filenames:
                    filename = utils.normalize_path(
                        os.path.dirname(config_filename),
                        filename
                    )

                    set_param_values_from_file(filename, params, prefix=prefix)

            else:

                # Just take the first value in case the file contain a value and
                # an uncertainty

                value = value.split()[0]

                full_key = ','.join([prefix, key])
                set_params(params, full_key, value)

    return params


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
    if data.ndim <= 2:
        data = data.reshape(1, -1)

    for line in data:

        # Don't take into account sparky peak file header
        if line[0] in ['Assignment']:
            continue

        assignments = parsing.parse_assignment(line[0])

        if len(assignments) == len(line) - 1:

            for index, assignment in enumerate(assignments, 1):
                subname = parsing.assignment_name([assignment])
                full_name = ','.join([prefix, subname])
                set_params(parameters, full_name, value=line[index])

        elif len(line) == 2:
            subname = parsing.assignment_name(assignments)
            full_name = ','.join([prefix, subname])
            parameters.append((full_name, line[1]))

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
    key_split = set([_.strip().lower() for _ in key.split(',')])

    for name in parameters:

        name_split = set([_.strip().lower() for _ in name.split('__')])

        if key_split.issubset(name_split):

            if value is not None:
                parameters[name].set(value=float(value))

            if vary is not None:
                parameters[name].set(vary=vary)


def main():
    pass


if __name__ == '__main__':
    main()

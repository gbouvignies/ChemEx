""" Reads the "experiment" files.
"""

import os
import os.path
import sys
import pkgutil
import ConfigParser


def read_cfg_file(input_file, res_incl=None, res_excl=None):
    """Reads the "experiment" file containing the experimental parameters and the location of the data files"""

    # Get the directory of the input file
    working_dir = os.path.dirname(input_file)

    # Parse the config file
    cfg = ConfigParser.ConfigParser()

    if os.path.isfile(input_file):
        try:
            cfg.read(input_file)
        except ConfigParser.MissingSectionHeaderError:
            exit('You are missing a section heading (default?) in {:s}\n'.format(input_file))
        except ConfigParser.ParsingError:
            exit('Having trouble reading your parameter file, have you forgotten \'=\' signs?\n{:s}'
                 .format(sys.exc_info()[1]))
    else:
        exit("The file \'{}\' is empty or does not exist!\n".format(input_file))

    try:
        # Reads experimental and estimated parameters shared by all residues
        experimental_parameters = read_experimental_paramaters(cfg)

        # Reads experimental measurements
        data = read_data(cfg, working_dir, experimental_parameters, res_incl, res_excl)

    except ConfigParser.NoSectionError:
        exit("\nIn {:s}, {:s}!\n".format(input_file, sys.exc_info()[1]))

    except (KeyboardInterrupt):
        exit("\n -- ChemEx killed while reading experiment and data files\n")

    return data


def read_experimental_paramaters(cfg):
    """Reads experimental and estimated parameters shared by all residues"""

    # This for legacy ("experimental_parameters" was previously called "global_parameters")
    if cfg.has_section('experimental_parameters'):
        experimental_paramaters = dict(cfg.items('experimental_parameters'))
    elif cfg.has_section('global_parameters'):
        experimental_paramaters = dict(cfg.items('global_parameters'))
    else:
        experimental_paramaters = dict()

    experimental_paramaters['experiment_type'] = cfg.get('experiment', 'type')

    if cfg.has_option('experiment', 'name'):
        experimental_paramaters['experiment_name'] = cfg.get('experiment', 'name').lower()

    return experimental_paramaters


def read_data(cfg, working_dir, global_parameters, res_incl=None, res_excl=None):
    """Reads experimental measurements"""

    exp_type = global_parameters['experiment_type']

    path = os.path.dirname(__file__)
    pkgs = [
        modname
        for _, modname, ispkg in pkgutil.iter_modules([path])
        if ispkg and modname in exp_type
    ]

    try:
        pkg = max(pkgs)
    except:
        exit("\nUnknown data type {:s}"
             "\nDid you forget _cpmg, _cest, etc?"
             "\n".format(global_parameters['experiment_type']))

    reading = __import__(pkg + '.reading', globals(), locals(), ['read_data'], -1)

    data = reading.read_data(cfg, working_dir, global_parameters, res_incl, res_excl)

    return data

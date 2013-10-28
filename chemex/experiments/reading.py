"""
Created on Mar 30, 2011

@author: guillaume
"""

# Standard libraries
import os
import sys
import pkgutil
import ConfigParser


def read_cfg_file(input_file, res_incl=None, res_excl=None):
    """Read the r2 data file containing the experimental parameters and the location of the fuda files"""

    # Get the directory of the input file
    working_dir = os.path.dirname(input_file)

    # Parse the config file
    cfg = ConfigParser.ConfigParser()
    cfg.optionxform = str
    cfg.read(input_file)

    try:
        # Reads experimental and estimated parameters shared by all residues
        global_parameters = read_global_paramaters(cfg)

        # Reads experimental measurements
        data = read_data(cfg, working_dir, global_parameters, res_incl, res_excl)

    except ConfigParser.NoSectionError:
        exit("\nIn {:s}, {:s}!\n".format(input_file, sys.exc_info()[1]))
    except (KeyboardInterrupt):
        exit("\n -- ChemEx killed while reading experiment and data files\n")

    return data


def read_global_paramaters(cfg):
    """Reads experimental and estimated parameters shared by all residues"""

    global_paramaters = dict(cfg.items('global_parameters'))
    global_paramaters['experiment_type'] = cfg.get('experiment', 'type')

    if cfg.has_option('experiment', 'name'):
        global_paramaters['experiment_name'] = cfg.get('experiment', 'name')

    return global_paramaters


def read_data(cfg, working_dir, global_parameters, res_incl=None, res_excl=None):
    """Reads experimental measurements"""

    exp_type = global_parameters['experiment_type']

    path = os.path.dirname(__file__)
    pkgs = [modname for _, modname, ispkg in pkgutil.iter_modules([path])
            if ispkg and modname in exp_type]

    try:
        pkg = max(pkgs)
    except:
        exit("\nUnknown data type {:s}"
             "\nDid you forget _cpmg, _cest, etc?"
             "\n".format(global_parameters['experiment_type']))

    reading = __import__(pkg + '.reading', globals(), locals(), ['read_data'], -1)

    data = reading.read_data(cfg, working_dir, global_parameters, res_incl, res_excl)

    return data

"""The util module contains a variety of utility functions."""

import configparser
import os
import sys


def make_dir(path=None):
    """Ensure existence of the directory."""
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def read_cfg_file(filename):
    """Read and parse the experiment configuration file with configparser."""
    config = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    config.optionxform = str

    if filename:
        try:
            out = config.read(filename)

            if not out:
                exit("\nERROR: The file \'{}\' is empty or does not exist!\n".format(filename))

        except configparser.MissingSectionHeaderError:
            exit("\nERROR: You are missing a section heading in {:s}\n".format(filename))

        except configparser.ParsingError:
            exit("\nERROR: Having trouble reading your parameter file, did you"
                 " forget '=' signs?\n{:s}".format(sys.exc_info()[1]))

    return config


def normalize_path(working_dir, filename):
    """Normalize the path of a filename relative to a specific directory."""
    if not os.path.isabs(filename):
        path = os.path.join(working_dir, filename)
    else:
        path = filename

    return path


def include_selection(data, selection):
    """Make a new dataset including points whose 'id' is in the selection."""
    new_data = [
        a_data_point for a_data_point in data
        if a_data_point.par.get('resonance_id', None) in selection
    ]

    return new_data


def exclude_selection(data, selection):
    """Make a new dataset excluding points whose 'id' is in the selection."""
    new_data = [
        a_data_point for a_data_point in data
        if a_data_point.par.get('resonance_id', None) not in selection
    ]

    if new_data == data:
        sys.stdout.write("\n No data removed! Aborting ...\n")
        exit(1)

    return new_data


def header1(string):
    """Print a formatted heading."""
    print(("\n".join(["", "", string, "=" * len(string), ""])))


def header2(string):
    """Print a formatted subheading."""
    print(("\n".join(["", string, "-" * len(string), ""])))

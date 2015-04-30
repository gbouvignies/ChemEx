import ConfigParser
import os
import sys
import math


def make_dir(path=None):
    """Ensure existence of the directory"""

    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def read_cfg_file(filename):
    """Parse config files with ConfigParser"""

    # Parse the config file
    config = ConfigParser.SafeConfigParser()

    try:
        config.read(filename)

    except ConfigParser.MissingSectionHeaderError:
        exit("You are missing a section heading in {:s}\n"
             .format(filename))

    except ConfigParser.ParsingError:
        exit("Having trouble reading your parameter file, have you "
             "forgotten '=' signs?\n{:s}".format(sys.exc_info()[1]))

    if not config:
        exit("The file \'{}\' is empty or does not exist!\n".format(filename))

    return config


class AutoDict(dict):
    """Implementation of perl's autovivification feature."""

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def normalize_path(working_dir, filename):
    """Normalizes the path of a file name relative to a specific directory."""

    if not os.path.isabs(filename):
        path = os.path.join(working_dir, filename)

    return path


def include_selection(data, selection):
    """Makes a new dataset including points whose 'id' is in selection."""

    new_data = [
        a_data_point for a_data_point in data
        if a_data_point.par.get('resonance_id', None) in selection
    ]

    return new_data


def exclude_selection(data, selection):
    """Makes a new dataset excluding points whose id is in selection."""

    new_data = [
        a_data_point for a_data_point in data
        if a_data_point.par.get('resonance_id', None) not in selection
    ]

    if new_data == data:
        sys.stdout.write("\n No Data removed! Aborting ...\n")
        exit(1)

    return new_data


def get_digit_number(integer):
    return int(math.log10(integer)) + 1


def header1(string):
    print("\n".join(["", "", string, "=" * len(string)]))


def header2(string):
    print("\n".join(["", string, "-" * len(string)]))

"""
Created on 2012-02-21

@author: guillaume
"""

import sys
import os


class autodict(dict):
    """
    Implementation of perl's autovivification feature.
    """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def normalize_path(working_dir, filename):
    """
    Checks whether the filename is a relative path and adjust it,
    so that it can be read from the active directory.
    """

    if not os.path.isabs(filename):
        filename = os.path.join(working_dir, filename)

    return filename


def include_selection(data, selection):
    """
    Makes a new dataset including points whose resonance_id is in selection.
    """

    new_data = list()

    for a_data_point in data:

        if ('resonance_id' in a_data_point.par and a_data_point.par[
            'resonance_id'] in selection):
            new_data.append(a_data_point)

    return new_data


def exclude_selection(data, selection):
    """
    Makes a new dataset excluding points whose resonance_id is in selection.
    """

    new_data = list()

    for a_data_point in data:

        if ('resonance_id' in a_data_point.par and a_data_point.par[
            'resonance_id'] not in selection):
            new_data.append(a_data_point)

    if new_data == data:
        sys.stdout.write("\n No Data removed! Aborting ...\n")
        exit(1)

    return new_data


def make_dir(path=None):
    """Make the directory if needed"""

    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError:
            exit("\nOSError: You can not use that directory!\n")


def header1(string):
    print("\n".join(["",
                     "",
                     string,
                     "=" * len(string)]))


def header2(string):
    print("\n".join(["",
                     string,
                     "-" * len(string)]))

"""
Created on 2012-02-21

@author: guillaume
"""

# Standar Libraies
import sys
import os
import re


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


def parse_assignment(assignment):
    """
    Parse assignment of form g1a1-g2a2 to get ((g1, a1), (g2, a2))
    Or g1a1-a2 to get ((g1, a1), (g1, a2))
    A '?' component is translated to ('', '')
    """

    res = assignment.lower().split('-')
    assignment = []
    last_group = None
    for s in res:
        ga = split_group_atom(s)
        if ga == None:
            if last_group:
                ga = (last_group, s)
            else:
                return None
        assignment.append((parse_group_name(ga[0]) + ga[1:]))
        last_group = ga[0]

    return tuple(assignment)


def split_group_atom(group_atom):
    """
    Skip to first digit, then skip to first H|C|N to find start of atom name.
    """

    s = re.search('[0-9]', group_atom)
    if s:
        first_digit = s.start()
        s = re.search('[hHcCnNqQmM]', group_atom[first_digit:])
        if s:
            HCNQM_offset = s.start()
            d = first_digit + HCNQM_offset
            return (group_atom[:d], group_atom[d:])
    if group_atom == '?':
        return ('', '')
    return None


def parse_group_name(g):
    s = re.search('[0-9]+', g)
    if s:
        return int(s.group()), g[:s.start()]
    return g, None


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

import contextlib
import os
import sys


@contextlib.contextmanager
def suppress(*exceptions):
    try:
        yield
    except exceptions:
        pass


def make_dir(path=None):
    """Make the directory if needed"""

    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError:
            exit("\nOSError: You can not use that directory!\n")


class autodict(dict):
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
        filename = os.path.join(working_dir, filename)

    return filename


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


def header1(string):
    print("\n".join(["", "", string, "=" * len(string)]))


def header2(string):
    print("\n".join(["", string, "-" * len(string)]))

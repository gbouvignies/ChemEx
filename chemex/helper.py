"""The util module contains a variety of utility functions."""
import re
import sys

import toml


RE_GROUPNAME = re.compile(r"^[A-Za-z0-9_-]+$")


def read_toml(filename=None):
    """Read and parse the experiment configuration file with configparser."""

    try:
        config = toml.load(filename)

    except FileNotFoundError:
        sys.exit(f"\nERROR: The file '{filename}' is empty or does not exist!\n")

    except (toml.TomlDecodeError, TypeError) as e:
        sys.exit(f"\nERROR: '{filename}': {e}\n")

    else:
        return config


def listfloat(text):
    return [float(val) for val in text.replace("\n", "").strip("[]").split(",")]


def liststr(text):
    return [val for val in text.replace("\n", "").strip("[],").split(",")]


def normalize_path(working_dir, filename):
    """Normalize the path of a filename relative to a specific directory."""

    path = filename

    if not path.is_absolute():
        path = working_dir / path

    return path.resolve()


def header1(string_):
    """Print a formatted heading."""
    if string_:
        print(("\n".join(["", "", string_, "=" * len(string_)])))


def header2(string_):
    """Print a formatted subheading."""
    if string_:
        print(("\n".join(["", string_, "-" * len(string_)])))

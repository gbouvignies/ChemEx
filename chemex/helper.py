"""The util module contains a variety of utility functions."""
import re

import toml


RE_GROUPNAME = re.compile(r"^[A-Za-z0-9_-]+$")


def read_toml(filename=None):
    """Read and parse the experiment configuration file with configparser."""

    try:
        config = toml.load(filename)

    except FileNotFoundError:
        exit(f"\nERROR: The file '{filename}' is empty or does not exist!\n")

    except (toml.TomlDecodeError, TypeError) as e:
        exit(f"\nERROR: '{filename}': {e}\n")

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


def header1(strings):
    """Print a formatted heading."""
    print(("\n".join(["", "", strings, "=" * len(strings), ""])))


def header2(strings):
    """Print a formatted subheading."""
    print(("\n".join(["", strings, "-" * len(strings), ""])))

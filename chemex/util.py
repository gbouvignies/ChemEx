"""The util module contains a variety of utility functions."""
import configparser
import sys


def read_cfg_file(filename):
    """Read and parse the experiment configuration file with configparser."""

    config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    config.optionxform = str

    try:
        out = config.read(str(filename))

        if not out and filename is not None:
            exit(f"\nERROR: The file '{filename}' is empty or does not exist!\n")

    except configparser.MissingSectionHeaderError:
        exit(f"\nERROR: You are missing a section heading in {filename:s}\n")

    except configparser.ParsingError:
        exit(
            "\nERROR: Having trouble reading your parameter file, did you"
            " forget '=' signs?\n{:s}".format(sys.exc_info()[1])
        )

    return config


def normalize_path(working_dir, filename):
    """Normalize the path of a filename relative to a specific directory."""

    path = filename

    if not path.is_absolute():
        path = working_dir / path

    return path.resolve()


def header1(string):
    """Print a formatted heading."""
    print(("\n".join(["", "", string, "=" * len(string), ""])))


def header2(string):
    """Print a formatted subheading."""
    print(("\n".join(["", string, "-" * len(string), ""])))

"""The parsing module contains the code for the parsing of command-line
arguments."""

import argparse
import importlib
import pkgutil
import sys

from chemex import __version__, experiments, util

FITMETHODS = {
    "cobyla",
    "slsqp",
    "tnc",
    "nelder",
    "cg",
    "bfgs",
    "powell",
    "l-bfgsb",
    "least_squares",
    "lbfgsb",
    "differential_evolution",
    "leastsq",
    "brute",
}


class MyParser(argparse.ArgumentParser):
    """Subclass of ArgumentParser to override the error method."""

    def error(self, message):
        sys.stderr.write("error: {}\n\n".format(message))
        self.print_help()
        sys.exit(2)


def build_parser(func_fit):
    """Parse the command-line arguments."""
    description = (
        "ChemEx is an analysis program for chemical exchange detected by "
        "NMR. It is designed to take almost any kind of NMR data to aid the "
        "analysis, but the principle techniques are CPMG relaxation "
        "dispersion and Chemical Exchange Saturation Transfer."
    )

    parser = MyParser(description=description, prog="chemex")

    parser.add_argument(
        "--version", action="version", version="{} {}".format(parser.prog, __version__)
    )

    commands = parser.add_subparsers(dest="commands")

    # parser for the positional argument "info"
    info_parser = commands.add_parser(
        "info",
        help="Shows experiments that can be fit",
        description="Enter the name of an experiment to obtain more info about it.",
    )

    info_parser.set_defaults(func=get_info)

    experiments_parser = info_parser.add_subparsers(dest="experiments")
    experiments_parser.required = True

    experiments_names = list_experiments()

    for name in experiments_names:
        experiment_module = importlib.import_module(name)
        experiment_shortname = compress_name(name)

        experiments_parser.add_parser(
            experiment_shortname,
            help=get_description_from_module(experiment_module),
            add_help=False,
        )

    # parser for the positional argument "fit"
    fit = commands.add_parser("fit", help="Starts a fit", prefix_chars="+-")

    fit.set_defaults(func=func_fit)

    fit.add_argument(
        "-e",
        dest="experiments",
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input files containing experimental setup and data location",
    )

    fit.add_argument(
        "-d",
        dest="model",
        metavar="MODEL",
        default="2st.pb_kex",
        help="Exchange model used to fit the data",
    )

    fit.add_argument(
        "-p",
        dest="parameters",
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input file containing the initial values of fitting parameters",
    )

    fit.add_argument(
        "-m",
        dest="method",
        metavar="FILE",
        help="Input file containing the fitting method",
    )

    fit.add_argument(
        "-o",
        dest="out_dir",
        metavar="DIR",
        default="./Output",
        help="Directory for output files",
    )

    fit.add_argument("--noplot", action="store_true", help="No plots of the fits")

    fit.add_argument(
        "-f",
        dest="fitmethod",
        metavar="FITMETHOD",
        default="leastsq",
        choices=sorted(FITMETHODS),
        help="Specify the fitting method",
    )

    selection = fit.add_mutually_exclusive_group()
    selection.add_argument(
        "+r",
        dest="res_incl",
        metavar="ID",
        nargs="+",
        help="residue(s) to include in the fit",
    )
    selection.add_argument(
        "-r",
        dest="res_excl",
        metavar="ID",
        nargs="+",
        help="residue(s) to exclude from the fit",
    )

    simulation = fit.add_mutually_exclusive_group()
    simulation.add_argument(
        "--mc", metavar="N", type=int, help="Run N Monte-Carlo simulations"
    )
    simulation.add_argument(
        "--bs", metavar="N", type=int, help="Run N Bootstrap simulations"
    )

    return parser


def get_description_from_module(module):
    return module.__doc__.strip().splitlines()[0]


def list_experiments():
    path = experiments.__path__
    prefix = experiments.__name__ + "."
    names = []
    types_modules = pkgutil.walk_packages(path, prefix)
    for __, name, ispkg in types_modules:
        module = importlib.import_module(name)
        if not ispkg and hasattr(module, "Profile"):
            names.append(name)
    return names


def get_info(args):
    name = expand_name(args.experiments)
    module = importlib.import_module(name)
    util.header1('Description of the "' "{}" '" experiment'.format(args.experiments))
    print(module.__doc__)


def compress_name(name):
    tree = name.split(".")
    shortname = ".".join((tree[-3], tree[-1]))
    return shortname


def expand_name(name):
    type1, type2 = name.split(".")
    fullname = ".".join(["chemex.experiments", type1, "profiles", type2])
    return fullname

"""The parsing module contains the code for the parsing of command-line arguments."""
import sys
from argparse import ArgumentParser
from argparse import Namespace
from pathlib import Path

from chemex import __version__
from chemex import chemex
from chemex import tools
from chemex.experiments.configurations import configurations
from chemex.experiments.descriptions import descriptions
from chemex.messages import console
from chemex.parameters.spin_system import SpinSystem


class MyParser(ArgumentParser):
    """Subclass of ArgumentParser to override the error method."""

    def error(self, message: str):
        console.print(f"[red]\nError: {message}\n")
        self.print_help()
        sys.exit(1)


def print_info(args: Namespace):
    """Print information for experiment."""
    descriptions.print(args.experiments)


def write_config(args: Namespace):
    """Write default configuration file for experiment."""
    configurations.write(args.experiments, args.output)
    print(
        f"\nSample configuration file for '{args.experiments}' written to "
        f"'{args.output}'.\n"
    )


def build_parser():
    """Parse the command-line arguments."""
    description = (
        "ChemEx is an analysis program for chemical exchange detected by "
        "NMR. It is designed to take almost any kind of NMR data to aid the "
        "analysis, but the principle techniques are CPMG relaxation "
        "dispersion and Chemical Exchange Saturation Transfer."
    )

    parser = MyParser(description=description, prog="chemex")

    parser.add_argument(
        "--version", action="version", version=f"{parser.prog} {__version__}"
    )

    commands = parser.add_subparsers(dest="commands")

    # parser for the positional argument "info"
    info_parser = commands.add_parser(
        "info",
        help="Show experiments that can be fit/simulated",
        description="Enter the name of an experiment to obtain more info about it.",
    )

    info_parser.set_defaults(func=print_info)

    experiments_parser = info_parser.add_subparsers(
        dest="experiments", metavar="experiment_type", required=True
    )

    for exp_name in descriptions.experiment_names():
        experiments_parser.add_parser(
            exp_name, help=descriptions.first_line(exp_name), add_help=False
        )

    # parser for the positional argument "info"
    config_parser = commands.add_parser(
        "config",
        help="Write sample configuration file for an experiment",
        description="Enter the name of an experiment to obtain the associated sample "
        "configuration file.",
    )

    config_parser.set_defaults(func=write_config)

    config_exp_parser = config_parser.add_subparsers(dest="experiments", required=True)

    for exp_name in descriptions.experiment_names():
        parser_ = config_exp_parser.add_parser(exp_name)
        parser_.add_argument(
            "-o",
            dest="output",
            type=Path,
            metavar="FILE",
            default="./experiment.toml",
            help="Name of the output file (default: './experiment.toml')",
        )

    # parser for the positional argument "fit"
    fit_parser = commands.add_parser("fit", help="Start a fit")

    fit_parser.set_defaults(func=chemex.run)

    fit_parser.add_argument(
        "-e",
        dest="experiments",
        type=Path,
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input file(s) containing experimental setup and data location",
    )

    fit_parser.add_argument(
        "-p",
        dest="parameters",
        type=Path,
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input file(s) containing the initial values of fitting parameters",
    )

    fit_parser.add_argument(
        "-m",
        dest="method",
        type=Path,
        metavar="FILE",
        nargs="+",
        help="Input file(s) containing the fitting method",
    )

    fit_parser.add_argument(
        "-o",
        dest="out_dir",
        type=Path,
        metavar="DIR",
        default="./Output",
        help="Directory for output files (default: './Output')",
    )

    fit_parser.add_argument(
        "-d",
        dest="model",
        metavar="MODEL",
        default="2st",
        help="Exchange model used to fit the data (default: '2st')",
    )

    fit_parser.add_argument(
        "--plot",
        type=str,
        choices=["nothing", "normal", "all"],
        default="normal",
        help="Plotting level (default: 'normal')",
    )

    fit_parser.add_argument(
        "--include",
        dest="include",
        metavar="ID",
        nargs="+",
        help="Residue(s) to include in the fit",
        type=SpinSystem,
    )

    fit_parser.add_argument(
        "--exclude",
        dest="exclude",
        metavar="ID",
        nargs="+",
        help="Residue(s) to exclude from the fit",
        type=SpinSystem,
    )

    # parser for the positional argument "simulate"
    simulate_parser = commands.add_parser("simulate", help="Start a simulation")

    simulate_parser.set_defaults(func=chemex.run)

    simulate_parser.add_argument(
        "-e",
        dest="experiments",
        type=Path,
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input file(s) containing experimental setup and data location",
    )

    simulate_parser.add_argument(
        "-p",
        dest="parameters",
        type=Path,
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input file(s) containing the values of parameters",
    )

    simulate_parser.add_argument(
        "-o",
        dest="out_dir",
        type=Path,
        metavar="DIR",
        default="./OutputSim",
        help="Directory for output files (default: './OutputSim')",
    )

    simulate_parser.add_argument(
        "-d",
        dest="model",
        metavar="MODEL",
        default="2st",
        help="Exchange model used to simulate the data (default: '2st')",
    )

    simulate_parser.add_argument(
        "--plot",
        type=str,
        choices=["nothing", "normal"],
        default="normal",
        help="Plotting level (default: 'normal')",
    )

    simulate_parser.add_argument(
        "--include",
        dest="include",
        metavar="ID",
        nargs="+",
        help="Residue(s) to include in the simulation",
        type=SpinSystem,
    )

    simulate_parser.add_argument(
        "--exclude",
        dest="exclude",
        metavar="ID",
        nargs="+",
        help="Residue(s) to exclude from the simulation",
        type=SpinSystem,
    )

    # parser for the positional argument "pick_cest"
    pick_cest_parser = commands.add_parser(
        "pick_cest", help="Plot CEST profiles for dip picking"
    )

    pick_cest_parser.set_defaults(func=tools.pick_cest.pick_cest)

    pick_cest_parser.add_argument(
        "-e",
        dest="experiments",
        type=Path,
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input file(s) containing experimental setup and data location",
    )

    pick_cest_parser.add_argument(
        "-o",
        dest="out_dir",
        type=Path,
        metavar="DIR",
        default="./Sandbox",
        help="Directory for output files (default: './Sandbox')",
    )

    # parser for the positional argument "pick_cest"
    plot_param_parser = commands.add_parser(
        "plot_param", help="Plot one selected parameter from a 'parameters.fit' file"
    )

    plot_param_parser.set_defaults(func=tools.plot_param.plot_param)

    plot_param_parser.add_argument(
        "-p",
        dest="parameters",
        type=Path,
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input file containing the fitted parameters to be plotted",
    )

    plot_param_parser.add_argument(
        "-n",
        dest="parname",
        metavar="NAME",
        required=True,
        help="Name of the parameter to plot",
    )

    return parser

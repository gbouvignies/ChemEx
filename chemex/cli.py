"""The parsing module contains the code for the parsing of command-line
arguments."""
import argparse
import pathlib
import sys

import chemex.chemex as cc
import chemex.experiments as ce
import chemex.helper as ch
import chemex.tools as ct


class MyParser(argparse.ArgumentParser):
    """Subclass of ArgumentParser to override the error method."""

    def error(self, message):
        sys.stderr.write(f"error: {message}\n\n")
        self.print_help()
        sys.exit(2)


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
        "--version", action="version", version=f"{parser.prog} {'__version__'}"
    )

    commands = parser.add_subparsers(dest="commands")

    # parser for the positional argument "info"
    info_parser = commands.add_parser(
        "info",
        help="Shows experiments that can be fit",
        description="Enter the name of an experiment to obtain more info about it.",
    )

    info_parser.set_defaults(func=print_info)

    experiments_parser = info_parser.add_subparsers(dest="experiments")
    experiments_parser.required = True

    docs = ce.get_info()

    for exp_name, doc in docs.items():
        experiments_parser.add_parser(
            exp_name, help=get_description_from_doc(doc), add_help=False
        )

    # parser for the positional argument "info"
    config_parser = commands.add_parser(
        "config",
        help="Shows sample configuration files of the modules",
        description="Enter the name of an experiment to obtain the associated sample "
        "configuration file.",
    )

    config_parser.set_defaults(func=write_config)

    config_exp_parser = config_parser.add_subparsers(dest="experiments")
    config_exp_parser.required = True

    infos = ce.get_config()

    for exp_name in infos:
        parser_ = config_exp_parser.add_parser(exp_name)
        parser_.add_argument(
            "-o",
            dest="output",
            type=pathlib.Path,
            metavar="FILE",
            default="./experiment.toml",
            help="Name of the output file",
        )

    # parser for the positional argument "fit"
    fit_parser = commands.add_parser("fit", help="Starts a fit", prefix_chars="+-")

    fit_parser.set_defaults(func=cc.fit)

    fit_parser.add_argument(
        "-e",
        dest="experiments",
        type=pathlib.Path,
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input files containing experimental setup and data location",
    )

    fit_parser.add_argument(
        "-d",
        dest="model",
        metavar="MODEL",
        default="2st.pb_kex",
        help="Exchange model used to fit the data",
    )

    fit_parser.add_argument(
        "-p",
        dest="parameters",
        type=pathlib.Path,
        metavar="FILE",
        nargs="+",
        required=True,
        help="Input file containing the initial values of fitting parameters",
    )

    fit_parser.add_argument(
        "-m",
        dest="method",
        type=pathlib.Path,
        metavar="FILE",
        help="Input file containing the fitting method",
    )

    fit_parser.add_argument(
        "-o",
        dest="out_dir",
        type=pathlib.Path,
        metavar="DIR",
        default="./Output",
        help="Directory for output files",
    )

    fit_parser.add_argument(
        "--noplot", action="store_true", help="No plots of the fits"
    )

    fit_parser.add_argument(
        "-f",
        dest="fitmethod",
        metavar="FITMETHOD",
        default="leastsq",
        help="Specify the fitting method",
    )

    selection = fit_parser.add_mutually_exclusive_group()
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

    simulation = fit_parser.add_mutually_exclusive_group()
    simulation.add_argument(
        "--mc", metavar="N", type=int, help="Run N Monte-Carlo simulations"
    )
    simulation.add_argument(
        "--bs", metavar="N", type=int, help="Run N Bootstrap simulations"
    )

    # parser for the positional argument "pick_cest"
    pick_cest_parser = commands.add_parser(
        "pick_cest", help="Plot CEST profiles for dip picking"
    )

    pick_cest_parser.set_defaults(func=ct.pick_cest.pick_cest)

    pick_cest_parser.add_argument(
        "-e",
        dest="experiments",
        type=pathlib.Path,
        metavar="FILE",
        nargs=1,
        required=True,
        help="Input files containing experimental setup and data location",
    )

    pick_cest_parser.add_argument(
        "-o",
        dest="out_dir",
        type=pathlib.Path,
        metavar="DIR",
        default="./Output",
        help="Directory for output files",
    )

    # parser for the positional argument "pick_cest"
    plot_param_parser = commands.add_parser(
        "plot_param", help="Plot one selected parameter from a 'parameters.fit' file"
    )

    plot_param_parser.set_defaults(func=ct.plot_param.plot_param)

    plot_param_parser.add_argument(
        "-p",
        dest="parameters",
        type=pathlib.Path,
        metavar="FILE",
        required=True,
        help="Input file containing the initial values of fitting parameters",
    )

    plot_param_parser.add_argument(
        "-n",
        dest="parname",
        metavar="NAME",
        required=True,
        help="Name of the parameter to plot",
    )

    return parser


def get_description_from_doc(doc):
    return doc.strip().splitlines()[0]


def print_info(args):
    docs = ce.get_info()
    ch.header1(f'Description of the "{args.experiments}" experiment')
    print(docs[args.experiments])


def write_config(args):
    config = ce.get_config()
    args.output.write_text(config[args.experiments])
    print(
        f"\nSample configuration file for '{args.experiments}' written to "
        f"'{args.output}'.\n"
    )

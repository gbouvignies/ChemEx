"""The chemex module provides the entry point for the chemex script."""
from __future__ import annotations

import sys
from argparse import Namespace

from chemex import model
from chemex.cli import build_parser
from chemex.configuration.methods import Method
from chemex.configuration.methods import read_methods
from chemex.configuration.methods import Selection
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.experiments.builder import build_experiments
from chemex.experiments.loader import register_experiments
from chemex.messages import print_logo
from chemex.messages import print_no_data
from chemex.messages import print_reading_defaults
from chemex.messages import print_reading_methods
from chemex.messages import print_running_simulations
from chemex.messages import print_start_fit
from chemex.optimize.fitting import run_methods
from chemex.optimize.helper import execute_post_fit
from chemex.parameters import database


def run_fit(args: Namespace, experiments: Experiments):
    # Filter datapoints out if necessary (e.g., on-resonance filter CEST)
    experiments.filter()

    if args.method is not None:
        print_reading_methods()
        methods = read_methods(args.method)
    else:
        methods = {"": Method()}

    print_start_fit()
    run_methods(experiments, methods, args.output, args.plot)


def run_sim(args: Namespace, experiments: Experiments):
    print_running_simulations()

    path = args.output
    plot = args.plot == "normal"

    database.fix_all_parameters()

    execute_post_fit(experiments, path, plot=plot)


def run(args: Namespace):
    """Run the fit or simulation."""

    # Parse kinetics model
    model.set_model(args.model)

    # Read experimental setup and data
    selection = Selection(args.include, args.exclude)
    experiments = build_experiments(args.experiments, selection)

    if not experiments:
        print_no_data()
        sys.exit()

    # Read initial values of fitting/fixed parameters
    print_reading_defaults()
    defaults = read_defaults(args.parameters)
    database.set_param_defaults(defaults)

    if args.commands == "simulate":
        run_sim(args, experiments)
    else:
        run_fit(args, experiments)


def main():
    """Do all the magic."""
    print_logo()

    register_experiments()

    parser = build_parser()
    args = parser.parse_args()
    args.func(args)

"""The chemex module provides the entry point for the chemex script."""

import sys
from argparse import Namespace

from chemex.cli import build_parser
from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.containers.experiments import Experiments
from chemex.experiments.builder import build_experiments
from chemex.messages import (
    print_logo,
    print_no_data,
    print_reading_defaults,
    print_reading_methods,
    print_running_simulations,
    print_start_fit,
)
from chemex.optimize.fitting import run_methods
from chemex.optimize.helper import execute_simulation
from chemex.runtime import AnalysisSession, ensure_plugins_registered


def run_fit(
    args: Namespace,
    experiments: Experiments,
    session: AnalysisSession,
) -> None:
    # Filter datapoints out if necessary (e.g., on-resonance filter CEST)
    experiments.filter()

    if args.method is not None:
        print_reading_methods()
        methods = read_methods(args.method)
    else:
        methods = {"": Method()}

    print_start_fit()
    run_methods(experiments, methods, args.output, args.plot, session=session)


def run_sim(
    args: Namespace,
    experiments: Experiments,
    session: AnalysisSession,
) -> None:
    print_running_simulations()

    path = args.output
    plot = args.plot == "normal"

    session.parameters.fix_all_parameters()

    execute_simulation(experiments, path, plot=plot, session=session)


def run(args: Namespace, session: AnalysisSession | None = None) -> None:
    """Run the fit or simulation."""
    if session is None:
        session = AnalysisSession.create()

    # Parse kinetics model
    session.set_model(args.model)

    # Read experimental setup and data
    selection = Selection(args.include, args.exclude)
    experiments = build_experiments(args.experiments, selection, session=session)

    if not experiments:
        print_no_data()
        sys.exit()

    # Read initial values of fitting/fixed parameters
    print_reading_defaults()
    defaults = read_defaults(args.parameters)
    session.parameters.set_param_defaults(defaults)

    if args.commands == "simulate":
        run_sim(args, experiments, session)
    else:
        run_fit(args, experiments, session)


def main() -> None:
    """Do all the magic."""
    print_logo()
    ensure_plugins_registered()

    parser = build_parser()
    args = parser.parse_args()
    args.func(args)

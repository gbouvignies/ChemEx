"""The chemex module provides the entry point for the chemex script."""
import sys

import chemex
import chemex.cli as cc
import chemex.containers.experiment as cce
import chemex.helper as ch
import chemex.optimize.fitting as cf
import chemex.optimize.helper as coh
import chemex.parameters.helper as cph
import chemex.parameters.kinetics as cpk
import chemex.parameters.settings as cps


LOGO = fr"""
* * * * * * * * * * * * * * * * * * * * * * * * *
*      ________                   ______        *
*     / ____/ /_  ___  ____ ___  / ____/  __    *
*    / /   / __ \/ _ \/ __ `__ \/ __/ | |/_/    *
*   / /___/ / / /  __/ / / / / / /____>  <      *
*   \____/_/ /_/\___/_/ /_/ /_/_____/_/|_|      *
*                                               *
*   Analysis of NMR Chemical Exchange data      *
*                                               *
*   Version: {chemex.__version__:<34s} *
*                                               *
* * * * * * * * * * * * * * * * * * * * * * * * *
"""


def main():
    """Do all the magic."""
    print(LOGO)

    parser = cc.build_parser()
    args = parser.parse_args()

    if args.commands is None:
        parser.print_help()
    else:
        args.func(args)


def fit(args):

    # Parse kinetics model
    model = cpk.parse_model(args.model)

    # Read initial values of fitting/fixed parameters
    defaults = cps.read_defaults(args.parameters)

    # Read experimental setup and data
    selection = {"include": args.include, "exclude": args.exclude}
    experiments = cce.read(args.experiments, model, selection, defaults)

    if not experiments:
        sys.exit("\nerror: No data to fit")

    # Create parameters
    params = cph.create_params(experiments, defaults)

    # Filter datapoints out if necessary (e.g., on-resonance filter CEST)
    experiments.filter(params)

    if args.commands == "simulate":
        _run_sim(args, experiments, params)
    else:
        _run_fit(args, defaults, experiments, params)


def _run_fit(args, defaults, experiments, params):
    ch.header1("Running the main fit")

    fitter = cf.Fit(experiments, args.out_dir, args.plot, defaults)
    fitter.read_methods(args.method)
    fitter.run_methods(params)


def _run_sim(args, experiments, params):

    ch.header1("Running the simulation")

    for param in params.values():
        if param.vary:
            param.vary = False

    path = args.out_dir
    plot = args.plot == "normal"

    coh.post_fit(experiments, params, path, plot=plot, simulation=True)

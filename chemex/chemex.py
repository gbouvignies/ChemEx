"""The chemex module provides the entry point for the chemex script."""
import sys

import chemex
import chemex.cli as cc
import chemex.containers.experiment as cce
import chemex.fitting as cf
import chemex.helper as ch
import chemex.parameters.kinetics as cpk
import chemex.parameters.settings as cpp


LOGO = r"""
* * * * * * * * * * * * * * * * * * * * * * * * *
*      ________                   ______        *
*     / ____/ /_  ___  ____ ___  / ____/  __    *
*    / /   / __ \/ _ \/ __ `__ \/ __/ | |/_/    *
*   / /___/ / / /  __/ / / / / / /____>  <      *
*   \____/_/ /_/\___/_/ /_/ /_/_____/_/|_|      *
*                                               *
*   Analysis of NMR Chemical Exchange data      *
*                                               *
*   Version: {:<34s} *
*                                               *
* * * * * * * * * * * * * * * * * * * * * * * * *
""".format(
    chemex.__version__
)


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
    model = cpk.parse_model(name=args.model)

    # Read experimental setup and data
    experiments = cce.read(
        filenames=args.experiments, model=model, selection=args.selection
    )

    # Create parameters
    params = experiments.params_default

    if not experiments:
        sys.exit("\nerror: No data to fit")

    # Update initial values of fitting/fixed parameters
    ch.header1("Reading default parameters")
    cpp.set_values(params, experiments, args.parameters)

    # Filter datapoints out if necessary (e.g., on-resonance filter CEST)
    experiments.filter(params)

    # Run the fit
    ch.header1("Running the main fit")
    fitter = cf.Fit(experiments, args.out_dir, args.plot)
    fitter.read_method(args.method)
    params_ = fitter.fit(params)
    if args.bs:
        fitter.bootstrap(params_, args.bs)
    elif args.mc:
        fitter.monte_carlo(params_, args.mc)

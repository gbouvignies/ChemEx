"""The chemex module provides the entry point for the chemex script."""
import sys

import chemex
import chemex.cli as cc
import chemex.containers.experiment as cce
import chemex.containers.plot as ccp
import chemex.fitting as cf
import chemex.helper as ch
import chemex.parameters.helper as cph
import chemex.parameters.kinetics as cpk
import chemex.parameters.settings as cps


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

    # Run the fit
    ch.header1("Running the main fit")
    fitter = cf.Fit(experiments, args.out_dir, args.plot)
    fitter.read_methods(args.method)
    params_fitted = fitter.run_methods(params)

    # Run Monte-Carlo simulations
    for mc_name in ("mc", "bs", "bsn"):
        n_iter = vars(args).get(mc_name)
        fitter.mc_simulations(params_fitted, n_iter, mc_name)


def simulate(args):

    # Parse kinetics model
    model = cpk.parse_model(name=args.model)

    # Read initial values of fitting/fixed parameters
    defaults = cps.read_defaults(args.parameters)

    # Read experimental setup and data
    selection = {"include": args.include, "exclude": args.exclude}
    experiments = cce.read(args.experiments, model, selection, defaults)

    if not experiments:
        sys.exit("\nerror: No data to simulate")

    # Create parameters
    params = experiments.params
    params.update_constraints()

    # Set to "fix" all parameters that are set to "fit"
    for param in params.values():
        if param.vary:
            param.vary = False

    # Filter datapoints out if necessary (e.g., on-resonance filter CEST)
    experiments.filter(params)

    # Run the simulation
    ch.header1("Running the simulation")
    path = args.out_dir
    print(f'\nWriting results -> "{path}/"')
    cps.write_par(params, path)
    ccp.write_plots(experiments, params, path, simulation=True)

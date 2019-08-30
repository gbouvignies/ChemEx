"""The chemex module provides the entry point for the chemex script."""
import copy
import shutil

import numpy as np

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
    experiments = cce.read(filenames=args.experiments, model=model)
    experiments.select(args.res_incl, args.res_excl)

    # Create and update initial values of fitting/fixed parameters
    ch.header1("Reading Default Parameters")
    params = experiments.params_default
    for filename in args.parameters:
        cpp.set_params_from_config_file(params, filename)

    # Filter datapoints out if necessary (e.g., on-resonance filter CEST)
    experiments.filter(params)

    # Customize the output directory
    out_dir = args.out_dir
    result = fit_write_plot(experiments, params, out_dir, args)

    if args.bs or args.mc:
        if args.bs:
            nmb = args.bs
        else:
            nmb = args.mc
        formatter_output_dir = "".join(["{:0", str(int(np.log10(nmb)) + 1), "d}"])
        for index in range(1, nmb + 1):
            if args.bs:
                experiments_index = experiments.bootstrap()
            else:
                experiments_index = experiments.monte_carlo(result.params)
            output_dir_ = out_dir / formatter_output_dir.format(index)
            params_mc = copy.deepcopy(result.params)
            fit_write_plot(experiments_index, params_mc, output_dir_, args)


def fit_write_plot(experiments, params, output_dir, args):
    """Perform the fit, write the output files and plot the results."""
    result = cf.run_fit(experiments, params, args.method, args.fitmethod)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_results(result, experiments, args.method, output_dir)
    if not args.noplot:
        plot_results(result, experiments, output_dir)
    return result


def write_results(result, experiments, method, output_dir):
    """Write the results of the fit to output files.

    The files below are created and contain the following information:
      - parameters.fit: fitting parameters and their uncertainties
      - contstraints.fit: expression used for constraining parameters
      - *.dat: experimental and fitted data
      - statistics.fit: statistics for the fit

    """
    ch.header1("Writing Results")

    print("\nFile(s):")

    if method:
        shutil.copyfile(method, output_dir / "fitting-method.cfg")

    cpp.write_par(result.params, path=output_dir)
    cpp.write_constraints(result.params, path=output_dir)
    experiments.write(path=output_dir, params=result.params)
    cf.write_statistics(result, path=output_dir)


def plot_results(result, experiments, path):
    """Plot the experimental and fitted data."""
    ch.header1("Plotting Data")
    print("\nFile(s):")
    path_plots = path / "Plots"
    path_plots.mkdir(parents=True, exist_ok=True)
    try:
        experiments.plot(path=path_plots, params=result.params)
    except KeyboardInterrupt:
        print("  - Plotting cancelled")

    # if result.method == "brute":
    #     labels = [
    #         cpn.ParamName.from_full_name(var).name.upper()
    #         for var in result.var_names
    #     ]
    #     outfile = path / "results_brute.pdf"
    #     plotting.plot_results_brute(result, varlabels=labels, output=outfile)
    #     print((f"  * {outfile}"))

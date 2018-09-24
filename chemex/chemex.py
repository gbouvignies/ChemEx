"""The chemex module provides the entry point for the chemex script."""

import copy
import os
import shutil

import numpy as np

from chemex import __version__, cli, datasets, fitting, parameters, util

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
    __version__
)


def main():
    """Do all the magic."""
    print(LOGO)

    parser = cli.build_parser(fit)
    args = parser.parse_args()

    if args.commands is None:
        parser.print_help()
    else:
        args.func(args)


def get_info(args):
    # cli.format_experiment_help(args.types, args.experiments)
    print(args.experiments)


def fit(args):
    # Read experimental setup and data
    data = datasets.read_data(args)

    # Create and update initial values of fitting/fixed parameters
    util.header1("Reading Default Parameters")
    params = parameters.create_params(data)

    for name in args.parameters:
        parameters.set_params_from_config_file(params, name)

    # Filter datapoints out if necessary (e.g., on-resonance filter CEST)
    for profile in data:
        profile.filter_points(params)

    data.ndata = sum([len(profile.val) for profile in data])

    # Customize the output directory
    output_dir = args.out_dir if args.out_dir else "./Output"

    if args.res_incl and len(args.res_incl) == 1:
        output_dir = os.path.join(output_dir, args.res_incl.pop().upper())

    result = fit_write_plot(args, params, data, output_dir)

    if args.bs or args.mc:
        if args.bs:
            nmb = args.bs
        else:
            nmb = args.mc

        formatter_output_dir = "".join(["{:0", str(int(np.log10(nmb)) + 1), "d}"])

        for index in range(1, nmb + 1):
            if args.bs:
                data_index = data.make_bs_dataset()
            else:
                data_index = data.make_mc_dataset(result.params)

            output_dir_ = os.path.join(output_dir, formatter_output_dir.format(index))

            params_mc = copy.deepcopy(result.params)

            fit_write_plot(args, params_mc, data_index, output_dir_)


def fit_write_plot(args, params, data, output_dir):
    """Perform the fit, write the output files and plot the results."""
    result = fitting.run_fit(args.method, params, data, args.fitmethod)

    util.make_dir(output_dir)

    write_results(result, data, args.method, output_dir)

    if not args.noplot:
        plot_results(result, data, output_dir)

    return result


def write_results(result, data, method, output_dir):
    """Write the results of the fit to output files.

    The files below are created and contain the following information:
      - parameters.fit: fitting parameters and their uncertainties
      - contstraints.fit: expression used for constraining parameters
      - *.dat: experimental and fitted data
      - statistics.fit: statistics for the fit

    """
    util.header1("Writing Results")

    print("\nFile(s):")

    if method:
        shutil.copyfile(method, os.path.join(output_dir, "fitting-method.cfg"))

    parameters.write_par(result.params, output_dir=output_dir)
    parameters.write_constraints(result.params, output_dir=output_dir)
    data.write_to(result.params, output_dir=output_dir)
    fitting.write_statistics(result, path=output_dir)


def plot_results(result, data, output_dir):
    """Plot the experimental and fitted data."""
    from chemex.experiments import plotting

    util.header1("Plotting Data")

    print("\nFile(s):")

    output_dir_plot = os.path.join(output_dir, "Plots")
    util.make_dir(output_dir_plot)

    try:
        plotting.plot_data(data, result.params, output_dir=output_dir_plot)
    except KeyboardInterrupt:
        print(" - Plotting cancelled")

    if result.method == "brute":
        labels = [
            parameters.ParameterName.from_full_name(var).name.upper()
            for var in result.var_names
        ]
        outfile = os.path.join(output_dir, "results_brute.pdf")
        plotting.plot_results_brute(result, varlabels=labels, output=outfile)
        print(("  * {}".format(outfile)))

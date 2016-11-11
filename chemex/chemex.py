"""The chemex module provides the entry point for the chemex script."""

import copy
import os
import shutil
import sys

import numpy as np

from chemex import datasets, fitting, parameters, parsing, util


def print_logo():
    """Print the ChemEx logo."""
    from chemex import __version__

    print((
        "\n"
        "* * * * * * * * * * * * * * * * * * * * * * * * *\n"
        "*      ________                   ______        *\n"
        "*     / ____/ /_  ___  ____ ___  / ____/  __    *\n"
        "*    / /   / __ \/ _ \/ __ `__ \/ __/ | |/_/    *\n"
        "*   / /___/ / / /  __/ / / / / / /____>  <      *\n"
        "*   \____/_/ /_/\___/_/ /_/ /_/_____/_/|_|      *\n"
        "*                                               *\n"
        "*   Analysis of NMR Chemical Exchange data      *\n"
        "*                                               *\n"
        "*   Version: {:<34s} *\n"
        "*                                               *\n"
        "* * * * * * * * * * * * * * * * * * * * * * * * *\n"
            .format(__version__)
    ))


def make_bootstrap_dataset(data):
    """Create a new dataset to run a bootstrap simulation."""
    # TODO: fix bootstrap error estimation, issue with copy.deepcoopy()
    data_bs = datasets.DataSet()

    for profile in data:
        data_bs.append(profile.make_bs_profile())

    return data_bs


def make_montecarlo_dataset(data, params):
    """Create a new dataset to run a Monte-Carlo simulation."""
    # TODO: fix Monte-Carlo error estimation, issue with copy.deepcoopy()
    data_mc = copy.deepcopy(data)

    for profile in data_mc:
        profile.val = profile.calculate_profile(params) + np.random.randn(len(profile.val)) * profile.err

    return data_mc


def read_data(args):
    """Read experimental setup and data."""
    util.header1("Reading Experimental Data")

    data = datasets.DataSet()

    if args.experiments:
        print(("{:<45s} {:<25s} {:<25s}".format("File Name", "Experiment", "Profiles")))
        print(("{:<45s} {:<25s} {:<25s}".format("---------", "----------", "--------")))
        for filename in args.experiments:
            data.add_dataset_from_file(filename, args.model, args.res_incl, args.res_excl)

    if not data.data:
        sys.exit("\nNo data to fit!\n")

    return data


def write_results(params, data, method, output_dir):
    """Write the results of the fit to output files.

    The files below are created and contain the following information:
      - parameters.fit: fitting parameters and their uncertainties
      - contstraints.fit: expression used for constraining parameters
      - *.dat: experimental and fitted data
      - chi2.fit: statistics for the fit
    """
    util.header1("Writing Results")

    print("\nFile(s):")

    if method:
        shutil.copyfile(method, os.path.join(output_dir, 'fitting-method.cfg'))

    # FIXME: change function name to "write_statistics_to" ?
    data.write_chi2_to(params, path=output_dir)
    parameters.write_par(params, output_dir=output_dir)
    parameters.write_constraints(params, output_dir=output_dir)
    data.write_to(params, output_dir=output_dir)


def plot_results(params, data, output_dir):
    """Plot the experimental and fitted data."""
    from chemex.experiments import plotting

    util.header1("Plotting Data")

    print("\nFile(s):")

    output_dir_plot = os.path.join(output_dir, 'Plots')
    util.make_dir(output_dir_plot)

    try:
        plotting.plot_data(data, params, output_dir=output_dir_plot)
    except KeyboardInterrupt:
        print(" - Plotting cancelled")


def fit_write_plot(args, params, data, output_dir):
    """Perform the fit, write the output files and plot the results."""
    params = fitting.run_fit(args.method, params, data)

    util.make_dir(output_dir)

    write_results(
        params,
        data,
        args.method,
        output_dir
    )

    if not args.noplot:
        plot_results(params, data, output_dir)

    return params


def main():
    """Do all the magic."""
    print_logo()

    args = parsing.arg_parse()

    if args.commands == 'info':
        parsing.format_experiment_help(args.types, args.experiments)

    elif args.commands == 'fit':

        # Read experimental setup and data
        data = read_data(args)

        # Create and update initial values of fitting/fixed parameters
        util.header1("Reading Default Parameters")
        params = parameters.create_params(data)
        parameters.set_params_from_config_file(params, args.parameters)

        # Filter datapoints out if necessary (e.g., on-resonance filter CEST)
        for profile in data:
            profile.filter_points(params)

        # Customize the output directory
        output_dir = args.out_dir if args.out_dir else './Output'
        if args.res_incl:
            if len(args.res_incl) == 1:
                output_dir = os.path.join(
                    output_dir, args.res_incl.pop().upper()
                )

        if not args.bs:
            params = fit_write_plot(
                args,
                params,
                data,
                output_dir
            )

        if args.bs or args.mc:

            if args.bs:
                n = int(args.bs)
            else:
                n = int(args.mc)

            formatter_output_dir = ''.join(
                ['{:0', str(int(np.log10(n)) + 1), 'd}']
            )

            for index in range(1, n + 1):

                if args.bs:
                    data_index = make_bootstrap_dataset(data)
                else:
                    data_index = make_montecarlo_dataset(data, params)

                output_dir_ = \
                    os.path.join(output_dir, formatter_output_dir.format(index))

                params_mc = copy.deepcopy(params)

                fit_write_plot(
                    args,
                    params_mc,
                    data_index,
                    output_dir_
                )

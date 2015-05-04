#!/usr/bin/env python

import os
import shutil
import random
import copy

from chemex import fitting
from chemex import writing
from chemex import parsing
from chemex import parameters
from chemex import utils
from chemex.experiments import reading as reading_exp
from chemex.experiments import utils as utils_exp


def print_logo():
    """Prints ChemEx logo"""

    from .version import __version__

    print(
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
        "*   Version: {:<15s}                    *\n"
        "*                                               *\n"
        "* * * * * * * * * * * * * * * * * * * * * * * * *\n"
        .format(__version__)
    )


def make_bootstrap_dataset(data):
    """Creates a new dataset to run a bootstrap simulation"""

    from random import choice

    profiles = {}

    for data_point in data:
        # The reference attribute is added to the profile id to separate the
        # reference points from the rest and make sure they are always present
        # in the bootstrapped sample
        reference = data_point.par.get('reference', False)
        profile_id = (data_point.par['profile_id'], reference)

        profiles.setdefault(profile_id, []).append(data_point)

    data_bs = []

    for profile in profiles.values():
        data_bs.extend(
            [choice(profile) for _ in xrange(len(profile))]
        )

    return data_bs


def make_montecarlo_dataset(data):
    """Creates a new dataset to run a Monte-Carlo simulation"""

    data_mc = copy.deepcopy(data)

    for data_pt in data_mc:
        data_pt.val = random.gauss(data_pt.cal, data_pt.err)

    return data_mc


def read_data(args):
    """Reads the files containing the experimental data point location and
    setup"""

    utils.header1("Reading Experimental Data")

    data = list()

    if args.experiments:
        print("\nFile(s):")
        for index, filename in enumerate(args.experiments, 1):
            print("  {}. {}".format(index, filename))
            data.extend(
                reading_exp.read_file_exp(
                    filename, args.res_incl, args.res_excl
                )
            )

    if not data:
        exit("\nNo Data to fit!\n")

    return data


def write_results(params, data, method, output_dir):
    """Writes the the chi2 of the fit, fitted parameters and the
    back-calculated points"""

    utils.header1("Writing Results")

    print("\nFile(s):")

    if method:
        shutil.copyfile(method, os.path.join(output_dir, 'fitting-method.cfg'))

    writing.write_chi2(params, data, output_dir=output_dir)
    writing.write_par(params, output_dir=output_dir)
    writing.write_dat(data, output_dir=output_dir)


def plot_results(params, data, output_dir):
    """Plots the the experimental and fitted points"""

    from chemex.experiments import plotting

    utils.header1("Plotting Data")

    print("\nFile(s):")

    output_dir_plot = os.path.join(output_dir, 'plots')
    utils.make_dir(output_dir_plot)

    try:
        plotting.plot_data(data, params, output_dir=output_dir_plot)
    except KeyboardInterrupt:
        print(" - Plotting cancelled")


def fit_write_plot(args, params, data, output_dir):
    params = fitting.run_fit(args.method, params, data)

    utils.make_dir(output_dir)

    write_results(
        params,
        data,
        args.method,
        output_dir
    )

    # Plot results
    if not args.noplot:
        plot_results(params, data, output_dir)

    return params


def main():
    """All the magic"""

    print_logo()

    args = parsing.arg_parse()

    if args.commands == 'info':

        utils_exp.format_experiment_help(args.types, args.experiments)

    elif args.commands == 'fit':

        # Read experimental points
        data = read_data(args)

        # Create the lists of both fitting and fixed parameters
        utils.header1("Reading Default Parameters")
        params = parameters.create_params(data)
        parameters.set_params_from_config_file(params, args.parameters)

        # Custom output directory
        output_dir = args.out_dir if args.out_dir else './output'
        if args.res_incl:
            if len(args.res_incl) == 1:
                output_dir = os.path.join(output_dir, args.res_incl[0].upper())

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
                ['{:0', utils.get_digit_number(n), 'd}']
            )

            for index in range(1, n + 1):

                if args.bs:
                    data_index = make_bootstrap_dataset(data)
                else:
                    data_index = make_montecarlo_dataset(data)

                output_dir_ = \
                    os.path.join(output_dir, formatter_output_dir.format(index))

                fit_write_plot(
                    args,
                    params,
                    data_index,
                    output_dir_
                )


if __name__ == '__main__':
    main()

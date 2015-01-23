#!/usr/bin/env python
"""
Created on Mar 31, 2011

@author: guillaume
"""

import os
import shutil
import random
from copy import deepcopy
from math import log10

from chemex import fitting, writing, parsing, reading, tools
from chemex.experiments.reading import read_file_exp
from chemex.experiments.misc import format_experiment_help


def print_logo():
    """ Prints ChemEx logo"""

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
        "* * * * * * * * * * * * * * * * * * * * * * * * *\n"
        "\n"
    )


def make_bootstrap_dataset(data):
    """Creates a new dataset to run a bootstrap simulation"""

    profiles = {}
    reference_points = {}

    for data_point in data:
        if data_point.par['reference']:
            reference_points.setdefault(data_point.par['profile_id'],
                []).append(data_point)
        else:
            profiles.setdefault(data_point.par['profile_id'], []).append(
                data_point)

    bootstrap_data = []

    for profile_id, profile in profiles.items():

        if profile_id in reference_points:
            bootstrap_data.extend(
                [random.choice(reference_points[profile_id]) for _ in
                 reference_points[profile_id]])

        bootstrap_data.extend([random.choice(profile) for _ in profile])

    return bootstrap_data


def make_montecarlo_dataset(data):
    """Creates a new dataset to run a Monte-Carlo simulation"""

    data_mc = deepcopy(data)

    for data_pt in data_mc:
        data_pt.val = random.gauss(data_pt.cal, data_pt.err)

    return data_mc


def read_data(args):
    """Reads the files containing the experimental data point location and
    setup"""

    tools.header1("Reading Experimental Data")

    data = list()

    if args.experiments:
        print("\nFile(s):")
        for index, filename in enumerate(args.experiments, 1):
            print("  {}. {}".format(index, filename))
            data.extend(read_file_exp(filename, args.res_incl, args.res_excl))

    if not data:
        exit("\nNo Data to fit!\n")

    return data


def write_results(par, par_err, par_indexes, par_fixed, data, method,
                  output_dir):
    """Writes the the chi2 of the fit, fitted parameters and the
    back-calculated points"""

    tools.header1("Writing Results")

    print("\nFile(s):")

    if method:
        shutil.copyfile(method, os.path.join(output_dir, 'fitting-method.cfg'))

    writing.write_chi2(par, par_indexes, par_fixed, data,
                       output_dir=output_dir)
    writing.write_par(par, par_err, par_indexes, par_fixed,
                      output_dir=output_dir)
    writing.write_dat(data, output_dir=output_dir)


def plot_results(par, par_indexes, par_fixed, data, output_dir):
    """Plots the the experimental and fitted points"""

    from chemex import plotting

    tools.header1("Plotting Data")

    print("\nFile(s):")

    output_dir_plot = os.path.join(output_dir, 'plots')
    tools.make_dir(output_dir_plot)

    try:
        plotting.plot_data(data, par, par_indexes, par_fixed,
                           output_dir=output_dir_plot)
    except KeyboardInterrupt:
        print(" - Plotting cancelled")


def fit_write_plot(args, par, par_indexes, par_fixed, data, output_dir):
    # Fit the data to the model
    par_fit, par_err, par_indexes, par_fixed = \
        fitting.run_fit(args.method, par, par_indexes, par_fixed, data)

    tools.make_dir(output_dir)

    write_results(par_fit, par_err, par_indexes, par_fixed, data, args.method,
                  output_dir)

    # Plot results
    if not args.noplot:
        plot_results(par_fit, par_indexes, par_fixed, data, output_dir)

    return par_fit, par_err, par_indexes, par_fixed


def main():
    """All the magic"""

    print_logo()

    args = parsing.arg_parse()

    if args.commands == 'info':

        format_experiment_help(args.types, args.experiments)

    elif args.commands == 'fit':

        # Read experimental points
        data = read_data(args)

        # Create the lists of both fitting and fixed parameters
        tools.header1("Reading Default Parameters")
        par, par_indexes, par_fixed, data = reading.create_par_list_to_fit(
            args.parameters, data)

        # Custom output directory
        output_dir = args.out_dir if args.out_dir else './output'
        if args.res_incl:
            if len(args.res_incl) == 1:
                output_dir = os.path.join(output_dir, args.res_incl[0].upper())

        if not args.bs:
            par, par_err, par_indexes, par_fixed = \
                fit_write_plot(args, par, par_indexes, par_fixed, data,
                               output_dir)

        if args.bs or args.mc:

            n = int(args.bs) if args.bs else int(args.mc)
            formatter_output_dir = ''.join(
                ['{:0', str(int(log10(n)) + 1), 'd}'])

            for index in range(1, n + 1):

                if args.bs:
                    data_index = make_bootstrap_dataset(data)
                elif args.mc:
                    data_index = make_montecarlo_dataset(data)

                output_dir_ = os.path.join(output_dir,
                                           formatter_output_dir.format(index))

                fit_write_plot(args, par, par_indexes, par_fixed, data_index,
                               output_dir_)


if __name__ == '__main__':
    main()

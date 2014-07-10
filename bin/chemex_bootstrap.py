#!/usr/bin/env python
"""
Created on Mar 31, 2011

@author: guillaume
"""

import os
import shutil
import random


# Local Libraries
try:
    from chemex import fitting, plotting, writing, parsing, reading, tools
    from chemex.experiments.reading import read_cfg_file as read_cfg_file_data

except (KeyboardInterrupt):
    exit("\n -- ChemEx killed before it could begin\n")


def make_bootstrap_dataset(data):
    profiles = {}
    reference_points = {}

    for data_point in data:
        if data_point.par['reference']:
            reference_points.setdefault(data_point.par['profile_id'], []).append(data_point)
        else:
            profiles.setdefault(data_point.par['profile_id'], []).append(data_point)

    bootstrap_data = []

    for profile_id, profile in profiles.iteritems():

        if profile_id in reference_points:
            bootstrap_data.extend([random.choice(reference_points[profile_id]) for _ in reference_points[profile_id]])

        bootstrap_data.extend([random.choice(profile) for _ in profile])

    return bootstrap_data


def main():
    """All the magic"""

    writing.print_logo()

    args = parsing.arg_parse()

    # Don't allow simultaneous include and exclude flags
    if args.res_incl and args.res_excl:
        exit('\nCan not simultaneously include and exclude residues!\n')
    elif args.res_incl:
        args.res_incl = [res.lower() for res in args.res_incl]
    elif args.res_excl:
        args.res_excl = [res.lower() for res in args.res_excl]

    # Read experimental data
    data = list()

    if args.experiments:
        for filename in args.experiments:
            data.extend(read_cfg_file_data(filename, args.res_incl, args.res_excl))

    if not data:
        exit("\nNo Data to fit!\n")

    # Create the lists of both fitting and fixed parameters
    par, par_indexes, par_fixed, data = reading.create_par_list_to_fit(args.parameters, data)

    data_bs = make_bootstrap_dataset(data)

    # Fit the data to the model
    par, par_err, par_indexes, par_fixed, reduced_chi2 = \
        fitting.run_fit(args.method, par, par_indexes, par_fixed, data_bs)

    # Write outputs
    print("")
    print("Reduced chi2: {:.3e}".format(reduced_chi2))

    # Custom output directory
    output_dir = args.out_dir if args.out_dir else './output'
    if args.res_incl:
        if len(args.res_incl) == 1:
            output_dir = os.path.join(output_dir, args.res_incl[0].upper())

    tools.make_dir(output_dir)

    print("")
    print(" - Writing results:")
    if args.method:
        shutil.copyfile(args.method, os.path.join(output_dir, 'fitting-method.cfg'))

    writing.write_chi2(par, par_indexes, par_fixed, data_bs, output_dir=output_dir)
    writing.write_par(par, par_err, par_indexes, par_fixed, output_dir=output_dir)
    writing.write_dat(data_bs, output_dir=output_dir)

    if not args.noplot:

        print("")
        print(" - Plotting data:")

        output_dir_plot = os.path.join(output_dir, 'plots')
        tools.make_dir(output_dir_plot)

        try:
            plotting.plot_data(data_bs, par, par_indexes, par_fixed, output_dir=output_dir_plot)
        except (KeyboardInterrupt):
            print(" - Plotting cancelled")


if __name__ == '__main__':
    main()

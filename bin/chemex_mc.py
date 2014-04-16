#!/usr/bin/env python
"""
Created on Mar 31, 2011

@author: guillaume
"""

import os
import random

from shutil import copyfile
from copy import deepcopy

try:
    from chemex import fitting, plotting, writing, parsing, reading, tools
    from chemex.experiments.reading import read_cfg_file as read_cfg_file_data
except (KeyboardInterrupt):
    exit("\n -- ChemEx killed before it could begin\n")


def main():
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
        exit("\nNo data to fit!\n")

    # Custom output directory
    output_dir = args.out_dir if args.out_dir else './'
    if args.res_incl:
        if len(args.res_incl) == 1:
            output_dir = os.path.join(output_dir, args.res_incl[0].upper())

    tools.make_dir(output_dir)

    # Copy the method to the output directory, as a backup
    method_file = None
    if args.method:
        copyfile(args.method, output_dir + "/fitting-method.cfg")

    # Create the lists of both fitting and fixed parameters
    par, par_indexes, par_fixed = reading.create_par_list_to_fit(args.parameters, data)

    # Fit the data to the model
    par, par_err, par_indexes, par_fixed, reduced_chi2 = \
        fitting.run_fit(args.method, par, par_indexes, par_fixed, data)


    for _ in range(1000):

        output_dir_ = os.path.join(output_dir, '{:03d}'.format(_))

        tools.make_dir(output_dir_)

        data_mc = deepcopy(data)

        for data_pt in data_mc:
            data_pt.val = random.gauss(data_pt.cal, data_pt.err)

        par_mc, par_err_mc, par_indexes_mc, par_fixed_mc, _rchi2 = \
            fitting.run_fit(args.method, par, par_indexes, par_fixed, data_mc)

        # Write outputs
        writing.write_chi2(par_mc, par_indexes_mc, par_fixed_mc, data_mc, output_dir=output_dir_)
        writing.write_par(par_mc, par_err_mc, par_indexes_mc, par_fixed_mc, output_dir=output_dir_)
        writing.write_dat(data_mc, output_dir=output_dir_)

        if not args.noplot:

            output_dir_plot = os.path.join(output_dir_, 'plots')
            tools.make_dir(output_dir_plot)

            print(' -- plotting data')
            try:
                plotting.plot_data(data_mc, par_mc, par_indexes_mc, par_fixed_mc, output_dir=output_dir_plot)
            except (KeyboardInterrupt):
                print('\n -- plotting cancelled')


if __name__ == '__main__':
    main()

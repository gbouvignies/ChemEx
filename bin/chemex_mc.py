#!/usr/bin/env python
"""
Created on Mar 31, 2011

@author: guillaume
"""

# Import
from copy import deepcopy
import os
import random
from shutil import copyfile

from chemex.chi2 import write_chi2, calc_chi2

try:
    from chemex import fitting, writing, parsing, reading, plotting
    from chemex.experiments.reading import read_cfg_file as read_cfg_file_data
except (KeyboardInterrupt):
    exit("\n -- ChemEx killed before it could begin\n")


def main():
    writing.print_logo()

    args = parsing.arg_parse()

    # Don't allow simultaneous include and exclude flags
    if args.res_incl and args.res_excl:
        exit('\nCan not simultaneously include and exclude residues!\n')

    # Reads and stores the datasets
    data = list()

    # Read experimental data
    if args.experiments:
        for filename in args.experiments:
            data.extend(read_cfg_file_data(filename, args.res_incl, args.res_excl))

    # Custom output directory
    output_dir = args.out_dir if args.out_dir else './'

    # Residue Selection
    if args.res_incl:
        if len(args.res_incl) == 1:
            output_dir = "/".join([output_dir, args.res_incl[0]])

    if not data:
        exit("\nNo Data to fit!\n")

    # Make the output directory if you need to
    if output_dir != './':
        if not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except OSError:
                exit("\nOSError: You can not use that directory!\n")

    # Copy the method to the output directory, as a backup
    method_file = None
    if args.method:
        copyfile(args.method, output_dir + "/fitting-method.cfg")

    # Create the lists of both fitting and fixed parameters
    par, par_indexes, par_fixed = reading.create_par_list_to_fit(args.parameters, data)

    calc_chi2(par, par_indexes, par_fixed, data)

    for _ in range(25):

        output_dir_ = os.path.join(output_dir, '{:03d}'.format(_))

        if not os.path.exists(output_dir_):
            try:
                os.makedirs(output_dir_)
            except OSError:
                exit("\nOSError: You can not use that directory!\n")

        data_mc = deepcopy(data)

        for data_pt in data_mc:
            data_pt.val = random.gauss(data_pt.cal, data_pt.err)

        par_mc, par_err_mc, par_indexes_mc, par_fixed_mc, _rchi2 = \
            fitting.run_fit(args.method, par, par_indexes, par_fixed, data_mc)

        # Write outputs
        write_chi2(par_mc, par_indexes_mc, par_fixed_mc, data_mc, output_dir=output_dir_)
        writing.write_par(par_mc, par_err_mc, par_indexes_mc, output_dir=output_dir_)
        writing.write_fixed_par(par_fixed_mc, output_dir=output_dir_)
        writing.write_dat(data_mc, output_dir=output_dir_)

        if not args.noplot:

            print(' -- plotting data')
            try:
                plotting.plot_data(data_mc, par_mc, par_indexes_mc, par_fixed_mc, output_dir=output_dir_)
            except (KeyboardInterrupt):
                print('\n -- plotting cancelled')


if __name__ == '__main__':
    main()

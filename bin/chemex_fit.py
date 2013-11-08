#!/usr/bin/env python
"""
Created on Mar 31, 2011

@author: guillaume
"""

# Import
import os
from shutil import copyfile

from chemex.chi2 import write_chi2


# Local Libraries
try:
    from chemex import fitting, plotting, writing, parsing, reading
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
            output_dir = "/".join([output_dir, args.res_incl[0].upper()])

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
    if args.method:
        copyfile(args.method, output_dir + "/fitting-method.cfg")

    # Create the lists of both fitting and fixed parameters
    par, par_indexes, par_fixed = reading.create_par_list_to_fit(args.parameters, data)

    # Fit the data to the model
    par, par_err, par_indexes, par_fixed, reduced_chi2 = \
        fitting.run_fit(args.method, par, par_indexes, par_fixed, data)

    # Write outputs
    print("")
    print("Reduced chi2: {:.3e}".format(reduced_chi2))
    print(" -- writing results")
    write_chi2(par, par_indexes, par_fixed, data, output_dir=output_dir)
    writing.write_par(par, par_err, par_indexes, par_fixed, output_dir=output_dir)
    writing.write_dat(data, output_dir=output_dir)

    if not args.noplot:

        print(' -- plotting data')
        try:
            plotting.plot_data(data, par, par_indexes, par_fixed, output_dir=output_dir)
        except (KeyboardInterrupt):
            print('\n -- plotting cancelled')


if __name__ == '__main__':
    main()

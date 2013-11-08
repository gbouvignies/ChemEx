"""
Created on Mar 31, 2011

@author: guillaume
"""

import os
import sys

from chemex.plotting import plot_data


def print_logo():
    """ Prints ChemEx logo to stdout """

    sys.stdout.write("""
* * * * * * * * * * * * * * * * * * * * * * * * *
*      ________                   ______        *
*     / ____/ /_  ___  ____ ___  / ____/  __    *
*    / /   / __ \/ _ \/ __ `__ \/ __/ | |/_/    *
*   / /___/ / / /  __/ / / / / / /____>  <      *
*   \____/_/ /_/\___/_/ /_/ /_/_____/_/|_|      *
*                                               *
*   Analysis of NMR Chemical Exchange data      *
*                                               *
* * * * * * * * * * * * * * * * * * * * * * * * *

""")


def write_dat(data, output_dir='./'):
    """Write dispersion profiles into a file"""

    datasets = dict()

    for data_point in data:
        experiment_name = data_point.par['experiment_name']
        datasets.setdefault(experiment_name, list()).append(data_point)

    for experiment_name, data in datasets.iteritems():

        filename = ''.join([experiment_name, '.dat'])
        filename = os.path.join(output_dir, filename)

        with open(filename, 'w') as f:
            for data_point in data:
                f.write(''.join([str(data_point), '\n']))


def write_par(par, par_err, par_indexes, par_fixed, output_dir='./'):
    """Write fitted parameters int a file"""

    if not par_indexes:
        return None

    filename = os.path.join(output_dir, 'parameters.fit')

    par_names = set(par_indexes) | set(par_fixed)

    with open(filename, 'w') as f:
        for par_name in sorted(par_names):

            if par_name in par_indexes:
                index = par_indexes[par_name]
                f.write(' '.join(str(a).upper() for a in par_name))
                f.write(' {: .5e} {: .5e}\n'.format(par[index], par_err[index]))

            elif par_name in par_fixed:
                f.write(' '.join(str(a).upper() for a in par_name))
                f.write(' {: .5e} fixed\n'.format(par_fixed[par_name]))


def dump_parameters(par, par_indexes, par_fixed, data):
    """ The program has failed. Dump parameters to chemex_dump """

    i = 0
    while os.path.exists('chemex_dump.' + str(i)):
        i += 1

    dump = 'chemex_dump.' + str(i)
    try:
        os.makedirs(dump)
    except OSError:
        exit("\nOSError: Cannot create the dump. Ending now.\n")

    sys.stderr.write("\n -- Writing current state to {:s}. Please wait ...".format(dump))
    try:
        write_par(par, par, par_indexes, par_fixed, output_dir=dump)
        write_dat(data, output_dir=dump)
        plot_data(data, par, par_indexes, par_fixed, output_dir=dump)
    except (TypeError, ValueError):
        sys.stderr.write("\n -- Save state cancelled. Not all data could not be plotted")
    #    except (KeyboardInterrupt, SystemExit):
    except (KeyboardInterrupt, SystemExit):
        exit("\n -- Dump has received a kill signal. Stopping immediately.\n")

#    exit("\n")


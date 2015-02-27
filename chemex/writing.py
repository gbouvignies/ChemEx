import os
import sys

import scipy as sc
import scipy.stats as st

from chemex.experiments import plotting


def write_dat(data, output_dir='./'):
    """Write dispersion profiles into a file"""

    datasets = dict()

    for data_point in data:
        experiment_name = data_point.par['experiment_name']
        datasets.setdefault(experiment_name, list()).append(data_point)

    for experiment_name, data in datasets.items():

        filename = ''.join([experiment_name, '.dat'])
        filename = os.path.join(output_dir, filename)

        print("  * {}".format(filename))

        with open(filename, 'w') as f:

            for data_point in data:
                f.write(''.join([str(data_point), '\n']))


def write_par(par, par_err, par_indexes, par_fixed, output_dir='./'):
    """Write fitted parameters int a file"""

    from ConfigParser import SafeConfigParser, DuplicateSectionError

    filename = os.path.join(output_dir, 'parameters.fit')

    print("  * {}".format(filename))

    par_names = set(par_indexes) | set(par_fixed)

    par_name_global = set(['KEX', 'KEX_AB', 'KEX_BC', 'KEX_AC', 'PB', 'PC'])

    par_dict = {}

    for name in par_names:

        if name in par_indexes:

            index = par_indexes[name]
            val = par[index]
            err = par_err[index]
            par_dict[name] = '{: .5e} {: .5e}'.format(val, err)

        else:

            val = par_fixed[name]
            par_dict[name] = '{: .5e} fixed'.format(val)

    cfg = SafeConfigParser()
    cfg.optionxform = str

    for name, val in sorted(par_dict.items()):

        name_list = list(name)

        if name_list[0].upper() in par_name_global:
            name_str = ', '.join([str(_).upper() for _ in name_list])
            section = 'global'

        else:
            name_str = str(name_list.pop(1)).upper()
            section = ', '.join([str(_).upper() for _ in name_list])

        try:
            cfg.add_section(section)

        except DuplicateSectionError:
            pass

        cfg.set(section, name_str, val)

    with open(filename, 'w') as f:
        cfg.write(f)

def write_chi2(par, par_indexes, par_fixed, data, output_dir='./'):
    """
    Write reduced chi2
    """

    data_nb = len(data)
    par_nb = len(par)

    residuals = sc.asarray(
        [data_point.calc_residual(par, par_indexes, par_fixed)
         for data_point in data])

    _ks_value, ks_p_value = st.kstest(residuals, 'norm')

    chi2 = sum(residuals ** 2)
    dof = data_nb - par_nb
    reduced_chi2 = chi2 / dof

    chi2_p_value = 1.0 - st.chi2.cdf(chi2, dof)

    filename = os.path.join(output_dir, 'chi2.fit')

    with open(filename, 'w') as f:
        print("  * {}".format(filename))

        f.write(
            '# {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}\n'
            .format('chi2', 'ndata', 'npar', 'rchi2', 'chi2-test', 'ks-test')
        )

        f.write(
            '  {: 15.5e} {: 15d} {: 15d} {: 15.5e} {: 15.5e} {: 15.5e}\n'
            .format(chi2, data_nb, par_nb, reduced_chi2, chi2_p_value,
                    ks_p_value)
        )


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

    sys.stderr.write(
        "\n - Writing current state to {:s}. Please wait ...".format(dump))
    try:
        write_par(par, par, par_indexes, par_fixed, output_dir=dump)
        write_dat(data, output_dir=dump)
        plotting.plot_data(data, par, par_indexes, par_fixed, output_dir=dump)

    except (TypeError, ValueError):
        sys.stderr.write(
            "\n - Save state cancelled. Not all data could not be plotted")

    # except (KeyboardInterrupt, SystemExit):
    except (KeyboardInterrupt, SystemExit):
        exit("\n - Dump has received a kill signal. Stopping immediately.\n")

# exit("\n")

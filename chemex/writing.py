import os
import sys

import scipy as sc
import scipy.stats as st

from chemex.experiments import plotting


def write_dat(data, output_dir='./'):
    """Write dispersion profiles into a file"""

    datasets = dict()

    for profile in data:
        experiment_name = profile.experiment_name
        datasets.setdefault(experiment_name, list()).append(profile)

    for experiment_name, data in datasets.items():

        filename = ''.join([experiment_name, '.dat'])
        filename = os.path.join(output_dir, filename)

        print("  * {}".format(filename))

        with open(filename, 'w') as f:

            for profile in data:
                f.write(profile.print_profile())


def write_par(params, output_dir='./'):
    """Write fitted parameters int a file"""

    from ConfigParser import SafeConfigParser, DuplicateSectionError

    filename = os.path.join(output_dir, 'parameters.fit')

    print("  * {}".format(filename))

    par_name_global = set(['KEX', 'KEX_AB', 'KEX_BC', 'KEX_AC', 'PB', 'PC'])

    par_dict = {}

    for name in params:

        val = params[name].value
        err = params[name].stderr

        if params[name].vary and err is not None:
            par_dict[name] = '{: .5e} {: .5e}'.format(val, err)
        else:
            par_dict[name] = '{: .5e} fixed'.format(val)

    cfg = SafeConfigParser()
    cfg.optionxform = str

    for name, val in sorted(par_dict.items()):

        name_list = name.replace('_p_', '.').split('__')

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


def write_chi2(params, data, output_dir='./'):
    """
    Write reduced chi2
    """

    data_nb = sum(len(profile.val) for profile in data)
    par_nb = len(params)

    residuals = sc.asarray([
        residual
        for profile in data
        for residual in profile.calculate_residuals(params)
    ])

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


def dump_parameters(params, data):
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
        "\n - Writing current state to {:s}. Please wait ...\n".format(dump))

    try:
        write_par(params, output_dir=dump)
        write_dat(data, output_dir=dump)
        plotting.plot_data(data, params, output_dir=dump)

    except (TypeError, ValueError):
        sys.stderr.write(
            "\n - Save state cancelled. Not all data could not be plotted")

    # except (KeyboardInterrupt, SystemExit):
    except (KeyboardInterrupt, SystemExit):
        exit("\n - Dump has received a kill signal. Stopping immediately.\n")

# exit("\n")

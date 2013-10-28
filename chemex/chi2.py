"""
Created on Mar 30, 2011

@author: guillaume
"""

# Standard Libraries
import sys
import os

# Specialized Libraries
import scipy as sc
import scipy.stats as st

# ChemEx Libraries
from chemex.writing import dump_parameters


def make_calc_residuals(verbose=True, threshold=1e-3):

    def calc_residuals(par, par_indexes, par_fixed, data):
        """
        Calculate the residuals for all values knowing the parameters par
        """

        try:
            residuals = [a_data_point.calc_residual(par, par_indexes, par_fixed)
                         for a_data_point in data]

        except KeyboardInterrupt:
            sys.stderr.write("\n -- Keyboard Interrupt: calculation stopped")
            dump_parameters(par, par_indexes, par_fixed, data)
            sys.exit()

        if verbose:

            chi2 = sum(sc.asarray(residuals) ** 2)

            if (calc_residuals.old_chi2 - chi2) / calc_residuals.old_chi2 > threshold:
                sys.stdout.write('Chi2: {:8.2e}\n'.format(chi2))
                sys.stdout.flush()
                calc_residuals.old_chi2 = chi2

        return residuals

    calc_residuals.old_chi2 = sys.float_info.max

    return calc_residuals


def calc_chi2(par, par_indexes, par_fixed, data):
    """
    Calculate the residuals for all values knowing the parameters par
    """

    chi2 = sum(data_point.calc_residual(par, par_indexes, par_fixed) ** 2
                     for data_point in data)

    return chi2


def calc_reduced_chi2(par, par_indexes, par_fixed, data):

    data_nb = len(data)
    par_nb = len(par)

    return calc_chi2(par, par_indexes, par_fixed, data) / (data_nb - par_nb)


def write_chi2(par, par_indexes, par_fixed, data, output_dir='./'):
    """
    Write reduced chi2
    """

    data_nb = len(data)
    par_nb = len(par)

    residuals = sc.asarray([data_point.calc_residual(par, par_indexes, par_fixed)
                            for data_point in data])

    _ks_value, ks_p_value = st.kstest(residuals, 'norm')

    chi2 = sum(residuals ** 2)
    dof = data_nb - par_nb
    reduced_chi2 = chi2 / dof

    chi2_p_value = 1.0 - st.chi2.cdf(chi2, dof)

    filename = os.path.join(output_dir, 'chi2.fit')

    with open(filename, 'w') as f:
        f.write('# {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}\n'
                .format('chi2', 'ndata', 'npar', 'rchi2', 'chi2-test', 'ks-test'))
        f.write('  {: 15.5e} {: 15d} {: 15d} {: 15.5e} {: 15.5e} {: 15.5e}\n'
                .format(chi2, data_nb, par_nb, reduced_chi2, chi2_p_value, ks_p_value))

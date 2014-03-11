"""
Created on Mar 30, 2011

@author: guillaume
"""

import sys
import scipy as sc

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

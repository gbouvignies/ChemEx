import sys

import scipy as sc

from chemex import writing


def make_calc_residuals(verbose=True, threshold=1e-3):
    def calc_residuals(params, data):
        """Calculate the residuals for all values knowing the parameters par"""

        try:
            residuals = [
                a_data_point.calc_residual(params)
                for a_data_point in data
            ]

        except KeyboardInterrupt:
            sys.stderr.write("\n -- Keyboard Interrupt: calculation stopped")
            writing.dump_parameters(params, data)
            sys.exit()

        if verbose:

            chi2 = sum(sc.asarray(residuals) ** 2)

            if (calc_residuals.old_chi2 - chi2) / calc_residuals.old_chi2 > threshold:
                sys.stdout.write(
                    '  * {:.3e} / {:.3e}\n'
                        .format(chi2, calc_reduced_chi2(params, data))
                )
                sys.stdout.flush()
                calc_residuals.old_chi2 = chi2

        return residuals

    calc_residuals.old_chi2 = sys.float_info.max

    return calc_residuals


def calc_chi2(params, data):
    """
    Calculate the residuals for all values knowing the parameters par
    """

    chi2 = sum(data_point.calc_residual(params) ** 2
               for data_point in data)

    return chi2


def calc_reduced_chi2(params, data):
    data_nb = len(data)
    par_nb = len(params)

    return calc_chi2(params, data) / (data_nb - par_nb)

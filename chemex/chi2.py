import sys

import scipy as sp


def make_calc_residuals(verbose=True, threshold=1e-3):
    def calc_residuals(params, data):
        """Calculate the residuals for all values knowing the parameters par"""

        try:
            residuals = sp.array([
                residual
                for profile in data
                for residual in profile.calculate_residuals(params)
            ])
        except KeyboardInterrupt:
            raise

        if verbose:

            chi2 = sum(residuals ** 2)

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

    chi2 = sum(
        point
        for profile in data
        for point in profile.calculate_residuals(params) ** 2
    )

    return chi2


def calc_reduced_chi2(params, data):
    data_nb = sum(len(profile.val) for profile in data)
    par_nb = len([param for param in params.values() if param.vary])
    return calc_chi2(params, data) / (data_nb - par_nb)

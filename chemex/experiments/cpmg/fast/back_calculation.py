from numpy.linalg import matrix_power
from scipy.linalg import expm

from ....bases.two_states.fast import P_180Y
from ....caching import lru_cache
from .liouvillian import compute_iy_eq, compute_liouvillians, get_iy


@lru_cache()
def make_calc_observable(time_t2=0.0, ppm_to_rads=1.0, _id=None):
    """
    Factory to make "calc_observable" function to calculate the intensity in
    presence of exchange after a CPMG block.

    Parameters
    ----------
    time_t2 : float
        Time of the CPMG block
    ncyc : integer
        Number of cycles, [t-180-2t-180-t]*n
    id : tuple
        Some type of identification for caching optimization

    Returns
    -------
    out : function
        Calculate intensity after the CEST block

    """

    @lru_cache(100)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, r_ixy=5.0, dr_ixy=0.0,
                         ncyc=0):
        """
        Calculate the intensity in presence of exchange during a cpmg-type pulse
        train.

        Parameters
        ----------
        i0 : float
            Initial intensity.
        pb : float
            Fractional population of state B,
            0.0 for 0%, 1.0 for 100%
        kex : float
            Exchange rate between state A and B in /s.
        dw : float
            Chemical shift difference between states A and B in rad/s.
        r_ixy : float
            Transverse relaxation rate of state a in /s.
        dr_ixy : float
            Transverse relaxation rate difference between states a and b in /s.

        Returns
        -------
        out : float
            Intensity after the CPMG block

        """

        dw *= ppm_to_rads

        mag_eq = compute_iy_eq(pb)

        if ncyc == 0:

            mag = mag_eq

        else:

            l_free = compute_liouvillians(
                pb=pb,
                kex=kex,
                dw=dw,
                r_ixy=r_ixy,
                dr_ixy=dr_ixy
            )

            t_cp = time_t2 / (4.0 * ncyc)
            p_free = expm(l_free * t_cp)

            mag = matrix_power(
                p_free
                .dot(P_180Y)
                .dot(p_free),
                2 * ncyc
            ).dot(mag_eq)

        magy_a, _ = get_iy(mag)

        return magy_a

    def calc_observable(i0=0.0, **kwargs):
        """
        Calculate the intensity in presence of exchange after a CPMG block.

        Parameters
        ----------
        i0 : float
            Initial intensity.

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        return i0 * _calc_observable(**kwargs)

    return calc_observable

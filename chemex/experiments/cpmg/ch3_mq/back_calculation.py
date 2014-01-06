'''
Created on Feb 23, 2012

@author: Mike Latham
@author: Guillaume Bouvignies
'''

from scipy import dot, diag
from scipy.linalg import expm
from numpy.linalg import matrix_power

from chemex.caching import lru_cache
from chemex.constants import gamma_ratio
from .liouvillian import \
    compute_2hxcy_eq, \
    get_2hxcy, \
    compute_liouvillian

RATIO = gamma_ratio['C']

# Constants
# Calculate the pulses
# Assuming all 180 degree pulses are perfect in the mq_ch3_8 basis set
P180_HX = diag([+1.0, +1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0], 0)
P180_HY = diag([-1.0, -1.0, +1.0, +1.0, -1.0, -1.0, +1.0, +1.0], 0)
P180_CX = diag([+1.0, -1.0, +1.0, -1.0, +1.0, -1.0, +1.0, -1.0], 0)
P180_CY = diag([-1.0, +1.0, -1.0, +1.0, -1.0, +1.0, -1.0, +1.0], 0)


@lru_cache()
def make_calc_observable(time_t2=0.0, ppm_to_rads_h=1.0, ppm_to_rads_c=1.0, smallflg='Y', _id=None):
    """
    Factory to make "calc_observable" function to calculate the intensity in presence
    of exchange after a CEST block.

    Parameters
    ----------
    pw : float
        Pulse width for a 90 degree pulse.
    time_T2 : float
        Time of the CPMG block.
    ncyc : integer
        Number of cycles, t-180-2t-180-t.
    id : tuple
        Some type of identification for caching optimization

    Returns
    -------
    out : function
        Calculate intensity after the CEST block

    """

    if smallflg == 'N':

        make_propagators = lru_cache(1)(compute_liouvillian)

        @lru_cache(100)
        def _calc_observable(pb=0.0, kex=0.0, dwh=0.0, dwc=0.0, r_mq=5.0, dr_mq=0.0, ncyc=0):
            """
            Calculate the intensity in presence of exchange during a cpmg-type pulse
            train. Based on the sequence "hmqc_CH3_exchange_bigprotein_600_lek_v2.c".
            Parameter names from the sequence are:

            chemex      procpar
            ------      -------
            time_T2     time_T2
            ncyc        ncyc_cp

            Parameters
            ----------
            i0 : float
                Initial intensity
            pb : float
                Fractional population of state b
            kex : float
                Exchange rate between state a and b in /s
            dwh : float
                Proton chemical shift difference between states A and B in ppm
            dwc : float
                Carbon chemical shift difference between states A and B in ppm
            r_mq : float
                Multiple quantum relaxation rate in /s
            dr_mq : float
                Multiple quantum relaxation rate difference between a and b in /s
            smallflg: string ['y', 'n']
               Flag to include small_protein_flg block in fitting

            Returns
            -------
            out : float
                Intensity after the CPMG block

            """

            dwh *= ppm_to_rads_h
            dwc *= ppm_to_rads_c

            l_free = make_propagators(
                pb=pb, kex=kex, dwh=dwh, dwc=dwc, r_mq=r_mq, dr_mq=dr_mq,
            )

            mag_eq = compute_2hxcy_eq(pb)

            if ncyc == 0:
                mag = mag_eq

            else:

                t_cp = time_t2 / (4.0 * ncyc)
                p_free = expm(l_free * t_cp)
                p_cpy = matrix_power(p_free.dot(P180_CY).dot(p_free), ncyc)

                mag = reduce(dot, [p_cpy, P180_HX, p_cpy, mag_eq])

            magz_a, _ = get_2hxcy(mag)

            return magz_a

    else:

        @lru_cache(1)
        def make_propagators(pb=0.0, kex=0.0, dwh=0.0, dwc=0.0, r_mq=5.0, dr_mq=0.0):

            l_free = compute_liouvillian(
                pb=pb, kex=kex, dwh=dwh, dwc=dwc, r_mq=r_mq, dr_mq=dr_mq,
            )

            p_zeta = expm(l_free * 1.0 / (8.0 * 125.3))

            return l_free, p_zeta

        @lru_cache(100)
        def _calc_observable(pb=0.0, kex=0.0, dwh=0.0, dwc=0.0, r_mq=5.0, dr_mq=0.0, ncyc=0):
            """
            Calculate the intensity in presence of exchange during a cpmg-type pulse
            train. Based on the sequence "hmqc_CH3_exchange_bigprotein_600_lek_v2.c".
            Parameter names from the sequence are:

            chemex      procpar
            ------      -------
            time_T2     time_T2
            ncyc        ncyc_cp

            Parameters
            ----------
            i0 : float
                Initial intensity
            pb : float
                Fractional population of state b
            kex : float
                Exchange rate between state a and b in /s
            dwh : float
                Proton chemical shift difference between states A and B in ppm
            dwc : float
                Carbon chemical shift difference between states A and B in ppm
            r_mq : float
                Multiple quantum relaxation rate in /s
            dr_mq : float
                Multiple quantum relaxation rate difference between a and b in /s
            smallflg: string ['y', 'n']
               Flag to include small_protein_flg block in fitting

            Returns
            -------
            out : float
                Intensity after the CPMG block

            """

            dwh *= ppm_to_rads_h
            dwc *= ppm_to_rads_c

            l_free, p_zeta = make_propagators(
                pb=pb, kex=kex, dwh=dwh, dwc=dwc, r_mq=r_mq, dr_mq=dr_mq,
            )

            mag_eq = compute_2hxcy_eq(pb)

            if ncyc == 0:

                mag = reduce(dot, [p_zeta, P180_HX, P180_CX, p_zeta, mag_eq])

            else:

                t_cp = time_t2 / (4.0 * ncyc)
                p_free = expm(l_free * t_cp)
                p_cpy = matrix_power(p_free.dot(P180_CY).dot(p_free), ncyc)

                mag = reduce(dot, [p_zeta, P180_HX, P180_CX, p_zeta, p_cpy, P180_HX, p_cpy, mag_eq])

            magz_a, _ = get_2hxcy(mag)

            return magz_a


    def calc_observable(i0=0.0, **kwargs):
        """
        Calculate the intensity in presence of exchange after a CEST block.

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
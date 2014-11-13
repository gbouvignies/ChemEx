import scipy as sc
from scipy.linalg import eigvals

from chemex.experiments.misc import correct_chemical_shift
from chemex.caching import lru_cache
from .liouvillian import compute_liouvillian


@lru_cache()
def make_calc_observable(time_t1=0.0, b1_offset=0.0, b1_frq=0.0, carrier=0.0, ppm_to_rads=0.0, _id=None):
    """
    Factory to make "calc_observable" function to calculate the intensity in presence
    of exchange after a CEST block.

    Parameters
    ----------
    time_t1 : float
        Duration of the CW block.
    b1_offset : float
        Frequency offset of the applied B1 field in Hz.
    b1_frq : float
        Strength of the applied B1 field in Hz.
    b1_inh : float
        B1 field inhomogeneity in Hz.
    b1_inh_res : int
        Resolution to model B1 field inhomogeneity.
    carrier : float
        Carrier position in rad/s.
    ppm_to_rads : float
        Conversion factor from ppm to rad/s
    id : tuple
        Some type of identification for caching optimization

    Returns
    -------
    out : function
        Calculate intensity after the CEST block

    """

    w1 = b1_frq * 2.0 * sc.pi
    w1_offset = b1_offset * 2.0 * sc.pi

    @lru_cache(5)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, r_nz=1.5, r_nxy=0.0, dr_nxy=0.0, cs=0.0):
        """
        Calculate the intensity in presence of exchange after a CEST block assuming
        initial intensity of 1.0.

        Parameters
        ----------
        pb : float
            Fractional population of state B,
            0.0 for 0%, 1.0 for 100%
        kex : float
            Exchange rate between state A and B in /s.
        dw : float
            Chemical shift difference between states A and B in rad/s.
        r_nz : float
            Longitudinal relaxation rate of state {a,b} in /s.
        r_nxy : float
            Transverse relaxation rate of state a in /s.
        dr_nxy : float
            Transverse relaxation rate difference between states a and b in /s.
        cs : float
            Resonance position in rad/s.

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        if abs(b1_offset) >= 10000.0:

            return (1.0 - pb)

        else:

            dw *= ppm_to_rads

            exchange_induced_shift, _ = correct_chemical_shift(pb=pb, kex=kex, dw=dw, r_ixy=r_nxy, dr_ixy=dr_nxy)

            wg = (cs - carrier) * ppm_to_rads - exchange_induced_shift - w1_offset

            liouvillian = compute_liouvillian(pb=pb, kex=kex, dw=dw, r_nxy=r_nxy, dr_nxy=dr_nxy, r_nz=r_nz,
                                              cs_offset=wg, w1=w1)

            # eval, evec = eig(liouvillian)
            #
            # index, r1 = max(enumerate(-sc.absolute(eval)), key=operator.itemgetter(1))
            # coeff = max(sc.absolute(evec[2:6:3, index])) ** 2
            #
            #
            # magz_a = (1.0 - pb) * sc.exp(r1 * time_t1) * coeff

            r1 = -sc.sort(sc.absolute(eigvals(liouvillian)))[0]

            magz_a = (1.0 - pb) * sc.exp(r1 * time_t1) * wg ** 2 / (wg ** 2 + w1 ** 2)

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

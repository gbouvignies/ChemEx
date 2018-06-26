import scipy as sc
from scipy.linalg import expm

from ....caching import lru_cache
from .liouvillian import compute_nz_eq, compute_base_liouvillians, \
    compute_free_liouvillian, get_nz


@lru_cache()
def make_calc_observable(time_t1=0.0, b1_offset=0.0, b1_frq=0.0, b1_inh=0.0,
                         b1_inh_res=5, carrier=0.0, ppm_to_rads=0.0, _id=None):
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

    base_liouvillians, weights = \
        compute_base_liouvillians(b1_offset, b1_frq, b1_inh, b1_inh_res)

    @lru_cache(5)
    def _calc_observable(pb=0.0, pc=0.0, kex_ab=0.0, kex_bc=0.0, kex_ac=0.0,
                         dw_ab=0.0, dw_ac=0.0, r_nz=1.5, r_nxy=0.0, dr_nxy_ab=0.0, dr_nxy_ac=0.0, cs=0.0):
        """
        Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        pb : float
            Fractional population of state B,
            0.0 for 0%, 1.0 for 100%
        pc : float
            Fractional population of state C,
            0.0 for 0%, 1.0 for 100%
        kex_ab : float
            Exchange rate between state A and B in /s.
        kex_bc : float
            Exchange rate between state B and C in /s.
        kex_ac : float
            Exchange rate between state A and C in /s.
        dw_ab : float
            Chemical shift difference between states A and B in rad/s.
        dw_ac : float
            Chemical shift difference between states A and C in rad/s.
        r_nz : float
            Longitudinal relaxation rate of state {a,b} in /s.
        r_nxy : float
            Transverse relaxation rate of state a in /s.
        dr_nxy_ab : float
            Transverse relaxation rate difference between states a and c in /s.
        dr_nxy_ac : float
            Transverse relaxation rate difference between states a and c in /s.
        cs : float
            Resonance position in rad/s.

        Returns
        -------
        out : float
            Intensity after the CEST block

        """

        if abs(b1_offset) >= 10000.0:

            return 1.0 - pb - pc

        else:

            dw_ab *= ppm_to_rads
            dw_ac *= ppm_to_rads

            mag_eq = compute_nz_eq(pb, pc)

            exchange_induced_shift = 0.0  # TODO
            wg = (cs - carrier) * ppm_to_rads - exchange_induced_shift

            liouvillians = \
                base_liouvillians + \
                compute_free_liouvillian(
                    pb=pb,
                    pc=pc,
                    kex_ab=kex_ab,
                    kex_bc=kex_bc,
                    kex_ac=kex_ac,
                    dw_ab=dw_ab,
                    dw_ac=dw_ac,
                    r_nxy=r_nxy,
                    r_nz=r_nz,
                    dr_nxy_ab=dr_nxy_ab,
                    dr_nxy_ac=dr_nxy_ac,
                    cs_offset=wg
            )

            propagator = sc.zeros_like(liouvillians[0])
            for liouvillian, weight in zip(liouvillians, weights):
                propagator += weight * expm(liouvillian * time_t1)
            propagator /= sum(weights)

            magz_a, _, _ = get_nz(sc.dot(propagator, mag_eq))

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

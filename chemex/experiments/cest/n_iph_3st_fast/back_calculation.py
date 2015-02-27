import scipy as sp
from scipy.linalg import eig, inv

from ....caching import lru_cache
from .liouvillian import compute_liouvillian


PI = sp.pi
dot = sp.dot
diag = sp.diag
exp = sp.exp
ix_ = sp.ix_


@lru_cache()
def make_calc_observable(time_t1=0.0, b1_offset=0.0, b1_frq=0.0, carrier=0.0,
                         ppm_to_rads=0.0, _id=None):
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

    w1 = b1_frq * 2.0 * PI
    w1_offset = b1_offset * 2.0 * PI

    @lru_cache(5)
    def _calc_observable(pb=0.0, pc=0.0, kex_ab=0.0, kex_bc=0.0, kex_ac=0.0,
                         dw_ab=0.0, dw_ac=0.0, r_nz=1.5, r_nxy=0.0,
                         dr_nxy_ab=0.0, dr_nxy_ac=0.0, cs=0.0):
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

            magz_a = 1.0 - pb - pc

        else:

            dw_ab *= ppm_to_rads
            dw_ac *= ppm_to_rads

            exchange_induced_shift = 0.0  # TODO

            wg = (
                (cs - carrier) * ppm_to_rads -
                exchange_induced_shift -
                w1_offset
            )

            liouvillian = compute_liouvillian(
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
                cs_offset=wg,
                w1=w1
            )

            s, vr = eig(liouvillian)
            vri = inv(vr)

            sl1 = [2, 5, 8]  # za, zb, zc
            sl2 = [i for i, w in enumerate(s.imag) if abs(w) < 1.0e-6]
            sl3 = [2]  # za

            vri = vri[ix_(sl2, sl1)].real
            t = diag(exp(s[sl2].real * time_t1))
            vr = vr[ix_(sl3, sl2)].real
            magz_eq = sp.asarray([[1 - pb - pc], [pb], [pc]])

            magz_a = dot(dot(dot(vr, t), vri), magz_eq)[0, 0]

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

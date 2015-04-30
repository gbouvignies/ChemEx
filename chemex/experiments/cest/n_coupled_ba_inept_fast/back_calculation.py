import scipy as sp
from scipy.linalg import eig, inv, expm
from chemex.experiments.misc import correct_chemical_shift, correct_intensities

from chemex.caching import lru_cache
from .liouvillian import compute_liouvillian


dot = sp.dot
diag = sp.diag
exp = sp.exp


@lru_cache()
def make_calc_observable(time_t1=0.0, b1_offset=0.0, b1_frq=0.0, carrier=0.0,
                         ppm_to_rads=0.0, ppm_to_rads_h=0.0, multiplet=None, _id=None):
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

    w1 = b1_frq * 2.0 * sp.pi
    w1_offset = b1_offset * 2.0 * sp.pi


    @lru_cache(5)
    def _calc_observable(pb=0.0, kex=0.0, dw=0.0, dw_h=0.0, r_nz=1.5, r_nxy=0.0,
                         dr_nxy=0.0, cs=0.0):
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

        dw_h *= ppm_to_rads_h
        dw *= ppm_to_rads

        kab = pb * kex
        kba = kex - kab

        liouv_h1 = sp.array([[-kab, +kba],
                             [+kab, -kba - dr_nxy + 1j * dw_h]])

        liouv_n1 = sp.array([[-kab, +kba],
                             [+kab, -kba - dr_nxy + 1j * dw]])

        magz = sp.array([[1.0 - pb], [pb]])

        t1 = 2.25e-3
        t2 = 2.75e-3

        magz = (expm(liouv_h1.conjugate() * t1).dot(expm(liouv_h1 * t1)).dot(magz)).real
        magz = (expm(liouv_n1.conjugate() * t2).dot(expm(liouv_n1 * t2)).dot(magz)).real

        if abs(b1_offset) < 10000.0:

            exchange_induced_shift, _ = correct_chemical_shift(
                pb=pb,
                kex=kex,
                dw=dw,
                r_ixy=r_nxy,
                dr_ixy=dr_nxy
            )

            wg = (
                (cs - carrier) * ppm_to_rads -
                exchange_induced_shift -
                w1_offset
            )

            magz_ = sp.zeros_like(magz)

            for j, weight in multiplet:
                liouvillian = compute_liouvillian(
                    pb=pb,
                    kex=kex,
                    dw=dw,
                    r_nxy=r_nxy,
                    dr_nxy=dr_nxy,
                    r_nz=r_nz,
                    cs_offset=(wg + j),
                    w1=w1
                )

                s, vr = eig(liouvillian)
                vri = inv(vr)

                sl1 = [2, 5]
                sl2 = [i for i, w in enumerate(s.imag) if abs(w) < 1.0e-6]
                sl3 = [2, 5]

                vri = vri[sp.ix_(sl2, sl1)].real
                t = diag(exp(s[sl2].real * time_t1))
                vr = vr[sp.ix_(sl3, sl2)].real

                magz_ += weight * dot(dot(dot(vr, t), vri), magz)

            magz = sp.array(magz_)

        magz_a, magz_b = magz[:, 0]

        magz_a, magz_b = correct_intensities(magz_a, magz_b, pb, kex, dw, r_nxy, dr_nxy)

        magz = sp.array([[magz_a], [magz_b]])

        magz = (expm(liouv_n1.conjugate() * t2).dot(expm(liouv_n1 * t2)).dot(magz)).real
        magz = (expm(liouv_h1.conjugate() * t1).dot(expm(liouv_h1 * t1)).dot(magz)).real

        magz_a, magz_b = magz[:, 0]

        magz_a, magz_b = correct_intensities(magz_a, magz_b, pb, kex, dw, 0.0, dr_nxy)

        return magz_b


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



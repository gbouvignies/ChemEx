import scipy as sp
import scipy.linalg as sl

import chemex.caching as caching
from chemex.experiments.utils import correct_chemical_shift
from chemex.constants import scalar_couplings
from .liouvillian import compute_liouvillian


inv = sl.inv
eig = sl.eig
ix_ = sp.ix_
diag = sp.diag
asarray = sp.asarray
dot = sp.dot
exp = sp.exp

TWO_PI = 2.0 * sp.pi

JHN = scalar_couplings['amide_HN']


@caching.lru_cache()
def make_calc_observable(time_t1=0.0, b1_offset=0.0, b1_frq=0.0, b1_frq_h=0.0,
                         carrier=0.0, carrier_h=0.0, ppm_to_rads=0.0,
                         ppm_to_rads_h=0.0, _id=None):
    @caching.lru_cache(5)
    def _calc_observable(pb=0.0, kex=0.0, dw_h=0.0, dw_n=0.0, r_nxy=5.0,
                         dr_nxy=None, r_nz=1.5, r_2hznz=None, r_2hxynxy=0.0,
                         r_hxy=10.0, r_hz=1.0, etaxy=0.0, etaz=0.0, j_hn=JHN,
                         cs_n=0.0, cs_h=0.0):
        """
        Calculate the intensity in presence of exchange during a cpmg-type pulse train.

        Keyword arguments:
        I0 -- Initial intensity,
              0.0 (default),
        pb -- population of state B,
              0.0 for 0% (default),
              1.0 for 100%
        kex -- exchange rate between state A and B in /s,
               0.0 (default)
        dw -- chemical shift difference between states A and B in /s,
              0.0 /s (default)
        r_Nxy -- transverse relaxation rate in /s,
                 5.0 (default)
        dr_Nxy -- transverse relaxation rate difference in /s between states a and b,
                  0.0 (default)
        r_Nz -- longitudinal relaxation rate in /s,
                1.5 (default)
        cs_offset -- chemical shift from the carrier in rad/s,
                     0.0 (default)
        B1_offset -- frequency offset of the applied B1 field in Hz,
                     0.0 Hz (default)
        B1_frq -- strength of the applied B1 field in Hz,
        time_t1 -- time of the CW block,
                   0 ms (default)

        Returns: float
        """

        wg_h = (cs_h - carrier_h) * ppm_to_rads_h
        wg_n = (cs_n - carrier) * ppm_to_rads
        dw_h_rads = dw_h * ppm_to_rads_h
        dw_n_rads = dw_n * ppm_to_rads

        if abs(b1_offset) > 9999.0:

            magz_a = 1.0 - pb

        else:

            exchange_induced_shift_n, _ = correct_chemical_shift(
                pb=pb,
                kex=kex,
                dw=dw_n_rads,
                r_ixy=r_nxy,
                dr_ixy=dr_nxy
            )

            offset_n = wg_n - exchange_induced_shift_n - b1_offset * TWO_PI

            liouvillian = compute_liouvillian(
                pb=pb,
                kex=kex,
                dw_h=dw_h_rads,
                dw_n=dw_n_rads,
                r_nxy=r_nxy,
                dr_nxy=dr_nxy,
                r_nz=r_nz,
                r_2hznz=r_2hznz,
                r_2hxynxy=r_2hxynxy,
                r_hxy=r_hxy,
                r_hz=r_hz,
                etaxy=etaxy,
                etaz=etaz,
                cs_offset_h=wg_h,
                cs_offset_n=offset_n,
                j_hn=j_hn,
                w1_h=TWO_PI * b1_frq_h,
                w1_n=TWO_PI * b1_frq,
            )

            s, vr = eig(liouvillian)
            vri = inv(vr)

            sl1 = [5, 20]
            sl2 = [i for i, w in enumerate(s.imag) if abs(w) < 1.0e-6]
            sl3 = [5]

            vri = vri[ix_(sl2, sl1)].real
            t = diag(exp(s[sl2].real * time_t1))
            vr = vr[ix_(sl3, sl2)].real
            magz_eq = asarray([[1 - pb], [pb]])

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


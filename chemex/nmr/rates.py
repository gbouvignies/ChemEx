import functools as ft

import numpy as np
import scipy.constants as sc

import chemex.nmr.constants as cnc


def _calculate_jw(tauc, s2, w):
    return 2.0 / 5.0 * tauc * s2 / (1.0 + (w * tauc) ** 2)


class RatesIS:
    gi: float
    gs: float
    ris3: float
    rih3: float
    rsh3: float
    csa_i: np.ndarray
    csa_s: np.ndarray
    phi_i: np.ndarray
    phi_s: np.ndarray
    j_is: float

    def __init__(self):
        self.gh = cnc.GAMMA["h"]

        # Dipolar factors
        self.dis = -sc.mu_0 * sc.hbar * self.gi * self.gs / (8.0 * np.pi * self.ris3)
        dih = -sc.mu_0 * sc.hbar * self.gi * self.gh / (8.0 * np.pi * self.rih3)
        dsh = -sc.mu_0 * sc.hbar * self.gs * self.gh / (8.0 * np.pi * self.rsh3)
        self.dis2 = self.dis * self.dis
        self.dih2 = dih * dih
        self.dsh2 = dsh * dsh

        # CSA factors
        self.delta_i = self.csa_i[0:2] - self.csa_i[2]
        self.delta_s = self.csa_s[0:2] - self.csa_s[2]

    def __call__(self, h_frq, tauc, s2):

        # B0 in Tesla
        b0 = 2.0 * np.pi * 1e6 * h_frq / self.gh

        # Larmor angular frequencies
        wi, ws, wh = -np.array([self.gi, self.gs, self.gh]) * b0

        # CSA factors
        ci_ = -wi * self.delta_i / 3.0
        cs_ = -ws * self.delta_s / 3.0
        ci2 = ci_[0] ** 2 + ci_[1] ** 2 - ci_[0] * ci_[1]
        cs2 = cs_[0] ** 2 + cs_[1] ** 2 - cs_[0] * cs_[1]
        p2i = 0.5 * (3.0 * np.cos(self.phi_i[0:2]) ** 2 - 1.0)
        p2s = 0.5 * (3.0 * np.cos(self.phi_s[0:2]) ** 2 - 1.0)
        ci = (ci_ * p2i).sum()
        cs = (cs_ * p2s).sum()

        j0, ji, js, jh, jmis, jpis, jmih, jpih, jmsh, jpsh = _calculate_jw(
            tauc,
            s2,
            np.array(
                [0.0, wi, ws, wh, wi - ws, wi + ws, wi - wh, wi + wh, ws - wh, ws + wh]
            ),
        )

        return {
            "r2_i": 0.5
            * (
                ci2 * (4 * j0 + 3 * ji)
                + self.dis2 * (4 * j0 + 3 * ji + jmis + 6 * js + 6 * jpis)
                + self.dih2 * (4 * j0 + 3 * ji + jmih + 6 * jh + 6 * jpih)
            ),
            "r2_s": 0.5
            * (
                cs2 * (4 * j0 + 3 * js)
                + self.dis2 * (4 * j0 + 3 * js + jmis + 6 * ji + 6 * jpis)
                + self.dsh2 * (4 * j0 + 3 * js + jmsh + 6 * jh + 6 * jpsh)
            ),
            "r1_i": (
                ci2 * 3 * ji
                + self.dis2 * (3 * ji + jmis + 6 * jpis)
                + self.dih2 * (3 * ji + jmih + 6 * jpih)
            ),
            "r1_s": (
                cs2 * 3 * js
                + self.dis2 * (3 * js + jmis + 6 * jpis)
                + self.dsh2 * (3 * js + jmsh + 6 * jpsh)
            ),
            "r2a_i": (
                0.5 * ci2 * (4 * j0 + 3 * ji)
                + cs2 * 3 * js
                + 0.5 * self.dis2 * (4 * j0 + 3 * ji + jmis + 6 * jpis)
                + 0.5 * self.dih2 * (4 * j0 + 3 * ji + jmih + 6 * jh + 6 * jpih)
                + self.dsh2 * (3 * js + jmsh + 6 * jpsh)
            ),
            "r2a_s": (
                0.5 * cs2 * (4 * j0 + 3 * js)
                + ci2 * 3 * ji
                + 0.5 * self.dis2 * (4 * j0 + 3 * js + jmis + 6 * jpis)
                + 0.5 * self.dsh2 * (4 * j0 + 3 * js + jmsh + 6 * jh + 6 * jpsh)
                + self.dih2 * (3 * ji + jmih + 6 * jpih)
            ),
            "r2mq_is": 0.5
            * (
                ci2 * (4 * j0 + 3 * ji)
                + cs2 * (4 * j0 + 3 * js)
                + self.dis2 * (3 * ji + jmis + 3 * js + 6 * jpis)
                + self.dih2 * (4 * j0 + 3 * ji + jmih + 6 * jh + 6 * jpih)
                + self.dsh2 * (4 * j0 + 3 * js + jmsh + 6 * jh + 6 * jpsh)
            ),
            "r1a_is": (
                ci2 * 3 * ji
                + cs2 * 3 * js
                + self.dis2 * (3 * ji + 3 * js)
                + self.dih2 * (3 * ji + jmih + 6 * jpih)
                + self.dsh2 * (3 * js + jmsh + 6 * jpsh)
            ),
            "etaxy_i": self.dis * ci * (4 * j0 + 3 * ji),
            "etaxy_s": self.dis * cs * (4 * j0 + 3 * js),
            "etaz_i": self.dis * ci * 6 * ji,
            "etaz_s": self.dis * cs * 6 * js,
            "sigma_is": self.dis2 * (-jmis + 6.0 * jpis),
            "mu_is": 0.5 * self.dis2 * (-jmis + 6.0 * jpis),
            "j_is": self.j_is,
        }


class RateNH(RatesIS):
    gi = cnc.GAMMA["n"]
    gs = cnc.GAMMA["h"]
    ris3 = 1.04e-10 ** 3
    rih3 = 1.79e-10 ** 3
    rsh3 = 1.85e-10 ** 3
    csa_i = np.array([69.0, 42.0, -111.0]) * 1e-6
    csa_s = np.array([5.7, 0.5, -6.2]) * 1e-6
    phi_i = np.deg2rad([109.6, 90.0, 19.6])
    phi_s = np.deg2rad([10.0, 80.0, 90.0])
    j_is = -93.0


class RateNH_D(RateNH):
    rih3 = 2.50e-10 ** 3
    rsh3 = 2.48e-10 ** 3


class RateHN(RatesIS):
    gi = cnc.GAMMA["h"]
    gs = cnc.GAMMA["n"]
    ris3 = 1.04e-10 ** 3
    rih3 = 1.85e-10 ** 3
    rsh3 = 1.79e-10 ** 3
    csa_i = np.array([5.7, 0.5, -6.2]) * 1e-6
    csa_s = np.array([69.0, 42.0, -111.0]) * 1e-6
    phi_i = np.deg2rad([10.0, 80.0, 90.0])
    phi_s = np.deg2rad([109.6, 90.0, 19.6])
    j_is = -93.0


class RateHN_D(RateHN):
    rih3 = 2.48e-10 ** 3
    rsh3 = 2.50e-10 ** 3


class RateCH(RatesIS):
    gi = cnc.GAMMA["c"]
    gs = cnc.GAMMA["h"]
    ris3 = 1.09e-10 ** 3
    rih3 = 1.70e-10 ** 3
    rsh3 = 1.85e-10 ** 3
    csa_i = np.array([-5.3, -5.3, 10.7]) / 3.0 * 1e-6
    csa_s = np.array([0.0, 0.0, 0.0]) * 1e-6
    phi_i = np.deg2rad([109.6, 90.0, 19.6])
    phi_s = np.deg2rad([0.0, 0.0, 0.0])
    j_is = 140.0


class RateCH_D(RateCH):
    rih3 = 2.03e-10 ** 3
    rsh3 = 2.52e-10 ** 3


rate_functions = {
    "nh": ft.lru_cache(32)(RateNH()),
    "nh_d": ft.lru_cache(32)(RateNH_D()),
    "hn": ft.lru_cache(32)(RateHN()),
    "hn_d": ft.lru_cache(32)(RateHN_D()),
}

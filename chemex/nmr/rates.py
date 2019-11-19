import dataclasses as dc

import numpy as np
import scipy.constants as sc

import chemex.nmr.constants as cnc


is_params = {
    "hn": {
        "gamma_i": cnc.GAMMA["h"],
        "gamma_s": cnc.GAMMA["n"],
        "ris3": 1.02e-10 ** 3,
        "rext3_p": 1.85e-10 ** 3,
        "rext3_d": 2.48e-10 ** 3,
        "csa_i": 10e-6,
        "csa_s": -164e-6,
        "phi_i": np.deg2rad(8.0),
        "phi_s": np.deg2rad(17.0),
        "j_is": -93.0,
    },
    "nh": {
        "gamma_i": cnc.GAMMA["n"],
        "gamma_s": cnc.GAMMA["h"],
        "ris3": 1.02e-10 ** 3,
        "rext3_p": 1.85e-10 ** 3,
        "rext3_d": 2.48e-10 ** 3,
        "csa_i": -164e-6,
        "csa_s": 10e-6,
        "phi_i": np.deg2rad(17.0),
        "phi_s": np.deg2rad(8.0),
        "j_is": -93.0,
    },
    "ch": {
        "gamma_i": cnc.GAMMA["c"],
        "gamma_s": cnc.GAMMA["h"],
        "ris3": 1.09e-10 ** 3,
        "rext3_p": 1.85e-10 ** 3,
        "rext3_d": 2.48e-10 ** 3,
        "csa_i": 16e-6,
        "csa_s": 0e-6,
        "phi_i": np.deg2rad(22.0),
        "phi_s": np.deg2rad(0.0),
        "j_is": 140.0,
    },
}


@dc.dataclass
class ModelFree:
    tauc: float = 5.0e-9
    taui: float = 50.0e-12
    s2: float = 0.9
    deuterated: bool = False

    def jw(self, w):
        tauc = self.tauc
        taui = self.taui
        taup = taui + tauc
        arg1 = self.s2 / (1.0 + (w * tauc) ** 2)
        arg2 = (1.0 - self.s2) * taup * taui / (taup ** 2 + (w * taui * tauc) ** 2)
        jw = 2.0 / 5.0 * tauc * (arg1 + arg2)
        return jw

    def jw_isx(self, w0_i, w0_s, w0_x):
        ws = np.array(
            [
                0.0,
                w0_i,
                w0_s,
                w0_x,
                w0_i - w0_s,
                w0_i + w0_s,
                w0_i - w0_x,
                w0_i + w0_x,
                w0_s - w0_x,
                w0_s + w0_x,
            ]
        )
        return self.jw(ws)


@dc.dataclass
class RatesIS:
    h_larmor_frq: float
    spins: str

    def __post_init__(self):
        is_param = is_params[self.spins]
        gi = is_param["gamma_i"]
        gs = is_param["gamma_s"]
        gh = cnc.GAMMA["h"]
        b0 = 2.0 * np.pi * 1e6 * self.h_larmor_frq / gh
        self.wi = -gi * b0
        self.ws = -gs * b0
        self.wh = -gh * b0
        self.dis = -sc.mu_0 * sc.hbar / (8.0 * np.pi) * gi * gs
        self.dih = -sc.mu_0 * sc.hbar / (8.0 * np.pi) * gi * gh
        self.dsh = -sc.mu_0 * sc.hbar / (8.0 * np.pi) * gs * gh
        self.ci = 1.0 / 3.0 * is_param["csa_i"] * gi * b0
        self.cs = 1.0 / 3.0 * is_param["csa_s"] * gs * b0
        self.p2i = 0.5 * (3.0 * np.cos(is_param["phi_i"]) ** 2 - 1.0)
        self.p2s = 0.5 * (3.0 * np.cos(is_param["phi_s"]) ** 2 - 1.0)

    def calculate(self, model_free):
        is_param = is_params[self.spins]
        ris3 = is_param["ris3"]
        rih3 = is_param["rext3_d"] if model_free.deuterated else is_param["rext3_p"]
        rsh3 = rih3
        dis = self.dis / ris3
        dih = self.dih / rih3
        dsh = self.dsh / rsh3
        ci = self.ci
        cs = self.cs
        p2i = self.p2i
        p2s = self.p2s
        js = model_free.jw_isx(self.wi, self.ws, self.wh)
        j0, ji, js, jh, jmis, jpis, jmih, jpih, jmsh, jpsh = js
        rates = {
            "r2_i": 0.5
            * (
                ci ** 2 * (4 * j0 + 3 * ji)
                + dis ** 2 * (4 * j0 + 3 * ji + jmis + 6 * js + 6 * jpis)
                + dih ** 2 * (4 * j0 + 3 * ji + jmih + 6 * jh + 6 * jpih)
            ),
            "r2_s": 0.5
            * (
                cs ** 2 * (4 * j0 + 3 * js)
                + dis ** 2 * (4 * j0 + 3 * js + jmis + 6 * ji + 6 * jpis)
                + dsh ** 2 * (4 * j0 + 3 * js + jmsh + 6 * jh + 6 * jpsh)
            ),
            "r1_i": (
                ci ** 2 * 3 * ji
                + dis ** 2 * (3 * ji + jmis + 6 * jpis)
                + dih ** 2 * (3 * ji + jmih + 6 * jpih)
            ),
            "r1_s": (
                cs ** 2 * 3 * js
                + dis ** 2 * (3 * js + jmis + 6 * jpis)
                + dsh ** 2 * (3 * js + jmsh + 6 * jpsh)
            ),
            "r2a_i": (
                0.5 * ci ** 2 * (4 * j0 + 3 * ji)
                + cs ** 2 * 3 * js
                + 0.5 * dis ** 2 * (4 * j0 + 3 * ji + jmis + 6 * jpis)
                + 0.5 * dih ** 2 * (4 * j0 + 3 * ji + jmih + 6 * jh + 6 * jpih)
                + dsh ** 2 * (3 * js + jmsh + 6 * jpsh)
            ),
            "r2a_s": (
                0.5 * cs ** 2 * (4 * j0 + 3 * js)
                + ci ** 2 * 3 * ji
                + 0.5 * dis ** 2 * (4 * j0 + 3 * js + jmis + 6 * jpis)
                + 0.5 * dsh ** 2 * (4 * j0 + 3 * js + jmsh + 6 * jh + 6 * jpsh)
                + dih ** 2 * (3 * ji + jmih + 6 * jpih)
            ),
            "r2mq_is": 0.5
            * (
                ci ** 2 * (4 * j0 + 3 * ji)
                + cs ** 2 * (4 * j0 + 3 * js)
                + dis ** 2 * (3 * ji + jmis + 3 * js + 6 * jpis)
                + dih ** 2 * (4 * j0 + 3 * ji + jmih + 6 * jh + 6 * jpih)
                + dsh ** 2 * (4 * j0 + 3 * js + jmsh + 6 * jh + 6 * jpsh)
            ),
            "r1a_is": (
                ci ** 2 * 3 * ji
                + cs ** 2 * 3 * js
                + dis ** 2 * (3 * ji + 3 * js)
                + dih ** 2 * (3 * ji + jmih + 6 * jpih)
                + dsh ** 2 * (3 * js + jmsh + 6 * jpsh)
            ),
            "etaxy_i": dis * ci * p2i * (4 * j0 + 3 * ji),
            "etaxy_s": dis * cs * p2s * (4 * j0 + 3 * js),
            "etaz_i": dis * ci * p2i * 6 * ji,
            "etaz_s": dis * cs * p2s * 6 * js,
            "sigma_is": dis ** 2 * (-jmis + 6.0 * jpis),
            "mu_is": 0.5 * dis ** 2 * (-jmis + 6.0 * jpis),
            "j_is": is_param["j_is"],
        }
        return rates

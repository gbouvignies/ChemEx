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
        "csa_i": -10e-6,
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
        "csa_s": -10e-6,
        "phi_i": np.deg2rad(17.0),
        "phi_s": np.deg2rad(8.0),
        "j_is": -93.0,
    },
}


@dc.dataclass
class ModelFree:
    tauc: float = 5.0e-9
    taui: float = 50.0e-12
    s2: float = 0.9
    deuterated: bool = False

    @property
    def tau_(self):
        return self.tauc * self.taui / (self.tauc + self.taui)

    def jw(self, w):
        tauc = self.tauc
        tau_ = self.tau_
        lc = tauc / (1.0 + (w * tauc) ** 2)
        l_ = tau_ / (1.0 + (w * tau_) ** 2)
        jw = 2.0 / 5.0 * (lc * self.s2 + l_ * (1.0 - self.s2))
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
        gamma_i = is_param["gamma_i"]
        gamma_s = is_param["gamma_s"]
        gamma_h = cnc.GAMMA["h"]
        b0 = 2.0 * np.pi * self.h_larmor_frq * 1e6 / gamma_h
        self.w_i = b0 * gamma_i
        self.w_s = b0 * gamma_s
        self.w_h = b0 * gamma_h
        ris3 = is_param["ris3"]
        scale = sc.mu_0 * sc.hbar / (8.0 * np.pi)
        self.ad_is = scale * gamma_i * gamma_s / ris3
        self.ad_ih = scale * gamma_i * gamma_h
        self.ad_sh = scale * gamma_s * gamma_h
        self.ac_i = -1.0 / 3.0 * is_param["csa_i"] * self.w_i
        self.ac_s = -1.0 / 3.0 * is_param["csa_s"] * self.w_s
        self.geo_i = 0.5 * (3.0 * np.cos(is_param["phi_i"]) ** 2 - 1.0)
        self.geo_s = 0.5 * (3.0 * np.cos(is_param["phi_s"]) ** 2 - 1.0)

    def calculate(self, model_free):
        is_param = is_params[self.spins]
        rext3 = is_param["rext3_d"] if model_free.deuterated else is_param["rext3_p"]
        ad_is = self.ad_is
        ad_ih = self.ad_ih / rext3
        ad_sh = self.ad_sh / rext3
        ac_i = self.ac_i
        ac_s = self.ac_s
        geo_i = self.geo_i
        geo_s = self.geo_s
        jws = model_free.jw_isx(self.w_i, self.w_s, self.w_h)
        j0, jwi, jws, jwh, jwm_is, jwp_is, jwm_ih, jwp_ih, jwm_sh, jwp_sh = jws
        rates = {
            "r2_i": (
                0.5 * ac_i ** 2 * (4 * j0 + 3 * jwi)
                + 0.5 * ad_is ** 2 * (4 * j0 + 3 * jwi + jwm_is + 6 * jws + 6 * jwp_is)
                + 0.5 * ad_ih ** 2 * (4 * j0 + 3 * jwi + jwm_ih + 6 * jwh + 6 * jwp_ih)
            ),
            "r2_s": (
                0.5 * ac_s ** 2 * (4 * j0 + 3 * jws)
                + 0.5 * ad_is ** 2 * (4 * j0 + 3 * jws + jwm_is + 6 * jwi + 6 * jwp_is)
                + 0.5 * ad_sh ** 2 * (4 * j0 + 3 * jws + jwm_sh + 6 * jwh + 6 * jwp_sh)
            ),
            "r1_i": (
                ac_i ** 2 * 3 * jwi
                + ad_is ** 2 * (3 * jwi + jwm_is + 6 * jwp_is)
                + ad_ih ** 2 * (3 * jwi + jwm_ih + 6 * jwp_ih)
            ),
            "r1_s": (
                ac_s ** 2 * 3 * jws
                + ad_is ** 2 * (3 * jws + jwm_is + 6 * jwp_is)
                + ad_sh ** 2 * (3 * jws + jwm_sh + 6 * jwp_sh)
            ),
            "r2a_i": (
                0.5 * ac_i ** 2 * (4 * j0 + 3 * jwi)
                + ac_s ** 2 * 3 * jws
                + 0.5 * ad_is ** 2 * (4 * j0 + 3 * jwi + jwm_is + 6 * jwp_is)
                + 0.5 * ad_ih ** 2 * (4 * j0 + 3 * jwi + jwm_ih + 6 * jwh + 6 * jwp_ih)
                + ad_sh ** 2 * (3 * jws + jwm_sh + 6 * jwp_sh)
            ),
            "r2a_s": (
                0.5 * ac_s ** 2 * (4 * j0 + 3 * jws)
                + ac_i ** 2 * 3 * jwi
                + 0.5 * ad_is ** 2 * (4 * j0 + 3 * jws + jwm_is + 6 * jwp_is)
                + 0.5 * ad_sh ** 2 * (4 * j0 + 3 * jws + jwm_sh + 6 * jwh + 6 * jwp_sh)
                + ad_ih ** 2 * (3 * jwi + jwm_ih + 6 * jwp_ih)
            ),
            "r2mq_is": (
                0.5 * ac_i ** 2 * (4 * j0 + 3 * jwi)
                + 0.5 * ac_s ** 2 * (4 * j0 + 3 * jws)
                + 0.5 * ad_is ** 2 * (3 * jwi + jwm_is + 3 * jws + 6 * jwp_is)
                + 0.5 * ad_ih ** 2 * (4 * j0 + 3 * jwi + jwm_ih + 6 * jwh + 6 * jwp_ih)
                + 0.5 * ad_sh ** 2 * (4 * j0 + 3 * jws + jwm_sh + 6 * jwh + 6 * jwp_sh)
            ),
            "r1a_is": (
                +ac_i ** 2 * 3 * jwi
                + ac_s ** 2 * 3 * jws
                + ad_is ** 2 * (3 * jwi + 3 * jws)
                + ad_ih ** 2 * (3 * jwi + jwm_ih + 6 * jwp_ih)
                + ad_sh ** 2 * (3 * jws + jwm_sh + 6 * jwp_sh)
            ),
            "etaxy_i": ad_is * ac_i * geo_i * (4 * j0 + 3 * jwi),
            "etaxy_s": ad_is * ac_s * geo_s * (4 * j0 + 3 * jws),
            "etaz_i": ad_is * ac_i * geo_i * 6 * jwi,
            "etaz_s": ad_is * ac_s * geo_s * 6 * jws,
            "sigma_is": ad_is ** 2 * (-jwm_is + 6.0 * jwp_is),
            "mu_is": 0.5 * ad_is ** 2 * (-jwm_is + 6.0 * jwp_is),
            "j_is": is_param["j_is"],
        }
        return rates

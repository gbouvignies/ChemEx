import numpy as np
import scipy.constants as sc

import chemex.nmr.constants as cnc


_PARAMS_IS = {
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


def _calculate_jw(tauc, s2, w):
    return 2.0 / 5.0 * tauc * s2 / (1.0 + (w * tauc) ** 2)


def calculate_rates(spin_system, h_frq, tauc, s2, kh2o, deuterated=False):
    ss_params = _PARAMS_IS[spin_system]

    gi, gs, gh = ss_params["gamma_i"], ss_params["gamma_s"], cnc.GAMMA["h"]
    ris3 = ss_params["ris3"]
    rih3 = rsh3 = ss_params["rext3_d"] if deuterated else ss_params["rext3_p"]
    csa_i, csa_s = ss_params["csa_i"], ss_params["csa_s"]
    phi_i, phi_s = ss_params["phi_i"], ss_params["phi_s"]

    b0 = 2.0 * np.pi * 1e6 * h_frq / gh
    wi, ws, wh = -np.array([gi, gs, gh]) * b0
    scale = -sc.mu_0 * sc.hbar / (8.0 * np.pi)
    dis, dih, dsh = scale * np.array([gi * gs / ris3, gi * gh / rih3, gs * gh / rsh3])
    ci, cs = -np.array([csa_i * wi, csa_s * ws]) / 3.0
    p2i, p2s = 0.5 * (3.0 * np.cos(np.array([phi_i, phi_s])) ** 2 - 1.0)

    ws = np.array(
        [0.0, wi, ws, wh, wi - ws, wi + ws, wi - wh, wi + wh, ws - wh, ws + wh]
    )
    jw = _calculate_jw(tauc, s2, ws)
    j0, ji, js, jh, jmis, jpis, jmih, jpih, jmsh, jpsh = jw

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
        "j_is": ss_params["j_is"],
    }
    return rates

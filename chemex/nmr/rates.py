from __future__ import annotations

from itertools import product
from typing import TypeVar

import numpy as np
from scipy.constants import hbar, mu_0

from chemex.configuration.conditions import Conditions
from chemex.models.model import model
from chemex.nmr.basis import Basis
from chemex.nmr.constants import GAMMA
from chemex.typing import ArrayFloat

# Type definition
T = TypeVar("T", float, ArrayFloat)


def _calculate_jw(tauc: float, s2: float, w: T) -> T:
    """Calculates J(w) for given tau_c, s2, and angular frequency w.

    Args:
        tauc (float): Correlation time in seconds.
        s2 (float): Order parameter squared, amplitude of motion measure.
        w (T): Angular frequency/frequencies for J(w) calculation. Can be float or
               ArrayFloat.

    Returns:
        T: Spectral density function(s) J(w).
    """
    tauc_ = tauc * 1e-9
    return 2.0 / 5.0 * tauc_ * s2 / (1.0 + (w * tauc_) ** 2)


class RatesIS:
    """Base class for calculating relaxation rates in IS spin systems.

    Attributes describe gyromagnetic ratios, distances, CSA parameters, and Euler angles
    for nuclei I and S.
    """

    gi: float
    gs: float
    ris3: float
    rih3: float
    rsh3: float
    csa_i: ArrayFloat
    csa_s: ArrayFloat
    phi_i: ArrayFloat
    phi_s: ArrayFloat

    def __init__(self) -> None:
        """Initializes RatesIS object with default gyromagnetic ratios and distances."""
        self.gh = GAMMA["h"]

        # Dipolar factors
        self.dis = -mu_0 * hbar * self.gi * self.gs / (8.0 * np.pi * self.ris3)
        dih = -mu_0 * hbar * self.gi * self.gh / (8.0 * np.pi * self.rih3)
        dsh = -mu_0 * hbar * self.gs * self.gh / (8.0 * np.pi * self.rsh3)
        self.dis2 = self.dis * self.dis
        self.dih2 = dih * dih
        self.dsh2 = dsh * dsh

        # CSA factors
        self.delta_i = self.csa_i[:2] - self.csa_i[2]
        self.delta_s = self.csa_s[:2] - self.csa_s[2]

    def __call__(self, h_frq: float, tauc: float, s2: float) -> dict[str, float]:
        """Calculates relaxation rates for given Larmor frequency, tau_c, and s2.

        Args:
            h_frq (float): Proton Larmor frequency in MHz.
            tauc (float): Correlation time in seconds.
            s2 (float): Order parameter squared.

        Returns:
            dict[str, float]: Dictionary of relaxation rates (R2, R1 for I and S).
        """
        # B0 in Tesla
        b0 = 2.0 * np.pi * 1e6 * h_frq / self.gh

        # Larmor angular frequencies
        wi, ws, wh = -np.array([self.gi, self.gs, self.gh]) * b0

        # CSA factors
        ci_ = -wi * self.delta_i / 3.0
        cs_ = -ws * self.delta_s / 3.0
        ci2 = ci_[0] ** 2 + ci_[1] ** 2 - ci_[0] * ci_[1]
        cs2 = cs_[0] ** 2 + cs_[1] ** 2 - cs_[0] * cs_[1]
        p2i = 0.5 * (3.0 * np.cos(self.phi_i[:2]) ** 2 - 1.0)
        p2s = 0.5 * (3.0 * np.cos(self.phi_s[:2]) ** 2 - 1.0)
        ci = (ci_ * p2i).sum()
        cs = (cs_ * p2s).sum()

        j0, ji, js, jh, jmis, jpis, jmih, jpih, jmsh, jpsh = _calculate_jw(
            tauc,
            s2,
            np.array(
                [0.0, wi, ws, wh, wi - ws, wi + ws, wi - wh, wi + wh, ws - wh, ws + wh],
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
        }


class RateNH(RatesIS):
    """Class for calculating relaxation rates in NH systems using model-free approach."""

    gi = GAMMA["n"]
    gs = GAMMA["h"]
    ris3 = 1.04e-10**3
    rih3 = 1.79e-10**3
    rsh3 = 1.85e-10**3
    csa_i = np.array([69.0, 42.0, -111.0]) * 1e-6
    csa_s = np.array([5.7, 0.5, -6.2]) * 1e-6
    phi_i = np.deg2rad([109.6, 90.0, 19.6])
    phi_s = np.deg2rad([10.0, 80.0, 90.0])

    def __call__(
        self,
        h_frq: float,
        tauc: float,
        s2: float,
        khh: float = 0.0,
    ) -> dict[str, float]:
        """Calculates rates for NH systems, including exchange contributions if any.

        Args:
            h_frq (float): Proton Larmor frequency in MHz.
            tauc (float): Correlation time in seconds.
            s2 (float): Order parameter squared.
            khh (float, optional): Exchange rate between H atoms. Defaults to 0.

        Returns:
            dict[str, float]: Relaxation rates for NH systems, including exchanges.
        """
        rates = super().__call__(h_frq, tauc, s2)
        if khh == 0:
            return rates
        # Make a copy of rates before adding 'khh' due to lru_cache"
        # on 'super().__call__'
        rates = rates.copy()
        rates["r2_s"] += khh
        rates["r1_s"] += khh
        rates["r2a_i"] += khh
        rates["r2a_s"] += khh
        rates["r2mq_is"] += khh
        rates["r1a_is"] += khh
        return rates


class RateNH_D(RateNH):
    rih3 = 2.50e-10**3
    rsh3 = 2.48e-10**3


class RateHN(RatesIS):
    gi = GAMMA["h"]
    gs = GAMMA["n"]
    ris3 = 1.04e-10**3
    rih3 = 1.85e-10**3
    rsh3 = 1.79e-10**3
    csa_i = np.array([5.7, 0.5, -6.2]) * 1e-6
    csa_s = np.array([69.0, 42.0, -111.0]) * 1e-6
    phi_i = np.deg2rad([10.0, 80.0, 90.0])
    phi_s = np.deg2rad([109.6, 90.0, 19.6])

    def __call__(
        self,
        h_frq: float,
        tauc: float,
        s2: float,
        khh: float = 0.0,
    ) -> dict[str, float]:
        rates = super().__call__(h_frq, tauc, s2)
        if khh == 0:
            return rates
        # Make a copy of rates before adding 'khh' due to lru_cache
        # on 'super().__call__'
        rates = rates.copy()
        rates["r2_i"] += khh
        rates["r1_i"] += khh
        rates["r2a_i"] += khh
        rates["r2a_s"] += khh
        rates["r2mq_is"] += khh
        rates["r1a_is"] += khh
        return rates


class RateHN_D(RateHN):
    rih3 = 2.48e-10**3
    rsh3 = 2.50e-10**3


class RateCH(RatesIS):
    gi = GAMMA["c"]
    gs = GAMMA["h"]
    ris3 = 1.09e-10**3
    rih3 = 1.70e-10**3
    rsh3 = 1.85e-10**3
    csa_i = np.array([-5.3, -5.3, 10.7]) / 3.0 * 1e-6
    csa_s = np.array([0.0, 0.0, 0.0]) * 1e-6
    phi_i = np.deg2rad([109.6, 90.0, 19.6])
    phi_s = np.deg2rad([0.0, 0.0, 0.0])


class RateCH_D(RateCH):
    rih3 = 2.03e-10**3
    rsh3 = 2.52e-10**3


class RateHC(RatesIS):
    gi = GAMMA["h"]
    gs = GAMMA["c"]
    ris3 = 1.09e-10**3
    rih3 = 1.85e-10**3
    rsh3 = 1.70e-10**3
    csa_i = np.array([0.0, 0.0, 0.0]) * 1e-6
    csa_s = np.array([-5.3, -5.3, 10.7]) / 3.0 * 1e-6
    phi_i = np.deg2rad([0.0, 0.0, 0.0])
    phi_s = np.deg2rad([109.6, 90.0, 19.6])


class RateHC_D(RateHC):
    rih3 = 2.52e-10**3
    rsh3 = 2.03e-10**3


rate_functions: dict[str, RatesIS] = {
    "nh": RateNH(),
    "nh_d": RateNH_D(),
    "hn": RateHN(),
    "hn_d": RateHN_D(),
    "ch": RateCH(),
    "ch_d": RateCH_D(),
    "hc": RateHC(),
    "hc_d": RateHC_D(),
    "cn": RateCH(),
    "cn_d": RateCH_D(),
}

_RATE_NAMES = [
    "r2_i",
    "r2_s",
    "r1_i",
    "r1_s",
    "r2a_i",
    "r2a_s",
    "r2mq_is",
    "r1a_is",
    "etaxy_i",
    "etaxy_s",
    "etaz_i",
    "etaz_s",
    "sigma_is",
    "mu_is",
]


def get_model_free_expressions(basis: Basis, conditions: Conditions) -> dict[str, str]:
    """Generates expressions for model-free analysis based on basis and conditions.

    Args:
        basis (Basis): Basis set for model-free analysis.
        conditions (Conditions): Conditions including Larmor frequency and deuteration.

    Returns:
        dict[str, str]: Mapping of rate names to expressions for model-free analysis.
    """
    deuterated_extension = "_d" if conditions.is_deuterated else ""
    rate_function_name = f"{basis.spin_system}{deuterated_extension}"

    h_frq_str = f"{conditions.h_larmor_frq}"
    has_h_exchange = basis.spin_system in {"nh", "hn"}

    model_free_expr: dict[str, str] = {}
    for state, name in product(model.states, _RATE_NAMES):
        rate_name = f"{name}_{state}"
        khh = f", {{khh_{state}}}" if has_h_exchange else ""
        arguments = f"{h_frq_str}, {{tauc_{state}}}, {{s2_{state}}}{khh}"
        expr = f"{rate_function_name}({arguments})['{name}']"
        model_free_expr[rate_name] = expr

    return model_free_expr

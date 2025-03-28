"""Module for NMR spectroscopy constants and multiplet pattern calculations.

This module defines essential constants and functions for Nuclear Magnetic Resonance
(NMR) spectroscopy, focusing on gyromagnetic ratios and scalar coupling constants.
"""

from collections import Counter
from dataclasses import dataclass

import numpy as np

from chemex.parameters.spin_system.nucleus import Nucleus
from chemex.typing import ArrayFloat

GAMMA = {
    Nucleus.H1: 26.752_212_8e07,
    Nucleus.N15: -2.712_618_04e07,
    Nucleus.C13: 6.728_284e07,
    Nucleus.F19: 25.18148e07,
    Nucleus.P31: 10.8394e07,
}

# Define nuclide frequency ratios wrt proton
# IUPAC values for bio NMR: Markley et al. Pure Appl. Chem. (1998) 70, 117.
# 19F value comes from: Harris et al. Pure Appl. Chem. (2001) 73, 1795.
XI_RATIO = {
    Nucleus.H1: 100.000_000_0e-02,
    Nucleus.N15: 10.132_911_8e-02,
    Nucleus.C13: 25.144_953_0e-02,
    Nucleus.F19: 94.094_011e-02,
    Nucleus.P31: 40.480_863_6e-02,
}

SIGNED_XI_RATIO: dict[Nucleus, float] = {
    key: np.sign(GAMMA[key]) * val for key, val in XI_RATIO.items()
}

# Generic scalar couplings used as starting values in IS systems (in Hz)
J_COUPLINGS = {
    "nh": -93.0,
    "hn": -93.0,
    "cn": -10.7,
    "ch": 130.0,
    "hc": 130.0,
}


# Residue-specific scalar coupling values with neighboring carbons (in Hz)
J_EFF: dict[str, dict[str, tuple[float, ...]]] = {
    "a": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0,),
        "c1'": (45.0, 11.0),
        "c2": (11.2,),
        "c8": (9.5, 9.5, 6.0, 10.0),
    },
    "c": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0,),
        "c1'": (45.0, 11.0),
        "c6": (66.0, 13.0, 3.0),
    },
    "d": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 52.5),
        "cg": (52.5,),
    },
    "e": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0),
        "cg": (35.0, 52.5),
        "cd": (52.5,),
    },
    "f": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 52.5),
    },
    "g": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, -10.7, -7.7),
        "c1'": (45.0, 11.0),
        "c8": (9.5, 9.5, 6.0, 10.0),
        "n1": (7.5, 15.0),
    },
    "h": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 52.5),
        "cg": (52.5, 72.0, -14.4),
        "cd2": (72.0, -14.4),
        "ce1": (-14.4, -14.4),
        "nd1": (-14.4, -14.4),
        "ne2": (-14.4, -14.4),
    },
    "i": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0, 35.0),
        "cg1": (35.0, 35.0),
        "cg2": (35.0,),
        "cd1": (35.0,),
    },
    "k": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0),
        "cg": (35.0, 35.0),
        "cd": (35.0, 35.0),
        "ce": (35.0, -10.7),
        "nz": (-10.7,),
    },
    "l": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0),
        "cg": (35.0, 35.0, 35.0),
        "cd1": (35.0,),
        "cd2": (35.0,),
    },
    "m": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0),
        "cg": (35.0,),
        "ce": (),
    },
    "n": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 52.5, -7.7),
        "cg": (52.5, -14.4),
        "nd2": (-14.4, -7.7),
    },
    "p": {
        "n": (-7.7, -10.7, -14.4, -10.7),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0),
        "cg": (35.0, 35.0),
        "cd": (35.0, -10.7),
    },
    "q": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0),
        "cg": (35.0, 52.5, -7.7),
        "cd": (52.5, -14.4),
        "ne2": (-14.4, -7.7),
    },
    "r": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0),
        "cg": (35.0, 35.0),
        "cd": (35.0, -10.7),
        "ne": (-10.7, -14.4),
        "cz": (-14.4, -14.4, -14.4),
        "nh1": (-14.4,),
        "nh2": (-14.4,),
    },
    "s": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0,),
    },
    "t": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0),
        "cg2": (35.0,),
        "c1'": (45.0, 11.0),
    },
    "v": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 35.0, 35.0),
        "cg1": (35.0,),
        "cg2": (35.0,),
    },
    "w": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 52.5),
        "ne1": (-14.4, -14.4),
    },
    "y": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0, 52.5),
    },
    "u": {
        "c1'": (45.0, 11.0),
        "c6": (66.0, 13.0, 3.0),
        "n3": (19.2, 10.7),
    },
    "": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
    },
}


@dataclass
class Distribution:
    values: ArrayFloat
    weights: ArrayFloat


def get_multiplet(symbol: str, nucleus: str) -> Distribution:
    """Calculate the multiplet pattern."""
    multiplet = [0.0]
    for coupling in J_EFF[symbol.lower()][nucleus.lower()]:
        doublet = coupling * 0.5 * np.array([[-1.0], [1.0]], dtype=np.float64)
        multiplet = (multiplet + doublet).flatten()
    counter = Counter(multiplet)
    values, weights = zip(*counter.items(), strict=True)
    return Distribution(
        np.array(values, dtype=np.float64), np.array(weights, dtype=np.float64)
    )

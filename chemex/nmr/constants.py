"""The constants module defines the gyromagnetic ratios and scalar coupling
constants.

Define the gyromagnetic ratios in rad/s/T
IUPAC values: Harris et al, Concepts in Magn. Reson., (2002) 14, p326
"""
import collections

import numpy as np


GAMMA = {
    "h": 26.752_212_8e07,
    "q": 26.752_212_8e07,
    "n": -2.712_618_04e07,
    "c": 6.728_284e07,
    "f": 25.18148e07,
    "p": 10.8394e07,
}

G_RATIO = {key: val / GAMMA["h"] for key, val in GAMMA.items()}

# Define nuclide frequency ratios wrt proton
# IUPAC values for bio NMR: Markley et al, Pure & Appl. Chem., (1998) 70, p117
XI_RATIO = {
    "h": 100.000_000_0e-02,
    "q": 100.000_000_0e-02,
    "n": 10.132_911_8e-02,
    "c": 25.144_953_0e-02,
    "f": 40.480_863_6e-02,
}

SIGNED_XI_RATIO = {key: np.sign(GAMMA[key]) * val for key, val in XI_RATIO.items()}

# Residue-specific scalar coupling values with neighbouring carbons (in Hz)
J_COUPLING = {
    "a": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0,),
    },
    "c": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
        "cb": (35.0,),
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
    "g": {"n": (-7.7, -10.7, -14.4), "c": (52.5, -14.4), "ca": (52.5, -10.7, -7.7)},
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
    "": {
        "n": (-7.7, -10.7, -14.4),
        "c": (52.5, -14.4),
        "ca": (52.5, 35.0, -10.7, -7.7),
    },
}

Distribution = collections.namedtuple("Distribution", ["values", "weights"])


def get_multiplet(symbol, nucleus):
    """Calculate the multiplet pattern."""
    multiplet = np.array([0.0])
    for coupling in J_COUPLING[symbol.lower()][nucleus.lower()]:
        doublet = coupling * 0.5 * np.array([-1.0, 1.0]).reshape(-1, 1)
        multiplet = (multiplet + doublet).reshape(-1)
    counter = collections.Counter(multiplet)
    values, weights = zip(*counter.items())
    return Distribution(np.array(values), np.array(weights))

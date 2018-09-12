"""The constants module defines the gyromagnetic ratios and scalar coupling
constants.

Define the gyromagnetic ratios in rad/s/T
IUPAC values: Harris et al, Concepts in Magn. Reson., (2002) 14, p326
"""

import collections

import numpy as np

gamma = {
    "h": 26.7522128e+07,
    "q": 26.7522128e+07,
    "n": -2.71261804e+07,
    "c": 6.728284e+07,
    "f": 25.18148e+07,
    "p": 10.8394e+07,
}

g_ratio = {key: val / gamma["h"] for key, val in gamma.items()}

# Define nuclide frequency ratios wrt proton
# IUPAC values for bio NMR: Markley et al, Pure & Appl. Chem., (1998) 70, p117
xi_ratio = {
    "h": 100.0000000e-02,
    "q": 100.0000000e-02,
    "n": 10.1329118e-02,
    "c": 25.1449530e-02,
    "f": 40.4808636e-02,
}

# Residue-specific scalar coupling values with neighbouring carbons (in Hz)
j_couplings = {
    "a": {"n": (7.7, 10.7, 14.4), "c": (52.0, 14.4), "ca": (52.0, 35.0), "cb": (35.0,)},
    "c": {"n": (7.7, 10.7, 14.4), "c": (52.0, 14.4), "ca": (52.0, 35.0), "cb": (35.0,)},
    "d": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 51.0),
    },
    "e": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0),
        "cg": (35.0, 51.0),
    },
    "f": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 51.0),
    },
    "g": {"n": (7.7, 10.7, 14.4), "c": (52.0, 14.4), "ca": (52.0,)},
    "h": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 51.0),
        "cg": (51.0, 72.0),
        "cd2": (72.0,),
        "ce1": (),
    },
    "i": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0, 35.0),
        "cg1": (35.0, 35.0),
        "cg2": (35.0,),
        "cd1": (35.0,),
    },
    "k": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0),
        "cg": (35.0, 35.0),
        "cd": (35.0, 35.0),
        "ce": (35.0,),
    },
    "l": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0),
        "cg": (35.0, 35.0, 35.0),
        "cd1": (35.0,),
        "cd2": (35.0,),
    },
    "m": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0),
        "cg": (35.0,),
        "ce": (),
    },
    "n": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 51.0),
    },
    "p": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0),
        "cg": (35.0, 35.0),
        "cd": (35.0,),
    },
    "q": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0),
        "cg": (35.0, 51.0),
    },
    "r": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0),
        "cg": (35.0, 35.0),
        "cd": (35.0,),
    },
    "s": {"n": (7.7, 10.7, 14.4), "c": (52.0, 14.4), "ca": (52.0, 35.0), "cb": (35.0,)},
    "t": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0),
        "cg2": (35.0,),
    },
    "v": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 35.0, 35.0),
        "cg1": (35.0,),
        "cg2": (35.0,),
    },
    "w": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 51.0),
    },
    "y": {
        "n": (7.7, 10.7, 14.4),
        "c": (52.0, 14.4),
        "ca": (52.0, 35.0),
        "cb": (35.0, 51.0),
    },
    "": {"n": (7.7, 10.7, 14.4), "c": (52.0, 14.4)},
}

Distribution = collections.namedtuple("Distribution", ["values", "weights"])


def get_multiplet(symbol, nucleus):
    """Calculate the multiplet pattern."""
    multiplet = np.array([0.0])
    for coupling in j_couplings[symbol][nucleus]:
        doublet = coupling * 0.5 * np.array([-1.0, 1.0]).reshape(-1, 1)
        multiplet = (multiplet + doublet).reshape(-1)
    counter = collections.Counter(multiplet)
    values, weights = zip(*counter.items())
    distribution = Distribution(np.array(values), np.array(weights))
    return distribution

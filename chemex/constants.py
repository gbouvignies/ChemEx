"""
Created on Aug 5, 2011

@author: guillaume
"""

from math import pi

# Define the gyromagnetic ratios in rad/s/T
# IUPAC values: Harris et al, Concepts in Magn. Reson., (2002) 14, p326
gamma = {
    'H': 26.7522128e+07,
    'N': -2.71261804e+07,
    'C': 6.728284e+07,
    'F': 25.18148e+07,
    'P': 10.8394e+07
}

# Define the gyromagnetic ratios in MHz/T
gamma_in_MHz_T = dict((key, val / (2.0 * pi) * 1.0e-06) for key, val in gamma.iteritems())

# Define nuclide frequency ratios wrt proton
# IUPAC values for bio NMR: Markley et al, Pure & Appl. Chem., (1998) 70, p117
xi_ratio = {
    'H': 100.0000000e-02,
    'N': 10.1329118e-02,
    'C': 25.1449530e-02,
    'F': 40.4808636e-02,
}

scalar_couplings = {
    'amide_HN': -92.0,
    'methyl_HC': 125.0,
    'methylene_HC': 130.0,
    'aromatic_HC': 180.0,
}

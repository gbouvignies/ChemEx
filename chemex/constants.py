# Define the gyromagnetic ratios in rad/s/T
# IUPAC values: Harris et al, Concepts in Magn. Reson., (2002) 14, p326

gamma = {
    '1H': 26.7522128e+07,
    '15N': -2.71261804e+07,
    '13C': 6.728284e+07,
    'F': 25.18148e+07,
    'P': 10.8394e+07
}

# Define nuclide frequency ratios wrt proton
# IUPAC values for bio NMR: Markley et al, Pure & Appl. Chem., (1998) 70, p117

xi_ratio = {
    '1H': 100.0000000e-02,
    '15N': 10.1329118e-02,
    '13C': 25.1449530e-02,
    'F': 40.4808636e-02,
}

scalar_couplings = {
    'amide_hn': -92.0,
    'methyl_ch': 125.0,
    'methylene_ch': 130.0,
    'aromatic_ch': 180.0,
}

"""
Created on Feb 29, 2012
Edited on March 12, 2012

@author: Alex Hansen
@author: Mike Latham
@author: Guillaume Bouvignies
"""

parse_line = "15N - Pure In-phase Nitrogen CPMG (3 states)"

description = \
    """    Analyzes 15N chemical exchange in the presence of high power 1H CW
        decoupling during the CPMG block. This keeps the spin system purely in-phase
        throughout, and is calculated using the 9x9, single spin matrix:

        [ Nx(a), Ny(a), Nz(a), Nx(b), Ny(b), Nz(b), Nx(c), Ny(c), Nz(c) ]

        Off resonance effects are taken into account. The calculation is designed
        specifically to analyze the experiment found in the reference."""

reference = {
    'journal': '',
    'year': 0,
    'volume': 0,
    'pages': ''
}

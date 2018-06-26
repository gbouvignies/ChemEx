"""
Created on May 20, 2017

@author: Tairan Yuwen
"""

parse_line = "15N - Pure In-phase Nitrogen CPMG"

description = """\
Analyzes 15N chemical exchange in the presence of high power 1H CW decoupling 
during the CPMG block. This keeps the spin system purely in-phase throughout, 
and is calculated using the 6x6, single spin matrix:

[ Nx(a), Ny(a), Nz(a), Nx(b), Ny(b), Nz(b) ]

This version is modified such that CPMG pulses are applied with [0013] phase cycle
as shown in Daiwen's paper. The cw decoupling on 1H is assumed to be 
strong enough (> 15 kHz) such that perfect 1H decoupling can be achieved.

Off resonance effects are taken into account. The calculation is designed
specifically to analyze the experiment found in the reference:"""

reference = {
    'journal': 'Journal of Magnetic Resonance',
    'year': 2015,
    'volume': 257,
    'pages': '1-7'
}

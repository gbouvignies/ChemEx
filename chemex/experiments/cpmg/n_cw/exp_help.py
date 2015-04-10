"""
Created on Feb 29, 2012
Edited on March 12, 2012

@author: Alex Hansen
@author: Mike Latham
@author: Guillaume Bouvignies
"""

parse_line = "15N - Pure In-phase Nitrogen CPMG"

description = """\
Analyzes 15N chemical exchange in the presence of high power 1H CW decoupling 
during the CPMG block. This keeps the spin system purely in-phase throughout, 
and is calculated using the 6x6, single spin matrix:

[ Nx(a), Ny(a), Nz(a), Nx(b), Ny(b), Nz(b) ]

Off resonance effects are taken into account. The calculation is designed
specifically to analyze the experiment found in the reference:"""

reference = {
    'journal': 'Journal of Physical Chemistry B',
    'year': 2008,
    'volume': 112,
    'pages': '5898-5904'
}

"""
Created on Mar 2, 2012

@author: Alex Hansen
"""

parse_line = "15N - Pure In-phase Nitrogen CEST"

description = """\
Analyzes 15N chemical exchange in the presence of 1H composite decoupling
during the CEST block. This keeps the spin system purely in-phase throughout,
and is calculated using the 6x6, single spin matrix:

[ Nx(a), Ny(a), Nz(a), Nx(b), Ny(b), Nz(b) ]

The calculation is designed specifically to analyze the experiment found in
the reference:"""

reference = {
    'journal': 'J Am Chem Soc',
    'year': 2012,
    'volume': 134,
    'pages': '8148-61'
}

"""
Created on Mar 2, 2012

@author: Alex Hansen
"""

parse_line = "13C - Pure In-phase Carbon CEST"

description = \
    """    Analyzes 13C chemical exchange in the presence of 1H composite
        decoupling during the CEST block. This keeps the spin system purely in-phase
        throughout, and is calculated using the 6x6, single spin matrix:

        [ Cx(a), Cy(a), Cz(a), Cx(b), Cy(b), Cz(b) ]

        The calculation is designed specifically to analyze the experiment
        found in the reference."""

reference = {'journal': 'J Biomol NMR',
             'year': 2012,
             'volume': 53,
             'pages': '303-10'
}

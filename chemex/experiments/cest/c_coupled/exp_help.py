"""
Created on May 1, 2013

@author: Guillaume Bouvignies
"""

parse_line = "13C - Coupled Aliphatic Carbon CEST"

description = \
    """    Analyzes aliphatic 13C chemical exchange in a uniformly 13C-labeled sample
        in the presence of 1H composite decoupling during the CEST block. This keeps
        the spin system only couple to the neighbouring carbons, and is calculated
        using the 6x6, single spin matrix for each component of the multiplet:

        [ Cx(a), Cy(a), Cz(a), Cx(b), Cy(b), Cz(b) ]

        The calculation is designed specifically to analyze the experiment
        found in the reference."""

reference = {'journal': 'In preparation',
             'year': 2013,
             'volume': 0,
             'pages': ''
}

"""
Created on May 1, 2013

@author: Guillaume Bouvignies
"""

parse_line = "15N - Coupled Nitrogen CEST"

description = \
    """    Analyzes 15N chemical exchange in a uniformly 13C-labeled sample
        in the presence of 1H composite decoupling during the CEST block.
        This keeps the spin system only coupled to the neighbouring carbons,
        and is calculated using the 6x6, single spin matrix for each component
        of the multiplet:

        [ Nx(a), Ny(a), Nz(a), Nx(b), Ny(b), Nz(b) ]

        The calculation is designed specifically to analyze the experiment
        found in the reference."""

reference = {'journal': 'In preparation',
             'year': 2013,
             'volume': 0,
             'pages': ''
}

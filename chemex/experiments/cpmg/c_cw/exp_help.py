"""
Created on Feb 29, 2012
Edited on March 12, 2012

@author: Alex Hansen
"""

parse_line = "13C - Pure In-phase Carbon CPMG"

description = \
    """    Analyzes C13 chemical exchange in the presence of high power 1H CW
        decoupling during the CPMG block. This keeps the spin system purely in-phase
        throughout, and is calculated using the 6x6, single spin matrix:

        [ Cx(a), Cy(a), Cz(a), Cx(b), Cy(b), Cz(b) ]

        Off resonance effects are taken into account and is still quite fast to
        calculated. The calculation is designed specifically to analyze the experiment
        found in the reference."""

reference = {
    'journal': 'Journal of the American Chemical Society',
    'year': 2008,
    'volume': 130,
    'pages': '2667-2675'
}

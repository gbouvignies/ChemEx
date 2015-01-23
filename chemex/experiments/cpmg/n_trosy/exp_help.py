# coding=utf-8
"""
Created on Oct 10, 2013

@author: Guillaume Bouvignies
"""

parse_line = "15N - N-H TROSY CPMG"

description = """\
Analyzes 15N constant-time TROSY CPMG relaxation dispersion experiments for
measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the 12x12, two spin matrix:

[ Nx(a), Ny(a), Nz(a), HzNx(a), HzNy(a), HzNz(a),
  Nx(b), Ny(b), Nz(b), HzNx(b), HzNy(b), HzNz(b) ]

Off resonance effects are taken into account. The calculation is designed
specifically to analyze the experiment found in the reference:"""

reference = {
    'journal': 'Journal',
    'year': 2007,
    'volume': 0,
    'pages': '1-10'
}

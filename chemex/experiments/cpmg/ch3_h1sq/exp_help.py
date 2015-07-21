"""
Created on July 21, 2015

@author: Alex Hansen
"""

# local import

parse_line = "1H(methyl - 13CH3) - Single-Quantum Proton CPMG "

description = """\
Measures methyl proton chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization is 
initally anti-phase and is read out as in-phase prior to 1H detection. The 
calculation uses a 13x13 basis set:

[ E/2, Hx(a), Hy(a), Hz(a), 2HxCz(a), 2HyCz(a), 2HzCz(a),
       Hx(b), Hy(b), Hz(b), 2HxCz(b), 2HyCz(b), 2HzCz(b)]

Off resonance effects are taken into account. The calculation is designed 
explicitly for analyzing the pulse sequence written by Alex Hansen:

ch3_1hsqcpmg_vXX.alh"""

reference = {
    'journal': 'Journal of the American Chemical Society',
    'year': 2007,
    'volume': 129,
    'pages': '9514-9521'
}

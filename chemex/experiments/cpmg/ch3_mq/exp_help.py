'''
Created on Feb 29, 2012
Edited on March 12, 2012

@author: Alex Hansen
@author: Mike Latham
@author: Guillaume Bouvignies
'''

# local import

parse_line = "1H-13C(methyl) - Multiple Quantum CPMG"

description = \
    ''' Analyzes HyCx methyl group multiple quantum CPMG measured on site-specific
        13CH3-labeled methyl groups in a highly deuterated background.  This is a
        simplified basis set, which assumes you are on-resonance for 13C (ie, off-
        resonance effects are not taken into account) as described in the reference:

        [HxCx(a), HyCx(a), HxCy(a), HyCy(b), HxCx(b), HyCx(b), HxCy(b), HyCy(b)]

        This calculation is designed specifically to analyze data from the
        experiment found in the reference and can be run with either
        small_protein_flag='y' or 'n'.

        Lewis Kay experiment: hmqc_CH3_exchange_bigprotein_*00_lek_v2
    '''

reference = {'journal': 'Journal of the American Chemical Society',
             'year': 2004,
             'volume': 126,
             'pages': '3964-3973'
}

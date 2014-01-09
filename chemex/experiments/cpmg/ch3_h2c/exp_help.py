'''
Created on Mar 14, 2012

@author: Mike Latham
'''

# local import

parse_line = "13C(methyl) - H to C CPMG "

description = \
    ''' Measures methyl carbon chemical exchange recorded on site-specifically
        13CH3-labeled proteins in a highly deuterated background.  Magnetization
        is initally anti-phase and is read out as in-phase.  Because of the
        P-element only even ncyc should be recorded.  The calculation uses a 13x13
        basis set:

        [ E/2, Hx(a), Hy(a), Hz(a), 2HxCz(a), 2HyCz(a), 2HzCz(a),
               Hx(b), Hy(b), Hz(b), 2HxCz(b), 2HyCz(b), 2HzCz(b)]

        Off resonance effects are taken into account. The calculation is designed
        explicitly for analyzing the Lewis Kay pulse sequence:

        HtoC_CH3_exchange_*00_lek_ILV
    '''

reference = {'journal': 'Journal of Biomolecular NMR',
             'year': 2007,
             'volume': 38,
             'pages': '79-88'
}
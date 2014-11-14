'''
Created on Nov 11, 2014

@author: Mike Latham
'''

parse_line = "13CO - Pure Anti-phase Carbonyl 13C CPMG"

description = \
    '''    Analyzes carbonyl chemical exchange that is maintained as anti-phase
        magnetization throughout the CPMG block. This results in lower intrinsic
        relaxation rates and therefore better sensitivity. The calculations use a
        12x12, 2-spin exchange matrix:

        [ COx(a), COy(a), COz(a), 2COxNz(a), 2COyNz(a), 2COzNz(a),
          COx(b), COy(b), COz(b), 2COxNz(b), 2COyNz(b), 2COzNz(b)]

        Because of the length of the shaped pulses used during the CPMG blocks,
        off resonance effects are taken into account only for the 90-degree pulses
        that create COxNz before the CPMG and COzNz after the CPMG.
        The calculation is designed explicitly for analyzing the Kay
        laboratory pulse sequence:

        CO_CPMG_SCFilter_x00_dfh1

        For uniformly 13C-labeled proteins.
    '''

reference = {'journal': 'Journal of Biomolecular NMR',
             'year': 2008,
             'volume': 42,
             'pages': '35-47'
}

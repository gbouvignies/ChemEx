'''
Created on Feb 12, 2014 

@author: mike latham
'''

# local import

parse_line = "1H(methyl - 13CHD2) - Pure Anti-Phase Proton CPMG "

description = \
    ''' Measures methyl proton chemical exchange recorded on site-specifically
        13CHD2-labeled proteins in a highly deuterated background.  
        Magnetization is initally anti-phase and is read out as anti-phase
        prior to 13C evolution. The calculation uses a 13x13 basis set:

        [ E/2, Hx(a), Hy(a), Hz(a), 2HxCz(a), 2HyCz(a), 2HzCz(a),
               Hx(b), Hy(b), Hz(b), 2HxCz(b), 2HyCz(b), 2HzCz(b)]

        Off resonance effects are taken into account. The calculation is 
        designed explicitly for analyzing the Lewis Kay pulse sequence:

        CHD2_H_SQ_exchange_hsqc_lek_*00_enhanced
    '''

reference = {'journal': 'Journal of the American Chemical Society',
             'year': 2010,
             'volume': 132,
             'pages': '10992-10995'
}

"""
Created on June 26, 2014

@author: Mike Latham
"""

parse_line = "1H - Pure Anti-phase Proton CPMG"

description = """\
Analyzes amide proton chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a 13x13,
2-spin exchange matrix:
    
[ E/2, Hx(a), Hy(a), Hz(a), 2HxNz(a), 2HyNz(a), 2HzNz(a),
       Hx(b), Hy(b), Hz(b), 2HxNz(b), 2HyNz(b), 2HzNz(b)]
    
Off resonance effects are taken into account. The calculation is designed
explicitly for analyzing the Lewis Kay pulse sequence:
    
H1_CPMG_Rex_hsqc_lek_x00
    
with antiphase_flg set to 'y'"""

reference = {'journal': 'Journal of Biomolecular NMR',
             'year': 2011,
             'volume': 50,
             'pages': '13-18'
            }

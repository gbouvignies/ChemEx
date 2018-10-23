"""1H(methyl - 13CHD2) - Pure Anti-Phase Proton CPMG

Measures methyl proton chemical exchange recorded on site-specifically
13CHD2-labeled proteins in a highly deuterated background. Magnetization is
initally anti-phase and is read out as anti-phase prior to 13C evolution. The
calculation uses a 12x12 basis set:

[ Hx(a), Hy(a), Hz(a), 2HxCz(a), 2HyCz(a), 2HzCz(a),
  Hx(b), Hy(b), Hz(b), 2HxCz(b), 2HyCz(b), 2HzCz(b)]

Note
----

Off resonance effects are taken into account. The calculation is designed
explicitly for analyzing the Lewis Kay pulse sequence:

CHD2_H_SQ_exchange_hsqc_lek_*00_enhanced

Reference
---------

Journal of the American Chemical Society (2010) 132, 10992-5
"""
from chemex.experiments.cpmg import hn_ap


ProfileCPMGCHD2H1SQ = hn_ap.ProfileCPMGHNAP

"""ChemEx is an analysis program for chemical exchange detected by NMR.

It is designed to take almost any kind of NMR data to aid the analysis,
but the principle techniques are CPMG relaxation dispersion and Chemical
Exchange Saturation Transfer.

"""
import importlib.metadata as importlib_metadata

__version__ = importlib_metadata.version(__name__)

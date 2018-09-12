"""ChemEx is an analysis program for chemical exchange detected by NMR.

It is designed to take almost any kind of NMR data to aid the analysis,
but the principle techniques are CPMG relaxation dispersion and Chemical
Exchange Saturation Transfer.

"""

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

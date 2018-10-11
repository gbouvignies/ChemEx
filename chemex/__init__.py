"""ChemEx is an analysis program for chemical exchange detected by NMR.

It is designed to take almost any kind of NMR data to aid the analysis,
but the principle techniques are CPMG relaxation dispersion and Chemical
Exchange Saturation Transfer.

"""
from pkg_resources import DistributionNotFound
from pkg_resources import get_distribution

__all__ = ["__version__"]

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

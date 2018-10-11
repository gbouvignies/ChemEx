"""ChemEx is an analysis program for chemical exchange detected by NMR.

It is designed to take almost any kind of NMR data to aid the analysis,
but the principle techniques are CPMG relaxation dispersion and Chemical
Exchange Saturation Transfer.

"""

__all__ = ["__version__"]

try:
    from chemex._version import version as __version__
except ImportError:
    # broken installation, we don't even try
    # unknown only works because we do poor mans version compare
    __version__ = "unknown"

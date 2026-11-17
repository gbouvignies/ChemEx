"""B1 field inhomogeneity distribution models and numerical quadrature.

Available Distributions
-----------------------
ChemEx supports five explicit B1 distribution models in CEST-family experiment
configurations:

- **homogeneous**: Ideal single-point B1 at the nominal value.
- **dephasing**: Asymptotic transverse-dephasing mode for uninterrupted
  single-block continuous-wave CEST (valid when major-dip offsets are masked).
- **gaussian**: Positive-support Gaussian distribution (parent mean ``nominal``,
  parent standard deviation ``scale * nominal``, conditioned on ``B1 >= 0``)
  integrated via Gauss-Legendre quadrature with automatic or explicit ``res``.
- **skewed**: Positive-support Azzalini skew-normal distribution parameterized
  by ``scale``, third standardized moment ``skewness``
  (``|skewness| < 0.9952717``), and explicit ``center`` (``"mode"`` or
  ``"mean"``), reducing smoothly and identically to ``gaussian`` at
  ``skewness = 0.0``.
- **custom**: User-defined relative ``scales`` and ``weights`` (e.g. from an
  experimental FT-nutation profile).
"""

from __future__ import annotations

from chemex.nmr.distributions import (
    custom,
    dephasing,
    gaussian,
    homogeneous,
    skewed,
)

__all__ = [
    "custom",
    "dephasing",
    "gaussian",
    "homogeneous",
    "skewed",
]

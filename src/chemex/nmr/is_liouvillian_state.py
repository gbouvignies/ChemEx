from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from chemex.nmr.b1 import B1Profile
from chemex.nmr.constants import Distribution


def _default_b1_profile() -> B1Profile:
    return B1Profile.gaussian(1e32, scale=0.0, res=11)


def _default_jeff_distribution() -> Distribution:
    return Distribution(np.array([0.0]), np.array([1.0]))


@dataclass(slots=True)
class ISLiouvillianState:
    """Mutable control state for the internal IS Liouvillian engine."""

    par_values: dict[str, float] = field(default_factory=dict)
    ppm_i: float = 0.0
    ppm_s: float = 0.0
    carrier_i: float = 0.0
    carrier_s: float = 0.0
    offset_i: float = 0.0
    offset_s: float = 0.0
    b1_i_profile: B1Profile = field(default_factory=_default_b1_profile)
    b1_s: float = 1e32
    jeff_i: Distribution = field(default_factory=_default_jeff_distribution)
    gradient_dephasing: float = 0.0

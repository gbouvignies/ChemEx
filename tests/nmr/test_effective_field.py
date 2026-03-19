from __future__ import annotations

import numpy as np

from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.effective_field import (
    build_i_effective_field_tilts,
    calculate_i_effective_field_angle,
    tilt_magnetization_along_i_effective_field,
)


def test_calculate_i_effective_field_angle_is_pi_over_two_on_resonance() -> None:
    angle = calculate_i_effective_field_angle(
        b1_i=2_000.0,
        cs_i=0.0,
        ppm_i=1.0,
        carrier_i=0.0,
        offset_i=0.0,
    )

    np.testing.assert_allclose(angle, np.pi * 0.5)


def test_tilt_magnetization_along_i_effective_field_round_trips() -> None:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    tilts = build_i_effective_field_tilts(
        basis,
        ["a", "b"],
        {"cs_i_a": 0.0, "cs_i_b": 0.0},
        b1_i=2_000.0,
        ppm_i=1.0,
        carrier_i=0.0,
        offset_i=0.0,
    )
    original = basis.vectors["ix_a"] + 2.0 * basis.vectors["iz_b"]

    tilted = tilt_magnetization_along_i_effective_field(original.copy(), tilts)
    restored = tilt_magnetization_along_i_effective_field(tilted.copy(), tilts, back=True)

    np.testing.assert_allclose(restored, original)

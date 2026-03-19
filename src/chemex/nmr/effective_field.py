"""Helpers for effective-field rotations in NMR Liouvillian calculations."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np

from chemex.nmr.basis import Basis
from chemex.typing import Array


@dataclass(frozen=True, slots=True)
class EffectiveFieldTilt:
    """Indices and angle for tilting one state's Ix/Iz pair."""

    index_x: int
    index_z: int
    angle: float


def calculate_i_effective_field_angle(
    *,
    b1_i: float,
    cs_i: float,
    ppm_i: float,
    carrier_i: float,
    offset_i: float,
) -> float:
    """Calculate the tilt angle between the z-axis and the effective field."""
    w1 = b1_i * 2.0 * np.pi
    wi = -(cs_i * ppm_i - carrier_i * ppm_i - offset_i * 2.0 * np.pi * np.sign(ppm_i))
    return float(np.arctan2(w1, wi))


def build_i_effective_field_tilts(
    basis: Basis,
    states: Iterable[str],
    par_values: dict[str, float],
    *,
    b1_i: float,
    ppm_i: float,
    carrier_i: float,
    offset_i: float,
) -> tuple[EffectiveFieldTilt, ...]:
    """Build per-state Ix/Iz tilts for rotation along the effective field."""
    component_indices = {
        component_state: index
        for index, component_state in enumerate(basis.components_states)
    }
    return tuple(
        EffectiveFieldTilt(
            index_x=component_indices[f"ix_{state}"],
            index_z=component_indices[f"iz_{state}"],
            angle=calculate_i_effective_field_angle(
                b1_i=b1_i,
                cs_i=par_values[f"cs_i_{state}"],
                ppm_i=ppm_i,
                carrier_i=carrier_i,
                offset_i=offset_i,
            ),
        )
        for state in states
    )


def tilt_magnetization_along_i_effective_field(
    magnetization: Array,
    tilts: Iterable[EffectiveFieldTilt],
    *,
    back: bool = False,
) -> Array:
    """Rotate Ix/Iz components for each state along the effective field."""
    for tilt in tilts:
        angle = -tilt.angle if back else tilt.angle
        rotation_matrix = np.array(
            [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        )
        components = magnetization[..., [tilt.index_x, tilt.index_z], :]
        magnetization[..., [tilt.index_x, tilt.index_z], :] = rotation_matrix @ components
    return magnetization

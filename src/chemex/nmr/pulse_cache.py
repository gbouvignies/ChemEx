"""Helpers for spectrometer pulse cache state."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from chemex.typing import Array


@dataclass(slots=True)
class _BasePulseCache:
    pw90: float = 0.0
    dirty: bool = True

    def invalidate(self) -> None:
        self.dirty = True

    def set_b1(self, value: float) -> None:
        self.pw90 = 1.0 / (4.0 * value) if value else 0.0
        self.dirty = True


@dataclass(slots=True)
class ISpinPulseCache(_BasePulseCache):
    p90: Array = field(default_factory=lambda: np.array(0.0))
    p180: Array = field(default_factory=lambda: np.array(0.0))
    p240: Array = field(default_factory=lambda: np.array(0.0))

    def pulse_widths(self) -> Array:
        return np.array([1.0, 2.0, 8.0 / 3.0]) * self.pw90

    def store(self, *, p90: Array, p180: Array, p240: Array) -> None:
        self.p90 = p90
        self.p180 = p180
        self.p240 = p240
        self.dirty = False


@dataclass(slots=True)
class SSpinPulseCache(_BasePulseCache):
    p180: Array = field(default_factory=lambda: np.array(0.0))

    def pulse_widths(self) -> Array:
        return np.array([2.0]) * self.pw90

    def store(self, *, p180: Array) -> None:
        self.p180 = p180
        self.dirty = False

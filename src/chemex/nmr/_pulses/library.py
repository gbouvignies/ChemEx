from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from chemex.nmr._pulses.kernel import PulseKernel
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
class _ISpinPulseCache(_BasePulseCache):
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
class _SSpinPulseCache(_BasePulseCache):
    p180: Array = field(default_factory=lambda: np.array(0.0))

    def pulse_widths(self) -> Array:
        return np.array([2.0]) * self.pw90

    def store(self, *, p180: Array) -> None:
        self.p180 = p180
        self.dirty = False


class PulseLibrary:
    """Cached base and composite pulse operators built from a pulse kernel."""

    def __init__(self, kernel: PulseKernel) -> None:
        self._kernel = kernel
        self._i_pulses = _ISpinPulseCache()
        self._s_pulses = _SSpinPulseCache()

    def invalidate(self, *spins: str) -> None:
        for spin in spins or ("i", "s"):
            if spin == "i":
                self._i_pulses.invalidate()
            elif spin == "s":
                self._s_pulses.invalidate()
            else:
                msg = f"Unknown spin cache {spin!r}"
                raise ValueError(msg)

    def set_b1_i(self, value: float) -> None:
        self._i_pulses.set_b1(value)

    def set_b1_s(self, value: float) -> None:
        self._s_pulses.set_b1(value)

    def _ensure_base_pulses_i(self) -> None:
        if self._i_pulses.dirty:
            base = self._kernel.pulse_i(self._i_pulses.pulse_widths(), 0.0)
            p90, p180, p240 = self._kernel.add_phases(base, "i").swapaxes(0, 1)
            self._i_pulses.store(p90=p90, p180=p180, p240=p240)

    @property
    def p90_i(self) -> Array:
        self._ensure_base_pulses_i()
        return self._i_pulses.p90

    @property
    def p180_i(self) -> Array:
        self._ensure_base_pulses_i()
        return self._i_pulses.p180

    @property
    def p240_i(self) -> Array:
        self._ensure_base_pulses_i()
        return self._i_pulses.p240

    @property
    def p9018090_i_1(self) -> Array:
        return self.p90_i[[3, 0, 1, 2]] @ self.p180_i @ self.p90_i[[3, 0, 1, 2]]

    @property
    def p9018090_i_2(self) -> Array:
        return self.p90_i[[1, 2, 3, 0]] @ self.p180_i @ self.p90_i[[1, 2, 3, 0]]

    @property
    def p9024090_i_1(self) -> Array:
        return self.p90_i[[3, 0, 1, 2]] @ self.p240_i @ self.p90_i[[3, 0, 1, 2]]

    @property
    def p9024090_i_2(self) -> Array:
        return self.p90_i[[1, 2, 3, 0]] @ self.p240_i @ self.p90_i[[1, 2, 3, 0]]

    def _ensure_base_pulses_s(self) -> None:
        if self._s_pulses.dirty:
            # S-spin caching only has a single base pulse width, so there is no
            # width axis to preserve here. Keep the broadcast distribution axes
            # intact and add only the phase axis.
            p180 = self._kernel.add_phases(
                self._kernel.pulse_s(self._s_pulses.pulse_widths(), 0.0),
                "s",
            )
            self._s_pulses.store(p180=p180)

    @property
    def p180_s(self) -> Array:
        self._ensure_base_pulses_s()
        return self._s_pulses.p180

    def p9024090_nh(self, *, reverse: bool = False) -> Array:
        ph_n = 1 if reverse else 3
        ph_h = 3 if reverse else 1
        pw240i, pw9024090i = np.array([8.0, 14.0]) * self._i_pulses.pw90 / 3.0
        pw240s, pw9024090s = np.array([8.0, 14.0]) * self._s_pulses.pw90 / 3.0
        t0, t1, t2, t3 = 0.5 * np.diff(
            np.sort([pw240i, pw240s, pw9024090i, pw9024090s]),
            prepend=0.0,
        )
        p0 = self._kernel.pulse_is(2.0 * t0, 0, 0)
        if pw9024090i <= pw9024090s:
            p1 = self._kernel.pulse_is(t1, ph_n, 0)
            p2 = (
                self._kernel.pulse_is(t2, ph_n, ph_h)
                if pw9024090i > pw240s
                else self._kernel.pulse_s(t2, 0)
            )
            p3 = self._kernel.pulse_s(t3, ph_h)
        else:
            p1 = self._kernel.pulse_is(t1, 0, ph_h)
            p2 = (
                self._kernel.pulse_is(t2, ph_n, ph_h)
                if pw9024090s > pw240i
                else self._kernel.pulse_i(t2, 0)
            )
            p3 = self._kernel.pulse_i(t3, ph_n)
        pw9024090is_xx = p3 @ p2 @ p1 @ p0 @ p1 @ p2 @ p3
        return self._kernel.add_phases(
            self._kernel.add_phases(pw9024090is_xx, "s"),
            "i",
        )

    @property
    def p9024090_nh_1(self) -> Array:
        return self.p9024090_nh()

    @property
    def p9024090_nh_2(self) -> Array:
        return self.p9024090_nh(reverse=True)

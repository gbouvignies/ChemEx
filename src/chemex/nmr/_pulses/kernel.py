from __future__ import annotations

from collections.abc import Iterable, Sequence
from functools import reduce
from typing import TYPE_CHECKING

import numpy as np

from chemex.nmr._pulses.propagators import (
    calculate_propagators,
    get_phases,
    make_perfect90,
    make_perfect180,
)
from chemex.typing import Array

if TYPE_CHECKING:
    from chemex.nmr._engine.engine import ISLiouvillianEngine


class PulseKernel:
    """Low-level pulse and delay propagator construction for a spectrometer."""

    def __init__(self, engine: ISLiouvillianEngine) -> None:
        self._engine = engine
        self._phases = get_phases(engine)
        self.perfect90_i = self.add_phases(make_perfect90(engine, "i"))
        self.perfect180_i = make_perfect180(engine, "i")
        self.perfect180_s = make_perfect180(engine, "s")

    def add_phases(self, propagator: Array, spin: str = "i") -> Array:
        phases = self._phases[spin]
        return np.array([phases[i] @ propagator @ phases[-i] for i in range(4)])

    def pulse_generation(self, spin: str) -> tuple[int, int]:
        return self._engine.pulse_generation(spin)

    @staticmethod
    def _phase_mixed_term(
        l_x: Array,
        l_y: Array,
        phase: float,
        scale: float = 1.0,
    ) -> Array:
        """Combine a spin's x/y B1 Liouvillian terms for an RF pulse phase."""
        rad = phase * np.pi * 0.5
        return scale * np.cos(rad) * l_x + scale * np.sin(rad) * l_y

    def delays(self, times: float | Iterable[float]) -> Array:
        return calculate_propagators(self._engine.l_free, times)

    def pulse_i(
        self,
        times: float | Iterable[float],
        phase: float,
        scale: float = 1.0,
    ) -> Array:
        dephased = self._engine.b1_i_dist.dephasing
        liouv = self._engine.l_free + self._phase_mixed_term(
            self._engine.l_b1x_i,
            self._engine.l_b1y_i,
            phase,
            scale,
        )
        return calculate_propagators(liouv, times, dephasing=dephased)

    def pulse_s(
        self,
        times: float | Iterable[float],
        phase: float,
        scale: float = 1.0,
    ) -> Array:
        liouv = self._engine.l_free + self._phase_mixed_term(
            self._engine.l_b1x_s,
            self._engine.l_b1y_s,
            phase,
            scale,
        )
        return calculate_propagators(liouv, times)

    def pulse_is(
        self,
        times: float | Iterable[float],
        phase_i: float,
        phase_s: float,
    ) -> Array:
        dephased = self._engine.b1_i_dist.dephasing
        liouv = (
            self._engine.l_free
            + self._phase_mixed_term(
                self._engine.l_b1x_i, self._engine.l_b1y_i, phase_i
            )
            + self._phase_mixed_term(
                self._engine.l_b1x_s, self._engine.l_b1y_s, phase_s
            )
        )
        return calculate_propagators(liouv, times, dephasing=dephased)

    def shaped_pulse_i(
        self,
        pw: float,
        amplitudes: Sequence[float],
        phases: Iterable[float],
    ) -> Array:
        time = pw / len(amplitudes)
        pairs = list(zip(amplitudes, phases, strict=True))
        pulses = {
            (amp, ph): self.pulse_i(time, ph, scale=amp) for amp, ph in set(pairs)
        }
        base = reduce(np.matmul, (pulses[pair] for pair in reversed(pairs)))
        return self.add_phases(base, "i")

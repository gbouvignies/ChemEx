from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Self

import numpy as np

from chemex.configuration.conditions import Conditions
from chemex.nmr import propagators as _propagators
from chemex.nmr.b1 import B1DistributionModel
from chemex.nmr.basis import Basis
from chemex.nmr.constants import Distribution
from chemex.nmr.is_liouvillian_engine import ISLiouvillianEngine
from chemex.nmr.pulse_kernel import PulseKernel
from chemex.nmr.pulse_library import PulseLibrary
from chemex.nmr.spectrometer_analysis import SpectrometerAnalysis
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spin_system.nucleus import Nucleus
from chemex.typing import Array

DictArray = dict[str, Array]

# Keep the historic module-level export for external callers and benchmarks that
# import `chemex.nmr.spectrometer.calculate_propagators`.
calculate_propagators = _propagators.calculate_propagators


class Spectrometer:
    @classmethod
    def from_spin_system(
        cls,
        spin_system: SpinSystem,
        basis: Basis,
        conditions: Conditions,
    ) -> Self:
        return cls(ISLiouvillianEngine(spin_system, basis, conditions))

    def _invalidate_base_pulses(self, *spins: str) -> None:
        self._pulse_library.invalidate(*spins)

    def _set_engine_attr(
        self,
        name: str,
        value: float | str | Distribution,
        *,
        invalidate: tuple[str, ...] = ("i", "s"),
    ) -> None:
        setattr(self._engine, name, value)
        self._invalidate_base_pulses(*invalidate)

    def __init__(self, engine: ISLiouvillianEngine) -> None:
        self._engine = engine
        self._pulse_kernel = PulseKernel(engine)
        self._pulse_library = PulseLibrary(self._pulse_kernel)
        self.analysis = SpectrometerAnalysis(self._engine)
        size = self._engine.size
        self.identity = np.eye(self._engine.size).reshape((1, 1, size, size))
        self.perfect90_i = self._pulse_kernel.perfect90_i
        self.perfect180_i = self._pulse_kernel.perfect180_i
        self.perfect180_s = self._pulse_kernel.perfect180_s
        self.zfilter = self.keep(
            self.identity,
            components={"ie", "se", "iz", "sz", "2izsz"}
            & set(engine.basis.components),
        )

    def keep(self, magnetization: Array, components: Iterable[str]) -> Array:
        return self._engine.keep(magnetization, components)

    def update(self, par_values: dict[str, float]) -> None:
        self._engine.update(par_values)
        self._invalidate_base_pulses()

    @property
    def par_values(self) -> dict[str, float]:
        return self._engine.par_values

    @property
    def basis(self) -> Basis:
        return self._engine.basis

    @property
    def spin_system(self) -> SpinSystem:
        return self._engine.spin_system

    @property
    def ppm_i(self) -> float:
        return self._engine.ppm_i

    @property
    def ppm_s(self) -> float:
        return self._engine.ppm_s

    @property
    def carrier_i(self) -> float:
        return self._engine.carrier_i

    @carrier_i.setter
    def carrier_i(self, value: float) -> None:
        self._set_engine_attr("carrier_i", value)

    @property
    def carrier_s(self) -> float:
        return self._engine.carrier_s

    @carrier_s.setter
    def carrier_s(self, value: float) -> None:
        self._set_engine_attr("carrier_s", value)

    @property
    def offset_i(self) -> float:
        return self._engine.offset_i

    @offset_i.setter
    def offset_i(self, value: float) -> None:
        self._set_engine_attr("offset_i", value)

    @property
    def offset_s(self) -> float:
        return self._engine.offset_s

    @offset_s.setter
    def offset_s(self, value: float) -> None:
        self._set_engine_attr("offset_s", value)

    @property
    def b1_i(self) -> float:
        return self._engine.b1_i

    @b1_i.setter
    def b1_i(self, value: float) -> None:
        self._pulse_library.set_b1_i(value)
        self._engine.b1_i = value
        self._invalidate_base_pulses("i")

    def set_b1_i_inhomogeneity(
        self,
        nominal: float,
        distribution: B1DistributionModel = None,
    ) -> None:
        self._pulse_library.set_b1_i(nominal)
        self._engine.set_b1_i_inhomogeneity(nominal, distribution)
        self._invalidate_base_pulses("i")

    def set_b1_i_distribution(
        self,
        distribution: Distribution,
        *,
        nominal: float | None = None,
    ) -> None:
        if nominal is None:
            nominal = float(
                np.average(distribution.values, weights=distribution.weights),
            )
        self._pulse_library.set_b1_i(nominal)
        self._engine.set_b1_i_distribution(distribution, nominal=nominal)
        self._invalidate_base_pulses("i")

    @property
    def b1_s(self) -> float:
        return self._engine.b1_s

    @b1_s.setter
    def b1_s(self, value: float) -> None:
        self._pulse_library.set_b1_s(value)
        self._engine.b1_s = value
        self._invalidate_base_pulses("s")

    @property
    def jeff_i(self) -> Distribution:
        return self._engine.jeff_i

    @jeff_i.setter
    def jeff_i(self, value: Distribution) -> None:
        self._set_engine_attr("jeff_i", value)

    def get_equilibrium(self) -> Array:
        return self._engine.get_equilibrium()

    def get_start_magnetization(
        self,
        terms: Iterable[str],
        atom: Nucleus = Nucleus.H1,
    ) -> Array:
        return self._engine.get_start_magnetization(terms=terms, atom=atom)

    def tilt_mag_along_weff_i(
        self, magnetization: Array, *, back: bool = False
    ) -> Array:
        return self._engine.tilt_mag_along_weff_i(magnetization, back=back)

    @property
    def detection(self) -> str:
        return self._engine.detection

    @detection.setter
    def detection(self, value: str) -> None:
        self._engine.detection = value

    def detect(self, magnetization: Array) -> float:
        return self._engine.detect(magnetization)

    @property
    def gradient_dephasing(self) -> float:
        return self._engine.gradient_dephasing

    @gradient_dephasing.setter
    def gradient_dephasing(self, value: float) -> None:
        self._set_engine_attr("gradient_dephasing", value)

    def delays(self, times: float | Iterable[float]) -> Array:
        return self._pulse_kernel.delays(times)

    def pulse_i(
        self,
        times: float | Iterable[float],
        phase: float,
        scale: float = 1.0,
    ) -> Array:
        return self._pulse_kernel.pulse_i(times, phase, scale)

    def pulse_s(
        self,
        times: float | Iterable[float],
        phase: float,
        scale: float = 1.0,
    ) -> Array:
        return self._pulse_kernel.pulse_s(times, phase, scale)

    def pulse_is(
        self,
        times: float | Iterable[float],
        phase_i: float,
        phase_s: float,
    ) -> Array:
        return self._pulse_kernel.pulse_is(times, phase_i, phase_s)

    def shaped_pulse_i(
        self,
        pw: float,
        amplitudes: Sequence[float],
        phases: Iterable[float],
    ) -> Array:
        return self._pulse_kernel.shaped_pulse_i(pw, amplitudes, phases)

    @property
    def p90_i(self) -> Array:
        return self._pulse_library.p90_i

    @property
    def p180_i(self) -> Array:
        return self._pulse_library.p180_i

    @property
    def p240_i(self) -> Array:
        return self._pulse_library.p240_i

    @property
    def p9018090_i_1(self) -> Array:
        return self._pulse_library.p9018090_i_1

    @property
    def p9018090_i_2(self) -> Array:
        return self._pulse_library.p9018090_i_2

    @property
    def p9024090_i_1(self) -> Array:
        return self._pulse_library.p9024090_i_1

    @property
    def p9024090_i_2(self) -> Array:
        return self._pulse_library.p9024090_i_2

    @property
    def p180_s(self) -> Array:
        return self._pulse_library.p180_s

    def p9024090_nh(self, *, reverse: bool = False) -> Array:
        return self._pulse_library.p9024090_nh(reverse=reverse)

    @property
    def p9024090_nh_1(self) -> Array:
        return self._pulse_library.p9024090_nh_1

    @property
    def p9024090_nh_2(self) -> Array:
        return self._pulse_library.p9024090_nh_2

    def calculate_shifts(self) -> Array:
        return self.analysis.calculate_shifts()

    def calculate_r1rho(self) -> float:
        return self.analysis.calculate_r1rho()

    def offsets_to_ppms(self, offsets: Array) -> Array:
        return self._engine.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms: Array | float) -> Array | float:
        return self._engine.ppms_to_offsets(ppms)

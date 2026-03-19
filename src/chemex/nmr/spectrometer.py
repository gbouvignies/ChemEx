from __future__ import annotations

from collections.abc import Iterable, Sequence
from functools import reduce

import numpy as np

from chemex.nmr.b1 import B1DistributionModel
from chemex.nmr.constants import Distribution
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.liouvillian_views import reshape_single_liouvillian
from chemex.nmr.propagators import (
    calculate_propagators,
    get_phases,
    make_perfect90,
    make_perfect180,
)
from chemex.nmr.pulse_cache import ISpinPulseCache, SSpinPulseCache
from chemex.parameters.spin_system.nucleus import Nucleus
from chemex.typing import Array

DictArray = dict[str, Array]


class Spectrometer:
    def _add_phases(self, propagator: Array, spin: str = "i") -> Array:
        phases = self._phases[spin]
        return np.array([phases[i] @ propagator @ phases[-i] for i in range(4)])

    def _invalidate_base_pulses(self, *spins: str) -> None:
        for spin in spins or ("i", "s"):
            if spin == "i":
                self._i_pulses.invalidate()
            elif spin == "s":
                self._s_pulses.invalidate()
            else:
                msg = f"Unknown spin cache {spin!r}"
                raise ValueError(msg)

    def _set_liouvillian_attr(
        self,
        name: str,
        value: float | str | Distribution,
        *,
        invalidate: tuple[str, ...] = ("i", "s"),
    ) -> None:
        setattr(self.liouvillian, name, value)
        self._invalidate_base_pulses(*invalidate)

    def __init__(self, liouvillian: LiouvillianIS) -> None:
        self.liouvillian = liouvillian
        size = self.liouvillian.size
        self.identity = np.eye(self.liouvillian.size).reshape((1, 1, size, size))
        self._phases = get_phases(liouvillian)
        self.perfect90_i = self._add_phases(make_perfect90(liouvillian, "i"))
        self.perfect180_i = make_perfect180(liouvillian, "i")
        self.perfect180_s = make_perfect180(liouvillian, "s")
        self.zfilter = self.keep(
            self.identity,
            components={"ie", "se", "iz", "sz", "2izsz"}
            & set(liouvillian.basis.components),
        )
        self._i_pulses = ISpinPulseCache()
        self._s_pulses = SSpinPulseCache()

    def keep(self, magnetization: Array, components: Iterable[str]) -> Array:
        return self.liouvillian.keep(magnetization, components)

    def update(self, par_values: dict[str, float]) -> None:
        self.liouvillian.update(par_values)
        self._invalidate_base_pulses()

    @property
    def par_values(self) -> dict[str, float]:
        return self.liouvillian.par_values

    @property
    def carrier_i(self) -> float:
        return self.liouvillian.carrier_i

    @carrier_i.setter
    def carrier_i(self, value: float) -> None:
        self._set_liouvillian_attr("carrier_i", value)

    @property
    def carrier_s(self) -> float:
        return self.liouvillian.carrier_s

    @carrier_s.setter
    def carrier_s(self, value: float) -> None:
        self._set_liouvillian_attr("carrier_s", value)

    @property
    def offset_i(self) -> float:
        return self.liouvillian.offset_i

    @offset_i.setter
    def offset_i(self, value: float) -> None:
        self._set_liouvillian_attr("offset_i", value)

    @property
    def offset_s(self) -> float:
        return self.liouvillian.offset_s

    @offset_s.setter
    def offset_s(self, value: float) -> None:
        self._set_liouvillian_attr("offset_s", value)

    @property
    def b1_i(self) -> float:
        return self.liouvillian.b1_i

    @b1_i.setter
    def b1_i(self, value: float) -> None:
        self._i_pulses.set_b1(value)
        self.liouvillian.b1_i = value
        self._invalidate_base_pulses("i")

    def set_b1_i_inhomogeneity(
        self,
        nominal: float,
        distribution: B1DistributionModel = None,
    ) -> None:
        self._i_pulses.set_b1(nominal)
        self.liouvillian.set_b1_i_inhomogeneity(nominal, distribution)
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
        self._i_pulses.set_b1(nominal)
        self.liouvillian.set_b1_i_distribution(distribution, nominal=nominal)
        self._invalidate_base_pulses("i")

    @property
    def b1_s(self) -> float:
        return self.liouvillian.b1_s

    @b1_s.setter
    def b1_s(self, value: float) -> None:
        self._s_pulses.set_b1(value)
        self.liouvillian.b1_s = value
        self._invalidate_base_pulses("s")

    @property
    def b1_i_inh_scale(self) -> float:
        return self.liouvillian.b1_i_inh_scale

    @b1_i_inh_scale.setter
    def b1_i_inh_scale(self, value: float) -> None:
        self._set_liouvillian_attr("b1_i_inh_scale", value, invalidate=("i",))

    @property
    def b1_i_inh_res(self) -> int:
        return self.liouvillian.b1_i_inh_res

    @b1_i_inh_res.setter
    def b1_i_inh_res(self, value: int) -> None:
        self.liouvillian.b1_i_inh_res = value
        self._invalidate_base_pulses("i")

    @property
    def jeff_i(self) -> Distribution:
        return self.liouvillian.jeff_i

    @jeff_i.setter
    def jeff_i(self, value: Distribution) -> None:
        self._set_liouvillian_attr("jeff_i", value)

    def get_equilibrium(self) -> Array:
        return self.liouvillian.get_equilibrium()

    def get_start_magnetization(
        self,
        terms: Iterable[str],
        atom: Nucleus = Nucleus.H1,
    ) -> Array:
        return self.liouvillian.get_start_magnetization(terms=terms, atom=atom)

    def tilt_mag_along_weff_i(
        self, magnetization: Array, *, back: bool = False
    ) -> Array:
        return self.liouvillian.tilt_mag_along_weff_i(magnetization, back=back)

    @property
    def detection(self) -> str:
        return self.liouvillian.detection

    @detection.setter
    def detection(self, value: str) -> None:
        self.liouvillian.detection = value

    def detect(self, magnetization: Array) -> float:
        return self.liouvillian.detect(magnetization)

    @property
    def gradient_dephasing(self) -> float:
        return self.liouvillian.gradient_dephasing

    @gradient_dephasing.setter
    def gradient_dephasing(self, value: float) -> None:
        self._set_liouvillian_attr("gradient_dephasing", value)

    def delays(self, times: float | Iterable[float]) -> Array:
        return calculate_propagators(self.liouvillian.l_free, times)

    def pulse_i(
        self,
        times: float | Iterable[float],
        phase: float,
        scale: float = 1.0,
    ) -> Array:
        dephased = self.liouvillian.b1_i_dist.dephasing
        rad = phase * np.pi * 0.5
        liouv = (
            self.liouvillian.l_free
            + scale * np.cos(rad) * self.liouvillian.l_b1x_i
            + scale * np.sin(rad) * self.liouvillian.l_b1y_i
        )
        return calculate_propagators(liouv, times, dephasing=dephased)

    def pulse_s(
        self,
        times: float | Iterable[float],
        phase: float,
        scale: float = 1.0,
    ) -> Array:
        rad = phase * np.pi * 0.5
        liouv = (
            self.liouvillian.l_free
            + scale * np.cos(rad) * self.liouvillian.l_b1x_s
            + scale * np.sin(rad) * self.liouvillian.l_b1y_s
        )
        return calculate_propagators(liouv, times)

    def pulse_is(
        self,
        times: float | Iterable[float],
        phase_i: float,
        phase_s: float,
    ) -> Array:
        dephased = self.liouvillian.b1_i_dist.dephasing
        liouv = (
            self.liouvillian.l_free
            + np.cos(phase_i * np.pi * 0.5) * self.liouvillian.l_b1x_i
            + np.sin(phase_i * np.pi * 0.5) * self.liouvillian.l_b1y_i
            + np.cos(phase_s * np.pi * 0.5) * self.liouvillian.l_b1x_s
            + np.sin(phase_s * np.pi * 0.5) * self.liouvillian.l_b1y_s
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
        return self._add_phases(base, "i")

    def _ensure_base_pulses_i(self) -> None:
        if self._i_pulses.dirty:
            base = self.pulse_i(self._i_pulses.pulse_widths(), 0.0)
            p90, p180, p240 = self._add_phases(base, "i").swapaxes(0, 1)
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
            base = self.pulse_s(self._s_pulses.pulse_widths(), 0.0)
            (p180,) = self._add_phases(base, "s").swapaxes(0, 1)
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
        p0 = self.pulse_is(2.0 * t0, 0, 0)
        if pw9024090i <= pw9024090s:
            p1 = self.pulse_is(t1, ph_n, 0)
            p2 = (
                self.pulse_is(t2, ph_n, ph_h)
                if pw9024090i > pw240s
                else self.pulse_s(t2, 0)
            )
            p3 = self.pulse_s(t3, ph_h)
        else:
            p1 = self.pulse_is(t1, 0, ph_h)
            p2 = (
                self.pulse_is(t2, ph_n, ph_h)
                if pw9024090s > pw240i
                else self.pulse_i(t2, 0)
            )
            p3 = self.pulse_i(t3, ph_n)
        pw9024090is_xx = p3 @ p2 @ p1 @ p0 @ p1 @ p2 @ p3
        return self._add_phases(self._add_phases(pw9024090is_xx, "s"), "i")

    @property
    def p9024090_nh_1(self) -> Array:
        return self.p9024090_nh()

    @property
    def p9024090_nh_2(self) -> Array:
        return self.p9024090_nh(reverse=True)

    def calculate_shifts(self) -> Array:
        liouv = reshape_single_liouvillian(
            self.liouvillian.l_free,
            self.liouvillian.size,
            purpose="Shift eigenvalue calculation",
        )
        return np.linalg.eigvals(liouv.astype(np.complex128)).imag

    def offsets_to_ppms(self, offsets: Array) -> Array:
        return self.liouvillian.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms: Array | float) -> Array | float:
        return self.liouvillian.ppms_to_offsets(ppms)

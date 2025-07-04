from __future__ import annotations

from collections.abc import Hashable, Iterable, Sequence
from functools import reduce
from itertools import product

import numpy as np
from cachetools import cached
from cachetools.keys import hashkey
from scipy.linalg import expm

from chemex.nmr.constants import Distribution
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.typing import ArrayFloat, ArrayNumber

DictArrayFloat = dict[str, ArrayFloat]

# A small value used for numerical stability
SMALL_VALUE = 1e-6


def calculate_propagators(
    liouv: ArrayNumber,
    delays: float | Iterable[float],
    *,
    dephasing: bool = False,
) -> ArrayNumber:
    """Calculate the propagators for the given delays.

    Parameters
    ----------
    liouv : ArrayNumber
        The Liouvillians of the system
    delays : float | Iterable[float]
        The delays to calculate the propagators for
    dephasing : bool, optional
        Whether to apply dephasing, by default False

    Returns
    -------
    ArrayNumber
        The propagators

    """
    # Ensure delays is a 1D NumPy array
    delays_array = np.asarray(delays).flatten()

    if delays_array.size == 1 and not dephasing:
        # If there's only one delay and no dephasing, calculate the propagator directly
        # using scipy expm function
        return expm(liouv * delays_array[0]).astype(np.float64)

    # Calculate eigenvalues and eigenvectors of the Liouvillian matrices
    eigenvalues, eigenvectors = np.linalg.eig(liouv)
    eigenvalues = eigenvalues.astype(
        np.complex128
    )  # Ensure eigenvalues are complex, mainly for type checking

    if dephasing:
        # If dephasing is requested, adjust eigenvalues
        eigenvalues = np.where(
            np.abs(eigenvalues.imag) < SMALL_VALUE,
            eigenvalues,
            eigenvalues * 1e9,
        )

    # Calculate the exponentiated eigenvalues for each delay
    evt = np.exp(np.multiply.outer(delays_array, eigenvalues))

    # Create an empty array for diagonal matrices
    diag = np.empty((*delays_array.shape, *liouv.shape), dtype=np.complex128)

    # Calculate diagonal matrices for each delay
    for index in product(*(range(x) for x in evt.shape[:-1])):
        diag[index] = np.diag(evt[index])

    # Calculate the propagators using matrix multiplication
    propagators = eigenvectors @ diag @ np.linalg.inv(eigenvectors)

    if propagators.shape[0] == 1:
        # If there's only one propagator, return it with no extra time dimension
        propagators = propagators[0]

    # Return the real part of the propagators
    return propagators.real.astype(np.float64)


def _get_key(
    liouvillian: LiouvillianIS, *args: Hashable, **kwargs: Hashable
) -> tuple[Hashable, ...]:
    return tuple(hashkey(liouvillian.basis, *args, **kwargs))


@cached(cache={}, key=_get_key)
def _make_perfect180(liouvillian: LiouvillianIS, spin: str) -> ArrayNumber:
    size = liouvillian.size
    identity = np.eye(size).reshape((1, 1, size, size))
    compx, compy, compz = (f"{spin}{axis}" for axis in "xyz")
    perfect180: dict[str, ArrayNumber] = {
        comp: identity.copy() for comp in (compx, compy)
    }
    for comp in liouvillian.basis.components:
        vect = liouvillian.basis.vectors[comp].ravel()
        if compx in comp or compz in comp:
            perfect180[compy] -= 2 * np.diag(vect)
        if compy in comp or compz in comp:
            perfect180[compx] -= 2 * np.diag(vect)
    p180 = [perfect180[comp] for comp in (compx, compy, compx, compy)]
    return np.array(p180)


@cached(cache={}, key=_get_key)
def _make_perfect90(liouvillian: LiouvillianIS, spin: str) -> ArrayFloat:
    size = liouvillian.size
    zeros = np.zeros((size, size))
    rot = liouvillian.basis.matrices.get(f"b1x_{spin}", zeros)
    return expm(0.25 * rot).reshape(1, 1, size, size).astype(np.float64)


@cached(cache={}, key=_get_key)
def _get_phases(liouvillian: LiouvillianIS) -> DictArrayFloat:
    phases = {}
    size = liouvillian.size
    zeros = np.zeros((size, size))
    for spin in "is":
        l_rotz = liouvillian.basis.matrices.get(f"rotz_{spin}", zeros)
        phases[spin] = np.array([expm(n * 0.5 * np.pi * l_rotz) for n in range(4)])
    return phases


class Spectrometer:
    def _add_phases(self, propagator: ArrayNumber, spin: str = "i") -> ArrayNumber:
        phases = self._phases[spin]
        return np.array([phases[i] @ propagator @ phases[-i] for i in range(4)])

    def __init__(self, liouvillian: LiouvillianIS) -> None:
        self.liouvillian = liouvillian
        size = self.liouvillian.size
        self.identity = np.eye(self.liouvillian.size).reshape((1, 1, size, size))
        self._phases = _get_phases(liouvillian)
        self.perfect90_i = self._add_phases(_make_perfect90(liouvillian, "i"))
        self.perfect180_i = _make_perfect180(liouvillian, "i")
        self.perfect180_s = _make_perfect180(liouvillian, "s")
        self.zfilter = self.keep(
            self.identity,
            components={"ie", "se", "iz", "sz", "2izsz"}
            & set(liouvillian.basis.components),
        )
        self.calculate_i_flag = True
        self.calculate_s_flag = True
        self._pw90_i = 0.0
        self._pw90_s = 0.0
        self._p90_i = np.array(0.0)
        self._p180_i = np.array(0.0)
        self._p240_i = np.array(0.0)
        self._p180_s = np.array(0.0)

    def keep(self, magnetization: ArrayFloat, components: Iterable[str]) -> ArrayNumber:
        return self.liouvillian.keep(magnetization, components)

    def update(self, par_values: dict[str, float]) -> None:
        self.liouvillian.update(par_values)
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    @property
    def par_values(self) -> dict[str, float]:
        return self.liouvillian.par_values

    @property
    def carrier_i(self) -> float:
        return self.liouvillian.carrier_i

    @carrier_i.setter
    def carrier_i(self, value: float) -> None:
        self.liouvillian.carrier_i = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    @property
    def carrier_s(self) -> float:
        return self.liouvillian.carrier_s

    @carrier_s.setter
    def carrier_s(self, value: float) -> None:
        self.liouvillian.carrier_s = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    @property
    def offset_i(self) -> float:
        return self.liouvillian.offset_i

    @offset_i.setter
    def offset_i(self, value: float) -> None:
        self.liouvillian.offset_i = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    @property
    def offset_s(self) -> float:
        return self.liouvillian.offset_s

    @offset_s.setter
    def offset_s(self, value: float) -> None:
        self.liouvillian.offset_s = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    @property
    def b1_i(self) -> float:
        return self.liouvillian.b1_i

    @b1_i.setter
    def b1_i(self, value: float) -> None:
        self._pw90_i = 1.0 / (4.0 * value) if value else 0.0
        self.liouvillian.b1_i = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    @property
    def b1_s(self) -> float:
        return self.liouvillian.b1_s

    @b1_s.setter
    def b1_s(self, value: float) -> None:
        self._pw90_s = 1.0 / (4.0 * value) if value else 0.0
        self.liouvillian.b1_s = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    @property
    def b1_i_inh_scale(self) -> float:
        return self.liouvillian.b1_i_inh_scale

    @b1_i_inh_scale.setter
    def b1_i_inh_scale(self, value: float) -> None:
        self.liouvillian.b1_i_inh_scale = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    @property
    def b1_i_inh_res(self) -> int:
        return self.liouvillian.b1_i_inh_res

    @b1_i_inh_res.setter
    def b1_i_inh_res(self, value: int) -> None:
        self.liouvillian.b1_i_inh_res = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    @property
    def jeff_i(self) -> Distribution:
        return self.liouvillian.jeff_i

    @jeff_i.setter
    def jeff_i(self, value: Distribution) -> None:
        self.liouvillian.jeff_i = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    def get_equilibrium(self) -> ArrayFloat:
        return self.liouvillian.get_equilibrium()

    def get_start_magnetization(
        self,
        terms: Iterable[str],
        atom: str = "h",
    ) -> ArrayFloat:
        return self.liouvillian.get_start_magnetization(terms=terms, atom=atom)

    def tilt_mag_along_weff_i(
        self, magnetization: ArrayFloat, *, back: bool = False
    ) -> ArrayNumber:
        return self.liouvillian.tilt_mag_along_weff_i(magnetization, back=back)

    @property
    def detection(self) -> str:
        return self.liouvillian.detection

    @detection.setter
    def detection(self, value: str) -> None:
        self.liouvillian.detection = value

    def detect(self, magnetization: ArrayFloat) -> float:
        return self.liouvillian.detect(magnetization)

    @property
    def gradient_dephasing(self) -> float:
        return self.liouvillian.gradient_dephasing

    @gradient_dephasing.setter
    def gradient_dephasing(self, value: float) -> None:
        self.liouvillian.gradient_dephasing = value
        self.calculate_i_flag = True
        self.calculate_s_flag = True

    def delays(self, times: float | Iterable[float]) -> ArrayNumber:
        return calculate_propagators(self.liouvillian.l_free, times)

    def pulse_i(
        self,
        times: float | Iterable[float],
        phase: float,
        scale: float = 1.0,
    ) -> ArrayNumber:
        dephased = self.b1_i_inh_scale == np.inf
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
    ) -> ArrayNumber:
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
    ) -> ArrayNumber:
        dephased = self.b1_i_inh_scale == np.inf
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
    ) -> ArrayNumber:
        time = pw / len(amplitudes)
        pairs = list(zip(amplitudes, phases, strict=True))
        pulses = {
            (amp, ph): self.pulse_i(time, ph, scale=amp) for amp, ph in set(pairs)
        }
        base = reduce(np.matmul, (pulses[pair] for pair in reversed(pairs)))
        return self._add_phases(base, "i")

    def _calculate_base_pulses_i(self) -> None:
        if self.calculate_i_flag:
            pws = np.array([1.0, 2.0, 8.0 / 3.0]) * self._pw90_i
            base = self.pulse_i(pws, 0.0)
            pulses = self._add_phases(base, "i")
            self._p90_i, self._p180_i, self._p240_i = pulses.swapaxes(0, 1)
            self.calculate_i_flag = False

    @property
    def p90_i(self) -> ArrayFloat:
        self._calculate_base_pulses_i()
        return self._p90_i

    @property
    def p180_i(self) -> ArrayFloat:
        self._calculate_base_pulses_i()
        return self._p180_i

    @property
    def p240_i(self) -> ArrayFloat:
        self._calculate_base_pulses_i()
        return self._p240_i

    @property
    def p9018090_i_1(self) -> ArrayFloat:
        return self.p90_i[[3, 0, 1, 2]] @ self.p180_i @ self.p90_i[[3, 0, 1, 2]]

    @property
    def p9018090_i_2(self) -> ArrayFloat:
        return self.p90_i[[1, 2, 3, 0]] @ self.p180_i @ self.p90_i[[1, 2, 3, 0]]

    @property
    def p9024090_i_1(self) -> ArrayFloat:
        return self.p90_i[[3, 0, 1, 2]] @ self.p240_i @ self.p90_i[[3, 0, 1, 2]]

    @property
    def p9024090_i_2(self) -> ArrayFloat:
        return self.p90_i[[1, 2, 3, 0]] @ self.p240_i @ self.p90_i[[1, 2, 3, 0]]

    def _calculate_base_pulses_s(self) -> None:
        if self.calculate_s_flag:
            pws = np.array([2.0]) * self._pw90_s
            base = self.pulse_s(pws, 0.0)
            pulses = self._add_phases(base, "s")
            self._p180_s = pulses.swapaxes(0, 1)
            self.calculate_s_flag = False

    @property
    def p180_s(self) -> ArrayNumber:
        self._calculate_base_pulses_s()
        return self._p180_s

    def p9024090_nh(self, *, reverse: bool = False) -> ArrayNumber:
        ph_n = 1 if reverse else 3
        ph_h = 3 if reverse else 1
        pw240i, pw9024090i = np.array([8.0, 14.0]) * self._pw90_i / 3.0
        pw240s, pw9024090s = np.array([8.0, 14.0]) * self._pw90_s / 3.0
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
    def p9024090_nh_1(self) -> ArrayNumber:
        return self.p9024090_nh()

    @property
    def p9024090_nh_2(self) -> ArrayNumber:
        return self.p9024090_nh(reverse=True)

    def calculate_shifts(self) -> ArrayNumber:
        liouv = self.liouvillian.l_free.reshape(
            (self.liouvillian.size, self.liouvillian.size),
        )
        return np.linalg.eigvals(liouv.astype(np.complex128)).imag

    def offsets_to_ppms(self, offsets: ArrayFloat) -> ArrayFloat:
        return self.liouvillian.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms: ArrayFloat | float) -> ArrayFloat | float:
        return self.liouvillian.ppms_to_offsets(ppms)

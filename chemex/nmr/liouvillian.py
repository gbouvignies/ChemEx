""" The ref module contains the reference matrices and code for calculating the
    Liouvillian.

    Operator basis::

       { Eq,
         Ix, Iy, Iz, Sx, Sy, Sz,
         2IxSz, 2IySz, 2IzSx, 2IzSy,
         2IxSx, 2IxSy, 2IySx, 2IySy,
         2IzSz }

"""
from __future__ import annotations

import re
from collections.abc import Iterable
from itertools import product
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats

from chemex.configuration.conditions import Conditions
from chemex.model import model
from chemex.nmr.basis import Basis
from chemex.nmr.constants import Distribution
from chemex.nmr.constants import SIGNED_XI_RATIO
from chemex.nmr.constants import XI_RATIO

if TYPE_CHECKING:
    from chemex.parameters.spin_system import SpinSystem


_RE_COMP = re.compile(r"\[(.+?)\]")

_Q_ORDER_I = {"sq": 1.0, "dq": 2.0, "tq": 3.0}


def _make_gaussian(
    value: float, scale: float, res: int
) -> tuple[np.ndarray, np.ndarray]:
    if scale not in (0.0, np.inf) and res > 1:
        grid = np.linspace(-2.0, 2.0, res)
        dist = grid * scale + 1.0
    else:
        grid = np.array([0.0])
        dist = np.array([1.0])
    weights = stats.norm.pdf(grid)  # type: ignore
    weights /= weights.sum()
    return dist * value, weights


class LiouvillianIS:
    def _scale_matrix(self, name: str, value: float):
        state_names = {name.format(state=state) for state in model.states}
        for state_name in state_names & set(self.matrices):
            self.matrices[state_name] = np.sign(self._matrices_ref[state_name]) * value

    def __init__(self, spin_system: SpinSystem, basis: Basis, conditions: Conditions):
        self.spin_system = spin_system.correct(basis)
        self.basis = basis
        self.h_frq = 0.0 if conditions.h_larmor_frq is None else conditions.h_larmor_frq
        self.par_values: dict[str, float] = {}
        self.size = len(model.states) * len(basis)
        self.vectors = basis.vectors.copy()
        self.matrices = basis.matrices.copy()
        self._matrices_ref = self.matrices.copy()
        self._detection: str = ""
        self._detect_vector: np.ndarray = np.array([])
        self._q_order_i = _Q_ORDER_I.get(self.basis.extension, 1.0)
        scale = -2.0 * np.pi * self.h_frq
        self.ppm_i = scale * SIGNED_XI_RATIO.get(basis.atoms.get("i", ""), 1.0)
        self.ppm_s = scale * SIGNED_XI_RATIO.get(basis.atoms.get("s", ""), 1.0)
        self.carrier_i = 0.0
        self.carrier_s = 0.0
        self.offset_i = 0.0
        self.offset_s = 0.0
        self._l_base = np.array(0.0)
        self._b1_i_inh_scale = 0.0
        self._b1_i_inh_res = 11
        self.b1_i = 1e32
        self.b1_s = 1e32
        self.jeff_i = Distribution(np.array([0.0]), np.array([1.0]))
        self.gradient_dephasing = 0.0
        self._scale_matrix("j_is_{state}", self._q_order_i * np.pi)

    def _build_base_liouvillian(self) -> None:
        self._l_base = sum(
            (
                self.matrices[name] * self.par_values.get(name, 0.0)
                for name in self.matrices
            ),
            start=np.array(0.0),
        )

    def update(self, par_values: dict[str, float]) -> None:
        self.par_values = par_values
        self._build_base_liouvillian()

    @property
    def ppm_i(self) -> float:
        return self._ppm_i

    @ppm_i.setter
    def ppm_i(self, value: float):
        self._ppm_i = value
        self._scale_matrix("cs_i_{state}", self._q_order_i * value)
        self._scale_matrix("carrier_i", self._q_order_i * value)

    @property
    def ppm_s(self) -> float:
        return self._ppm_s

    @ppm_s.setter
    def ppm_s(self, value: float):
        self._ppm_s = value
        self._scale_matrix("cs_s_{state}", value)
        self._scale_matrix("carrier_s", value)

    @property
    def carrier_i(self) -> float:
        return float(self._carrier_i)

    @carrier_i.setter
    def carrier_i(self, value: float):
        self._carrier_i = np.array(value)
        self._l_carrier_i = self.matrices.get("carrier_i", np.array(0.0)) * value

    @property
    def carrier_s(self) -> float:
        return float(self._carrier_s)

    @carrier_s.setter
    def carrier_s(self, value: float):
        self._carrier_s = np.array(value)
        self._l_carrier_s = self.matrices.get("carrier_s", np.array(0.0)) * value

    @property
    def offset_i(self) -> float:
        return float(self._offset_i)

    @offset_i.setter
    def offset_i(self, value: float):
        self._offset_i = np.array(value)
        self._l_offset_i = (
            self.matrices.get("offset_i", np.array(0.0)) * value * np.sign(self.ppm_i)
        )

    @property
    def offset_s(self) -> float:
        return float(self._offset_s)

    @offset_s.setter
    def offset_s(self, value: float):
        self._offset_s = np.array(value)
        self._l_offset_s = (
            self.matrices.get("offset_s", 0.0) * value * np.sign(self.ppm_s)
        )

    @property
    def b1_i_inh_scale(self) -> float:
        return self._b1_i_inh_scale

    @b1_i_inh_scale.setter
    def b1_i_inh_scale(self, value: float):
        self._b1_i_inh_scale = value
        self.b1_i = self._b1_i

    @property
    def b1_i_inh_res(self) -> int:
        return self._b1_i_inh_res

    @b1_i_inh_res.setter
    def b1_i_inh_res(self, value: int):
        self._b1_i_inh_res = value
        self.b1_i = self._b1_i

    @property
    def b1_i(self) -> float:
        return self._b1_i

    @b1_i.setter
    def b1_i(self, value: float):
        self._b1_i = value
        self._b1_i_dist, self._b1_i_weights = _make_gaussian(
            self.b1_i, self.b1_i_inh_scale, self.b1_i_inh_res
        )
        self._b1_i_dist = self._b1_i_dist.reshape((-1, 1, 1))
        self._b1_i_weights = self._b1_i_weights.reshape((-1, 1, 1))
        self.l_b1x_i = self.matrices.get("b1x_i", 0.0) * self._b1_i_dist
        self.l_b1y_i = self.matrices.get("b1y_i", 0.0) * self._b1_i_dist

    @property
    def b1_s(self) -> float:
        return self._b1_s

    @b1_s.setter
    def b1_s(self, value: float):
        self._b1_s = value
        self.l_b1x_s = self.matrices.get("b1x_s", 0.0) * value
        self.l_b1y_s = self.matrices.get("b1y_s", 0.0) * value

    @property
    def jeff_i(self) -> Distribution:
        return self._jeff_i

    @jeff_i.setter
    def jeff_i(self, value: Distribution):
        self._jeff_i = value
        values = value.values.reshape((-1, 1, 1, 1))
        self._l_jeff_i = self.matrices.get("jeff_i", np.array(0.0)) * values
        self._jeff_i_weights = value.weights.reshape((-1, 1, 1, 1))

    @property
    def gradient_dephasing(self) -> float:
        return self._gradient_dephasing

    @gradient_dephasing.setter
    def gradient_dephasing(self, value: float):
        self._gradient_dephasing = value
        self._scale_matrix("d_{state}", value * 1e-12)
        self._build_base_liouvillian()

    @property
    def l_free(self) -> np.ndarray:
        return sum(
            (
                self._l_base,
                self._l_offset_i,
                self._l_offset_s,
                self._l_carrier_i,
                self._l_carrier_s,
                self._l_jeff_i,
            ),
            start=np.array(0.0),
        )

    @property
    def weights(self) -> np.ndarray:
        return self._b1_i_weights * self._jeff_i_weights

    @property
    def detection(self):
        return self._detection

    @detection.setter
    def detection(self, value: str):
        self._detection = value
        expr = _RE_COMP.sub(r'self.vectors.get("\g<1>")', value)
        vector: np.ndarray = eval(expr)
        self._detect_vector = vector.transpose()

    def detect(self, magnetization: np.ndarray) -> float:
        shape = -1, *magnetization.shape[-2:]
        mag_weighted = self.weights * magnetization
        mag = mag_weighted.reshape(shape).sum(axis=0)
        detected = self._detect_vector @ mag
        if np.iscomplexobj(detected):
            detected = np.sign(detected.real) * np.abs(detected)
        return float(detected)

    def get_equilibrium(self) -> np.ndarray:
        mag = np.zeros((self.size, 1))
        for state, (name, atom) in product(model.states, self.basis.atoms.items()):
            scale = self.par_values.get(f"p{state}", 0.0) * XI_RATIO.get(atom, 1.0)
            mag += self.vectors.get(f"{name}e_{state}", 0.0) * scale
            mag += self.vectors.get(f"{name}z_{state}", 0.0) * scale
        return mag

    def get_start_magnetization(
        self, terms: Iterable[str], atom: str = "h"
    ) -> np.ndarray:
        ratio = XI_RATIO.get(atom, 1.0)
        mag = np.zeros((self.size, 1))
        for term, state, (comp, vector) in product(
            terms, model.states, self.vectors.items()
        ):
            if comp.startswith(term) and comp.endswith(f"_{state}"):
                mag += self.par_values[f"p{state}"] * ratio * vector
        return mag

    def keep(self, magnetization: np.ndarray, components: Iterable[str]) -> np.ndarray:
        keep = sum(
            (self.vectors[name] for name in components), start=np.zeros((self.size, 1))
        )
        keep[keep > 0.0] = 1.0
        return keep * magnetization

    def offsets_to_ppms(self, offsets: np.ndarray) -> np.ndarray:
        return self.carrier_i + 2.0 * np.pi * offsets / abs(self.ppm_i)

    def ppms_to_offsets(self, ppms: np.ndarray | float) -> np.ndarray | float:
        return (ppms - self.carrier_i) * abs(self.ppm_i) / (2.0 * np.pi)

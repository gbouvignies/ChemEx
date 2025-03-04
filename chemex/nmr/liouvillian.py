"""The ref module contains the code for calculating the Liouvillians.

Operator basis:

{ Eq,
Ix, Iy, Iz, Sx, Sy, Sz,
2IxSz, 2IySz, 2IzSx, 2IzSy,
2IxSx, 2IxSy, 2IySx, 2IySy,
2IzSz }

"""

import re
from collections.abc import Iterable
from itertools import product

import numpy as np
from scipy import stats

from chemex.configuration.conditions import Conditions
from chemex.models.model import model
from chemex.nmr.basis import Basis
from chemex.nmr.constants import SIGNED_XI_RATIO, XI_RATIO, Distribution
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spin_system.nucleus import Nucleus
from chemex.typing import ArrayFloat, ArrayNumber

_RE_COMP = re.compile(r"\[(.+?)\]")

_Q_ORDER_I = {"sq": 1.0, "dq": 2.0, "tq": 3.0}

# A small value used for numerical stability
SMALL_VALUE = 1e-6


def _generate_gaussian_distribution(
    value: float,
    scale: float,
    res: int,
) -> Distribution:
    if scale not in (0.0, np.inf) and res > 1:
        grid = np.linspace(-2.0, 2.0, res)
        distribution = grid * scale + 1.0
    else:
        grid = np.array([0.0])
        distribution = np.array([1.0])

    weights = stats.norm.pdf(grid)
    weights /= weights.sum()

    return Distribution(distribution * value, weights)


class LiouvillianIS:
    def __init__(
        self,
        spin_system: SpinSystem,
        basis: Basis,
        conditions: Conditions,
    ) -> None:
        # Initialize attributes
        self.spin_system = spin_system.correct(basis)
        self.basis = basis
        self.h_frq = 0.0 if conditions.h_larmor_frq is None else conditions.h_larmor_frq
        self.par_values: dict[str, float] = {}
        self.size = len(model.states) * len(basis)
        self._detection: str = ""
        self._detect_vector: ArrayNumber = np.array([])
        self._q_order_i = _Q_ORDER_I.get(self.basis.extension, 1.0)
        scale = -2.0 * np.pi * self.h_frq
        self.ppm_i = scale * SIGNED_XI_RATIO.get(basis.nuclei.get("i", Nucleus.X), 1.0)
        self.ppm_s = scale * SIGNED_XI_RATIO.get(basis.nuclei.get("s", Nucleus.X), 1.0)
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
        self.basis.scale_matrices("j_is_{state}", self._q_order_i * np.pi)

    def _build_base_liouvillian(self) -> None:
        self._l_base = sum(
            (
                self.basis.matrices[name] * self.par_values.get(name, 0.0)
                for name in self.basis.matrices
            ),
            start=np.array(0.0),
        )

    def _collapse(self, magnetization: ArrayNumber) -> ArrayNumber:
        """Collapse the distribution of magnetization into an average."""
        if (ndim := magnetization.ndim) < 3:
            return magnetization
        magnetization_weighted = self.weights * magnetization
        sum_axes = tuple(range(ndim - 2))
        return magnetization_weighted.sum(axis=sum_axes)

    @property
    def ppm_i(self) -> float:
        return self._ppm_i

    @ppm_i.setter
    def ppm_i(self, value: float) -> None:
        self._ppm_i = value
        self.basis.scale_matrices(
            ["cs_i_{state}", "carrier_i"], self._q_order_i * value
        )

    @property
    def ppm_s(self) -> float:
        return self._ppm_s

    @ppm_s.setter
    def ppm_s(self, value: float) -> None:
        self._ppm_s = value
        self.basis.scale_matrices(["cs_s_{state}", "carrier_s"], value)

    @property
    def carrier_i(self) -> float:
        return float(self._carrier_i)

    @carrier_i.setter
    def carrier_i(self, value: float) -> None:
        self._carrier_i = np.array(value)
        self._l_carrier_i = self.basis.matrices.get("carrier_i", np.array(0.0)) * value

    @property
    def carrier_s(self) -> float:
        return float(self._carrier_s)

    @carrier_s.setter
    def carrier_s(self, value: float) -> None:
        self._carrier_s = np.array(value)
        self._l_carrier_s = self.basis.matrices.get("carrier_s", np.array(0.0)) * value

    @property
    def offset_i(self) -> float:
        return float(self._offset_i)

    @offset_i.setter
    def offset_i(self, value: float) -> None:
        self._offset_i = np.array(value)
        self._l_offset_i = (
            self.basis.matrices.get("offset_i", np.array(0.0))
            * value
            * np.sign(self.ppm_i)
        )

    @property
    def offset_s(self) -> float:
        return float(self._offset_s)

    @offset_s.setter
    def offset_s(self, value: float) -> None:
        self._offset_s = np.array(value)
        self._l_offset_s = (
            self.basis.matrices.get("offset_s", np.array(0.0))
            * value
            * np.sign(self.ppm_s)
        )

    @property
    def b1_i_inh_scale(self) -> float:
        return self._b1_i_inh_scale

    @b1_i_inh_scale.setter
    def b1_i_inh_scale(self, value: float) -> None:
        self._b1_i_inh_scale = value
        self.b1_i = self._b1_i

    @property
    def b1_i_inh_res(self) -> int:
        return self._b1_i_inh_res

    @b1_i_inh_res.setter
    def b1_i_inh_res(self, value: int) -> None:
        self._b1_i_inh_res = value
        self.b1_i = self._b1_i

    @property
    def b1_i(self) -> float:
        return self._b1_i

    @b1_i.setter
    def b1_i(self, value: float) -> None:
        self._b1_i = value
        gaussian_dist = _generate_gaussian_distribution(
            self.b1_i,
            self.b1_i_inh_scale,
            self.b1_i_inh_res,
        )
        gaussian_dist.values = gaussian_dist.values.reshape((-1, 1, 1))
        gaussian_dist.weights = gaussian_dist.weights.reshape((-1, 1, 1))
        self._b1_i_dist = gaussian_dist
        self.l_b1x_i = self.basis.matrices.get("b1x_i", 0.0) * gaussian_dist.values
        self.l_b1y_i = self.basis.matrices.get("b1y_i", 0.0) * gaussian_dist.values

    @property
    def b1_s(self) -> float:
        return self._b1_s

    @b1_s.setter
    def b1_s(self, value: float) -> None:
        self._b1_s = value
        self.l_b1x_s = self.basis.matrices.get("b1x_s", 0.0) * value
        self.l_b1y_s = self.basis.matrices.get("b1y_s", 0.0) * value

    @property
    def jeff_i(self) -> Distribution:
        return self._jeff_i

    @jeff_i.setter
    def jeff_i(self, distribution: Distribution) -> None:
        self._jeff_i = distribution
        values = distribution.values.reshape((-1, 1, 1, 1))
        self._l_jeff_i = self.basis.matrices.get("jeff_i", np.array(0.0)) * values
        self._jeff_i_weights = distribution.weights.reshape((-1, 1, 1, 1))

    @property
    def gradient_dephasing(self) -> float:
        return self._gradient_dephasing

    @gradient_dephasing.setter
    def gradient_dephasing(self, value: float) -> None:
        self._gradient_dephasing = value
        self.basis.scale_matrices("d_{state}", value * 1e-12)
        self._build_base_liouvillian()

    @property
    def l_free(self) -> ArrayNumber:
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
    def weights(self) -> ArrayNumber:
        return self._b1_i_dist.weights * self._jeff_i_weights

    @property
    def detection(self) -> str:
        return self._detection

    @detection.setter
    def detection(self, value: str) -> None:
        self._detection = value
        expr = _RE_COMP.sub(r'self.basis.vectors.get("\1")', value)
        vector: ArrayNumber = eval(expr)
        self._detect_vector = vector.transpose()

    def update(self, par_values: dict[str, float]) -> None:
        self.par_values = par_values
        self._build_base_liouvillian()

    def detect(self, magnetization: ArrayNumber) -> float:
        collapsed_magnetization = self._collapse(magnetization)
        detected = self._detect_vector @ collapsed_magnetization
        if np.iscomplexobj(detected):
            detected = np.sign(detected.real) * np.abs(detected)
        return float(detected)

    # def detect_spectrum(
    #     self, magnetization: ArrayNumber, observed_state: str = "a"
    # ) -> ArrayNumber:
    #     collapsed_magnetization = self._collapse(magnetization)
    #     component, _state = self.detection.split("_")
    #     self.basis.vectors.get("ix", 0.0) - 1j * self.basis.vectors.get(
    #         "iy", 0.0
    #     )

    #     # Getting the eigenvalues and eigenvectors of the Liouvillian
    #     eigen_values, eigen_vectors = np.linalg.eig(self._l_base)

    #     return self._detect_vector @ collapsed_magnetization

    def calculate_r1rho(self) -> float:
        liouv = self.l_free + self.l_b1x_i
        liouv = liouv.reshape((self.size, self.size))
        eigenvalues = np.linalg.eigvals(liouv)
        real_eigenvalues = eigenvalues[abs(eigenvalues.imag) < SMALL_VALUE].real
        return -float(np.max(real_eigenvalues))

    def tilt_mag_along_weff_i(
        self, magnetization: ArrayNumber, *, back: bool = False
    ) -> ArrayNumber:
        basis = self.basis

        w1 = self.b1_i * 2.0 * np.pi

        for state in model.states:
            index_x, index_z = (
                basis.components_states.index(f"ix_{state}"),
                basis.components_states.index(f"iz_{state}"),
            )
            components = magnetization[..., [index_x, index_z], :]

            cs = self.par_values[f"cs_i_{state}"]
            wi = -(
                cs * self.ppm_i
                - self.carrier_i * self.ppm_i
                - self.offset_i * 2.0 * np.pi * np.sign(self.ppm_i)
            )
            theta = np.arctan2(w1, wi)
            if back:
                theta = -theta
            rotation_matrix = np.array(
                [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
            )
            components_tilted = rotation_matrix @ components
            magnetization[..., [index_x, index_z], :] = components_tilted

        return magnetization

    def get_equilibrium(self) -> ArrayFloat:
        mag = np.zeros((self.size, 1))
        for state, (name, atom) in product(model.states, self.basis.atoms.items()):
            scale = self.par_values.get(f"p{state}", 0.0) * XI_RATIO.get(atom, 1.0)
            mag += self.basis.vectors.get(f"{name}e_{state}", 0.0) * scale
            mag += self.basis.vectors.get(f"{name}z_{state}", 0.0) * scale
        return mag

    def get_start_magnetization(
        self, terms: Iterable[str], atom: str = "h"
    ) -> ArrayFloat:
        ratio = XI_RATIO.get(atom, 1.0)
        terms_set = set(terms)
        mag = np.zeros((self.size, 1))

        for component, vector in self.basis.vectors.items():
            state_specific = "_" in component
            if not state_specific:
                continue
            state = component[-1]
            for term in terms_set:
                match = component.startswith(term.strip("+-"))
                if match:
                    sign = -1.0 if term.startswith("-") else 1.0
                    population = self.par_values.get(f"p{state}", 0.0)
                    mag += sign * population * ratio * vector

        return mag

    def keep(
        self, magnetization: ArrayNumber, components: Iterable[str]
    ) -> ArrayNumber:
        keep = sum(
            (self.basis.vectors[name] for name in components),
            start=np.zeros((self.size, 1)),
        )
        keep[keep > 0] = 1.0
        return keep * magnetization

    def offsets_to_ppms(self, offsets: ArrayFloat) -> ArrayFloat:
        return self.carrier_i + 2.0 * np.pi * offsets / abs(self.ppm_i)

    def ppms_to_offsets(self, ppms: ArrayFloat | float) -> ArrayFloat | float:
        return (ppms - self.carrier_i) * abs(self.ppm_i) / (2.0 * np.pi)

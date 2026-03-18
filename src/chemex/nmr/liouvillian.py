"""The ref module contains the code for calculating the Liouvillians.

Operator basis:

{ Eq,
Ix, Iy, Iz, Sx, Sy, Sz,
2IxSz, 2IySz, 2IzSx, 2IzSy,
2IxSx, 2IxSy, 2IySx, 2IySy,
2IzSz }

"""

from collections.abc import Iterable

import numpy as np

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.constants import SIGNED_XI_RATIO, Distribution
from chemex.nmr.detection import build_detection_vector
from chemex.nmr.distributions import get_b1_distribution
from chemex.nmr.magnetization import (
    build_equilibrium_magnetization,
    build_start_magnetization,
    detect_signal,
    keep_components,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spin_system.nucleus import Nucleus
from chemex.typing import Array

_Q_ORDER_I = {"sq": 1.0, "dq": 2.0, "tq": 3.0}

# A small value used for numerical stability
SMALL_VALUE = 1e-6


class LiouvillianIS:
    def __init__(
        self,
        spin_system: SpinSystem,
        basis: Basis,
        conditions: Conditions,
    ) -> None:
        # Initialize attributes
        self.model: ModelSpec = basis.model
        self.spin_system = spin_system.correct(basis)
        self.basis = basis
        self._matrices = basis.copy_matrices()
        self.h_frq = 0.0 if conditions.h_larmor_frq is None else conditions.h_larmor_frq
        self.par_values: dict[str, float] = {}
        self.size = len(self.model.states) * len(basis)
        self._detection: str = ""
        self._detect_vector: Array = np.array([])
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
        self._scale_matrices("j_is_{state}", self._q_order_i * np.pi)

    def _scale_matrices(self, names: list[str] | str, value: float) -> None:
        if isinstance(names, str):
            names = [names]
        for name in names:
            matrix_names = {name.format(state=state) for state in self.model.states}
            for matrix_name in matrix_names & self._matrices.keys():
                self._matrices[matrix_name] = (
                    np.sign(self.basis.matrices[matrix_name]) * value
                )

    def _build_base_liouvillian(self) -> None:
        self._l_base = sum(
            (
                self._matrices[name] * self.par_values.get(name, 0.0)
                for name in self._matrices
            ),
            start=np.array(0.0),
        )

    @property
    def ppm_i(self) -> float:
        return self._ppm_i

    @ppm_i.setter
    def ppm_i(self, value: float) -> None:
        self._ppm_i = value
        self._scale_matrices(["cs_i_{state}", "carrier_i"], self._q_order_i * value)

    @property
    def ppm_s(self) -> float:
        return self._ppm_s

    @ppm_s.setter
    def ppm_s(self, value: float) -> None:
        self._ppm_s = value
        self._scale_matrices(["cs_s_{state}", "carrier_s"], value)

    @property
    def carrier_i(self) -> float:
        return float(self._carrier_i)

    @carrier_i.setter
    def carrier_i(self, value: float) -> None:
        self._carrier_i = np.array(value)
        self._l_carrier_i = self._matrices.get("carrier_i", np.array(0.0)) * value

    @property
    def carrier_s(self) -> float:
        return float(self._carrier_s)

    @carrier_s.setter
    def carrier_s(self, value: float) -> None:
        self._carrier_s = np.array(value)
        self._l_carrier_s = self._matrices.get("carrier_s", np.array(0.0)) * value

    @property
    def offset_i(self) -> float:
        return float(self._offset_i)

    @offset_i.setter
    def offset_i(self, value: float) -> None:
        self._offset_i = np.array(value)
        self._l_offset_i = (
            self._matrices.get("offset_i", np.array(0.0))
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
            self._matrices.get("offset_s", np.array(0.0))
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
        gaussian_dist = get_b1_distribution(
            distribution_type="gaussian",
            value=self.b1_i,
            scale=self.b1_i_inh_scale,
            res=self.b1_i_inh_res,
        )
        gaussian_dist.values = gaussian_dist.values.reshape((-1, 1, 1))
        gaussian_dist.weights = gaussian_dist.weights.reshape((-1, 1, 1))
        self._b1_i_dist = gaussian_dist
        self.l_b1x_i = self._matrices.get("b1x_i", 0.0) * gaussian_dist.values
        self.l_b1y_i = self._matrices.get("b1y_i", 0.0) * gaussian_dist.values

    def set_b1_i_distribution(self, distribution: Distribution) -> None:
        """Set B1 inhomogeneity distribution directly.

        This method allows setting a custom B1 distribution, bypassing the
        default Gaussian distribution generated from b1_i_inh_scale and
        b1_i_inh_res.

        Parameters
        ----------
        distribution : Distribution
            B1 distribution with values and weights

        """
        distribution.values = distribution.values.reshape((-1, 1, 1))
        distribution.weights = distribution.weights.reshape((-1, 1, 1))
        self._b1_i_dist = distribution
        self._b1_i = float(distribution.values.mean())
        self.l_b1x_i = self._matrices.get("b1x_i", 0.0) * distribution.values
        self.l_b1y_i = self._matrices.get("b1y_i", 0.0) * distribution.values

    @property
    def b1_i_dist(self) -> Distribution:
        """Get the B1_i distribution object."""
        return self._b1_i_dist

    @property
    def b1_s(self) -> float:
        return self._b1_s

    @b1_s.setter
    def b1_s(self, value: float) -> None:
        self._b1_s = value
        self.l_b1x_s = self._matrices.get("b1x_s", 0.0) * value
        self.l_b1y_s = self._matrices.get("b1y_s", 0.0) * value

    @property
    def jeff_i(self) -> Distribution:
        return self._jeff_i

    @jeff_i.setter
    def jeff_i(self, distribution: Distribution) -> None:
        self._jeff_i = distribution
        values = distribution.values.reshape((-1, 1, 1, 1))
        self._l_jeff_i = self._matrices.get("jeff_i", np.array(0.0)) * values
        self._jeff_i_weights = distribution.weights.reshape((-1, 1, 1, 1))

    @property
    def gradient_dephasing(self) -> float:
        return self._gradient_dephasing

    @gradient_dephasing.setter
    def gradient_dephasing(self, value: float) -> None:
        self._gradient_dephasing = value
        self._scale_matrices("d_{state}", value * 1e-12)
        self._build_base_liouvillian()

    @property
    def l_free(self) -> Array:
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
    def weights(self) -> Array:
        return self._b1_i_dist.weights * self._jeff_i_weights

    @property
    def detection(self) -> str:
        return self._detection

    @detection.setter
    def detection(self, value: str) -> None:
        self._detection = value
        self._detect_vector = build_detection_vector(value, self.basis.vectors).transpose()

    def update(self, par_values: dict[str, float]) -> None:
        self.par_values = par_values
        self._build_base_liouvillian()

    def detect(self, magnetization: Array) -> float:
        return detect_signal(self._detect_vector, magnetization, self.weights)

    # def detect_spectrum(
    #     self, magnetization: Array, observed_state: str = "a"
    # ) -> Array:
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
        self, magnetization: Array, *, back: bool = False
    ) -> Array:
        basis = self.basis

        w1 = self.b1_i * 2.0 * np.pi

        for state in self.model.states:
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

    def get_equilibrium(self) -> Array:
        return build_equilibrium_magnetization(self.basis, self.par_values)

    def get_start_magnetization(
        self, terms: Iterable[str], atom: Nucleus = Nucleus.H1
    ) -> Array:
        return build_start_magnetization(self.basis, self.par_values, terms, atom)

    def keep(self, magnetization: Array, components: Iterable[str]) -> Array:
        return keep_components(self.basis, magnetization, components)

    def offsets_to_ppms(self, offsets: Array) -> Array:
        return self.carrier_i + 2.0 * np.pi * offsets / abs(self.ppm_i)

    def ppms_to_offsets(self, ppms: Array | float) -> Array | float:
        return (ppms - self.carrier_i) * abs(self.ppm_i) / (2.0 * np.pi)

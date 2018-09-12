""" The ref module contains the reference matrices and code for calculating the
    Liouvillian.

    Operator basis:
       {Ix, Iy, Iz, Sx, Sy, Sz,
        2IxSz, 2IySz, 2IzSx, 2IzSy,
        2IxSx, 2IxSy, 2IySx, 2IySy,
        2IzSz}
"""

import itertools

import numpy as np
from scipy import linalg, stats

from chemex.spindynamics import constants

COMPONENTS = {
    name: index
    for index, name in enumerate(
        (
            "ix",  # 0
            "iy",  # 1
            "iz",  # 2
            "sx",  # 3
            "sy",  # 4
            "sz",  # 5
            "2ixsz",  # 6
            "2iysz",  # 7
            "2izsx",  # 8
            "2izsy",  # 9
            "2ixsx",  # 10
            "2ixsy",  # 11
            "2iysx",  # 12
            "2iysy",  # 13
            "2izsz",  # 14
        )
    )
}

SPIN_SYSTEMS = {
    "iz": np.array([2]),
    "ixy": np.array([0, 1]),
    "ixyz": np.array([0, 1, 2]),
    "ixyzsz": np.array([0, 1, 2, 6, 7, 14]),
    "ixyzsxyz": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]),
    "ixysxy": np.array([10, 11, 12, 13]),
}

STATES_M = {
    "a": np.diag([1.0, 0.0, 0.0, 0.0]),
    "b": np.diag([0.0, 1.0, 0.0, 0.0]),
    "c": np.diag([0.0, 0.0, 1.0, 0.0]),
    "d": np.diag([0.0, 0.0, 0.0, 1.0]),
}

STATES_M["all"] = sum(STATES_M.values())
STATES_V = {state: np.diag(matrix).reshape(-1, 1) for state, matrix in STATES_M.items()}

POP_PAIRS = {
    component: [
        ("_".join((component, state)), "".join(("p", state)))
        for state in {"a", "b", "c", "d"}
    ]
    for component in COMPONENTS
}

N_SINGLE = len(COMPONENTS)
N_FULL = N_SINGLE * 4 + 1

ZEROS_M_SINGLE = np.zeros((N_SINGLE, N_SINGLE))
ZEROS_M_FULL = np.zeros((N_FULL, N_FULL))

ZEROS_V_SINGLE = np.zeros((N_SINGLE, 1))
ZEROS_V_FULL = np.zeros((N_FULL, 1))


def build_4st_is_spin_system():
    indexes = {
        "l": {
            "j": ((7, 9, 1, 4, 0, 3, 6, 8), (0, 3, 6, 8, 7, 9, 1, 4)),
            "w_i": ((1, 7, 12, 13, 0, 6, 10, 11), (0, 6, 10, 11, 1, 7, 12, 13)),
            "w_s": ((4, 9, 11, 13, 3, 8, 10, 12), (3, 8, 10, 12, 4, 9, 11, 13)),
            "w1x_i": ((2, 14, 8, 9, 1, 7, 12, 13), (1, 7, 12, 13, 2, 14, 8, 9)),
            "w1y_i": ((0, 10, 11, 6, 2, 8, 9, 14), (2, 8, 9, 14, 0, 10, 11, 6)),
            "w1x_s": ((5, 14, 6, 7, 4, 9, 11, 13), (4, 9, 11, 13, 5, 14, 6, 7)),
            "w1y_s": ((3, 10, 12, 8, 5, 6, 7, 14), (5, 6, 7, 14, 3, 10, 12, 8)),
        },
        "r": {
            "r2_i": ((0, 1), (0, 1)),
            "r1_i": ((2), (2)),
            "r2_s": ((3, 4), (3, 4)),
            "r1_s": ((5), (5)),
            "r2a_i": ((6, 7), (6, 7)),
            "r2a_s": ((8, 9), (8, 9)),
            "r2_mq": ((10, 11, 12, 13), (10, 11, 12, 13)),
            "r1a": ((14), (14)),
            "etaxy_i": ((6, 7, 0, 1), (0, 1, 6, 7)),
            "etaxy_s": ((8, 9, 3, 4), (8, 9, 3, 4)),
            "etaz_i": ((14, 2), (2, 14)),
            "etaz_s": ((14, 5), (5, 14)),
            "sigma": ((5, 0), (0, 5)),
            "mu_mq": ((10, 11, 12, 13), (13, 12, 11, 10)),
        },
        "k": {
            "kab": ((0, 1), (0)),
            "kba": ((1, 0), (1)),
            "kac": ((0, 2), (0)),
            "kca": ((2, 0), (2)),
            "kad": ((0, 3), (0)),
            "kda": ((3, 0), (3)),
            "kbc": ((1, 2), (1)),
            "kcb": ((2, 1), (2)),
            "kbd": ((1, 3), (1)),
            "kdb": ((3, 1), (3)),
            "kcd": ((2, 3), (2)),
            "kdc": ((3, 2), (3)),
        },
        "t": {"theta_i": (2), "theta_s": (5), "theta_is": (14)},
    }

    values = {
        "l": np.array([+1.0, +1.0, +1.0, +1.0, -1.0, -1.0, -1.0, -1.0]),
        "r": np.array([-1.0]),
        "k": np.array([-1.0, +1.0]),
        "t": np.array([+1.0]),
    }

    scales = {"j": np.pi}

    matrices = {}

    for kind in ("l", "r"):
        for name, ij in indexes[kind].items():
            matrix = ZEROS_M_SINGLE.copy()
            matrix[ij] = values[kind] * scales.get(name, 1.0)
            for state_name, state_matrix in STATES_M.items():
                if state_name != "all" or name.startswith(("w1", "w_")):
                    fname = "_".join([name, state_name])
                    fname = fname.replace("_all", "")
                    matrices[fname] = ZEROS_M_FULL.copy()
                    matrices[fname][:-1, :-1] = np.kron(state_matrix, matrix)

    for name, ij in indexes["k"].items():
        matrix = np.zeros((4, 4))
        matrix[ij] = values["k"]
        matrices[name] = ZEROS_M_FULL.copy()
        matrices[name][:-1, :-1] = np.kron(matrix, np.eye(15))

    for name, i in indexes["t"].items():
        vector = ZEROS_V_SINGLE.copy()
        vector[i] = values["t"]
        for state_name, state_vector in STATES_V.items():
            fname = "_".join([name, state_name])
            fname = fname.replace("_all", "")
            matrices[fname] = ZEROS_M_FULL.copy()
            matrices[fname][:-1, [-1]] = np.kron(state_vector, vector)

    for name in list(matrices):
        if name.startswith("w1") and name.endswith(("_a", "_b", "_c", "_d")):
            del matrices[name]

    vectors = {}

    for name, index in COMPONENTS.items():
        vector = ZEROS_V_SINGLE.copy()
        vector[index] = 1.0
        for state_name, state_vector in STATES_V.items():
            fname = "_".join([name, state_name])
            fname = fname.replace("_all", "")
            vectors[fname] = ZEROS_V_FULL.copy()
            vectors[fname][:-1, :] = np.kron(state_vector, vector)

    vectors["identity"] = ZEROS_V_FULL.copy()
    vectors["identity"][-1, 0] = +1.0

    return vectors, matrices


def build_basis(spin_system, state_nb, equilibrium=True):
    vectors, matrices = build_4st_is_spin_system()

    active_indexes = (
        SPIN_SYSTEMS[spin_system].reshape(1, -1)
        + N_SINGLE * np.arange(state_nb).reshape(-1, 1)
    ).reshape(-1)

    if equilibrium:
        active_indexes = np.concatenate((active_indexes, [-1]))

    vec_reduced = {}

    for name, vector in vectors.items():
        vector_ = vector[active_indexes]
        if vector_.any():
            vec_reduced[name] = vector_

    mesh = np.ix_(active_indexes, active_indexes)

    mat_reduced = {}

    for name, matrix in matrices.items():
        matrix_ = matrix[mesh]
        m_invalid = not matrix_.any()
        k_invalid = name.startswith("k") and matrix_.sum() != 0.0
        if not (m_invalid or k_invalid):
            mat_reduced[name] = matrix[mesh]

    return vec_reduced, mat_reduced


def add_cs_and_carrier(matrices, ppms):

    matrices_ = matrices.copy()

    for spin, state in itertools.product(("i", "s"), ("a", "b", "c", "d")):
        w_name = "_".join(("w", spin, state))
        cs_name = "_".join(("cs", spin, state))
        if w_name in matrices:
            matrices_[cs_name] = matrices[w_name] * ppms.get(spin, 1.0)

    for spin in ("i", "s"):
        w_name = "_".join(("w", spin))
        carrier_name = "_".join(("carrier", spin))
        j_eff_name = "_".join(("j_eff", spin))
        if w_name in matrices:
            matrices_[carrier_name] = (
                -1.0 * matrices_.get(w_name, 0.0) * ppms.get(spin, 1.0)
            )
            matrices_[j_eff_name] = np.pi * matrices_.get(w_name, 0.0)

    return matrices_


def calculate_propagators(liouvillian, delays, dephasing=False):
    """TODO: function docstring."""

    delays_ = np.asarray(delays).reshape(-1)
    shape = liouvillian.shape

    propagators = []

    for l in liouvillian.reshape(-1, *shape[-2:]):
        s, vr = linalg.eig(l)
        vri = linalg.inv(vr)

        if dephasing:
            sl = np.where(abs(s.imag) < 1e-6)[0]
            vr, s, vri = vr[:, sl], s[sl], vri[sl, :]

        d = np.asarray([np.diag(np.exp(s * t)) for t in delays_])

        propagators.append((vr @ d @ vri).real)

    propagators = np.asarray(propagators).swapaxes(0, 1).reshape(-1, *shape)

    return propagators


def make_perfect180(vectors):

    vect_size = list(vectors.values())[0].size

    perfect180 = {ptype: np.identity(vect_size) for ptype in ("ix", "iy", "sx", "sy")}

    for key in set(vectors) & set(COMPONENTS):
        if "ix" in key:
            perfect180["iy"] -= 2 * np.diag(vectors[key].reshape(-1))
        if "iy" in key:
            perfect180["ix"] -= 2 * np.diag(vectors[key].reshape(-1))
        if "iz" in key:
            perfect180["ix"] -= 2 * np.diag(vectors[key].reshape(-1))
            perfect180["iy"] -= 2 * np.diag(vectors[key].reshape(-1))
        if "sx" in key:
            perfect180["sy"] -= 2 * np.diag(vectors[key].reshape(-1))
        if "sy" in key:
            perfect180["sx"] -= 2 * np.diag(vectors[key].reshape(-1))
        if "sz" in key:
            perfect180["sx"] -= 2 * np.diag(vectors[key].reshape(-1))
            perfect180["sy"] -= 2 * np.diag(vectors[key].reshape(-1))

    return perfect180


class Liouvillian(object):
    """TODO"""

    def __init__(self, system, state_nb, atoms, h_larmor_frq, equilibrium=True):

        if atoms is None:
            atoms = {"i": "h"}

        self.ppms = {
            spin: 2.0 * np.pi * h_larmor_frq * constants.xi_ratio[atom]
            for spin, atom in atoms.items()
        }

        self._vectors, matrices = build_basis(system, state_nb, equilibrium)
        self._matrices = add_cs_and_carrier(matrices, self.ppms)

        self.keys = self._matrices.keys()

        self.detect = {name: vector.T for name, vector in self._vectors.items()}

        self.perfect180 = make_perfect180(self._vectors)

        self._w1_i_weights = 1.0
        self._w1_s_weights = 1.0
        self._w1_i_inh = 0.0
        self._w1_s_inh = 0.0
        self._w1_i_inh_res = 11
        self._w1_s_inh_res = 11

        self.update(cs_i_a=0.0)
        self.carrier_i = 0.0
        self.carrier_s = 0.0
        self.w1_i = 1e32
        self.w1_s = 1e32
        self.j_eff_i = 0.0
        self.j_eff_i_weights = 1.0

        zeros = np.zeros_like(self._matrices["kab"])
        self._rot90zp_i = linalg.expm(+0.5 * np.pi * self._matrices.get("w_i", zeros))
        self._rot90zm_i = linalg.expm(-0.5 * np.pi * self._matrices.get("w_i", zeros))
        self._rot90zp_s = linalg.expm(+0.5 * np.pi * self._matrices.get("w_s", zeros))
        self._rot90zm_s = linalg.expm(-0.5 * np.pi * self._matrices.get("w_s", zeros))

    @property
    def carrier_i(self):
        return self._carrier_i

    @carrier_i.setter
    def carrier_i(self, value):
        self._carrier_i = np.asarray(value)
        self._l_carrier_i = self._matrices.get("carrier_i", 0.0) * self._carrier_i

    @property
    def carrier_s(self):
        return self._carrier_s

    @carrier_s.setter
    def carrier_s(self, value):
        self._carrier_s = np.asarray(value)
        self._l_carrier_s = self._matrices.get("carrier_s", 0.0) * self._carrier_s

    @property
    def w1_i(self):
        return self._w1_i

    @w1_i.setter
    def w1_i(self, value):

        self._w1_i = value

        if self._w1_i_inh not in (0.0, np.inf) and self._w1_i_inh_res > 1:
            dist = np.linspace(-2.0, 2.0, self._w1_i_inh_res)
            w1_i_inh = self._w1_i_inh
        else:
            dist = np.array(0.0)
            w1_i_inh = 0.0

        dist = dist.reshape(-1, 1, 1)

        w1_i_dist = (dist * w1_i_inh + 1.0) * self._w1_i

        self._w1_i_weights = stats.norm.pdf(dist)
        self._w1_i_weights /= self._w1_i_weights.sum()

        self._l_w1x_i = self._matrices.get("w1x_i", 0.0) * w1_i_dist
        self._l_w1y_i = self._matrices.get("w1y_i", 0.0) * w1_i_dist

    @property
    def w1_i_inh(self):
        return self._w1_i_inh

    @w1_i_inh.setter
    def w1_i_inh(self, value):
        self._w1_i_inh = value
        self.w1_i = self._w1_i

    @property
    def w1_i_inh_res(self):
        return self._w1_i_inh_res

    @w1_i_inh_res.setter
    def w1_i_inh_res(self, value):
        self._w1_i_inh_res = value
        self.w1_i = self._w1_i

    @property
    def w1_s(self):
        return self._w1_s

    @w1_s.setter
    def w1_s(self, value):

        self._w1_s = value

        if self._w1_s_inh not in (0.0, np.inf) and self._w1_s_inh_res > 1:
            dist = np.linspace(-2.0, 2.0, self._w1_s_inh_res)
            w1_s_inh = self._w1_s_inh
        else:
            dist = np.array(0.0)
            w1_s_inh = 0.0

        dist = dist.reshape(-1, 1, 1)

        w1_s_dist = (dist * w1_s_inh + 1.0) * self._w1_s

        self._w1_s_weights = stats.norm.pdf(dist)
        self._w1_s_weights /= self._w1_s_weights.sum()

        self._l_w1x_s = self._matrices.get("w1x_s", 0.0) * w1_s_dist
        self._l_w1y_s = self._matrices.get("w1y_s", 0.0) * w1_s_dist

    @property
    def w1_s_inh(self):
        return self._w1_s_inh

    @w1_s_inh.setter
    def w1_s_inh(self, value):
        self._w1_s_inh = value
        self.w1_s = self._w1_s

    @property
    def w1_s_inh_res(self):
        return self._w1_s_inh_res

    @w1_s_inh_res.setter
    def w1_s_inh_res(self, value):
        self._w1_s_inh_res = value
        self.w1_s = self._w1_s

    @property
    def j_eff_i(self):
        return self._j_eff_i.reshape(-1)

    @j_eff_i.setter
    def j_eff_i(self, value):
        self._j_eff_i = np.asarray(value).reshape(-1, 1, 1, 1, 1)
        self._l_j_eff_i = self._matrices.get("j_eff_i", 0.0) * self._j_eff_i

    @property
    def j_eff_i_weights(self):
        return self._j_eff_i_weights.reshape(-1)

    @j_eff_i_weights.setter
    def j_eff_i_weights(self, value):
        self._j_eff_i_weights = np.asarray(value).reshape(-1, 1, 1, 1, 1)

    def compute_mag_eq(self, term="iz", **parvals):
        return sum(
            self._vectors.get(name1, 0.0) * parvals.get(name2, 0.0)
            for name1, name2 in POP_PAIRS[term]
        )

    def update(self, **parvals):
        self._l_free = sum(
            self._matrices[name] * parvals[name]
            for name in set(parvals) & set(self._matrices)
        )

    def collapse(self, vector):
        weights = (
            self._w1_i_weights * self._w1_s_weights * self._j_eff_i_weights
        ).reshape(-1, 1, 1)
        vector_ = vector.reshape(-1, *vector.shape[-2:])
        return (vector_ * weights).sum()

    def delays(self, times):
        liouv = self._l_free + self._l_carrier_i + self._l_carrier_s + self._l_j_eff_i
        return calculate_propagators(liouv, times)

    def pulse_i(self, times, phase, dephasing=False):
        l_w1_i = self._l_w1x_i * np.cos(phase * np.pi * 0.5) + self._l_w1y_i * np.sin(
            phase * np.pi * 0.5
        )

        liouv = (
            self._l_free
            + self._l_carrier_i
            + self._l_carrier_s
            + self._l_j_eff_i
            + l_w1_i
        )

        return calculate_propagators(liouv, times, dephasing)

    def pulses_90_180_i(self):
        pulses = {}
        liouv = (
            self._l_free
            + self._l_carrier_i
            + self._l_carrier_s
            + self._l_j_eff_i
            + self._l_w1x_i
        )
        t90 = 0.5 * np.pi / self._w1_i
        pulses["90px"] = calculate_propagators(liouv, t90)
        rot90zp, rot90zm = self._rot90zp_i, self._rot90zm_i
        pulses["90py"] = rot90zp @ pulses["90px"] @ rot90zm
        pulses["90mx"] = rot90zp @ pulses["90py"] @ rot90zm
        pulses["90my"] = rot90zp @ pulses["90mx"] @ rot90zm
        pulses["180px"] = pulses["90px"] @ pulses["90px"]
        pulses["180py"] = pulses["90py"] @ pulses["90py"]
        pulses["180mx"] = pulses["90mx"] @ pulses["90mx"]
        pulses["180my"] = pulses["90my"] @ pulses["90my"]
        return pulses

    def pulse_s(self, times, phase, dephasing=False):
        l_w1_s = self._l_w1x_s * np.cos(phase * np.pi * 0.5) + self._l_w1y_s * np.sin(
            phase * np.pi * 0.5
        )
        liouv = (
            self._l_free
            + self._l_carrier_i
            + self._l_carrier_s
            + self._l_j_eff_i
            + l_w1_s
        )
        return calculate_propagators(liouv, times, dephasing)

    def pulses_90_180_s(self):
        pulses = {}
        liouv = self._l_free + self._l_j_eff_i + self._l_w1x_s
        t90 = 0.5 * np.pi / self._w1_s
        rot90zp, rot90zm = self._rot90zp_s, self._rot90zm_s
        pulses["90px"] = calculate_propagators(liouv, t90)
        pulses["90py"] = rot90zp @ pulses["90px"] @ rot90zm
        pulses["90mx"] = rot90zp @ pulses["90py"] @ rot90zm
        pulses["90my"] = rot90zp @ pulses["90mx"] @ rot90zm
        pulses["180px"] = pulses["90px"] @ pulses["90px"]
        pulses["180py"] = pulses["90py"] @ pulses["90py"]
        pulses["180mx"] = pulses["90mx"] @ pulses["90mx"]
        pulses["180my"] = pulses["90my"] @ pulses["90my"]
        return pulses

    def pulse_is(self, times, phase_i, phase_s, dephasing=False):
        l_w1_i = self._l_w1x_i * np.cos(phase_i * np.pi * 0.5) + self._l_w1y_i * np.sin(
            phase_i * np.pi * 0.5
        )
        l_w1_s = self._l_w1x_s * np.cos(phase_s * np.pi * 0.5) + self._l_w1y_s * np.sin(
            phase_s * np.pi * 0.5
        )
        liouv = (
            self._l_free
            + self._l_carrier_i
            + self._l_carrier_s
            + self._l_j_eff_i
            + l_w1_i
            + l_w1_s
        )
        return calculate_propagators(liouv, times, dephasing)

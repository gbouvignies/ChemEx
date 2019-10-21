""" The ref module contains the reference matrices and code for calculating the
    Liouvillian.

    Operator basis:
       {Eq,
        Ix, Iy, Iz, Sx, Sy, Sz,
        2IxSz, 2IySz, 2IzSx, 2IzSy,
        2IxSx, 2IxSy, 2IySx, 2IySy,
        2IzSz}
"""
import itertools as it
import re

import asteval
import numpy as np
import scipy.stats as ss

import chemex.nmr.constants as cnc
import chemex.nmr.helper as cnh


class LiouvillianIS:
    BASES = {
        "ixy": ["ix", "iy"],
        "iz": ["iz"],
        "iz_eq": ["ie", "iz"],
        "ixyz": ["ix", "iy", "iz"],
        "ixyz_eq": ["ie", "ix", "iy", "iz"],
        "ixysxy": ["2ixsx", "2ixsy", "2iysx", "2iysy"],
        "ixyzsz": ["ix", "iy", "iz", "2ixsz", "2iysz", "2izsz"],
        "ixyzsz_eq": ["ie", "ix", "iy", "iz", "2ixsz", "2iysz", "2izsz"],
        "ixyzsxyz": [
            "ix",
            "iy",
            "iz",
            "sx",
            "sy",
            "sz",
            "2ixsz",
            "2iysz",
            "2izsx",
            "2izsy",
            "2ixsx",
            "2ixsy",
            "2iysx",
            "2iysy",
            "2izsz",
        ],
        "ixyzsxyz_eq": [
            "ie",
            "se",
            "ix",
            "iy",
            "iz",
            "sx",
            "sy",
            "sz",
            "2ixsz",
            "2iysz",
            "2izsx",
            "2izsy",
            "2ixsx",
            "2ixsy",
            "2iysx",
            "2iysy",
            "2izsz",
        ],
        "ch3_1htq_grad": ["ix", "iy", "iz", "2ixsz", "2iysz", "2izsz"],
    }
    TRANSITIONS = {
        "r2_i_{state}": (("ix", "ix", -1.0), ("iy", "iy", -1.0)),
        "r2_s_{state}": (("sx", "sx", -1.0), ("sy", "sy", -1.0)),
        "r1_i_{state}": (("iz", "iz", -1.0), ("iz", "ie", +1.0)),
        "r1_s_{state}": (("sz", "sz", -1.0), ("sz", "se", +1.0)),
        "r2a_i_{state}": (("2ixsz", "2ixsz", -1.0), ("2iysz", "2iysz", -1.0)),
        "r2a_s_{state}": (("2izsx", "2izsx", -1.0), ("2izsy", "2izsy", -1.0)),
        "r2mq_is_{state}": (
            ("2ixsx", "2ixsx", -1.0),
            ("2ixsy", "2ixsy", -1.0),
            ("2iysx", "2iysx", -1.0),
            ("2iysy", "2iysy", -1.0),
        ),
        "r1a_is_{state}": (("2izsz", "2izsz", -1.0),),
        "etaxy_i_{state}": (
            ("ix", "2ixsz", -1.0),
            ("iy", "2iysz", -1.0),
            ("2ixsz", "ix", -1.0),
            ("2iysz", "iy", -1.0),
        ),
        "etaxy_s_{state}": (
            ("sx", "2izsx", -1.0),
            ("sy", "2izsy", -1.0),
            ("2izsx", "sx", -1.0),
            ("2izsy", "sy", -1.0),
        ),
        "etaz_i_{state}": (
            ("iz", "2izsz", -1.0),
            ("2izsz", "iz", -1.0),
            ("2izsz", "ie", +1.0),
        ),
        "etaz_s_{state}": (
            ("sz", "2izsz", -1.0),
            ("2izsz", "sz", -1.0),
            ("2izsz", "se", +1.0),
        ),
        "sigma_is_{state}": (
            ("iz", "sz", -1.0),
            ("sz", "iz", -1.0),
            ("sz", "ie", +1.0),
            ("iz", "se", +1.0),
        ),
        "mu_is_{state}": (
            ("2ixsx", "2iysy", +1.0),
            ("2ixsy", "2iysx", -1.0),
            ("2iysx", "2ixsy", -1.0),
            ("2iysy", "2ixsx", +1.0),
        ),
        "rotz_i": (
            ("ix", "iy", -1.0),
            ("iy", "ix", +1.0),
            ("2ixsx", "2iysx", -1.0),
            ("2iysx", "2ixsx", +1.0),
            ("2ixsy", "2iysy", -1.0),
            ("2iysy", "2ixsy", +1.0),
            ("2ixsz", "2iysz", -1.0),
            ("2iysz", "2ixsz", +1.0),
        ),
        "rotz_s": (
            ("sx", "sy", -1.0),
            ("sy", "sx", +1.0),
            ("2ixsx", "2ixsy", -1.0),
            ("2ixsy", "2ixsx", +1.0),
            ("2iysx", "2iysy", -1.0),
            ("2iysy", "2iysx", +1.0),
            ("2izsx", "2izsy", -1.0),
            ("2izsy", "2izsx", +1.0),
        ),
        "cs_i_{state}": (
            ("ix", "iy", -1.0),
            ("iy", "ix", +1.0),
            ("2ixsx", "2iysx", -1.0),
            ("2iysx", "2ixsx", +1.0),
            ("2ixsy", "2iysy", -1.0),
            ("2iysy", "2ixsy", +1.0),
            ("2ixsz", "2iysz", -1.0),
            ("2iysz", "2ixsz", +1.0),
        ),
        "cs_s_{state}": (
            ("sx", "sy", -1.0),
            ("sy", "sx", +1.0),
            ("2ixsx", "2ixsy", -1.0),
            ("2ixsy", "2ixsx", +1.0),
            ("2iysx", "2iysy", -1.0),
            ("2iysy", "2iysx", +1.0),
            ("2izsx", "2izsy", -1.0),
            ("2izsy", "2izsx", +1.0),
        ),
        "carrier_i": (
            ("ix", "iy", +1.0),
            ("iy", "ix", -1.0),
            ("2ixsx", "2iysx", +1.0),
            ("2iysx", "2ixsx", -1.0),
            ("2ixsy", "2iysy", +1.0),
            ("2iysy", "2ixsy", -1.0),
            ("2ixsz", "2iysz", +1.0),
            ("2iysz", "2ixsz", -1.0),
        ),
        "carrier_s": (
            ("sx", "sy", +1.0),
            ("sy", "sx", -1.0),
            ("2ixsx", "2ixsy", +1.0),
            ("2ixsy", "2ixsx", -1.0),
            ("2iysx", "2iysy", +1.0),
            ("2iysy", "2iysx", -1.0),
            ("2izsx", "2izsy", +1.0),
            ("2izsy", "2izsx", -1.0),
        ),
        "offset_i": (
            ("ix", "iy", +2.0 * np.pi),
            ("iy", "ix", -2.0 * np.pi),
            ("2ixsx", "2iysx", +2.0 * np.pi),
            ("2iysx", "2ixsx", -2.0 * np.pi),
            ("2ixsy", "2iysy", +2.0 * np.pi),
            ("2iysy", "2ixsy", -2.0 * np.pi),
            ("2ixsz", "2iysz", +2.0 * np.pi),
            ("2iysz", "2ixsz", -2.0 * np.pi),
        ),
        "offset_s": (
            ("sx", "sy", +2.0 * np.pi),
            ("sy", "sx", -2.0 * np.pi),
            ("2ixsx", "2ixsy", +2.0 * np.pi),
            ("2ixsy", "2ixsx", -2.0 * np.pi),
            ("2iysx", "2iysy", +2.0 * np.pi),
            ("2iysy", "2iysx", -2.0 * np.pi),
            ("2izsx", "2izsy", +2.0 * np.pi),
            ("2izsy", "2izsx", -2.0 * np.pi),
        ),
        "jeff_i": (
            ("ix", "iy", -2.0 * np.pi),
            ("iy", "ix", +2.0 * np.pi),
            ("2ixsx", "2iysx", -2.0 * np.pi),
            ("2iysx", "2ixsx", +2.0 * np.pi),
            ("2ixsy", "2iysy", -2.0 * np.pi),
            ("2iysy", "2ixsy", +2.0 * np.pi),
            ("2ixsz", "2iysz", -2.0 * np.pi),
            ("2iysz", "2ixsz", +2.0 * np.pi),
        ),
        "j_is_{state}": (
            ("ix", "2iysz", -np.pi),
            ("2iysz", "ix", +np.pi),
            ("2ixsz", "iy", -np.pi),
            ("iy", "2ixsz", +np.pi),
            ("sx", "2izsy", -np.pi),
            ("2izsy", "sx", +np.pi),
            ("2izsx", "sy", -np.pi),
            ("sy", "2izsx", +np.pi),
        ),
        "b1x_i": (
            ("iy", "iz", -2.0 * np.pi),
            ("iz", "iy", +2.0 * np.pi),
            ("2iysx", "2izsx", -2.0 * np.pi),
            ("2izsx", "2iysx", +2.0 * np.pi),
            ("2iysy", "2izsy", -2.0 * np.pi),
            ("2izsy", "2iysy", +2.0 * np.pi),
            ("2iysz", "2izsz", -2.0 * np.pi),
            ("2izsz", "2iysz", +2.0 * np.pi),
        ),
        "b1y_i": (
            ("iz", "ix", -2.0 * np.pi),
            ("ix", "iz", +2.0 * np.pi),
            ("2izsx", "2ixsx", -2.0 * np.pi),
            ("2ixsx", "2izsx", +2.0 * np.pi),
            ("2izsy", "2ixsy", -2.0 * np.pi),
            ("2ixsy", "2izsy", +2.0 * np.pi),
            ("2izsz", "2ixsz", -2.0 * np.pi),
            ("2ixsz", "2izsz", +2.0 * np.pi),
        ),
        "b1x_s": (
            ("sy", "sz", -2.0 * np.pi),
            ("sz", "sy", +2.0 * np.pi),
            ("2ixsy", "2ixsz", -2.0 * np.pi),
            ("2ixsz", "2ixsy", +2.0 * np.pi),
            ("2iysy", "2iysz", -2.0 * np.pi),
            ("2iysz", "2iysy", +2.0 * np.pi),
            ("2izsy", "2izsz", -2.0 * np.pi),
            ("2izsz", "2izsy", +2.0 * np.pi),
        ),
        "b1y_s": (
            ("sz", "sx", -2.0 * np.pi),
            ("sx", "sz", +2.0 * np.pi),
            ("2ixsz", "2ixsx", -2.0 * np.pi),
            ("2ixsx", "2ixsz", +2.0 * np.pi),
            ("2iysz", "2iysx", -2.0 * np.pi),
            ("2iysx", "2iysz", +2.0 * np.pi),
            ("2izsz", "2izsx", -2.0 * np.pi),
            ("2izsx", "2izsz", +2.0 * np.pi),
        ),
    }

    def __init__(self, name, model, atoms, h_frq):
        self._b1_i_inh_scale = 0.0
        self._b1_i_inh_res = 11
        self._ppm_i = self._ppm_s = 1.0
        self._parvals = None
        self._l_base = None
        self.name = name
        self.state_nb = model.state_nb
        self.atoms = atoms
        self.h_frq = h_frq
        self.states = cnh.get_state_names(model.state_nb)
        self.basis = self.BASES[name]
        self.size = len(self.basis) * model.state_nb
        self.matrices = {
            **self._build_transition_matrices(),
            **self._build_exchange_matrices(),
        }
        self._matrices_ref = self.matrices.copy()
        self.vectors = self._build_component_vectors()
        self._interpreter = asteval.Interpreter(
            symtable={"vectors": self.vectors}, minimal=True
        )
        self.ppm_i = -2.0 * np.pi * h_frq * cnc.SIGNED_XI_RATIO.get(atoms.get("i"), 1.0)
        self.ppm_s = -2.0 * np.pi * h_frq * cnc.SIGNED_XI_RATIO.get(atoms.get("s"), 1.0)
        self.carrier_i, self.carrier_s = 0.0, 0.0
        self.offset_i, self.offset_s = 0.0, 0.0
        self.b1_i, self.b1_s = 1e32, 1e32
        self.jeff_i = cnc.Distribution(0.0, 1.0)
        self.detection = ""
        self.update((("cs_i_a", 0.0), ("pa", 1.0)))

    def update(self, parvals):
        self._parvals = parvals
        self._l_base = sum(
            self.matrices.get(name, 0.0) * value for name, value in parvals
        )

    @property
    def ppm_i(self):
        return self._ppm_i

    @ppm_i.setter
    def ppm_i(self, value):
        self._ppm_i = value
        self._scale_matrix("cs_i_{state}", value)
        self._scale_matrix("carrier_i", value)

    @property
    def ppm_s(self):
        return self._ppm_s

    @ppm_s.setter
    def ppm_s(self, value):
        self._ppm_s = value
        self._scale_matrix("cs_s_{state}", value)
        self._scale_matrix("carrier_s", value)

    @property
    def carrier_i(self):
        return self._carrier_i

    @carrier_i.setter
    def carrier_i(self, value):
        self._carrier_i = np.asarray(value)
        self._l_carrier_i = self.matrices.get("carrier_i", 0.0) * value

    @property
    def carrier_s(self):
        return self._carrier_s

    @carrier_s.setter
    def carrier_s(self, value):
        self._carrier_s = np.asarray(value)
        self._l_carrier_s = self.matrices.get("carrier_s", 0.0) * value

    @property
    def offset_i(self):
        return self._offset_i

    @offset_i.setter
    def offset_i(self, value):
        self._offset_i = np.asarray(value)
        self._l_offset_i = (
            self.matrices.get("offset_i", 0.0) * value * np.sign(self.ppm_i)
        )

    @property
    def offset_s(self):
        return self._offset_s

    @offset_s.setter
    def offset_s(self, value):
        self._offset_s = np.asarray(value)
        self._l_offset_s = (
            self.matrices.get("offset_s", 0.0) * value * np.sign(self.ppm_s)
        )

    @property
    def b1_i(self):
        return self._b1_i

    @b1_i.setter
    def b1_i(self, value):
        self._b1_i = value
        self._b1_i_dist, self._b1_i_weights = make_gaussian(
            self.b1_i, self.b1_i_inh_scale, self.b1_i_inh_res
        )
        self._b1_i_dist = self._b1_i_dist.reshape((-1, 1, 1))
        self._b1_i_weights = self._b1_i_weights.reshape((-1, 1, 1))
        self.l_b1x_i = self.matrices.get("b1x_i", 0.0) * self._b1_i_dist
        self.l_b1y_i = self.matrices.get("b1y_i", 0.0) * self._b1_i_dist

    @property
    def b1_i_inh_scale(self):
        return self._b1_i_inh_scale

    @b1_i_inh_scale.setter
    def b1_i_inh_scale(self, value):
        self._b1_i_inh_scale = value
        self.b1_i = self._b1_i

    @property
    def b1_i_inh_res(self):
        return self._b1_i_inh_res

    @b1_i_inh_res.setter
    def b1_i_inh_res(self, value):
        self._b1_i_inh_res = value
        self.b1_i = self._b1_i

    @property
    def b1_s(self):
        return self._b1_s

    @b1_s.setter
    def b1_s(self, value):
        self._b1_s = value
        self.l_b1x_s = self.matrices.get("b1x_s", 0.0) * value
        self.l_b1y_s = self.matrices.get("b1y_s", 0.0) * value

    @property
    def jeff_i(self):
        return self._jeff_i.reshape(-1)

    @jeff_i.setter
    def jeff_i(self, value):
        self._jeff_i = value
        values = np.asarray(value.values).reshape((-1, 1, 1, 1))
        self._l_jeff_i = self.matrices.get("jeff_i", 0.0) * values
        self._jeff_i_weights = np.asarray(value.weights).reshape((-1, 1, 1, 1))

    @property
    def l_free(self):
        return sum(
            [
                self._l_base,
                self._l_offset_i,
                self._l_offset_s,
                self._l_carrier_i,
                self._l_carrier_s,
                self._l_jeff_i,
            ]
        )

    @property
    def weights(self):
        return self._b1_i_weights * self._jeff_i_weights

    @property
    def detection(self):
        return self._detection

    @detection.setter
    def detection(self, value):
        self._detection = value
        expr = re.sub(r"(\w+\b(?<!\b1j))", r'vectors["\1"].reshape(1, -1)', value)
        self._detect_vector = self._interpreter(expr)

    def detect(self, magnetization):
        shape = -1, *magnetization.shape[-2:]
        mag_weighted = self.weights * magnetization
        mag = mag_weighted.reshape(shape).sum(axis=0)
        detected = self._detect_vector @ mag
        if np.iscomplexobj(detected):
            detected = np.sign(detected.real) * np.abs(detected)
        return np.float64(detected)

    def get_equilibrium(self):
        parvals = dict(self._parvals)
        mag = np.zeros(self.size)
        for state, (name, atom) in it.product(self.states, self.atoms.items()):
            scale = parvals.get(f"p{state}", 0.0) * cnc.XI_RATIO.get(atom, 1.0)
            mag += self.vectors.get(f"{name}e_{state}", 0.0) * scale
            mag += self.vectors.get(f"{name}z_{state}", 0.0) * scale
        return mag.reshape(-1, 1)

    def get_start_magnetization(self, terms=None, atom=None):
        if terms is None:
            terms = []
        elif isinstance(terms, str):
            terms = [terms]
        parvals = dict(self._parvals)
        ratio = cnc.XI_RATIO.get(atom, 1.0)
        mag = np.zeros(self.size)
        for term in terms:
            for state, (comp, vector) in it.product(self.states, self.vectors.items()):
                if comp.startswith(term) and comp.endswith(f"_{state}"):
                    mag += parvals[f"p{state}"] * ratio * vector
        return mag.reshape(-1, 1)

    def keep_components(self, magnetization, terms=None):
        if terms is None:
            terms = []
        elif isinstance(terms, str):
            terms = [terms]
        keep = np.zeros((self.size, 1))
        for term in terms:
            for state, (comp, vector) in it.product(self.states, self.vectors.items()):
                if comp.startswith(term) and comp.endswith(f"_{state}"):
                    keep += vector.reshape(-1, 1)
        return keep * magnetization

    def offsets_to_ppms(self, offsets):
        return self.carrier_i + 2.0 * np.pi * offsets / abs(self.ppm_i)

    def ppms_to_offsets(self, ppms):
        return (ppms - self.carrier_i) * abs(self.ppm_i) / (2.0 * np.pi)

    def _build_transition_matrices(self):
        m_zeros = np.zeros((self.size, self.size))
        matrices = {}
        for transition_name, state in it.product(self.TRANSITIONS, self.states):
            name_ = transition_name.format(state=state)
            indices, values = self._get_indices(transition_name, state)
            if values:
                matrices.setdefault(name_, m_zeros.copy())
                matrices[name_][indices] = values
        return matrices

    def _build_exchange_matrices(self):
        matrices = {}
        for (i1, s1), (i2, s2) in it.permutations(enumerate(self.states), r=2):
            name = f"k{s1}{s2}"
            matrix = np.zeros((self.state_nb, self.state_nb))
            matrix[((i1, i2), i1)] = -1.0, 1.0
            matrices[name] = np.kron(matrix, np.eye(len(self.basis)))
        return matrices

    def _build_component_vectors(self):
        names = [f"{name}_{state}" for state in self.states for name in self.basis]
        vectors = dict(zip(names, np.eye(self.size)))
        for state, name in it.product(self.states, self.basis):
            vectors.setdefault(name, 0.0)
            vectors[name] += vectors.get(f"{name}_{state}", 0.0)
        return vectors

    def _get_indices(self, transition_name, state):
        rows = []
        cols = []
        vals = []
        offset = self.states.index(state) * len(self.basis)
        for start, end, value in self.TRANSITIONS[transition_name]:
            if {start, end}.issubset(self.basis):
                rows.append(self.basis.index(start) + offset)
                cols.append(self.basis.index(end) + offset)
                vals.append(value)
        return (rows, cols), vals

    def _scale_matrix(self, name, value):
        names = {name.format(state=state) for state in self.states}
        for name_ in names & set(self.matrices):
            self.matrices[name_] = np.sign(self._matrices_ref[name_]) * value


class Liouvillian1HTQDif(LiouvillianIS):
    TRANSITIONS = dict(
        **LiouvillianIS.TRANSITIONS,
        **{
            "d_{state}": (
                ("ix", "ix", -1.0),
                ("iy", "iy", -1.0),
                ("iz", "iz", -1.0),
                ("sx", "sx", -1.0),
                ("sy", "sy", -1.0),
                ("sz", "sz", -1.0),
                ("2ixsz", "2ixsz", -1.0),
                ("2iysz", "2iysz", -1.0),
                ("2izsx", "2izsx", -1.0),
                ("2izsy", "2izsy", -1.0),
                ("2ixsx", "2ixsx", -1.0),
                ("2ixsy", "2ixsy", -1.0),
                ("2iysx", "2iysx", -1.0),
                ("2iysy", "2iysy", -1.0),
                ("2izsz", "2izsz", -1.0),
            )
        },
    )

    def __init__(self, name, model, atoms, h_frq):
        if "ch3_1htq" not in name:
            print(
                "'Liouvillian1HTQDif' is only compatible with the 'ch3_1htq_grad' "
                "basis. Basis changed to 'ch3_1htq_grad'"
            )
            name = "ch3_1htq_grad"
        super().__init__(name, model, atoms, h_frq)
        self._scale_matrix("j_is_{state}", 3.0 * np.pi)
        self.gradient_dephasing = 0.0

    @property
    def ppm_i(self):
        return self._ppm_i

    @ppm_i.setter
    def ppm_i(self, value):
        self._ppm_i = value
        self._scale_matrix("cs_i_{state}", 3.0 * value)
        self._scale_matrix("carrier_i", 3.0 * value)

    @property
    def gradient_dephasing(self):
        return self._gradient_dephasing

    @gradient_dephasing.setter
    def gradient_dephasing(self, value):
        self._gradient_dephasing = value
        self._scale_matrix("d_{state}", value)
        self.update(self._parvals)


def make_gaussian(value, scale, res):
    if scale not in (0.0, np.inf) and res > 1:
        grid = np.linspace(-2.0, 2.0, res)
        dist = grid * scale + 1.0
    else:
        grid = np.array(0.0)
        dist = np.array(1.0)
    weights = ss.norm.pdf(grid)
    weights /= weights.sum()
    return dist * value, weights

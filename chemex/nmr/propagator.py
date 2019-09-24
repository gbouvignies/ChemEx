import numpy as np
from scipy import linalg

from chemex.nmr import liouvillian


class PropagatorIS:
    LIOUV = liouvillian.LiouvillianIS

    def __init__(self, basis, model, atoms, h_frq):
        self._pw90_i = self._pw90_s = None
        self._p90_i = self._p180_i = self._p240_i = None
        self._p90_s = self._p180_s = self._p240_s = None
        self.liouvillian = self.LIOUV(basis, model, atoms, h_frq)
        size = self.liouvillian.size
        self.identity = np.eye(self.liouvillian.size).reshape((1, 1, size, size))
        self.ppm_i = self.liouvillian.ppm_i
        self.ppm_s = self.liouvillian.ppm_s
        self._phases = self._get_phases()
        self.perfect90_i = self._make_perfect90("i")
        self.perfect90_s = self._make_perfect90("s")
        self.perfect180_i = self._make_perfect180("i")
        self.perfect180_s = self._make_perfect180("s")

    def update(self, parvals):
        self.liouvillian.update(parvals)
        self._p90_i = self._p90_s = None

    @property
    def carrier_i(self):
        return self.liouvillian.carrier_i

    @carrier_i.setter
    def carrier_i(self, value):
        self.liouvillian.carrier_i = value
        self._p90_i = self._p90_s = None

    @property
    def carrier_s(self):
        return self.liouvillian.carrier_s

    @carrier_s.setter
    def carrier_s(self, value):
        self.liouvillian.carrier_s = value
        self._p90_i = self._p90_s = None

    @property
    def offset_i(self):
        return self.liouvillian.offset_i

    @offset_i.setter
    def offset_i(self, value):
        self.liouvillian.offset_i = value
        self._p90_i = self._p90_s = None

    @property
    def offset_s(self):
        return self.liouvillian.offset_s

    @offset_s.setter
    def offset_s(self, value):
        self.liouvillian.offset_s = value
        self._p90_i = self._p90_s = None

    @property
    def b1_i(self):
        return self.liouvillian.b1_i

    @b1_i.setter
    def b1_i(self, value):
        self._pw90_i = 1.0 / (4.0 * value) if value else 0.0
        self.liouvillian.b1_i = value
        self._p90_i = self._p90_s = None

    @property
    def b1_s(self):
        return self.liouvillian.b1_s

    @b1_s.setter
    def b1_s(self, value):
        self._pw90_s = 1.0 / (4.0 * value) if value else 0.0
        self.liouvillian.b1_s = value
        self._p90_i = self._p90_s = None

    @property
    def b1_i_inh_scale(self):
        return self.liouvillian.b1_i_inh_scale

    @b1_i_inh_scale.setter
    def b1_i_inh_scale(self, value):
        self.liouvillian.b1_i_inh_scale = value
        self._p90_i = self._p90_s = None

    @property
    def b1_i_inh_res(self):
        return self.liouvillian.b1_i_inh_res

    @b1_i_inh_res.setter
    def b1_i_inh_res(self, value):
        self.liouvillian.b1_i_inh_res = value
        self._p90_i = self._p90_s = None

    @property
    def jeff_i(self):
        return self.liouvillian.jeff_i

    @jeff_i.setter
    def jeff_i(self, value):
        self.liouvillian.jeff_i = value
        self._p90_i = self._p90_s = None

    def get_equilibrium(self):
        return self.liouvillian.get_equilibrium()

    def get_start_magnetization(self, terms=None, atom=None):
        return self.liouvillian.get_start_magnetization(terms=terms, atom=atom)

    def keep_components(self, magnetization, terms=None):
        return self.liouvillian.keep_components(magnetization, terms)

    @property
    def detection(self):
        return self.liouvillian.detection

    @detection.setter
    def detection(self, value):
        self.liouvillian.detection = value

    def detect(self, magnetization):
        return self.liouvillian.detect(magnetization)

    def delays(self, times):
        return calculate_propagators(self.liouvillian.l_free, times)

    def pulse_i(self, times, phase, dephasing=False):
        liouv = (
            self.liouvillian.l_free
            + np.cos(phase * np.pi * 0.5) * self.liouvillian.l_b1x_i
            + np.sin(phase * np.pi * 0.5) * self.liouvillian.l_b1y_i
        )
        return calculate_propagators(liouv, times, dephasing)

    def pulse_s(self, times, phase, dephasing=False):
        liouv = (
            self.liouvillian.l_free
            + np.cos(phase * np.pi * 0.5) * self.liouvillian.l_b1x_s
            + np.sin(phase * np.pi * 0.5) * self.liouvillian.l_b1y_s
        )
        return calculate_propagators(liouv, times, dephasing)

    def pulse_is(self, times, phase_i, phase_s, dephasing=False):
        liouv = (
            self.liouvillian.l_free
            + np.cos(phase_i * np.pi * 0.5) * self.liouvillian.l_b1x_i
            + np.sin(phase_i * np.pi * 0.5) * self.liouvillian.l_b1y_i
            + np.cos(phase_s * np.pi * 0.5) * self.liouvillian.l_b1x_s
            + np.sin(phase_s * np.pi * 0.5) * self.liouvillian.l_b1y_s
        )
        return calculate_propagators(liouv, times, dephasing)

    @property
    def p90_i(self):
        self._calculate_base_pulses_i()
        return self._p90_i

    @property
    def p180_i(self):
        self._calculate_base_pulses_i()
        return self._p180_i

    @property
    def p240_i(self):
        self._calculate_base_pulses_i()
        return self._p240_i

    @property
    def p9018090_i_1(self):
        return self.p90_i[[3, 0, 1, 2]] @ self.p180_i @ self.p90_i[[3, 0, 1, 2]]

    @property
    def p9018090_i_2(self):
        return self.p90_i[[1, 2, 3, 0]] @ self.p180_i @ self.p90_i[[1, 2, 3, 0]]

    @property
    def p9024090_i_1(self):
        return self.p90_i[[3, 0, 1, 2]] @ self.p240_i @ self.p90_i[[3, 0, 1, 2]]

    @property
    def p9024090_i_2(self):
        return self.p90_i[[1, 2, 3, 0]] @ self.p240_i @ self.p90_i[[1, 2, 3, 0]]

    @property
    def p90_s(self):
        self._calculate_base_pulses_s()
        return self._p90_s

    @property
    def p180_s(self):
        self._calculate_base_pulses_s()
        return self._p180_s

    @property
    def p240_s(self):
        self._calculate_base_pulses_s()
        return self._p240_s

    @property
    def p9018090_s_1(self):
        return self.p90_s[3, 0, 1, 2] @ self.p180_s @ self.p90_s[3, 0, 1, 3]

    @property
    def p9018090_s_2(self):
        return self.p90_s[1, 2, 3, 0] @ self.p180_s @ self.p90_s[1, 2, 3, 0]

    @property
    def p9024090_s_1(self):
        return self.p90_s[3, 0, 1, 2] @ self.p240_s @ self.p90_s[3, 0, 1, 3]

    @property
    def p9024090_s_2(self):
        return self.p90_s[1, 2, 3, 0] @ self.p240_s @ self.p90_s[1, 2, 3, 0]

    @property
    def p9024090_is(self):
        pw240_i, pw9024090_i = np.array([8.0, 14.0]) * self._pw90_i / 3.0
        pw240_s, pw9024090_s = np.array([8.0, 14.0]) * self._pw90_s / 3.0
        times = np.sort([0.0, pw240_i, pw240_s, pw9024090_i, pw9024090_s])
        deltas = (times[1:] - times[:-1]) * np.array([1.0, 0.5, 0.5, 0.5])
        phase_i, phase_s = (0, 3), (0, 3)
        p0 = self.pulse_is(deltas[0], phase_i[0], phase_s[0])
        if pw9024090_i <= pw9024090_s:
            p1 = self.pulse_is(deltas[1], phase_i[1], phase_s[0])
            if pw9024090_i <= pw240_s:
                p2 = self.pulse_s(deltas[2], phase_s[0])
            else:
                p2 = self.pulse_is(deltas[2], phase_i[1], phase_s[1])
            p3 = self.pulse_s(deltas[3], phase_s[1])
        else:
            p1 = self.pulse_is(deltas[1], phase_i[0], phase_s[1])
            if pw9024090_s <= pw240_i:
                p2 = self.pulse_i(deltas[2], phase_i[0])
            else:
                p2 = self.pulse_is(deltas[2], phase_i[1], phase_s[1])
            p3 = self.pulse_i(deltas[3], phase_i[1])
        p_xx = p3 @ p2 @ p1 @ p0 @ p1 @ p2 @ p3
        phases_i, phases_s = self._phases["i"], self._phases["s"]
        p_comp_is = np.array(
            [
                [
                    phases_i[k] @ phases_s[j] @ p_xx @ phases_s[-j] @ phases_i[-k]
                    for j in range(4)
                ]
                for k in range(4)
            ]
        )
        return p_comp_is

    def offsets_to_ppms(self, offsets):
        return self.liouvillian.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.liouvillian.ppms_to_offsets(ppms)

    def _calculate_base_pulses_i(self):
        if self._p90_i is None:
            phases = self._phases["i"]
            pws = np.array([1.0, 2.0, 8.0 / 3.0]) / (4.0 * self.liouvillian.b1_i)
            pulses_ = self.pulse_i(pws, 0.0)
            pulses = np.array([phases[j] @ pulses_ @ phases[-j] for j in range(4)])
            self._p90_i, self._p180_i, self._p240_i = pulses.swapaxes(0, 1)

    def _calculate_base_pulses_s(self):
        if self._p90_s is None:
            phases = self._phases["s"]
            pws = np.array([1.0, 2.0, 8.0 / 3.0]) / (4.0 * self.liouvillian.b1_s)
            pulses_ = self.pulse_s(pws, 0.0)
            pulses = np.array([phases[j] @ pulses_ @ phases[-j] for j in range(4)])
            self._p90_s, self._p180_s, self._p240_s = pulses.swapaxes(0, 1)

    def _make_perfect180(self, spin):
        compx, compy, compz = (f"{spin}{axis}" for axis in "xyz")
        perfect180 = {comp: self.identity.copy() for comp in (compx, compy)}
        for comp in self.liouvillian.basis:
            vect = self.liouvillian.vectors[comp]
            if compx in comp or compz in comp:
                perfect180[compy] -= 2 * np.diag(vect)
            if compy in comp or compz in comp:
                perfect180[compx] -= 2 * np.diag(vect)
        p180 = [perfect180[comp] for comp in (compx, compy, compx, compy)]
        return np.array(p180)

    def _make_perfect90(self, spin):
        phases = self._phases[spin]
        zeros = np.zeros((self.liouvillian.size, self.liouvillian.size))
        rot = self.liouvillian.matrices.get(f"b1x_{spin}", zeros)
        base = linalg.expm(0.25 * rot)
        p90 = np.array([phases[i] @ base @ phases[-i] for i in range(4)])
        return p90

    def _get_phases(self):
        phases = {}
        zeros = np.zeros((self.liouvillian.size, self.liouvillian.size))
        for spin in "is":
            l_rotz = self.liouvillian.matrices.get(f"rotz_{spin}", zeros)
            phases[spin] = np.array(
                [linalg.expm(n * 0.5 * np.pi * l_rotz) for n in range(4)]
            )
        return phases


class Propagator1HTQDif(PropagatorIS):
    LIOUV = liouvillian.Liouvillian1HTQDif

    @property
    def gradient_dephasing(self):
        return self.liouvillian.gradient_dephasing

    @gradient_dephasing.setter
    def gradient_dephasing(self, value):
        self.liouvillian.gradient_dephasing = value
        self._p90_i = self._p90_s = None


def calculate_propagators(liouv, delays, dephasing=False):
    delays_ = np.asarray(delays).reshape(-1)
    shape = liouv.shape
    propagators = []
    for a_liouvillian in liouv.reshape(-1, *shape[-2:]):
        s, vr = linalg.eig(a_liouvillian)
        vri = linalg.inv(vr)
        if dephasing:
            sl = np.where(abs(s.imag) < 1e-6)[0]
            vr, s, vri = vr[:, sl], s[sl], vri[sl, :]
        d = np.asarray([np.diag(np.exp(s * t)) for t in delays_])
        propagators.append((vr @ d @ vri).real)
    propagators = np.asarray(propagators).swapaxes(0, 1).reshape(-1, *shape)
    if propagators.shape[0] == 1:
        propagators = propagators[0]
    return propagators

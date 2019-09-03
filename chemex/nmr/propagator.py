import numpy as np
from scipy import linalg

from chemex.nmr import liouvillian


class PropagatorIS:
    LIOUV = liouvillian.LiouvillianIS

    def __init__(self, basis, model, atoms, h_frq):
        self.liouvillian = self.LIOUV(basis, model, atoms, h_frq)
        self.identity = np.eye(self.liouvillian.size)
        self.ppm_i = self.liouvillian.ppm_i
        self.ppm_s = self.liouvillian.ppm_s
        self._phases = self._get_phases()
        self.perfect90_i = self._make_perfect90("i")
        self.perfect90_s = self._make_perfect90("s")
        self.perfect180_i = self._make_perfect180("i")
        self.perfect180_s = self._make_perfect180("s")

    def update(self, parvals):
        return self.liouvillian.update(parvals)

    @property
    def carrier_i(self):
        return self.liouvillian.carrier_i

    @carrier_i.setter
    def carrier_i(self, value):
        self.liouvillian.carrier_i = value

    @property
    def carrier_s(self):
        return self.liouvillian.carrier_s

    @carrier_s.setter
    def carrier_s(self, value):
        self.liouvillian.carrier_s = value

    @property
    def offset_i(self):
        return self.liouvillian.offset_i

    @offset_i.setter
    def offset_i(self, value):
        self.liouvillian.offset_i = value

    @property
    def offset_s(self):
        return self.liouvillian.offset_s

    @offset_s.setter
    def offset_s(self, value):
        self.liouvillian.offset_s = value

    @property
    def b1_i(self):
        return self.liouvillian.b1_i

    @b1_i.setter
    def b1_i(self, value):
        self.liouvillian.b1_i = value

    @property
    def b1_s(self):
        return self.liouvillian.b1_s

    @b1_s.setter
    def b1_s(self, value):
        self.liouvillian.b1_s = value

    @property
    def b1_i_inh_scale(self):
        return self.liouvillian.b1_i_inh_scale

    @b1_i_inh_scale.setter
    def b1_i_inh_scale(self, value):
        self.liouvillian.b1_i_inh_scale = value

    @property
    def b1_i_inh_res(self):
        return self.liouvillian.b1_i_inh_res

    @b1_i_inh_res.setter
    def b1_i_inh_res(self, value):
        self.liouvillian.b1_i_inh_res = value

    @property
    def jeff_i(self):
        return self.liouvillian.jeff_i

    @jeff_i.setter
    def jeff_i(self, value):
        self.liouvillian.jeff_i = value

    def get_equilibrium(self):
        return self.liouvillian.get_equilibrium()

    def get_start_magnetization(self, terms=None, atom=None):
        return self.liouvillian.get_start_magnetization(terms=terms, atom=atom)

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

    def pulses_90_180_i(self):
        pw = 1.0 / (4.0 * self.liouvillian.b1_i)
        phases = self._phases["i"]
        base = self.pulse_i(pw, 0.0)
        p90 = np.array([phases[i] @ base @ phases[-i] for i in range(4)])
        p180 = np.array([pulse @ pulse for pulse in p90])
        return p90, p180

    def pulses_90_180_s(self):
        pw = 1.0 / (4.0 * self.liouvillian.b1_s)
        phases = self._phases["s"]
        base = self.pulse_i(pw, 0.0)
        p90 = np.array([phases[i] @ base @ phases[-i] for i in range(4)])
        p180 = np.array([pulse @ pulse for pulse in p90])
        return p90, p180

    def offsets_to_ppms(self, offsets):
        return self.liouvillian.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.liouvillian.ppms_to_offsets(ppms)

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
        return p180

    def _make_perfect90(self, spin):
        phases = self._phases[spin]
        zeros = self.identity * 0.0
        rot = self.liouvillian.matrices.get(f"b1x_{spin}", zeros)
        base = linalg.expm(0.25 * rot)
        p90 = np.array([phases[i] @ base @ phases[-i] for i in range(4)])
        return p90

    def _get_phases(self):
        phases = {}
        zeros = self.identity * 0.0
        for spin in "is":
            l_rotz = self.liouvillian.matrices.get(f"rotz_{spin}", zeros)
            phases[spin] = [linalg.expm(n * 0.5 * np.pi * l_rotz) for n in range(4)]
        return phases


class Propagator1HTQDif(PropagatorIS):
    LIOUV = liouvillian.Liouvillian1HTQDif

    @property
    def gradient_dephasing(self):
        return self.liouvillian.gradient_dephasing

    @gradient_dephasing.setter
    def gradient_dephasing(self, value):
        self.liouvillian.gradient_dephasing = value


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
    return propagators

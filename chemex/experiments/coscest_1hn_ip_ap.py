"""
1HN In-phase/Anti-phase Proton COS-CEST
=======================================

Analyzes chemical exchange during the COS-CEST block. Magnetization evolution
is calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states::

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }


References
----------

Yuwen, Bouvignies and Kay. J Mag Reson (2018) 292:1-7


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config coscest_1hn_ip_ap

"""
import functools as ft

import numpy as np
import numpy.linalg as nl

import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.liouvillian as cnl


_SCHEMA = {
    "type": "object",
    "properties": {
        "experiment": {
            "type": "object",
            "properties": {
                "d1": {"type": "number"},
                "time_t1": {"type": "number"},
                "carrier": {"type": "number"},
                "sw": {"type": "number"},
                "cos_n": {"type": "integer"},
                "cos_res": {"type": "integer", "default": 10},
                "b1_frq": {"type": "number"},
                "b1_inh_scale": {"type": "number", "default": 0.1},
                "b1_inh_res": {"type": "integer", "default": 11},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
            },
            "required": ["d1", "time_t1", "carrier", "b1_frq", "sw", "cos_n"],
        }
    },
}


def read(config):
    ch.validate(config, _SCHEMA)
    config["basis"] = cnl.Basis(type="ixyzsz_eq", spin_system="hn")
    config["fit"] = _fit_this(config)
    config["data"]["filter_ref_planes"] = True
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this(config):
    return {
        "rates": [
            "r2_i_{states}",
            "r1_i_{observed_state}",
            "r1_s_{observed_state}",
            "etaxy_i_{observed_state}",
            "etaz_i_{observed_state}",
        ],
        "model_free": [
            "tauc_{observed_state}",
            "s2_{observed_state}",
            "khh_{observed_state}",
        ],
    }


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t1 = settings["time_t1"]
        self.d1 = settings["d1"]
        self.taua = 2.38e-3
        self.sw = settings["sw"]
        self.cos_n = settings["cos_n"]
        self.cos_res = settings["cos_res"]
        self.prop.carrier_i = settings["carrier"]
        self.prop.b1_i = settings["b1_frq"]
        self.prop.b1_i_inh_scale = settings["b1_inh_scale"]
        self.prop.b1_i_inh_res = settings["b1_inh_res"]
        self.observed_state = settings["observed_state"]
        self.prop.detection = f"2izsz_{self.observed_state}"
        self.dephased = settings["b1_inh_scale"] == np.inf
        self.taud = self.d1
        self.p90_i = self.prop.perfect90_i
        self.p180_sx = self.prop.perfect180_s[0]
        self.p180_isx = self.prop.perfect180_i[0] @ self.prop.perfect180_s[0]
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, offsets, params_local):
        self.prop.update(params_local)
        self.prop.offset_i = 0.0
        d_taud, d_taua = self.prop.delays([self.taud, self.taua])
        start = d_taud @ self.prop.get_start_magnetization(terms="ie")
        start = self.prop.keep_components(start, terms=["ie", "iz"])
        intst = {}
        for offset in set(offsets):
            self.prop.offset_i = offset
            cest = self._calc_cosine_shape()
            if abs(offset) >= 1e4:
                inept = self.p90_i[3] @ d_taua @ self.p180_isx @ d_taua @ self.p90_i[0]
                cest = inept @ cest
            intst[offset] = self.prop.detect(cest @ start)
        return np.array([intst[offset] for offset in offsets])

    def offsets_to_ppms(self, offsets):
        return self.prop.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.prop.ppms_to_offsets(ppms)

    def _calc_cosine_shape(self):
        dt = 1.0 / (self.cos_res * self.sw)
        n_periods = int(self.time_t1 * self.sw)
        n_left = int((self.time_t1 * self.sw - n_periods) * self.cos_res)
        grid = np.linspace(-np.pi, np.pi, self.cos_res, endpoint=False)
        n_values = (np.arange(self.cos_n) - 0.5 * (self.cos_n - 1)).reshape(-1, 1)
        amplitudes = np.cos(n_values * grid).sum(axis=0)
        phases = np.zeros(self.cos_res)
        pulse = self.prop.shaped_pulse_i(self.cos_res * dt, amplitudes, phases)
        pulse = nl.matrix_power(pulse, n_periods)
        if n_left:
            pulse_left = self.prop.shaped_pulse_i(n_left * dt, amplitudes[:n_left], 0.0)
            pulse = pulse_left @ pulse
        return pulse
"""
15N - N-H TROSY CEST
====================

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)x(6n), single spin matrix, where n is the number of
states:

[ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
  Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b),
  ... ]


Reference
---------

Dong, Bouvignies and Kay. PNAS (2014) 24:8820-5


Note
----

A sample configuration file for this module is available using the command:

    chemex config cest_15N_tr

"""
import functools as ft

import numpy as np

import chemex.containers.cest as ccc
import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.propagator as cnp
import chemex.nmr.rates as cnr


_SCHEMA = {
    "type": "object",
    "properties": {
        "experiment": {
            "type": "object",
            "properties": {
                "time_t1": {"type": "number"},
                "carrier": {"type": "number"},
                "b1_frq": {"type": "number"},
                "b1_inh_scale": {"type": "number", "default": 0.1},
                "b1_inh_res": {"type": "integer", "default": 11},
                "antitrosy": {"type": "boolean", "default": False},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
            },
            "required": ["time_t1", "carrier", "b1_frq"],
        }
    },
}
_FIT_SETTING = {
    "dw_ab": "fit",
    "r1_a": "fit",
    "r2_a": "fit",
    "r2_b": "fit",
    "r1a_a": "fit",
    "etaxy_a": "fit",
}


def read(config):
    config["spin_system"] = {
        "basis": "ixyzsz",
        "atoms": {"i": "n", "s": "h"},
        "constraints": ["nh"],
        "rates": "nh",
    }
    ch.validate(config, _SCHEMA)
    ch.validate(config, ccc.CEST_SCHEMA)
    if config["experiment"]["antitrosy"]:
        _FIT_SETTING["etaz_a"] = "fit"
    experiment = ceh.read(
        config=config,
        pulse_seq_cls=PulseSeq,
        propagator_cls=cnp.PropagatorIS,
        container_cls=ccc.CestProfile,
        rates_cls=cnr.RatesIS,
        fit_setting=_FIT_SETTING,
    )
    return experiment


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t1 = settings["time_t1"]
        self.prop.carrier_i = settings["carrier"]
        self.prop.b1_i = settings["b1_frq"]
        self.prop.b1_i_inh_scale = settings["b1_inh_scale"]
        self.prop.b1_i_inh_res = settings["b1_inh_res"]
        self.antitrosy = settings["antitrosy"]
        self.observed_state = settings["observed_state"]
        self.prop.detection = self._get_detection(settings["observed_state"])
        self.dephased = settings["b1_inh_scale"] == np.inf
        self.calculate = ft.lru_cache(maxsize=5)(self.calculate_)

    def calculate_(self, offsets, params_local):
        self.prop.update(params_local)
        start = self._get_start()
        ref = abs(offsets) < 1e4
        intst = {offset: self.prop.detect(start) for offset in offsets[ref]}
        for offset in set(offsets[~ref]):
            self.prop.offset_i = offset
            intst[offset] = self.prop.detect(
                self.prop.pulse_i(self.time_t1, 0.0, self.dephased) @ start
            )
        return np.array([intst[offset] for offset in offsets])

    def _get_start(self):
        start = self.prop.get_start_magnetization("2izsz")
        if self.antitrosy:
            start += self.prop.get_start_magnetization("iz")
        else:
            start -= self.prop.get_start_magnetization("iz")
        start *= 0.5
        return start

    def _get_detection(self, state):
        if self.antitrosy:
            detection = f"2izsz_{state} + iz_{state}"
        else:
            detection = f"2izsz_{state} - iz_{state}"
        return detection

    def offsets_to_ppms(self, offsets):
        return self.prop.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.prop.ppms_to_offsets(ppms)

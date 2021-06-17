"""
15N–1HN TROSY CEST
==================

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states::

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

References
----------

Long, Bouvignies and Kay. PNAS (2014) 111:8820-8825


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config cest_15N_tr

"""
import functools as ft

import numpy as np

import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.liouvillian as cnl


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


def read(config):
    ch.validate(config, _SCHEMA)
    config["basis"] = cnl.Basis(type="ixyzsz", spin_system="nh")
    config["fit"] = _fit_this(config)
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this(config):
    this = {
        "rates": ["r2_i_{states}", "r1_i_{observed_state}", "r1a_is_{states}"],
        "model_free": ["tauc_{observed_state}", "s2_{observed_state}", "khh_{states}"],
    }
    if config["experiment"]["antitrosy"]:
        this["rates"].extend(["etaxy_i_{observed_state}", "etaz_i_{observed_state}"])
    return this


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
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, offsets, params_local):
        self.prop.update(params_local)
        start = self._get_start()
        intst = {}
        for offset in set(offsets):
            if abs(offset) <= 1e4:
                self.prop.offset_i = offset
                intst[offset] = (
                    self.prop.pulse_i(self.time_t1, 0.0, self.dephased) @ start
                )
            else:
                intst[offset] = start
        return np.array([self.prop.detect(intst[offset]) for offset in offsets])

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
            detection = f"[2izsz_{state}] + [iz_{state}]"
        else:
            detection = f"[2izsz_{state}] - [iz_{state}]"
        return detection

    def offsets_to_ppms(self, offsets):
        return self.prop.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.prop.ppms_to_offsets(ppms)

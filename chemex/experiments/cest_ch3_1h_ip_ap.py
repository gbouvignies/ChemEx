"""
1H (Methyl-13CH3/13CHD2) In-phase/Anti-phase Proton CEST
========================================================

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states::

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

References
----------

| Yuwen and Kay. J Biomol NMR (2017) 68:215-224
| Yuwen and Kay. J Biomol NMR (2018) 70:93-102


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config cest_ch3_1h_ip_ap

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
                "d1": {"type": "number"},
                "time_t1": {"type": "number"},
                "carrier": {"type": "number"},
                "b1_frq": {"type": "number"},
                "b1_inh_scale": {"type": "number", "default": 0.1},
                "b1_inh_res": {"type": "integer", "default": 11},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
            },
            "required": ["d1", "time_t1", "carrier", "b1_frq"],
        }
    },
}


def read(config):
    ch.validate(config, _SCHEMA)
    config["basis"] = cnl.Basis(type="ixyzsz_eq", spin_system="hc")
    config["fit"] = _fit_this()
    config["data"]["filter_ref_planes"] = True
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this():
    return {
        "rates": [
            "r2_i_{states}",
            "r1_i_{observed_state}",
            "r1_s_{observed_state}",
            "etaxy_i_{observed_state}",
            "etaz_i_{observed_state}",
        ],
        "model_free": ["tauc_{observed_state}", "s2_{observed_state}"],
    }


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t1 = settings["time_t1"]
        self.d1 = settings["d1"]
        self.taua = 2.00e-3
        self.prop.carrier_i = settings["carrier"]
        self.prop.b1_i = settings["b1_frq"]
        self.prop.b1_i_inh_scale = settings["b1_inh_scale"]
        self.prop.b1_i_inh_res = settings["b1_inh_res"]
        self.observed_state = settings["observed_state"]
        self.prop.detection = f"[2izsz_{self.observed_state}]"
        self.dephased = settings["b1_inh_scale"] == np.inf
        self.p90_i = self.prop.perfect90_i
        self.p180_isx = self.prop.perfect180_i[0] @ self.prop.perfect180_s[0]

    @ft.lru_cache(maxsize=10000)
    def calculate(self, offsets, params_local):
        self.prop.update(params_local)
        self.prop.offset_i = 0.0
        d_d1, d_taua = self.prop.delays([self.d1, self.taua])
        start = d_d1 @ self.prop.get_start_magnetization(terms="ie")
        start = self.prop.keep_components(start, terms=["ie", "iz"])
        intst = {}
        for offset in set(offsets):
            self.prop.offset_i = offset
            mag = self.prop.pulse_i(self.time_t1, 0.0, self.dephased) @ start
            if abs(offset) > 1e4:
                inept = self.p90_i[3] @ d_taua @ self.p180_isx @ d_taua @ self.p90_i[0]
                mag = inept @ mag
            intst[offset] = self.prop.detect(mag)
        return np.array([intst[offset] for offset in offsets])

    def offsets_to_ppms(self, offsets):
        return self.prop.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.prop.ppms_to_offsets(ppms)

"""
Pure anti-phase 1HN CEST
========================

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)x(6n), single spin matrix, where n is the number of
states:

[ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
  Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b),
  ... ]


Reference
---------
Sekhar, Rosenzweig, Bouvignies and Kay. PNAS (2016) 113:E2794-E2801


Note
----

A sample configuration  file for this module is available using the command:

    chemex config cest_1hn_ip_ap

"""
import functools as ft

import numpy as np

import chemex.containers.cest as ccc
import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.propagator as cnp
import chemex.nmr.rates as cnr


TYPE = __name__.split(".")[-1]
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
    config["spin_system"] = {
        "basis": "ixyzsz_eq",
        "atoms": {"i": "h", "s": "n"},
        "constraints": ["hn"],
        "rates": "hn",
    }
    ch.validate(config, _SCHEMA)
    ch.validate(config, ccc.CEST_SCHEMA)
    experiment = ceh.read(
        config=config,
        pulse_seq_cls=PulseSeq,
        propagator_cls=cnp.PropagatorIS,
        container_cls=ccc.CestProfile,
        rates_cls=cnr.RatesIS,
    )
    return experiment


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t1 = settings["time_t1"]
        self.d1 = settings["d1"]
        self.taua = 1.0 / (4.0 * 105.0)
        self.prop.carrier_i = settings["carrier"]
        self.prop.b1_i = settings["b1_frq"]
        self.prop.b1_i_inh_scale = settings["b1_inh_scale"]
        self.prop.b1_i_inh_res = settings["b1_inh_res"]
        self.observed_state = settings["observed_state"]
        self.dephased = settings["b1_inh_scale"] == np.inf
        self.calculate = ft.lru_cache(maxsize=5)(self.calculate_)

    def calculate_(self, offsets, params_local):
        self.prop.update(params_local)
        d_d1, d_t1, d_taua = self.prop.delays([self.d1, self.time_t1, self.taua])
        start = d_d1 @ self.prop.get_start_magnetization(terms=f"ie")
        intst = {}
        for offset in set(offsets):
            if abs(offset) < 1e4:
                self.prop.detection = f"2izsz_{self.observed_state}"
                self.prop.offset_i = offset
                intst[offset] = self.prop.detect(
                    self.prop.pulse_i(self.time_t1, 0.0, self.dephased) @ start
                )
            else:
                p90 = self.prop.perfect90_i
                p180_ix = self.prop.perfect180_i[0]
                p180_sx = self.prop.perfect180_s[0]
                inept = p90[3] @ d_taua @ p180_sx @ p180_ix @ d_taua @ p90[0]
                intst[offset] = self.prop.detect(inept @ d_t1 @ start)
        return np.array([intst[offset] for offset in offsets])

    def offsets_to_ppms(self, offsets):
        return self.prop.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.prop.ppms_to_offsets(ppms)

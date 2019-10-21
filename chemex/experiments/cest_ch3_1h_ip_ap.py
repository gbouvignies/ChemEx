"""
In-phase/Anti-phase CH3/CHD2 1H CEST
========================

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)x(6n), single spin matrix, where n is the number of
states:

[ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
  Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b),
  ... ]


Reference
---------

Yuwen and Kay. J Biomol NMR (2017) 68:215-224
Yuwen and Kay. J Biomol NMR (2018) 70:93-102


Note
----

A sample configuration file for this module is available using the command:

    chemex config cest_ch3_1h_ip_ap

"""
import functools as ft

import numpy as np

import chemex.containers.cest as ccc
import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.propagator as cnp


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
_FIT_SETTING = {
    "dw_ab": "fit",
    "r1a_a": "fit",
    "r1_a, nuc->c": "fit",
    "r2_a": "fit",
    "r2_b": "fit",
    "etaxy_a": "fit",
    "etaz_a": "fit",
}


def read(config):
    config["spin_system"] = {
        "basis": "ixyzsz_eq",
        "atoms": {"i": "h", "s": "c"},
        "constraints": ["hn"],
    }
    config["data"]["filter_ref_planes"] = True
    ch.validate(config, _SCHEMA)
    ch.validate(config, ccc.CEST_SCHEMA)
    experiment = ceh.read(
        config=config,
        pulse_seq_cls=PulseSeq,
        propagator_cls=cnp.PropagatorIS,
        container_cls=ccc.CestProfile,
        fit_setting=_FIT_SETTING,
    )
    return experiment


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
        self.prop.detection = f"2izsz_{self.observed_state}"
        self.dephased = settings["b1_inh_scale"] == np.inf
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)
        self.p90_i = self.prop.perfect90_i
        self.p180_isx = self.prop.perfect180_i[0] @ self.prop.perfect180_s[0]

    def _calculate(self, offsets, params_local):
        self.prop.update(params_local)
        self.prop.offset_i = 0.0
        d_d1, d_taua = self.prop.delays([self.d1, self.taua])
        start = d_d1 @ self.prop.get_start_magnetization(terms="ie")
        start = self.prop.keep_components(start, terms=["ie", "iz"])
        intst = {}
        for offset in set(offsets):
            self.prop.offset_i = offset
            mag = self.prop.pulse_i(self.time_t1, 0.0, self.dephased) @ start
            if abs(offset) >= 1e4:
                inept = self.p90_i[3] @ d_taua @ self.p180_isx @ d_taua @ self.p90_i[0]
                mag = inept @ mag
            intst[offset] = self.prop.detect(mag)
        return np.array([intst[offset] for offset in offsets])

    def offsets_to_ppms(self, offsets):
        return self.prop.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.prop.ppms_to_offsets(ppms)

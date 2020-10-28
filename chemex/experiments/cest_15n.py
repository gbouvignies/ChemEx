"""
15N Pure In-phase CEST
======================

Analyzes chemical exchange in the presence of 1H composite decoupling during
the CEST block. This keeps the spin system purely in-phase throughout, and is
calculated using the (3n)×(3n), single-spin matrix, where n is the number of
states::

    { Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... }

References
----------

Vallurupalli, Bouvignies and Kay. J Am Chem Soc (2012) 134:8148-8161


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config cest_15n

"""
import functools as ft

import numpy as np

import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.constants as cnc
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
    config["basis"] = cnl.Basis(type="ixyz", spin_system="nh")
    config["fit"] = _fit_this(config)
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this(config):
    return {
        "rates": ["r2_i_{states}", "r1_i_{observed_state}"],
        "model_free": ["tauc_{observed_state}", "s2_{observed_state}"],
    }


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t1 = settings["time_t1"]
        self.prop.carrier_i = settings["carrier"]
        self.prop.b1_i = settings["b1_frq"]
        self.prop.b1_i_inh_scale = settings["b1_inh_scale"]
        self.prop.b1_i_inh_res = settings["b1_inh_res"]
        if "13C" in config["conditions"].label:
            self.prop.jeff_i = cnc.get_multiplet("", "n")
        self.prop.detection = f"iz_{settings['observed_state']}"
        self.dephased = settings["b1_inh_scale"] == np.inf
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, offsets, params_local):
        self.prop.update(params_local)
        start = self.prop.get_equilibrium()
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

    def offsets_to_ppms(self, offsets):
        return self.prop.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.prop.ppms_to_offsets(ppms)

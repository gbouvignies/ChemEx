"""
13C Pure In-phase D-CEST
========================

Analyzes chemical exchange in the presence of 1H composite decoupling during
the D-CEST block. This keeps the spin system purely in-phase throughout, and
is calculated using the (3n)×(3n), single-spin matrix, where n is the number
of states::

    { Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... }

References
----------

| Yuwen, Kay and Bouvignies. ChemPhysChem (2018) 19:1707-1710
| Yuwen, Bouvignies and Kay. J Mag Reson (2018) 292:1-7


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config dcest_13c

"""
import functools as ft

import numpy as np
import numpy.linalg as la

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
                "pw90": {"type": "number"},
                "sw": {"type": "number"},
                "time_equil": {"type": "number", "default": 0.0},
                "b1_inh_scale": {"type": "number", "default": 0.1},
                "b1_inh_res": {"type": "integer", "default": 11},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
            },
            "required": ["time_t1", "carrier", "b1_frq", "pw90", "sw"],
        }
    },
}


def read(config):
    ch.validate(config, _SCHEMA)
    config["basis"] = cnl.Basis(type="ixyz", spin_system="ch")
    config["fit"] = _fit_this()
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this():
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
        self.sw = settings["sw"]
        self.pw90 = 4.0 * settings["pw90"] * settings["b1_frq"] / settings["sw"]
        self.time_eq = settings["time_equil"]
        self.tau_dante = 1.0 / settings["sw"] - self.pw90
        self.ncyc = int(settings["time_t1"] * settings["sw"] + 0.1)
        self.prop.b1_i = 1.0 / (4.0 * settings["pw90"])
        self.prop.b1_i_inh_scale = settings["b1_inh_scale"]
        self.prop.b1_i_inh_res = settings["b1_inh_res"]
        if "13C" in config["conditions"].label:
            spin_system = config["spin_system"]
            symbol = spin_system.symbols["i"]
            nucleus = spin_system.nuclei["i"]
            self.prop.jeff_i = cnc.get_multiplet(symbol, nucleus)
        self.observed_state = settings["observed_state"]
        self.prop.detection = f"iz_{self.observed_state}"
        self.dephased = settings["b1_inh_scale"] == np.inf
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, offsets, params_local):
        self.prop.update(params_local)
        start = self.prop.get_equilibrium()
        intst = {}
        d_eq = (
            self.prop.delays(self.time_eq) if self.time_eq > 0.0 else self.prop.identity
        )
        for offset in set(offsets):
            if abs(offset) <= 1e4:
                self.prop.offset_i = offset
                p_delay = self.prop.delays(self.tau_dante)
                p_pulse = self.prop.pulse_i(self.pw90, 0.0)
                intst[offset] = (
                    d_eq @ la.matrix_power(p_delay @ p_pulse, self.ncyc) @ start
                )
            else:
                intst[offset] = d_eq @ start
        return np.array([self.prop.detect(intst[offset]) for offset in offsets])

    def offsets_to_ppms(self, offsets):
        return self.prop.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.prop.ppms_to_offsets(ppms)

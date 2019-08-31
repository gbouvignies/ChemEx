"""
Pure In-phase 15N D-CEST
========================

Analyzes chemical exchange in the presence of 1H composite decoupling during
the D-CEST block. This keeps the spin system purely in-phase throughout, and
is calculated using the (3n)x(3n), single spin matrix, where n is the number
of states:

[ Ix(a), Iy(a), Iz(a),
  Ix(b), Iy(b), Iz(b),
   ... ]

Reference
---------

Yuwen, Kay and Bouvignies. ChemPhysChem (2018) 19:1707-1710


Note
----

A sample configuration  file for this module is available using the command:

    chemex config dcest_15n

"""
import functools as ft

import numpy as np
import numpy.linalg as la

import chemex.containers.cest as ccc
import chemex.experiments.helper as ceh
import chemex.nmr.constants as cnc


TYPE = __name__.split(".")[-1]
_SCHEMA = {
    "type": "object",
    "properties": {
        "experiment": {
            "type": "object",
            "properties": {
                "time_t1": {"type": "number"},
                "carrier": {"type": "number"},
                "b1_frq": {"type": "number"},
                "pw90_dante": {"type": "number"},
                "sw_dante": {"type": "number"},
                "b1_inh_scale": {"type": "number", "default": 0.1},
                "b1_inh_res": {"type": "integer", "default": 11},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
                "hd_exchange": {"type": "boolean", "default": False},
                "cn_label": {"type": "boolean", "default": False},
            },
            "required": ["time_t1", "carrier", "b1_frq", "pw90_dante", "sw_dante"],
        }
    },
}


def read(config):
    config["spin_system"] = {"basis": "ixyz", "atoms": {"i": "n"}}
    ceh.validate(config, _SCHEMA)
    ceh.validate(config, ccc.CEST_SCHEMA)
    experiment = ceh.read(
        config=config,
        pulse_seq_cls=PulseSeq,
        propagator_cls=cnp.PropagatorIS,
        container_cls=ccc.CestProfile,
    )
    return experiment


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t1 = settings["time_t1"]
        self.prop.carrier_i = settings["carrier"]
        self.sw_dante = settings["sw_dante"]
        self.pw_dante = (
            4.0 * settings["pw90_dante"] * settings["b1_frq"] / settings["sw_dante"]
        )
        self.tau_dante = 1.0 / settings["sw_dante"] - self.pw_dante
        self.ncyc = int(settings["time_t1"] * settings["sw_dante"] + 0.1)
        self.prop.b1_i = 1.0 / (4.0 * settings["pw90_dante"])
        self.prop.b1_i_inh_scale = settings["b1_inh_scale"]
        self.prop.b1_i_inh_res = settings["b1_inh_res"]
        if settings["cn_label"]:
            self.prop.jeff_i = cnc.get_multiplet("", "n")
        self.prop.detection = f"iz_{settings['observed_state']}"
        self.hd_exchange = settings["hd_exchange"]
        self.dephased = settings["b1_inh_scale"] == np.inf
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, offsets, params_local):
        self.prop.update(params_local)
        if not self.hd_exchange:
            start = self.prop.get_equilibrium()
        else:
            start = self.prop.get_start_magnetization(["iz_a", "iz_c"])
        intst = {}
        for offset in set(offsets):
            if abs(offset) < 1e4:
                self.prop.offset_i = offset
                p_delay = self.prop.delays(self.tau_dante)
                p_pulse = self.prop.pulse_i(self.pw_dante, 0.0)
                intst[offset] = la.matrix_power(p_pulse @ p_delay, self.ncyc) @ start
            else:
                intst[offset] = start
        return np.array([self.prop.detect(intst[offset]) for offset in offsets])

    def offsets_to_ppms(self, offsets):
        return self.prop.offsets_to_ppms(offsets)

    def ppms_to_offsets(self, ppms):
        return self.prop.ppms_to_offsets(ppms)
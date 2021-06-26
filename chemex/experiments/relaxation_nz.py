"""
15N T1
======

Analyzes 15N T1 experiments. This keeps the spin system purely in-phase
throughout, and is calculated using the (1n)×(1n), single-spin matrix,
where n is the number of states::

    { Iz(a), Iz(b), ... }

References
----------

Kay, Nicholson, Delaglio, Bax, and Torchia. J Mag Reson (1992) 97:359-375


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config relaxation_nz

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
                "observed_state": {"type": "string", "pattern": "[a-z]", "default": "a"}
            },
        }
    },
}


def read(config):
    ch.validate(config, _SCHEMA)
    config["basis"] = cnl.Basis(type="iz", spin_system="nh")
    config["fit"] = _fit_this()
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this():
    return {
        "rates": ["r1_i_{observed_state}"],
        "model_free": ["tauc_{observed_state}", "s2_{observed_state}"],
    }


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.prop.detection = f"[iz_{settings['observed_state']}]"

    @ft.lru_cache(maxsize=10000)
    def calculate(self, times, params_local):
        self.prop.update(params_local)
        start = self.prop.get_equilibrium()
        delays = self.prop.delays(times)
        return np.array([self.prop.detect(delay @ start) for delay in delays])

"""
15N T1
======

Analyzes 15N T1 experiments. This keeps the spin system purely in-phase throughout, and
is calculated using the (1n)x(1n), single spin matrix, where n is the number
of states:

[ Iz(a), Iz(b), ... ]

Reference
---------

TODO

Note
----

A sample configuration file for this module is available using the command:

    chemex config relaxation_nz

"""
import functools as ft

import numpy as np

import chemex.experiments.helper as ceh
import chemex.helper as ch


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
    config["spin_system"] = {"basis": "iz", "atoms": {"i": "n"}, "rates": "nh"}
    ch.validate(config, _SCHEMA)
    fit_setting = {"r1_a": "fit"}
    experiment = ceh.load_experiment(
        config=config, pulse_seq_cls=PulseSeq, fit_setting=fit_setting
    )
    return experiment


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.prop.detection = f"iz_{settings['observed_state']}"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, times, params_local):
        self.prop.update(params_local)
        start = self.prop.get_equilibrium()
        delays = self.prop.delays(times)
        return np.array([self.prop.detect(delay @ start) for delay in delays])

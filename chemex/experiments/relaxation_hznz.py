"""
15N–1HN R1(2HzNz)
=================

Analyzes 1H-15N longitudinal two-spin order relaxation experiments. Decay
is calculated using the (2n)×(2n), two-spin matrix, where n is the number
of states::

    { Iz(a), 2IzSz(a),
      Iz(b), 2IzSz(b), ... }

References
----------

Hansen, Yang, Feng, Zhou, Wiesner, Bai and Kay. J Am Chem Soc (2007) 129:11468-11479


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config relaxation_hznz

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
    config["basis"] = cnl.Basis(type="izsz", spin_system="nh")
    config["fit"] = _fit_this()
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this():
    return {
        "rates": ["r1a_is_{observed_state}"],
        "model_free": [
            "tauc_{observed_state}",
            "s2_{observed_state}",
            "khh_{observed_state}",
        ],
    }


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.prop.detection = f"[2izsz_{settings['observed_state']}]"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, times, params_local):
        self.prop.update(params_local)
        start = self.prop.get_start_magnetization(["2izsz"])
        delays = self.prop.delays(0.25 * np.array(times))
        p180_i = self.prop.perfect180_i[0]
        p180_s = self.prop.perfect180_s[0]
        return np.array(
            [
                self.prop.detect(
                    delay @ p180_s @ delay @ p180_i @ delay @ p180_s @ delay @ start
                )
                for delay in delays
            ]
        )

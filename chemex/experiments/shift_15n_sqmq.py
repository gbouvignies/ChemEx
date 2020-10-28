"""
15N Exchange Induced Shifts with (15N–1HN) HSQC/HMQC
====================================================

Analyzes exchange induced 15N chemical shift changes measured in
(15N–1HN) HMQC and HSQC data sets.

References
----------

Vallurupalli, Bouvignies and Kay. J Phys Chem B (2011) 115:14891-14900


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config shift_15n_sqmq

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
    config["basis"] = cnl.Basis(type="ixy_ixysxy", spin_system="nh")
    config["fit"] = {"rates": {}, "model_free": {}}
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.ppm_i = self.prop.ppm_i
        self.ppm_s = self.prop.ppm_s
        self.cs_i_state = f"cs_i_{settings['observed_state']}"
        self.cs_s_state = f"cs_s_{settings['observed_state']}"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, params_local):
        self.prop.update(params_local)
        params_map = dict(params_local)
        ref_shift_i = params_map[self.cs_i_state] * self.ppm_i
        ref_shift_s = params_map[self.cs_s_state] * self.ppm_s
        ref_shift_dq = ref_shift_i + ref_shift_s
        ref_shift_zq = ref_shift_i - ref_shift_s
        shifts = self.prop.calculate_shifts()
        shift_sq = _find_nearest(shifts, ref_shift_i)
        shift_dq = _find_nearest(shifts, ref_shift_dq)
        shift_zq = _find_nearest(shifts, ref_shift_zq)
        return 1e3 * (shift_sq - 0.5 * (shift_dq + shift_zq)) / self.ppm_i


def _find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

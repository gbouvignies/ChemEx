"""
15N Pure In-phase CPMG
======================

Analyzes 15N chemical exchange in the presence of high power 1H CW decoupling
during the CPMG block. This keeps the spin system purely in-phase throughout,
and is calculated using the (3n)×(3n), single-spin matrix, where n is the
number of states::

    { Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... }

The CW decoupling on 1H is assumed to be strong enough (> 15 kHz) such that
perfect 1H decoupling can be achieved.

References
----------

Hansen, Vallurupalli and Kay. J Phys Chem B (2008) 112:5898-5904


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config cpmg_15n_ip

"""
import functools as ft

import numpy as np
import numpy.linalg as nl

import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.liouvillian as cnl


_SCHEMA = {
    "type": "object",
    "properties": {
        "experiment": {
            "type": "object",
            "properties": {
                "time_t2": {"type": "number"},
                "carrier": {"type": "number"},
                "pw90": {"type": "number"},
                "time_equil": {"type": "number", "default": 0.0},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
            },
            "required": ["time_t2", "carrier", "pw90"],
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
        "rates": ["r2_i_{observed_state}"],
        "model_free": ["tauc_{observed_state}"],
    }


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t2 = settings["time_t2"]
        self.time_eq = settings["time_equil"]
        self.prop.carrier_i = settings["carrier"]
        self.pw90 = settings["pw90"]
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.prop.detection = f"iz_{settings['observed_state']}"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p90 = self.prop.p90_i
        p180 = self.prop.p180_i
        p180pmx = 0.5 * (p180[0] + p180[2])  # +/- phase cycling

        # Getting the starting magnetization
        start = self.prop.get_equilibrium()

        # Calculating the instensities as a function of ncyc
        part1 = d_neg @ p90[0] @ start
        part2 = d_eq @ p90[0] @ d_neg
        intst = {
            0: self.prop.detect(d_eq @ p90[0] @ p180pmx @ p90[0] @ start),
            -1: self.prop.detect(part2 @ d_cp[-1] @ p180pmx @ d_cp[-1] @ part1),
        }
        for ncyc in set(ncycs) - {0, -1}:
            echo = d_cp[ncyc] @ p180[1] @ d_cp[ncyc]
            cpmg = nl.matrix_power(echo, int(ncyc))
            intst[ncyc] = self.prop.detect(part2 @ cpmg @ p180pmx @ cpmg @ part1)

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_) - self.pw90))
        tau_cps[-1] = 0.5 * self.time_t2
        delays = [self.t_neg, self.time_eq, *tau_cps.values()]
        return tau_cps, delays

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.array(ncycs, dtype=np.float)
        ncycs_[ncycs_ == -1.0] = 0.5
        return ncycs_[ncycs_ > 0.0] / self.time_t2

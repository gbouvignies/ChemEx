"""
15N - Pure In-phase Nitrogen CPMG with 0013 phase cycle
=======================================================

Analyzes 15N chemical exchange in the presence of high power 1H CW decoupling
during the CPMG block. This keeps the spin system purely in-phase throughout,
and is calculated using the (3n)x(3n), single spin matrix, where n is the number
of states:

[ Ix(a), Iy(a), Iz(a),
  Ix(b), Iy(b), Iz(b),
   ... ]

This version is modified such that CPMG pulses are applied with [0013] phase cycle
as shown in Daiwen's paper. The cw decoupling on 1H is assumed to be
strong enough (> 15 kHz) such that perfect 1H decoupling can be achieved.

References
----------

Jiang, Yu, Zhang, Liu and Yang. J Magn Reson (2015) 257:1-7
Hansen, Vallurupalli and Kay. J Phys Chem B (2008) 112:5898-5904


Note
----
A sample configuration file for this module is available using the command:

    chemex config cpmg_15n_ip_0013

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
                "time_t2": {"type": "number"},
                "carrier": {"type": "number"},
                "pw90": {"type": "number"},
                "time_equil": {"type": "number", "default": 0.0},
                "ncyc_max": {"type": "integer"},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
            },
            "required": ["time_t2", "carrier", "pw90", "ncyc_max"],
        }
    },
}
_FIT_SETTING = {"dw_ab": "fit", "r2_a": "fit"}


def read(config):
    config["spin_system"] = {"basis": "ixyz", "atoms": {"i": "n"}, "rates": "nh"}
    ch.validate(config, _SCHEMA)
    experiment = ceh.load_experiment(
        config=config, pulse_seq_cls=PulseSeq, fit_setting=_FIT_SETTING
    )
    return experiment


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t2 = settings["time_t2"]
        self.time_eq = settings["time_equil"]
        self.ncyc_max = settings["ncyc_max"]
        self.prop.carrier_i = settings["carrier"]
        self.pw90 = settings["pw90"]
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.t_pos2 = +4.0 * self.pw90 / np.pi
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.prop.detection = f"iz_{settings['observed_state']}"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps, deltas, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_neg = delays[self.t_neg]
        d_pos2 = delays[self.t_pos2]
        d_delta = {ncyc: delays[delay] for ncyc, delay in deltas.items()}
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p90 = self.prop.p90_i
        p180 = self.prop.p180_i

        # Getting the starting magnetization
        start = self.prop.get_equilibrium()

        # Calculating the instensities as a function of ncyc
        intst = {
            0: self.prop.detect(
                d_delta[0]
                @ p90[3]
                @ p180[[0, 3]]
                @ d_pos2
                @ p180[[0, 1]]
                @ p90[1]
                @ start
            )
        }
        for ncyc in set(ncycs) - {0}:
            phases = self._get_phases(ncyc)
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg = ft.reduce(np.matmul, echo[phases.T])
            intst[ncyc] = self.prop.detect(
                d_delta[ncyc] @ p90[3] @ d_neg @ cpmg @ d_neg @ p90[1] @ start
            )

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_) - 0.75 * self.pw90))
        deltas = dict(zip(ncycs_, self.pw90 * (self.ncyc_max - ncycs_) + self.time_eq))
        deltas[0] = self.pw90 * (self.ncyc_max - 1) + self.time_eq
        delays = [self.t_neg, self.t_pos2]
        delays.extend(tau_cps.values())
        delays.extend(deltas.values())

        return tau_cps, deltas, delays

    @staticmethod
    @ft.lru_cache()
    def _get_phases(ncyc):
        cp_phases = [
            [0, 0, 1, 3, 0, 0, 3, 1, 0, 0, 3, 1, 0, 0, 1, 3],
            [1, 3, 2, 2, 3, 1, 2, 2, 3, 1, 2, 2, 1, 3, 2, 2],
        ]
        phases = np.take(cp_phases, np.flip(np.arange(2 * ncyc)), mode="wrap", axis=1)
        return phases

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        return ncycs_ / self.time_t2

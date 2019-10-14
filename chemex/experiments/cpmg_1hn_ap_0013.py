"""
1HN - Pure Anti-phase Proton CPMG with 0013 phase cycle
=======================================================

Analyzes amide proton chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use
the (6n)x(6n), two-spin matrix, where n is the number of states:

[ Hx(a), Hy(a), Hz(a), 2HxNz(a), 2HyNz(a), 2HzNz(a),
  Hx(b), Hy(b), Hz(b), 2HxNz(b), 2HyNz(b), 2HzNz(b),
   ... ]

This version is modified such that CPMG pulses are applied with [0013] phase cycle
as shown in Daiwen's paper.

References
----------
Yuwen T and Kay LE. Journal of magnetic resonance (2019) in press

Note
----
A sample configuration  file for this module is available using the command:

    chemex config cpmg_1hn_ap_0013

"""
import functools as ft

import numpy as np
import numpy.linalg as nl

import chemex.containers.cpmg as ccc
import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.propagator as cnp
import chemex.nmr.rates as cnr


_SCHEMA = {
    "type": "object",
    "properties": {
        "experiment": {
            "type": "object",
            "properties": {
                "time_t2": {"type": "number"},
                "carrier": {"type": "number"},
                "pw90": {"type": "number"},
                "ncyc_max": {"type": "integer"},
                "eburp_flg": {"type": "boolean", "default": False},
                "reburp_flg": {"type": "boolean", "default": False},
                "pw_eburp": {"type": "number", "default": 1.4e-3},
                "pw_reburp": {"type": "number", "default": 1.52e-3},
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
    config["spin_system"] = {
        "basis": "ixyzsz",
        "atoms": {"i": "h", "s": "n"},
        "constraints": ["hn"],
        "rates": "hn",
    }
    ch.validate(config, _SCHEMA)
    ch.validate(config, ccc.CPMG_SCHEMA)
    experiment = ceh.read(
        config=config,
        pulse_seq_cls=PulseSeq,
        propagator_cls=cnp.PropagatorIS,
        container_cls=ccc.CpmgProfile,
        rates_cls=cnr.RatesIS,
        fit_setting=_FIT_SETTING,
    )
    return experiment


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t2 = settings["time_t2"]
        self.prop.carrier_i = settings["carrier"]
        self.pw90 = settings["pw90"]
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.ncyc_max = settings["ncyc_max"]
        self.eburp_flg = settings["eburp_flg"]
        self.reburp_flg = settings["reburp_flg"]
        self.pw_eburp = settings["pw_eburp"]
        self.pw_reburp = settings["pw_reburp"]
        self.observed_state = settings["observed_state"]
        self.prop.detection = f"2izsz_{self.observed_state}"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps, deltas, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_neg = delays[self.t_neg]
        d_delta = {ncyc: delays[delay] for ncyc, delay in deltas.items()}
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p90 = self.prop.p90_i
        p180 = self.prop.p180_i

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization(terms=f"2izsz_{self.observed_state}")

        # Calculating the intensities as a function of ncyc
        if self.eburp_flg:
            p180pmy = p180[[1, 3]]
            pp90pmy = self.prop.perfect90_i[[1, 3]]
            d_eburp = delays[self.pw_eburp]
            e180e_pmy = pp90pmy @ d_eburp @ p180pmy @ d_eburp @ pp90pmy
            centre = np.mean([p180pmy @ e180e_pmy, e180e_pmy @ p180pmy], axis=0)
        elif self.reburp_flg:
            pp180pmy = self.prop.perfect180_i[[1, 3]]
            d_reburp = delays[0.5 * self.pw_reburp]
            centre = d_reburp @ pp180pmy @ d_reburp
        else:
            centre = p180[[1, 3]]
        centre = np.mean(centre, axis=0)

        intst = {0: self.prop.detect(d_delta[0] @ p90[0] @ centre @ p90[0] @ start)}

        for ncyc in set(ncycs) - {0}:
            phases1, phases2 = self._get_phases(ncyc)
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg1 = ft.reduce(np.matmul, echo[phases1.T])
            cpmg2 = ft.reduce(np.matmul, echo[phases2.T])
            intst[ncyc] = self.prop.detect(
                d_delta[ncyc]
                @ p90[0]
                @ d_neg
                @ cpmg2
                @ centre
                @ cpmg1
                @ d_neg
                @ p90[0]
                @ start
            )

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_) - 0.75 * self.pw90))
        deltas = dict(zip(ncycs_, self.pw90 * (self.ncyc_max - ncycs_)))
        deltas[0] = self.pw90 * self.ncyc_max
        delays = [self.t_neg, self.pw_eburp, 0.5 * self.pw_reburp]
        delays.extend(tau_cps.values())
        delays.extend(deltas.values())
        return tau_cps, deltas, delays

    @staticmethod
    @ft.lru_cache()
    def _get_phases(ncyc):
        cp_phases1 = [
            [1, 1, 2, 0, 1, 1, 0, 2, 1, 1, 0, 2, 1, 1, 2, 0],
            [2, 0, 3, 3, 0, 2, 3, 3, 0, 2, 3, 3, 2, 0, 3, 3],
        ]
        cp_phases2 = [
            [3, 3, 2, 0, 3, 3, 0, 2, 3, 3, 0, 2, 3, 3, 2, 0],
            [2, 0, 1, 1, 0, 2, 1, 1, 0, 2, 1, 1, 2, 0, 1, 1],
        ]
        phases1 = np.take(cp_phases1, np.flip(np.arange(ncyc)), mode="wrap", axis=1)
        phases2 = np.take(cp_phases2, np.arange(ncyc), mode="wrap", axis=1)
        return phases1, phases2

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.array(ncycs, dtype=np.float)
        ncycs_[ncycs_ == -1.0] = 0.5
        return ncycs_[ncycs_ > 0.0] / self.time_t2

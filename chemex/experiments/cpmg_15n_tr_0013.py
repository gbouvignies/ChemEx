"""
15N - N-H TROSY CPMG with 0013 phase cycle
==========================================

Analyzes 15N constant-time TROSY CPMG relaxation dispersion experiments for
measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the (6n)x(6n), two spin matrix, where n is the number of
states:

[ Nx(a), Ny(a), Nz(a), 2HzNx(a), 2HzNy(a), 2HzNz(a),
  Nx(b), Ny(b), Nz(b), 2HzNx(b), 2HzNy(b), 2HzNz(b),
  ... ]

This version is modified such that CPMG pulses are applied with [0013] phase cycle
as shown in Daiwen's paper.

References
----------

Jiang, Yu, Zhang, Liu and Yang. J Magn Reson (2015) 257:1-7
Vallurupalli, Hansen, Stollar, Meirovitch and Kay. PNAS (2007) 104:18473-18477


Note
----

A sample configuration file for this module is available using the command:

    chemex config cpmg_15n_tr_0013

"""
import functools as ft

import numpy as np

import chemex.containers.cpmg as ccc
import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.propagator as cnp


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
                "taub": {"type": "number", "default": 2.68e-3},
                "antitrosy": {"type": "boolean", "default": False},
                "s3e": {"type": "boolean", "default": True},
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
_FIT_SETTING = {"dw_ab": "fit", "r2_a": "fit"}


def read(config):
    config["spin_system"] = {
        "basis": "ixyzsz",
        "atoms": {"i": "n", "s": "h"},
        "constraints": ["nh"],
        "rates": "nh",
    }
    ch.validate(config, _SCHEMA)
    ch.validate(config, ccc.CPMG_SCHEMA)
    if config["experiment"]["antitrosy"]:
        _FIT_SETTING["etaxy_a"] = "fit"
    experiment = ceh.read(
        config=config,
        pulse_seq_cls=PulseSeq,
        propagator_cls=cnp.PropagatorIS,
        container_cls=ccc.CpmgProfile,
        fit_setting=_FIT_SETTING,
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
        self.taub = settings["taub"] - 2.0 * self.pw90 - 2.0 * self.pw90 / np.pi
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.antitrosy = settings["antitrosy"]
        self.s3e = settings["s3e"]
        self.prop.detection = self._get_detection(settings["observed_state"])
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps, deltas, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]
        d_taub = delays[self.taub]
        d_delta = {ncyc: delays[delay] for ncyc, delay in deltas.items()}
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p90 = self.prop.p90_i
        p180 = self.prop.p180_i
        p180_sx = self.prop.perfect180_s[0]

        # Getting the starting magnetization
        start = self._get_start()

        # Calculating the p-element
        if self.antitrosy:
            p180_is = p90[1] @ p90[0] @ p180_sx @ p90[0] @ p90[1]
        else:
            p180_is = p90[0] @ p90[1] @ p180_sx @ p90[1] @ p90[0]
        palmer = np.mean(
            p90[[0, 2]] @ d_taub @ p180_is @ d_taub @ p180_sx @ p90[[1, 3]], axis=0
        )

        # Calculating the cpmg trains
        intst = {
            0: self.prop.detect(
                d_eq
                @ d_delta[0]
                @ p90[3]
                @ p180[[1, 0]]
                @ palmer
                @ p180[[0, 1]]
                @ p90[0]
                @ d_delta[0]
                @ start
            )
        }
        for ncyc in set(ncycs) - {0}:
            phases1, phases2 = self._get_phases(ncyc)
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg1 = d_neg @ ft.reduce(np.matmul, echo[phases1.T]) @ d_neg
            cpmg2 = d_neg @ ft.reduce(np.matmul, echo[phases2.T]) @ d_neg
            intst[ncyc] = self.prop.detect(
                d_eq
                @ d_delta[ncyc]
                @ p90[3]
                @ cpmg2
                @ palmer
                @ cpmg1
                @ p90[0]
                @ d_delta[ncyc]
                @ start
            )

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_) - 0.75 * self.pw90))
        deltas = dict(zip(ncycs_, 0.5 * self.pw90 * (self.ncyc_max - ncycs_)))
        deltas[0] = 0.5 * self.pw90 * (self.ncyc_max - 1)
        delays = [self.t_neg, self.taub, self.time_eq]
        delays.extend(tau_cps.values())
        delays.extend(deltas.values())

        return tau_cps, deltas, delays

    @staticmethod
    @ft.lru_cache()
    def _get_phases(ncyc):
        cp_phases1 = [
            [1, 1, 0, 2, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 0, 2],
            [0, 2, 3, 3, 2, 0, 3, 3, 2, 0, 3, 3, 0, 2, 3, 3],
        ]
        cp_phases2 = [
            [0, 0, 1, 3, 0, 0, 3, 1, 0, 0, 3, 1, 0, 0, 1, 3],
            [1, 3, 2, 2, 3, 1, 2, 2, 3, 1, 2, 2, 1, 3, 2, 2],
        ]
        phases1 = np.take(cp_phases1, np.flip(np.arange(ncyc)), mode="wrap", axis=1)
        phases2 = np.take(cp_phases2, np.arange(ncyc), mode="wrap", axis=1)
        return phases1, phases2

    def _get_start(self):
        start = self.prop.get_start_magnetization("2izsz")
        if self.s3e:
            if self.antitrosy:
                start += self.prop.get_start_magnetization("iz")
            else:
                start -= self.prop.get_start_magnetization("iz")
            start *= 0.5
        return start

    def _get_detection(self, state):
        if self.antitrosy:
            detection = f"2izsz_{state} + iz_{state}"
        else:
            detection = f"2izsz_{state} - iz_{state}"
        return detection

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        return ncycs_ / self.time_t2

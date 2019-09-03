"""
15N - N-H TROSY CPMG
====================

Analyzes 15N constant-time TROSY CPMG relaxation dispersion experiments for
measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the (6n)x(6n), two spin matrix, where n is the number of
states:

[ Nx(a), Ny(a), Nz(a), 2HzNx(a), 2HzNy(a), 2HzNz(a),
  Nx(b), Ny(b), Nz(b), 2HzNx(b), 2HzNy(b), 2HzNz(b),
  ... ]

References
----------
Vallurupalli et al. Proc Natl Acad Sci USA (2007) 104


Note
----

A sample configuration  file for this module is available using the command:

    chemex config cpmg_15n_tr

"""
import functools as ft

import numpy as np
import numpy.linalg as nl

import chemex.containers.cpmg as ccc
import chemex.experiments.helper as ceh
import chemex.nmr.propagator as cnp


TYPE = __name__.split(".")[-1]
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
    config["spin_system"] = {
        "basis": "ixyzsz",
        "atoms": {"i": "n", "s": "h"},
        "constraints": ["nh"],
    }
    ceh.validate(config, _SCHEMA)
    ceh.validate(config, ccc.CPMG_SCHEMA)
    experiment = ceh.read(
        config=config,
        pulse_seq_cls=PulseSeq,
        propagator_cls=cnp.PropagatorIS,
        container_cls=ccc.CpmgProfile,
    )
    return experiment


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t2 = settings["time_t2"]
        self.time_eq = settings["time_equil"]
        self.prop.carrier_i = settings["carrier"]
        self.pw90 = settings["pw90"]
        self.taub = settings["taub"] - 2.0 * self.pw90 - 2.0 * self.pw90 / np.pi
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.antitrosy = settings["antitrosy"]
        self.prop.detection = self._get_detection(settings["observed_state"])
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]
        d_taub = delays[self.taub]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p90, p180 = self.prop.pulses_90_180_i()
        p180_sx = self.prop.perfect180_s[0]

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization("2izsz")

        # Calculating the p-element
        if self.antitrosy:
            palmer0 = (
                p180_sx @ d_taub @ p90[2] @ p90[1] @ p180_sx @ p90[1] @ p90[2] @ d_taub
            )
        else:
            palmer0 = (
                p180_sx @ d_taub @ p90[1] @ p90[0] @ p180_sx @ p90[0] @ p90[1] @ d_taub
            )
        palmer = np.mean(p90[[0, 2]] @ palmer0 @ p90[[1, 3]], axis=0)

        # Calculating the cpmg trains
        part1 = p90[0] @ start
        part2 = d_eq @ p90[1]
        intst = {0: self.prop.detect(part2 @ palmer0 @ part1)}
        for ncyc in set(ncycs) - {0}:
            echo = d_cp[ncyc] @ p180[[1, 0]] @ d_cp[ncyc]
            cpmg1, cpmg2 = nl.matrix_power(echo, ncyc)
            intst[ncyc] = self.prop.detect(
                part2 @ d_neg @ cpmg2 @ palmer @ cpmg1 @ d_neg @ part1
            )

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_) - self.pw90))
        delays = [self.t_neg, self.taub, self.time_eq]
        delays.extend(tau_cps.values())

        return tau_cps, delays

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

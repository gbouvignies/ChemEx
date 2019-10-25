"""
1H-15N DQ/ZQ CPMG
=================

Analyzes 15N and 1H chemical exchange by applying CPMG pulses on 15N and 1H
simultaneously. The spin system is maintained as DQ or ZQ during Trelax, and 
is calculated using the (15n)x(15n), two spin matrix, where n is the number
of states:

[        Ix(a),   Iy(a),   Iz(a),   Sx(a), IxSx(a), IySx(a), IzSx(a), 
  Sy(a), IxSy(a), IySy(a), IzSy(a), Sz(a), IxSz(a), IySz(a), IzSz(a),
         Ix(b),   Iy(b),   Iz(b),   Sx(b), IxSx(b), IySx(b), IzSx(b), 
  Sy(b), IxSy(b), IySy(b), IzSy(b), Sz(b), IxSz(b), IySz(b), IzSz(b),
   ... ]

The phase cycle of CPMG pulses is chosen based on νCPMG as described in the 
reference, which is a mixture of constant-phase and XY-family phase cycles.

References
----------

Orekhov, Korzhnev and Kay. J Am Chem Soc (2004) 126:1886-1891


Note
----

A sample configuration file for this module is available using the command:

    chemex config cpmg_hn_dq_zq

"""
import functools as ft

import numpy as np

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
                "carrier_h": {"type": "number"},
                "carrier_n": {"type": "number"},
                "pw90_h": {"type": "number"},
                "pw90_n": {"type": "number"},
                "dq_flg": {"type": "boolean"},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
            },
            "required": [
                "time_t2",
                "carrier_h",
                "carrier_n",
                "pw90_h",
                "pw90_n",
                "dq_flg",
            ],
        }
    },
}
_FIT_SETTING = {"dw_ab": "fit", "r2mq_a": "fit", "mu_a": "fit"}


def read(config):
    config["spin_system"] = {
        "basis": "ixyzsxyz",
        "atoms": {"i": "n", "s": "h"},
        "constraints": ["nh"],
        "rates": "nh",
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
        self.prop.carrier_i = settings["carrier_n"]
        self.prop.carrier_s = settings["carrier_h"]
        self.pw90_i = settings["pw90_n"]
        self.pw90_s = settings["pw90_h"]
        self.prop.b1_i = 1 / (4.0 * self.pw90_i)
        self.prop.b1_s = 1 / (4.0 * self.pw90_s)
        self.dq_flg = settings["dq_flg"]
        self.prop.detection = self._get_detection(settings["observed_state"])
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps = self._get_tau_cps(ncycs)
        tau_cp_list = list(tau_cps.values())
        delays = dict(zip(tau_cp_list, self.prop.delays(tau_cp_list)))
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p9024090 = self.prop.p9024090_is[[0, 1], [0, 1]]

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization(["2ixsx"])

        # Calculating the cpmg trains
        intst = {0: self.prop.detect(start)}
        for ncyc in set(ncycs) - {0}:
            phases = self._get_phases(ncyc)
            echo = d_cp[ncyc] @ p9024090 @ d_cp[ncyc]
            cpmg = ft.reduce(np.matmul, echo[phases])
            intst[ncyc] = self.prop.detect(cpmg @ start)
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_tau_cps(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        tau_cps = dict(
            zip(ncycs_, self.time_t2 / (4.0 * ncycs_) - 7.0 / 3.0 * self.pw90_i)
        )
        return tau_cps

    @ft.lru_cache()
    def _get_phases(self, ncyc):
        nu_cpmg = self.ncycs_to_nu_cpmgs(ncyc)
        if nu_cpmg < 51.0:
            phases_ = [0, 1, 0, 1]
        elif nu_cpmg < 255.0:
            phases_ = [0]
        else:
            phases_ = [0, 1, 0, 1, 1, 0, 1, 0]
        phases = np.take(phases_, np.flip(np.arange(2 * ncyc)), mode="wrap")
        return phases

    def _get_detection(self, state):
        if self.dq_flg:
            detection = f"2ixsx_{state} - 2iysy_{state}"
        else:
            detection = f"2ixsx_{state} + 2iysy_{state}"
        return detection

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        return ncycs_ / self.time_t2

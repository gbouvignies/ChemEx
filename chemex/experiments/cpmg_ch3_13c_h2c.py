"""
13C (Methyl) H-to-C CPMG
========================

Measures methyl carbon chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization is
initially anti-phase and is read out as in-phase. Because of the P-element 
only even ncyc should be recorded. Resulting magnetization intensity after 
the CPMG block is calculated using the (6n)×(6n), two-spin matrix, where n 
is the number of states::

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

References
----------
Lundström, Vallurupalli, Religa, Dahlquist and Kay. J Biomol NMR (2007) 38, 79-88


Note
----
A sample configuration file for this module is available using the command::

    $ chemex config cpmg_ch3_13c_h2c

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
                "taub": {"type": "number", "default": 2.0e-3},
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
    config["basis"] = cnl.Basis(type="ixyzsz", spin_system="ch")
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
        self.taub = settings["taub"] - 2.0 * self.pw90
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
        d_taub = delays[self.taub]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p90 = self.prop.p90_i
        p180 = self.prop.p180_i
        p180_sx = self.prop.perfect180_s[0]

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization("2izsz")

        # Calculating the p-element
        palmer = d_taub @ p90[0] @ p180_sx @ p90[0] @ d_taub

        # Calculating the inensities as a function of ncyc
        part1 = d_neg @ p90[0] @ start
        part2 = d_eq @ p90[1] @ d_neg
        intst = {0: self.prop.detect(part2 @ palmer @ part1)}
        for ncyc in set(ncycs) - {0}:
            echo = d_cp[ncyc] @ p180[[1, 0]] @ d_cp[ncyc]
            cpmg1, cpmg2 = nl.matrix_power(echo, ncyc)
            end = part2 @ cpmg2 @ palmer @ cpmg1 @ part1
            intst[ncyc] = self.prop.detect(end)

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

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        return ncycs_ / self.time_t2

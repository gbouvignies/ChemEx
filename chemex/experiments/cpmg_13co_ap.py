"""
13C (Carbonyl) Pure Anti-phase CPMG
===================================

Analyzes carbonyl chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a
(6n)×(6n), two-spin matrix, where n is the number of states::

    { COx(a), COy(a), COz(a), 2COxNz(a), 2COyNz(a), 2COzNz(a),
      COx(b), COy(b), COz(b), 2COxNz(b), 2COyNz(b), 2COzNz(b), ... }

Because of the length of the shaped pulses used during the CPMG blocks,
off-resonance effects are taken into account only for the 90-degree pulses
that create COxNz before the CPMG and COzNz after the CPMG.

The calculation can be run with or without C–C J-coupling refocusing element
via the "refocusing" flag, with such option ncyc_cp should be set as even.

References
----------

Lundström, Hansen and Kay. J Biomol NMR (2008) 42:35-47


Note
----

A sample configuration file for this module is available using the command::

    $ chemex config cpmg_13co_ap

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
                "refocusing": {"type": "boolean", "default": False},
                "taucc": {"type": "number", "default": 9.09e-3},
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
    config["basis"] = cnl.Basis(type="ixyzsz", spin_system="cn")
    config["even_ncycs"] = config["experiment"]["refocusing"]
    config["fit"] = _fit_this()
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this():
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
        self.refocusing = settings["refocusing"]
        self.taucc = settings["taucc"] - self.pw90
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.observed_state = settings["observed_state"]
        self.prop.detection = f"[2izsz_{self.observed_state}]"

    @ft.lru_cache(maxsize=10000)
    def calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]
        d_taucc = delays[self.taucc]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p90 = self.prop.p90_i
        p180 = self.prop.p180_i
        p180pmy = 0.5 * (p180[1] + p180[3])  # +/- phase cycling
        perfect180x = self.prop.perfect180_i[0]

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization(["2izsz"])

        # Calculate the flip block
        if self.refocusing:
            p_flip0 = p90[3] @ d_taucc @ p180pmy @ d_taucc @ p90[1]
            p_flip = d_neg @ p_flip0 @ d_neg
        else:
            p_flip = p_flip0 = p180pmy

        # Calculating the intensities as a function of ncyc
        intst = {0: self.prop.detect(d_eq @ p90[1] @ p_flip0 @ p90[1] @ start)}

        part1 = d_neg @ p90[1] @ start
        part2 = d_eq @ p90[1] @ d_neg
        for ncyc in set(ncycs) - {0}:
            echo = d_cp[ncyc] @ perfect180x @ d_cp[ncyc]
            cpmg = nl.matrix_power(echo, int(ncyc))
            end = part2 @ cpmg @ p_flip @ cpmg @ part1
            intst[ncyc] = self.prop.detect(end)

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_)))
        delays = [self.t_neg, self.time_eq, self.taucc, *tau_cps.values()]
        return tau_cps, delays

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        return ncycs_ / self.time_t2

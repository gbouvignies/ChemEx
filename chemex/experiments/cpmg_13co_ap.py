"""
13CO - Pure Anti-phase Carbonyl 13C CPMG
========================================

Analyzes carbonyl chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a
(6n)x(6n), two spin matrix, where n is the number of states:

[ COx(a), COy(a), COz(a), 2COxNz(a), 2COyNz(a), 2COzNz(a),
  COx(b), COy(b), COz(b), 2COxNz(b), 2COyNz(b), 2COzNz(b),
  ...]

Because of the length of the shaped pulses used during the CPMG blocks, off-
resonance effects are taken into account only for the 90-degree pulses that
create COxNz before the CPMG and COzNz after the CPMG.

The calculation can be run with or without sidechain CO inversion via the
'sidechain' flag for uniformly 13C-labeled proteins.

References
----------
Journal of Biomolecular NMR (2008) 42, 35-47


Note
----

A sample configuration  file for this module is available using the command:

    chemex config cpmg_13co_ap

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
                "sidechain": {"type": "boolean", "default": False},
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
    config["spin_system"] = {
        "basis": "ixyzsz",
        "atoms": {"i": "c", "s": "n"},
        "constraints": ["hn"],
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
        self.t_neg = -2.0 * self.pw90 / np.pi
        self.sidechain = settings["sidechain"]
        self.taucc = settings["taucc"]
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.observed_state = settings["observed_state"]
        self.prop.detection = f"2izsz_{self.observed_state}"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_neg = delays[self.t_neg]
        d_eq = delays[self.time_eq]
        d_taucc = delays[self.taucc]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p90, p180 = self.prop.pulses_90_180_i()
        perfect180x = self.prop.perfect180["ix"]
        p180pmy = 0.5 * (p180[1] + p180[3])  # +/- phase cycling

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization(terms=f"2izsz_{self.observed_state}")

        # Calculate the flip block
        if self.sidechain:
            p_flip = p180pmy
        else:
            p_flip = p90[3] @ d_taucc @ p180pmy @ d_taucc @ p90[1]

        # Calculating the intensities as a function of ncyc
        part1 = p90[1] @ start
        part2 = d_eq @ p90[1]
        intst = {0: self.prop.detect(part2 @ p_flip @ part1)}

        for ncyc in set(ncycs) - {0, -1}:
            echo = d_cp[ncyc] @ perfect180x @ d_cp[ncyc]
            cpmg = nl.matrix_power(echo, int(ncyc))
            end = part2 @ d_neg @ cpmg @ p180pmy @ cpmg @ d_neg @ part1
            intst[ncyc] = self.prop.detect(end)

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_)))
        tau_cps[-1] = 0.5 * self.time_t2
        delays = [self.t_neg, self.time_eq, self.taucc]
        delays.extend(tau_cps.values())
        return tau_cps, delays

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.array(ncycs, dtype=np.float)
        ncycs_[ncycs_ == -1.0] = 0.5
        return ncycs_[ncycs_ > 0.0] / self.time_t2
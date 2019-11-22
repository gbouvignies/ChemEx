"""
1H-13C(methyl) - Multiple Quantum CPMG
======================================

Analyzes HyCx methyl group multiple quantum CPMG measured on site-specific
13CH3-labeled methyl groups in a highly deuterated background. Resulting
magnetization intensity after the CPMG block is calculated using the (4n)x(4n),
 where n is the number of states:

[HxCx(a), HyCx(a), HxCy(a), HyCy(a),
 HxCx(b), HyCx(b), HxCy(b), HyCy(b)]

This is a simplified basis set, which assumes you are on-resonance for 13C
(ie, off-resonance effects are not taken into account)

This calculation is designed specifically to analyze data from the experiment
found in the reference and can be run with either small_protein_flag='y' or 'n'.

Reference
---------
Korzhnev et al. Journal of the American Chemical Society (2004), 126, 3964-73

Note
----
A sample configuration  file for this module is available using the command:

    chemex config cpmg_ch3_mq

"""
import functools as ft

import numpy as np
import numpy.linalg as nl

import chemex.experiments.helper as ceh
import chemex.helper as ch


_SCHEMA = {
    "type": "object",
    "properties": {
        "experiment": {
            "type": "object",
            "properties": {
                "time_t2": {"type": "number"},
                "small_protein": {"type": "boolean", "default": False},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
            },
            "required": ["time_t2"],
        }
    },
}
_FIT_SETTING = {"dw_ab, nuc->c": "fit", "r2mq_a": "fit"}


def read(config):
    config["spin_system"] = {"basis": "ixysxy", "atoms": {"i": "c", "s": "h"}}
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
        self.small_protein = settings["small_protein"]
        self.t_zeta = 1.0 / (8.0 * 125.3)
        self.prop.detection = f"2iysx_{settings['observed_state']}"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_zeta = delays[self.t_zeta]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the propagators corresponding to all the pulses
        p180_sx = self.prop.perfect180_s[0]
        p180_ix = self.prop.perfect180_i[0]
        p180_iy = self.prop.perfect180_i[1]

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization(terms=["2iysx"])

        # Calculating the instensities as a function of ncyc
        if self.small_protein:
            part1 = d_zeta @ p180_sx @ p180_ix @ d_zeta @ start
        else:
            part1 = start
        intst = {0: self.prop.detect(part1)}
        for ncyc in set(ncycs) - {0}:
            echo = d_cp[ncyc] @ p180_iy @ d_cp[ncyc]
            cpmg = nl.matrix_power(echo, int(ncyc))
            intst[ncyc] = self.prop.detect(cpmg @ p180_sx @ cpmg @ part1)

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_)))
        delays = [self.t_zeta]
        delays.extend(tau_cps.values())
        return tau_cps, delays

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.array(ncycs, dtype=np.float)
        ncycs_[ncycs_ == -1.0] = 0.5
        return ncycs_[ncycs_ > 0.0] / self.time_t2

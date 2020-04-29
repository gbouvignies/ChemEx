"""
1H (Methyl-13CH3) Double-Quantum Proton CPMG
============================================

Measures methyl proton chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization
is initially anti-phase and is read out as anti-phase. The calculation
uses a simplified scheme that includes only (6n)×(6n) basis set, 
two-spin matrix, where n is the number of states::

    { 2HDQxCz(a), 2HDQyCz(a), 2HzCz(a), HDQx(a), HDQy(a), Hz(a),
      2HDQxCz(b), 2HDQyCz(b), 2HzCz(b), HDQx(b), HDQy(b), Hz(b), ... }

References
----------

Gopalan, Yuwen, Kay and Vallurupalli. J Biomol NMR (2018) 72:79-91


Note
----

A sample configuration  file for this module is available using the command::

    $ chemex config cpmg_ch3_1h_dq

"""
import functools as ft

import numpy as np

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
                "comp180_flg": {"type": "boolean", "default": True},
                "ipap_flg": {"type": "boolean", "default": False},
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
    config["basis"] = cnl.Basis(type="ixyzsz", extension="dq", spin_system="hc")
    config["fit"] = _fit_this(config)
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this(config):
    this = {
        "rates": ["r2_i_{observed_state}"],
        "model_free": ["tauc_{observed_state}"],
    }
    if config["experiment"]["ipap_flg"]:
        this["rates"].append("r2a_i_{observed_state}")
        this["model_free"].append("s2_{observed_state}")
    return this


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t2 = settings["time_t2"]
        self.tauc = 0.67e-3  # ~ 1/(12*J[HC])
        self.pw90 = settings["pw90"]
        self.ipap_flg = settings["ipap_flg"]
        self.comp180_flg = settings["comp180_flg"]
        self.prop.carrier_i = settings["carrier"]
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.prop.detection = f"2ixsz_{settings['observed_state']}"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization(terms="2ixsz", atom="h")

        # Calculation of the propagators corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}
        d_tauc = delays[self.tauc]

        # Calculation of the propagators corresponding to all the pulses
        p180 = self.prop.p180_i
        p180pmy = 0.5 * (p180[1] + p180[3])  # +/- phase cycling
        if self.comp180_flg:
            p180_cp1 = self.prop.p9018090_i_1
            p180_cp2 = self.prop.p9018090_i_2
        else:
            p180_cp1 = p180_cp2 = p180

        # Calculating the intensities as a function of ncyc
        if self.ipap_flg:
            p180_sx = self.prop.perfect180_s[0]
            part1 = d_tauc @ start
            part2 = d_tauc
            centre0 = 0.5 * (p180pmy + p180_sx @ p180pmy @ p180_sx)
        else:
            part1 = start
            part2 = self.prop.identity
            centre0 = p180pmy

        intst = {0: self.prop.detect(part2 @ centre0 @ part1)}

        for ncyc in set(ncycs) - {0}:
            phases1, phases2 = self._get_phases(ncyc)
            echo1 = d_cp[ncyc] @ p180_cp1 @ d_cp[ncyc]
            echo2 = d_cp[ncyc] @ p180_cp2 @ d_cp[ncyc]
            cpmg1 = ft.reduce(np.matmul, echo1[phases1])
            cpmg2 = ft.reduce(np.matmul, echo2[phases2])
            centre = cpmg2 @ p180pmy @ cpmg1
            if self.ipap_flg:
                centre = 0.5 * (centre + p180_sx @ centre @ p180_sx)
            intst[ncyc] = self.prop.detect(part2 @ centre @ part1)

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        if self.comp180_flg:
            tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_) - 2 * self.pw90))
        else:
            tau_cps = dict(zip(ncycs_, self.time_t2 / (4.0 * ncycs_) - self.pw90))
        delays = [self.tauc]
        delays.extend(tau_cps.values())
        return tau_cps, delays

    @staticmethod
    def _get_phases(ncyc):
        cp_phases1 = [0, 1, 0, 1, 1, 0, 1, 0, 2, 3, 2, 3, 3, 2, 3, 2]
        cp_phases2 = [0, 3, 0, 3, 3, 0, 3, 0, 2, 1, 2, 1, 1, 2, 1, 2]
        phases1 = np.take(cp_phases1, np.flip(np.arange(ncyc)), mode="wrap")
        phases2 = np.take(cp_phases2, np.arange(ncyc), mode="wrap")
        return phases1, phases2

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        return ncycs_ / self.time_t2
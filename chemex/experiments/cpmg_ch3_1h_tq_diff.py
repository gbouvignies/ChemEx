"""
1H (Methyl-13CH3) Triple-Quantum Proton CPMG with Gradients
===========================================================

Measures methyl proton chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. The experiment
is specifically designed for measuring translational diffusion constants
of invisible states using a pulsed-field gradient approach that exploits
methyl 1H triple-quantum coherences. Magnetization is initially anti-phase
and is read out as anti-phase. The calculation uses a simplified scheme
that includes only (6n)×(6n) basis set, two-spin matrix, where n is the
number of states::

    { 2HTQxCz(a), 2HTQyCz(a), 2HzCz(a), HTQx(a), HTQy(a), Hz(a),
      2HTQxCz(b), 2HTQyCz(b), 2HzCz(b), HTQx(b), HTQy(b), Hz(b), ... }

References
----------
Yuwen, Sekhar, Baldwin, Vallurupalli and Kay. Angew Chem Int Ed (2018) 57, 16777-16780


Note
----

A sample configuration  file for this module is available using the command::

    $ chemex config cpmg_ch3_1h_tq_diff

"""
import functools as ft

import numpy as np

import chemex.experiments.helper as ceh
import chemex.helper as ch
import chemex.nmr.constants as cnc
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
                "delta": {"type": "number"},
                "gradient": {"type": "number"},
                "tau": {"type": "number", "default": 500.0e-6},
                "comp180_flg": {"type": "boolean", "default": True},
                "ipap_flg": {"type": "boolean", "default": False},
                "observed_state": {
                    "type": "string",
                    "pattern": "[a-z]",
                    "default": "a",
                },
            },
            "required": ["time_t2", "carrier", "pw90", "delta", "gradient"],
        }
    },
}


def read(config):
    ch.validate(config, _SCHEMA)
    config["basis"] = cnl.Basis(type="ixyzsz_diff", extension="tq", spin_system="hc")
    config["fit"] = _fit_this(config)
    return ceh.load_experiment(config=config, pulse_seq_cls=PulseSeq)


def _fit_this(config):
    this = {
        "rates": ["r2_i_{observed_state}", "d_{states}"],
        "model_free": ["tauc_{observed_state}", "d_{states}"],
    }
    if not config["experiment"]["ipap_flg"]:
        this["rates"].append("r2a_i_{observed_state}")
        this["model_free"].append("s2_{observed_state}")
    return this


class PulseSeq:
    def __init__(self, config, propagator):
        self.prop = propagator
        settings = config["experiment"]
        self.time_t2 = settings["time_t2"]
        self.delta = settings["delta"]
        self.tau = settings["tau"]
        self.tauc = 0.67e-3  # ~ 1/(12*J[HC])
        self.pw90 = settings["pw90"]
        self.k2_factor = (
            3.0 * cnc.GAMMA["h"] * settings["gradient"] * settings["delta"]
        ) ** 2
        self.ipap_flg = settings["ipap_flg"]
        self.comp180_flg = settings["comp180_flg"]
        self.prop.carrier_i = settings["carrier"]
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.prop.detection = f"[2ixsz_{settings['observed_state']}]"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization(terms="2ixsz", atom="h")

        # Calculation of the propagators with gradient 0
        self.prop.gradient_dephasing = 0.0
        d_tau_0, d_tauc_0 = self.prop.delays([self.tau, self.tauc])
        p180_sx = self.prop.perfect180_s[0]

        # Calculation of the propagators with gradient 1
        self.prop.gradient_dephasing = 1.0 / 3.0 * self.k2_factor
        (d_delta_1,) = self.prop.delays(self.delta)

        # Calculation of the propagators with gradient 2
        self.prop.gradient_dephasing = self.k2_factor
        (d_tau_2,) = self.prop.delays(self.tau)
        p90_2 = self.prop.p90_i
        p180_2 = self.prop.p180_i

        # Calculation of the propagators with gradient 3
        self.prop.gradient_dephasing = 7.0 / 3.0 * self.k2_factor
        (d_delta_3,) = self.prop.delays(self.delta)

        # Calculation of the propagators with gradient 4
        self.prop.gradient_dephasing = 4.0 * self.k2_factor
        tau_cps, all_delays = self._get_delays(ncycs)
        delays_4 = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_cp_4 = {ncyc: delays_4[delay] for ncyc, delay in tau_cps.items()}
        d_tau_4 = delays_4[self.tau]
        p180_4 = self.prop.p180_i
        p180pmy_4 = 0.5 * (p180_4[1] + p180_4[3])  # +/- phase cycling

        if self.comp180_flg:
            p180_cp1 = self.prop.p9018090_i_1
            p180_cp2 = self.prop.p9018090_i_2
        else:
            p180_cp1 = p180_cp2 = p180_4

        # The liouvillian for evolution before and after the CPMG period
        p180_grad1, p180_grad2 = p90_2[[1, 3]] @ p180_2[0] @ p90_2[[1, 3]]
        grad1 = d_tau_4 @ d_delta_3 @ p180_grad1 @ d_tau_2 @ d_delta_1
        grad2 = d_tau_0 @ d_delta_1 @ p180_grad2 @ d_tau_2 @ d_delta_3
        if self.ipap_flg:
            part1 = d_tauc_0 @ start
            part2 = d_tauc_0
            centre0 = grad2 @ p180pmy_4 @ grad1
            centre0 = 0.5 * (centre0 + p180_sx @ centre0 @ p180_sx)
        else:
            part1 = start
            part2 = self.prop.identity
            centre0 = grad2 @ p180pmy_4 @ grad1

        # Calculating the intensities as a function of ncyc
        intst = {0: self.prop.detect(part2 @ centre0 @ part1)}

        for ncyc in set(ncycs) - {0}:
            phases1, phases2 = self._get_phases(ncyc)
            echo1 = d_cp_4[ncyc] @ p180_cp1 @ d_cp_4[ncyc]
            echo2 = d_cp_4[ncyc] @ p180_cp2 @ d_cp_4[ncyc]
            cpmg1 = ft.reduce(np.matmul, echo1[phases1])
            cpmg2 = ft.reduce(np.matmul, echo2[phases2])
            centre = grad2 @ cpmg2 @ p180pmy_4 @ cpmg1 @ grad1
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
        delays = [self.tau, *tau_cps.values()]
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

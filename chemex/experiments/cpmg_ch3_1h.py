"""
1H(methyl - 13CH3) - Single Quantum Proton CPMG
===============================================

Measures methyl proton chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization is
initally anti-phase and is read out as anti-phase prior to 1H detection.
Resulting magnetization intensity after the CPMG block is calculated using
the (6n)x(6n), two spin matrix, where n is the number of states:

[ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
  Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b),
  ... ]

Reference
---------
Yuwen, Sekhar, Baldwin, Vallurupalli and Kay. Angewandte Chemie (2018) 51, 16777-80

Note
----
A sample configuration  file for this module is available using the command:

    chemex config cpmg_ch3_1h

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
                "ncyc_max": {"type": "integer"},
                "taua": {"type": "number", "default": 2e-3},
                "comp180_flg": {"type": "boolean", "default": True},
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
        "atoms": {"i": "h", "s": "c"},
        "constraints": ["hn"],
    }
    ch.validate(config, _SCHEMA)
    ch.validate(config, ccc.CPMG_SCHEMA)
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
        self.pw90 = settings["pw90"]
        self.ncyc_max = settings["ncyc_max"]
        self.prop.carrier_i = settings["carrier"]
        self.prop.b1_i = 1 / (4.0 * self.pw90)
        self.taua = settings["taua"]
        self.comp180_flg = settings["comp180_flg"]
        self.observed_state = settings["observed_state"]
        self.prop.detection = f"iy_{self.observed_state}"
        self.calculate = ft.lru_cache(maxsize=5)(self._calculate)

    def _calculate(self, ncycs, params_local):
        self.prop.update(params_local)

        # Calculation of the propagators corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, self.prop.delays(all_delays)))
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}
        d_taua = delays[self.taua]

        # Calculation of the propagators corresponding to all the pulses
        perfect180y = self.prop.perfect180_i[1]
        p180 = self.prop.p180_i
        p180c_py = self.prop.p9018090_i_1[1]
        p180c_my = self.prop.p9018090_i_2[3]
        p180pmy = 0.5 * (p180[1] + p180[3])  # +/- phase cycling
        if self.comp180_flg:
            p180_cp1 = self.prop.p9024090_i_1
            p180_cp2 = self.prop.p9024090_i_2
        else:
            p180_cp1 = p180_cp2 = p180

        # Getting the starting magnetization
        start = self.prop.get_start_magnetization(terms=f"iy")
        start = perfect180y @ d_taua @ d_taua @ start
        start = self.prop.keep_components(start, ["2ixsz_a", "2iysz_a"])

        # Calculating the intensities as a function of ncyc
        intst = {
            0: self.prop.detect(
                d_taua @ (p180pmy @ p180c_py + p180c_my @ p180pmy) @ d_taua @ start
            )
        }

        phases1, phases2 = self._get_phases()
        for ncyc in set(ncycs) - {0}:
            echo1 = d_cp[ncyc] @ p180_cp1 @ d_cp[ncyc]
            echo2 = d_cp[ncyc] @ p180_cp2 @ d_cp[ncyc]
            cpmg1 = ft.reduce(np.matmul, echo1[phases1[-ncyc:]])
            cpmg2 = ft.reduce(np.matmul, echo2[phases2[:ncyc]])
            if ncyc < self.ncyc_max:
                cpmg1 = ft.reduce(np.matmul, p180_cp1[phases1[:-ncyc]]) @ cpmg1
                cpmg2 = cpmg2 @ ft.reduce(np.matmul, p180_cp2[phases2[ncyc:]])
            centre = cpmg2 @ p180pmy @ cpmg1
            intst[ncyc] = self.prop.detect(
                d_taua @ (centre @ p180c_py + p180c_my @ centre) @ d_taua @ start
            )

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @ft.lru_cache()
    def _get_delays(self, ncycs):
        ncycs_ = np.asarray(ncycs)
        ncycs_ = ncycs_[ncycs_ > 0]
        frac = 7.0 / 3.0 if self.comp180_flg else 1.0
        tau_cps = dict(
            zip(
                ncycs_,
                (self.time_t2 - frac * self.pw90 * 4.0 * self.ncyc_max)
                / (4.0 * ncycs_),
            )
        )
        delays = [self.taua]
        delays.extend(tau_cps.values())
        return tau_cps, delays

    def _get_phases(self):
        cp_phases1 = [0, 1]
        cp_phases2 = [0, 3]
        phases1 = np.take(cp_phases1, np.flip(np.arange(self.ncyc_max)), mode="wrap")
        phases2 = np.take(cp_phases2, np.arange(self.ncyc_max), mode="wrap")
        return phases1, phases2

    def ncycs_to_nu_cpmgs(self, ncycs):
        ncycs_ = np.array(ncycs, dtype=np.float)
        return ncycs_[ncycs_ > 0.0] / self.time_t2

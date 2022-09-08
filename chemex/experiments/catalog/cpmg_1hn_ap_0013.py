from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import CpmgSettings
from chemex.configuration.experiment import ExperimentConfig
from chemex.configuration.experiment import ToBeFitted
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators
from chemex.experiments.factories import factories
from chemex.filterers import PlanesFilterer
from chemex.nmr.liouvillian import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters import CpmgPlotter
from chemex.printers.data import CpmgPrinter

# Type definitions
NDArrayFloat = NDArray[np.float_]
NDArrayInt = NDArray[np.int_]
NDArrayBool = NDArray[np.bool_]


EXPERIMENT_NAME = "cpmg_1hn_ap_0013"


class Cpmg1HnAp0013Settings(CpmgSettings):
    name: Literal["cpmg_1hn_ap_0013"]
    time_t2: float
    carrier: float
    pw90: float
    ncyc_max: int
    taua: float = 2.38e-3
    ipap_flg: bool = False
    eburp_flg: bool = False
    reburp_flg: bool = False
    pw_eburp: float = 1.4e-3
    pw_reburp: float = 1.52e-3
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def t_neg(self) -> float:
        return -2.0 * self.pw90 / np.pi

    @property
    def start(self) -> list[str]:
        return [f"2izsz_{self.observed_state}"]

    @property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class Cpmg1HnAp0013Config(
    ExperimentConfig[Cpmg1HnAp0013Settings, RelaxationDataSettings]
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])


def build_spectrometer(
    config: Cpmg1HnAp0013Config, spin_system: SpinSystem
) -> Spectrometer:

    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyzsz", spin_system="hn")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class Cpmg1HnAp0013Sequence:
    settings: Cpmg1HnAp0013Settings

    def _get_delays(
        self, ncycs: NDArrayFloat
    ) -> tuple[dict[float, float], dict[float, float], list[float]]:
        ncycs_no_ref = ncycs[ncycs > 0.0]
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - 0.75 * self.settings.pw90
            for ncyc in ncycs_no_ref
        }
        deltas = {
            ncyc: self.settings.pw90 * (self.settings.ncyc_max - ncyc)
            for ncyc in ncycs_no_ref
        }
        deltas[0.0] = self.settings.pw90 * self.settings.ncyc_max
        delays = [
            self.settings.taua,
            self.settings.t_neg,
            self.settings.pw_eburp,
            0.5 * self.settings.pw_reburp,
            *tau_cps.values(),
            *deltas.values(),
        ]
        return tau_cps, deltas, delays

    def _get_phases(self, ncyc: NDArrayFloat) -> tuple[NDArrayInt, NDArrayInt]:
        cp_phases1 = np.array(
            [
                [1, 1, 2, 0, 1, 1, 0, 2, 1, 1, 0, 2, 1, 1, 2, 0],
                [2, 0, 3, 3, 0, 2, 3, 3, 0, 2, 3, 3, 2, 0, 3, 3],
            ]
        )
        cp_phases2 = np.array(
            [
                [3, 3, 2, 0, 3, 3, 0, 2, 3, 3, 0, 2, 3, 3, 2, 0],
                [2, 0, 1, 1, 0, 2, 1, 1, 0, 2, 1, 1, 2, 0, 1, 1],
            ]
        )
        indexes = np.arange(int(ncyc))
        phases1 = np.take(cp_phases1, np.flip(indexes), mode="wrap", axis=1)
        phases2 = np.take(cp_phases2, indexes, mode="wrap", axis=1)
        return phases1, phases2

    def calculate(self, spectrometer: Spectrometer, data: Data) -> np.ndarray:

        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, deltas, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays)))
        d_neg = delays[self.settings.t_neg]
        d_taua = delays[self.settings.taua]
        d_eburp = delays[self.settings.pw_eburp]
        d_reburp = delays[0.5 * self.settings.pw_reburp]
        d_delta = {ncyc: delays[delay] for ncyc, delay in deltas.items()}
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the spectrometers corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i
        pp90_i = spectrometer.perfect90_i
        pp180_isx = spectrometer.perfect180_i[0] @ spectrometer.perfect180_s[0]

        # Calculation of the spectrometers for INEPT and purge elements
        inept = pp90_i[1] @ d_taua @ pp180_isx @ d_taua @ pp90_i[0]
        zfilter = spectrometer.zfilter

        # Getting the starting magnetization
        start = spectrometer.get_start_magnetization(terms=self.settings.start)

        # Calculating the central refocusing block
        if self.settings.eburp_flg:
            p180pmy = p180[[1, 3]]
            pp90pmy = spectrometer.perfect90_i[[1, 3]]
            e180e_pmy = pp90pmy @ d_eburp @ p180pmy @ d_eburp @ pp90pmy
            middle = [p180pmy @ e180e_pmy, e180e_pmy @ p180pmy]
        elif self.settings.reburp_flg:
            pp180pmy = spectrometer.perfect180_i[[1, 3]]
            middle = d_reburp @ pp180pmy @ d_reburp
        else:
            middle = p180[[1, 3]]
        middle = np.mean(middle, axis=0)

        # Calculating the instensities as a function of ncyc
        centre = {0.0: d_delta[0] @ p90[0] @ middle @ p90[0]}

        for ncyc in set(ncycs) - {0.0}:
            phases1, phases2 = self._get_phases(ncyc)
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg1 = reduce(np.matmul, echo[phases1.T])
            cpmg2 = reduce(np.matmul, echo[phases2.T])
            centre[ncyc] = (
                d_delta[ncyc] @ p90[0] @ d_neg @ cpmg2 @ middle @ cpmg1 @ d_neg @ p90[0]
            )

        intst = {
            ncyc: spectrometer.detect(inept @ zfilter @ centre[ncyc] @ start)
            for ncyc in set(ncycs)
        }

        if self.settings.ipap_flg:
            intst = {
                ncyc: intst[ncyc]
                + spectrometer.detect(centre[ncyc] @ zfilter @ inept @ start)
                for ncyc in set(ncycs)
            }

        # Return profile
        return np.array([intst[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: NDArrayFloat) -> NDArrayBool:
        return metadata == 0.0


def register():
    creators = Creators(
        config_creator=Cpmg1HnAp0013Config,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cpmg1HnAp0013Sequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)

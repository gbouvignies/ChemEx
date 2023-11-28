from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.linalg import matrix_power

from chemex.configuration.data import RelaxationDataSettings
from chemex.configuration.experiment import CpmgSettings, ExperimentConfig, ToBeFitted
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import PlanesFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.plotters.cpmg import CpmgPlotter
from chemex.printers.data import CpmgPrinter

if TYPE_CHECKING:
    from chemex.containers.data import Data
    from chemex.parameters.spin_system import SpinSystem
    from chemex.typing import ArrayBool, ArrayFloat


EXPERIMENT_NAME = "cpmg_15n_ip"


class Cpmg15NIpSettings(CpmgSettings):
    name: Literal["cpmg_15n_ip"]
    time_t2: float
    carrier: float
    pw90: float
    time_equil: float = 0.0
    observed_state: Literal["a", "b", "c", "d"] = "a"

    @property
    def t_neg(self) -> float:
        return -2.0 * self.pw90 / np.pi

    @property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class Cpmg15NIpConfig(ExperimentConfig[Cpmg15NIpSettings, RelaxationDataSettings]):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(rates=[f"r2_i_{state}"], model_free=[f"tauc_{state}"])


def build_spectrometer(
    config: Cpmg15NIpConfig,
    spin_system: SpinSystem,
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = 1 / (4.0 * settings.pw90)
    spectrometer.detection = settings.detection

    return spectrometer


@dataclass
class Cpmg15NIpSequence:
    settings: Cpmg15NIpSettings

    def _get_delays(self, ncycs: ArrayFloat) -> tuple[dict[float, float], list[float]]:
        ncycs_no_ref = ncycs[ncycs > 0]
        tau_cps = {
            ncyc: self.settings.time_t2 / (4.0 * ncyc) - self.settings.pw90
            for ncyc in ncycs_no_ref
        }
        # An ncyc of -1 corresponds to the experiment without the CPMG
        # refocusisng pulses
        tau_cps[-1.0] = 0.5 * self.settings.time_t2
        delays = [self.settings.t_neg, self.settings.time_equil, *tau_cps.values()]
        return tau_cps, delays

    def calculate(self, spectrometer: Spectrometer, data: Data) -> ArrayFloat:
        ncycs = data.metadata

        # Calculation of the spectrometers corresponding to all the delays
        tau_cps, all_delays = self._get_delays(ncycs)
        delays = dict(zip(all_delays, spectrometer.delays(all_delays), strict=True))
        d_neg = delays[self.settings.t_neg]
        d_eq = delays[self.settings.time_equil]
        d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}

        # Calculation of the spectrometers corresponding to all the pulses
        p90 = spectrometer.p90_i
        p180 = spectrometer.p180_i
        p180pmx = 0.5 * (p180[0] + p180[2])  # +/- phase cycling

        # Getting the starting magnetization
        start = spectrometer.get_equilibrium()

        # Calculating the instensities as a function of ncyc
        part1 = d_neg @ p90[0] @ start
        part2 = d_eq @ p90[0] @ d_neg
        intensities = {
            0.0: spectrometer.detect(d_eq @ p90[0] @ p180pmx @ p90[0] @ start),
            -1.0: spectrometer.detect(
                part2 @ d_cp[-1.0] @ p180pmx @ d_cp[-1.0] @ part1,
            ),
        }
        for ncyc in set(ncycs) - {0.0, -1.0}:
            echo = d_cp[ncyc] @ p180[1] @ d_cp[ncyc]
            cpmg = matrix_power(echo, int(ncyc))
            intensities[ncyc] = spectrometer.detect(
                part2 @ cpmg @ p180pmx @ cpmg @ part1,
            )

        # Return profile
        return np.array([intensities[ncyc] for ncyc in ncycs])

    @staticmethod
    def is_reference(metadata: ArrayFloat) -> ArrayBool:
        return metadata == 0


def register() -> None:
    creators = Creators(
        config_creator=Cpmg15NIpConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cpmg15NIpSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=PlanesFilterer,
        printer_creator=CpmgPrinter,
        plotter_creator=CpmgPlotter,
    )
    factories.register(type=EXPERIMENT_NAME, creators=creators)

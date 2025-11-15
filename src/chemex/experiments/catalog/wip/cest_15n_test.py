from __future__ import annotations

# pyright: reportGeneralTypeIssues=false
from dataclasses import dataclass
from functools import cached_property
from typing import Any, Literal

import numpy as np
from scipy.optimize import OptimizeResult, minimize_scalar

from chemex.configuration.base import ExperimentConfiguration, ToBeFitted
from chemex.configuration.conditions import ConditionsWithValidations
from chemex.configuration.data import CestDataSettings
from chemex.configuration.experiment import CestSettings
from chemex.containers.data import Data
from chemex.containers.dataset import load_relaxation_dataset
from chemex.experiments.factories import Creators, factories
from chemex.filterers import CestFilterer
from chemex.nmr.basis import Basis
from chemex.nmr.constants import get_multiplet
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.cest import CestPlotter
from chemex.printers.data import CestPrinter
from chemex.typing import Array

EXPERIMENT_NAME = "cest_15n_test"

OFFSET_REF = 1e4


class Cest15NSettings(CestSettings):
    name: Literal["cest_15n_test"]
    time_t1: float
    carrier: float
    b1_frq: float
    b1_inh_scale: float = 0.1
    b1_inh_res: int = 11

    @cached_property
    def start_terms(self) -> list[str]:
        return [f"iz{self.suffix}"]

    @cached_property
    def detection(self) -> str:
        return f"[iz_{self.observed_state}]"


class Cest15NTestConfig(
    ExperimentConfiguration[
        Cest15NSettings, ConditionsWithValidations, CestDataSettings
    ],
):
    @property
    def to_be_fitted(self) -> ToBeFitted:
        state = self.experiment.observed_state
        return ToBeFitted(
            rates=["r2_i", f"r1_i_{state}"],
            model_free=[f"tauc_{state}", f"s2_{state}"],
        )


def build_spectrometer(
    config: Cest15NTestConfig, spin_system: SpinSystem
) -> Spectrometer:
    settings = config.experiment
    conditions = config.conditions

    basis = Basis(type="ixyz", spin_system="nh")
    liouvillian = LiouvillianIS(spin_system, basis, conditions)
    spectrometer = Spectrometer(liouvillian)

    spectrometer.carrier_i = settings.carrier
    spectrometer.b1_i = settings.b1_frq
    spectrometer.b1_i_inh_scale = settings.b1_inh_scale
    spectrometer.b1_i_inh_res = settings.b1_inh_res

    spectrometer.detection = settings.detection

    if "13c" in conditions.label:
        liouvillian.jeff_i = get_multiplet("", "n")

    return spectrometer


@dataclass
class Cest15NTestSequence:
    settings: Cest15NSettings

    @staticmethod
    def is_reference(metadata: Array) -> Array:
        return np.abs(metadata) > OFFSET_REF

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        offsets = data.metadata

        start = spectrometer.get_start_magnetization(self.settings.start_terms)

        intensities: dict[float, Array] = {}

        for offset in set(offsets):
            intensities[offset] = start

            if self.is_reference(offset):
                continue

            spectrometer.offset_i = offset

            intensities[offset] = (
                spectrometer.pulse_i(self.settings.time_t1, 0.0) @ intensities[offset]
            )

        params = ExchangeParameters(
            kab=spectrometer.par_values["kab"],
            kba=spectrometer.par_values["kba"],
            m0a=spectrometer.par_values["m0a"],
            m0b=spectrometer.par_values["m0b"],
            ra=spectrometer.par_values["ra"],
            rb=spectrometer.par_values["rb"],
        )

        return np.array(
            [spectrometer.detect_ls_2st(intensities[offset]) for offset in offsets],
        )


def find_lineshape_maximum(
    par_values: dict[str, float],
    intensities: Array,
    observed_state: str,
    ppm_i: float,
) -> dict[str, Any]:
    """Find the closest maximum of the lineshape from a starting frequency.

    This function uses numerical optimization to find the closest maximum
    of the two-state exchange lineshape from a given starting point.

    Parameters
    ----------
    params : ExchangeParameters
        Exchange system parameters including rate constants, magnetizations,
        relaxation rates, and chemical shift frequencies
    omega_start : float
        Starting frequency for the search (Hz)
    search_range : float, optional
        Search range around the starting point (Hz), by default 1000.0

    Returns
    -------
    dict[str, Any]
        Dictionary containing:
        - 'omega_max': frequency of the maximum (Hz)
        - 'lineshape_max': lineshape value at the maximum
        - 'success': whether the optimization was successful
        - 'message': optimization message

    """

    def negative_lineshape_real(omega: float) -> float:
        """Negative real part of lineshape for maximization."""
        # Calculate alpha values
        alpha_a_val = par_values["r2_i_a"] + 1j * (omega - par_values["cs_i_a"] * ppm_i)
        alpha_b_val = par_values["r2_i_b"] + 1j * (omega - par_values["cs_i_b"] * ppm_i)

        # Calculate beta values (diagonal elements)
        beta_a_val = alpha_a_val + par_values["kab"]
        beta_b_val = alpha_b_val + par_values["kba"]

        # Calculate lineshape using the analytical expression
        m0a = intensities[..., 3, 1]
        m0b = intensities[..., 6, 1]
        numerator = m0a * (beta_b_val + par_values["kab"]) + m0b * (
            beta_a_val + par_values["kba"]
        )
        denominator = beta_a_val * beta_b_val - par_values["kab"] * par_values["kba"]

        # Return negative real part for maximization
        lineshape = numerator / denominator
        return -lineshape.real

    # Define search bounds
    bracket = (params.omega_a, params.omega_b)

    try:
        # Use scipy's minimize_scalar with bounded method
        result: OptimizeResult = minimize_scalar(
            negative_lineshape_real,
            bracket=bracket,
            options={"xatol": 1e-8},
        )
    except (ValueError, RuntimeError, ArithmeticError) as e:
        return {
            "omega_max": omega_start,
            "lineshape_max": 0.0,
            "success": False,
            "message": f"Optimization failed: {e!s}",
        }
    else:
        omega_max = result.x
        lineshape_max = -result.fun  # Convert back from negative

        return {
            "omega_max": omega_max,
            "lineshape_max": lineshape_max,
            "success": result.success,
            "message": (
                "Optimization completed successfully"
                if result.success
                else result.message
            ),
        }


def register() -> None:
    creators = Creators(
        config_creator=Cest15NTestConfig,
        spectrometer_creator=build_spectrometer,
        sequence_creator=Cest15NTestSequence,
        dataset_creator=load_relaxation_dataset,
        filterer_creator=CestFilterer,
        printer_creator=CestPrinter,
        plotter_creator=CestPlotter,
    )
    factories.register(name=EXPERIMENT_NAME, creators=creators)

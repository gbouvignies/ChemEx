from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from functools import cached_property
from typing import Protocol, Self, runtime_checkable

from lmfit import Parameters as ParametersLF

from chemex.containers.data import Data
from chemex.evaluation import ProfileEvaluator
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.printers.data import Printer
from chemex.typing import Array


@runtime_checkable
class PulseSequence(Protocol):
    """Defines a protocol for pulse sequence classes."""

    def calculate(self, spectrometer: Spectrometer, data: Data) -> Array:
        """Calculate the Array given a Spectrometer and Data."""
        ...

    def is_reference(self, metadata: Array) -> Array:
        """Check if the metadata is a reference."""
        ...


class Filterer(Protocol):
    """Defines a protocol for filterer classes."""

    def filter(self, data: Data) -> None:
        """Filter the given Data."""
        ...


@dataclass(order=True)
class Profile:
    """Represents the profile for a specific experiment."""

    data: Data = field(compare=False)
    spectrometer: Spectrometer = field(compare=False)
    pulse_sequence: PulseSequence = field(compare=False)
    name_map: dict[str, str] = field(compare=False)
    printer: Printer = field(compare=False)
    filterer: Filterer | None = field(compare=False, default=None)
    is_scaled: bool = field(compare=False, default=True)
    evaluator: ProfileEvaluator = field(compare=False, default_factory=ProfileEvaluator)
    spin_system: SpinSystem = field(compare=True, init=False)

    def __post_init__(self) -> None:
        """Initialize derived attributes."""
        self.spin_system = self.spectrometer.liouvillian.spin_system
        self.data.refs = self.pulse_sequence.is_reference(self.data.metadata)

    @cached_property
    def param_ids(self) -> set[str]:
        """Get the set of parameter IDs."""
        return set(self.name_map.values())

    @cached_property
    def cache_param_ids(self) -> tuple[str, ...]:
        """Get the stable ordered parameter ids used in residual cache keys."""
        return tuple(sorted(self.param_ids))

    def clear_cache(self) -> None:
        """Drop cached residuals after the profile inputs change."""
        self.evaluator.clear()

    def _get_parameter_values(self, params: ParametersLF) -> dict[str, float]:
        """Get the parameter values from the provided parameters."""
        return {
            local_name: params[param_id].value
            for local_name, param_id in self.name_map.items()
        }

    def update_spectrometer(self, params: ParametersLF) -> None:
        """Update the spectrometer settings."""
        parameter_values = self._get_parameter_values(params)
        self.spectrometer.update(parameter_values)

    def calculate(self, params: ParametersLF) -> Array:
        """Calculate and return the Array."""
        self.update_spectrometer(params)
        self.data.calc_unscaled = self.pulse_sequence.calculate(
            self.spectrometer, self.data
        )
        if self.is_scaled:
            self.data.calc = self.data.scale * self.data.calc_unscaled
        else:
            self.data.calc = self.data.calc_unscaled
        return self.data.calc

    def residuals(self, params: ParametersLF) -> Array:
        """Calculate and return residuals."""
        return self.evaluator.residuals(
            params,
            param_ids=self.cache_param_ids,
            data_revision=self.data.revision,
            calculate=self.calculate,
            exp=self.data.exp,
            err=self.data.err,
            mask=self.data.mask,
        )

    def filter(self, params: ParametersLF) -> None:
        """Apply a filter to the data, if a filterer is available."""
        if self.filterer is not None:
            self.update_spectrometer(params)
            self.filterer.filter(self.data)
            self.data.mark_dirty()
            self.clear_cache()

    def set_noise(self, value: float) -> None:
        """Set the noise value."""
        self.data.err[:] = value
        self.data.mark_dirty()
        self.clear_cache()

    def prepare_for_simulation(self) -> None:
        """Prepare data for simulation."""
        self.data.exp = self.data.calc
        self.printer.simulation = True
        self.data.mark_dirty()
        self.clear_cache()

    def monte_carlo(self: Self) -> Self:
        """Generate a Monte Carlo variant of the profile."""
        profile = deepcopy(self)
        profile.data = profile.data.monte_carlo()
        profile.evaluator = ProfileEvaluator()
        return profile

    def bootstrap(self) -> Self:
        """Generate a bootstrap variant of the profile."""
        profile = deepcopy(self)
        profile.data = profile.data.bootstrap()
        profile.evaluator = ProfileEvaluator()
        return profile

    def any_duplicate(self) -> bool:
        """Check for duplicate data points."""
        return self.data.any_duplicate()

    def __add__(self, other: Self) -> Self:
        """Combine two profiles."""
        profile = deepcopy(self)
        profile.data = self.data + other.data
        profile.evaluator = ProfileEvaluator()
        return profile

    def __str__(self) -> str:
        """Return string representation of the Profile."""
        return self.printer.print(str(self.spin_system), self.data)

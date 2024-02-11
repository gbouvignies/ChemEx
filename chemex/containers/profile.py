from __future__ import annotations

from collections.abc import Hashable
from copy import deepcopy
from dataclasses import dataclass, field
from functools import cached_property
from operator import attrgetter
from typing import Protocol

from cachetools import LRUCache, cachedmethod
from lmfit import Parameters as ParametersLF
from typing import Self

from chemex.containers.data import Data
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.printers.data import Printer
from chemex.typing import ArrayBool, ArrayFloat


class PulseSequence(Protocol):
    """Defines a protocol for pulse sequence classes."""

    def calculate(self, spectrometer: Spectrometer, data: Data) -> ArrayFloat:
        """Calculate the ArrayFloat given a Spectrometer and Data."""
        ...

    def is_reference(self, metadata: ArrayFloat) -> ArrayBool:
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
    spin_system: SpinSystem = field(compare=True, init=False)
    cache: LRUCache = field(compare=False, init=False)

    def _cache_key(self, params: ParametersLF) -> tuple[Hashable, ...]:
        return (
            *(params[param_id].value for param_id in self.param_ids),
            self.data.metadata.tostring(),
        )

    def __post_init__(self):
        """Initialize derived attributes."""
        self.spin_system = self.spectrometer.liouvillian.spin_system
        self.data.refs = self.pulse_sequence.is_reference(self.data.metadata)
        self.cache = LRUCache(maxsize=5)

    @cached_property
    def param_ids(self) -> set[str]:
        """Get the set of parameter IDs."""
        return set(self.name_map.values())

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

    def calculate(self, params: ParametersLF) -> ArrayFloat:
        """Calculate and return the ArrayFloat."""
        self.update_spectrometer(params)
        self.data.calc_unscaled = self.pulse_sequence.calculate(
            self.spectrometer, self.data
        )
        if self.is_scaled:
            self.data.calc = self.data.scale * self.data.calc_unscaled
        else:
            self.data.calc = self.data.calc_unscaled
        return self.data.calc

    @cachedmethod(attrgetter("cache"), key=_cache_key)
    def residuals(self, params: ParametersLF) -> ArrayFloat:
        """Calculate and return residuals."""
        residuals = (self.calculate(params) - self.data.exp) / self.data.err
        return residuals[self.data.mask]

    def filter(self, params: ParametersLF) -> None:
        """Apply a filter to the data, if a filterer is available."""
        if self.filterer is not None:
            self.update_spectrometer(params)
            self.filterer.filter(self.data)

    def set_noise(self, value: float) -> None:
        """Set the noise value."""
        self.data.err[:] = value

    def prepare_for_simulation(self) -> None:
        """Prepare data for simulation."""
        self.data.exp = self.data.calc
        self.printer.simulation = True

    def monte_carlo(self: Self) -> Self:
        """Generate a Monte Carlo variant of the profile."""
        profile = deepcopy(self)
        profile.data = profile.data.monte_carlo()
        return profile

    def bootstrap(self) -> Self:
        """Generate a bootstrap variant of the profile."""
        profile = deepcopy(self)
        profile.data = profile.data.bootstrap()
        return profile

    def any_duplicate(self):
        """Check for duplicate data points."""
        return self.data.any_duplicate()

    def __add__(self, other: Self) -> Self:
        """Combine two profiles."""
        profile = deepcopy(self)
        profile.data = self.data + other.data
        return profile

    def __str__(self) -> str:
        """String representation of the Profile."""
        return self.printer.print(str(self.spin_system), self.data)

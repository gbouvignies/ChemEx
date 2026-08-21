from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from functools import cached_property
from typing import Protocol, Self, runtime_checkable

import numpy as np

from chemex.containers.data import Data
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
    spin_system: SpinSystem = field(compare=True, init=False)

    def __post_init__(self) -> None:
        """Initialize derived attributes."""
        self.spin_system = self.spectrometer.spin_system
        self.data.refs = self.pulse_sequence.is_reference(self.data.metadata)

    @cached_property
    def param_ids(self) -> set[str]:
        """Get the set of parameter IDs."""
        return set(self.name_map.values())

    def _local_parameter_values(
        self, parameter_values: Mapping[str, float]
    ) -> dict[str, float]:
        return {
            local_name: parameter_values[param_id]
            for local_name, param_id in self.name_map.items()
        }

    def update_spectrometer_from_values(
        self,
        parameter_values: Mapping[str, float],
    ) -> None:
        """Update the spectrometer directly from stable native parameter values."""
        self.spectrometer.update(self._local_parameter_values(parameter_values))

    def calculate_unscaled(self, parameter_values: Mapping[str, float]) -> Array:
        """Run only this profile's scientific calculation.

        Scaling, residual construction, masking, and publication remain owned by
        the caller. This is the narrow calculation seam used by native evaluation.
        """
        self.spectrometer.update(self._local_parameter_values(parameter_values))
        # Existing pulse implementations are typed against Data but consume only
        # metadata.  Give them a fresh zeroed carrier so kernels never see this
        # profile's observations, errors, masks, or calculated arrays.
        kernel_data = Data(
            exp=np.zeros_like(self.data.exp),
            err=np.zeros_like(self.data.err),
            # Kernel metadata is a per-call carrier, never an alias of the
            # authoritative profile data.  Native evaluation owns all
            # observations, uncertainties, masks, and residual semantics.
            metadata=np.array(self.data.metadata, copy=True),
        )
        return self.pulse_sequence.calculate(self.spectrometer, kernel_data)

    def calculate_from_values(
        self,
        parameter_values: Mapping[str, float],
    ) -> Array:
        """Calculate and publish data from resolved native parameter values."""
        self.data.calc_unscaled = self.calculate_unscaled(parameter_values)
        if self.is_scaled:
            self.data.calc = self.data.scale * self.data.calc_unscaled
        else:
            self.data.calc = self.data.calc_unscaled
        return self.data.calc

    def filter_from_values(self, parameter_values: Mapping[str, float]) -> None:
        """Apply the configured filter using resolved native parameter values."""
        if self.filterer is not None:
            self.update_spectrometer_from_values(parameter_values)
            self.filterer.filter(self.data)
            self.data.mark_dirty()

    def set_noise(self, value: float) -> None:
        """Set the noise value."""
        self.data.err[:] = value
        self.data.mark_dirty()

    def prepare_for_simulation(self) -> None:
        """Prepare data for simulation."""
        self.data.exp = self.data.calc
        self.printer.simulation = True
        self.data.mark_dirty()

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

    def any_duplicate(self) -> bool:
        """Check for duplicate data points."""
        return self.data.any_duplicate()

    def __add__(self, other: Self) -> Self:
        """Combine two profiles."""
        profile = deepcopy(self)
        profile.data = self.data + other.data
        return profile

    def __str__(self) -> str:
        """Return string representation of the Profile."""
        return self.printer.print(str(self.spin_system), self.data)

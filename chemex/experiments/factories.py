"""Factories for creating different parts of an experiment."""
from __future__ import annotations

from collections.abc import Callable, MutableMapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar

from chemex.configuration.base import ExperimentConfiguration
from chemex.containers.data import Data
from chemex.containers.profile import Filterer, PulseSequence
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.plotter import Plotter
from chemex.printers.data import Printer

Dataset = list[tuple[SpinSystem, Data]]
GenericConfig = ExperimentConfiguration[Any, Any, Any]
ConfigCreator = Callable[..., GenericConfig]
PropagatorCreator = Callable[..., Spectrometer]
SequenceCreator = Callable[..., PulseSequence]
DatasetCreator = Callable[..., Dataset]
FiltererCreator = Callable[..., Filterer]
PrinterCreator = Callable[[], Printer]
PlotterCreator = Callable[..., Plotter]


@dataclass
class Creators:
    config_creator: ConfigCreator
    spectrometer_creator: PropagatorCreator
    sequence_creator: SequenceCreator
    dataset_creator: DatasetCreator
    filterer_creator: FiltererCreator
    printer_creator: PrinterCreator
    plotter_creator: PlotterCreator

    def create_config(
        self,
        config_dict: MutableMapping[str, Any],
    ) -> GenericConfig:
        """Create a configuration object of a specific type."""
        return self.config_creator(**config_dict)

    def create_spectrometer(
        self,
        config: GenericConfig,
        spin_system: SpinSystem,
    ) -> Spectrometer:
        """Create and initialize a spectrometer of a specific type."""
        return self.spectrometer_creator(config, spin_system)

    def create_sequence(self, config: GenericConfig) -> PulseSequence:
        """Create a sequence of a specific type."""
        return self.sequence_creator(config.experiment)

    def create_dataset(self, base_path: Path, settings: GenericConfig) -> Dataset:
        """Create a dataset of a specific type."""
        return self.dataset_creator(base_path, settings)

    def create_filterer(
        self,
        config: GenericConfig,
        spectrometer: Spectrometer,
    ) -> Filterer:
        """Create a filterer used to remove undesired data points."""
        return self.filterer_creator(config=config, spectrometer=spectrometer)

    def create_printer(self) -> Printer:
        """Create a filterer used to remove undesired data points."""
        return self.printer_creator()

    def create_plotter(self, filename: Path, config: GenericConfig) -> Plotter:
        """Create a filterer used to remove undesired data points."""
        return self.plotter_creator(filename=filename, config=config)


class Factories:
    """Factory for creating all parts of an experiment."""

    creators_registry: ClassVar[dict[str, Creators]] = {}

    def register(self, type: str, creators: Creators):
        """Register a new propagtor type."""
        self.creators_registry[type] = creators

    def get(self, type: str) -> Creators:
        try:
            return self.creators_registry[type]
        except KeyError:
            msg = f"Unknown  type {type!r}"
            raise ValueError(msg) from None


factories = Factories()

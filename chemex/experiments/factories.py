"""Factories for creating different parts of an experiment."""
from __future__ import annotations

from collections.abc import Callable
from collections.abc import MutableMapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from chemex.configuration.experiment import ExperimentConfig
from chemex.containers.dataset import Dataset
from chemex.containers.profile import Filterer
from chemex.containers.profile import PulseSequence
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.plotter import Plotter
from chemex.printers.data import Printer

ConfigCreator = Callable[..., ExperimentConfig]
PropagatorCreator = Callable[..., Spectrometer]
SequenceCreator = Callable[..., PulseSequence]
DatasetCreator = Callable[[Path, ExperimentConfig], Dataset]
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

    def create_config(self, config_dict: MutableMapping[str, Any]) -> ExperimentConfig:
        """Create a configuration object of a specific type."""
        return self.config_creator(**config_dict)

    def create_spectrometer(
        self, config: ExperimentConfig, spin_system: SpinSystem
    ) -> Spectrometer:
        """Create and initialize a spectrometer of a specific type."""
        return self.spectrometer_creator(config, spin_system)

    def create_sequence(self, config: ExperimentConfig) -> PulseSequence:
        """Create a sequence of a specific type."""
        return self.sequence_creator(config.experiment)

    def create_dataset(self, base_path: Path, settings: ExperimentConfig) -> Dataset:
        """Create a dataset of a specific type."""
        return self.dataset_creator(base_path, settings)

    def create_filterer(
        self, config: ExperimentConfig, spectrometer: Spectrometer
    ) -> Filterer:
        """Create a filterer used to remove undesired data points."""
        return self.filterer_creator(config=config, spectrometer=spectrometer)

    def create_printer(self) -> Printer:
        """Create a filterer used to remove undesired data points."""
        return self.printer_creator()

    def create_plotter(self, filename: Path, config: ExperimentConfig) -> Plotter:
        """Create a filterer used to remove undesired data points."""
        return self.plotter_creator(filename=filename, config=config)


class Factories:
    """Factory for creating all parts of an experiment."""

    creators_registry: dict[str, Creators] = {}

    def register(self, type: str, creators: Creators):
        """Register a new propagtor type."""
        self.creators_registry[type] = creators

    def get(self, type: str) -> Creators:
        try:
            return self.creators_registry[type]
        except KeyError:
            raise ValueError(f"Unknown  type {type!r}") from None


factories = Factories()

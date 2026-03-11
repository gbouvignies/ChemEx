"""Factories for creating different parts of an experiment."""

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar

from chemex.configuration.base import ExperimentConfiguration
from chemex.containers.data import Data
from chemex.containers.profile import Filterer, PulseSequence
from chemex.models.model import ModelSpec
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.plotter import Plotter
from chemex.printers.data import Printer

Dataset = list[tuple[SpinSystem, Data]]
GenericConfig = ExperimentConfiguration[Any, Any, Any]
ConfigCreator = type[GenericConfig]
PropagatorCreator = Callable[..., Spectrometer]
SequenceCreator = Callable[[Any], object]
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
        config_dict: Mapping[str, object],
        *,
        model: ModelSpec,
    ) -> GenericConfig:
        """Create a configuration object of a specific type."""
        config = self.config_creator.model_validate(
            config_dict,
            context={"model": model},
        )
        config.model = model
        config.experiment.model_name = model.name
        return config

    def create_spectrometer(
        self,
        config: GenericConfig,
        spin_system: SpinSystem,
    ) -> Spectrometer:
        """Create and initialize a spectrometer of a specific type."""
        return self.spectrometer_creator(config, spin_system)

    def create_sequence(self, config: GenericConfig) -> PulseSequence:
        """Create a sequence of a specific type."""
        sequence = self.sequence_creator(config.experiment)
        if isinstance(sequence, PulseSequence):
            return sequence
        msg = "Sequence creator must return an object implementing PulseSequence"
        raise TypeError(msg)

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

    def register(self, name: str, creators: Creators) -> None:
        """Register a new propagtor type."""
        self.creators_registry[name] = creators

    def get(self, name: str) -> Creators:
        try:
            return self.creators_registry[name]
        except KeyError:
            msg = f"Unknown  type {name!r}"
            raise ValueError(msg) from None


factories = Factories()

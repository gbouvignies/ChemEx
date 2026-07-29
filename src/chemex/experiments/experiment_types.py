"""Deep interface for defining and constructing Experiment Types."""

from __future__ import annotations

import tomllib
from collections.abc import Callable, Mapping
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Protocol, Self, TypeVar
from weakref import WeakKeyDictionary

from pydantic import ValidationError

from chemex.configuration.base import (
    BaseSettings,
    ExperimentConfiguration,
    ToBeFitted,
)
from chemex.configuration.conditions import Conditions
from chemex.configuration.data import (
    CestDataSettings,
    RelaxationDataSettings,
    ShiftDataSettings,
)
from chemex.configuration.methods import Selection
from chemex.containers.dataset import (
    Dataset,
    DatasetLoadError,
    load_relaxation_dataset,
    load_shift_dataset,
)
from chemex.containers.experiment import Experiment, NoiseEstimateNotice
from chemex.containers.profile import Filterer, Profile, PulseSequence
from chemex.filterers import (
    CestExperimentSettings,
    CestFilterer,
    NoFilterer,
    PlanesFilterer,
)
from chemex.models.model import ModelSpec
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.factory import ParameterFactory
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.cest import CestPlotter
from chemex.plotters.cpmg import CpmgExperimentSettings, CpmgPlotter
from chemex.plotters.plotter import Plotter
from chemex.plotters.relaxation import RelaxationPlotter
from chemex.plotters.shift import ShiftPlotter
from chemex.printers.data import (
    CestPrinter,
    CpmgPrinter,
    Printer,
    RelaxationPrinter,
    ShiftPrinter,
)
from chemex.toml import load_toml

GenericConfig = ExperimentConfiguration[Any, Any, Any]
ConfigT = TypeVar("ConfigT")


@dataclass(frozen=True)
class ProfileCalculation:
    """Experiment-specific machinery used to calculate a Profile."""

    spectrometer: Spectrometer
    pulse_sequence: PulseSequence

    def __post_init__(self) -> None:
        if not isinstance(self.pulse_sequence, PulseSequence):
            message = "Profile calculation must provide a pulse sequence"
            raise TypeError(message)


class ProfileCalculationFactory(Protocol[ConfigT]):
    """Strongly typed catalog adapter for per-profile calculation construction."""

    def __call__(
        self,
        config: ConfigT,
        spin_system: SpinSystem,
    ) -> ProfileCalculation: ...


class DatasetLoader(Protocol[ConfigT]):
    def __call__(self, base_path: Path, settings: ConfigT) -> Dataset: ...


class FiltererFactory(Protocol[ConfigT]):
    def __call__(
        self,
        *,
        config: ConfigT,
        spectrometer: Spectrometer,
    ) -> Filterer: ...


class PrinterFactory(Protocol):
    def __call__(self) -> Printer: ...


class PlotterFactory(Protocol[ConfigT]):
    def __call__(
        self,
        *,
        filename: Path,
        config: ConfigT,
    ) -> Plotter: ...


@dataclass(frozen=True)
class ExperimentSupport[ConfigT]:
    """Repeated data and presentation policy for an Experiment Type."""

    load_dataset: DatasetLoader[ConfigT]
    create_filterer: FiltererFactory[ConfigT]
    create_printer: PrinterFactory
    create_plotter: PlotterFactory[ConfigT]


class BuildExperimentSettings(Protocol):
    """Experiment settings mutated by generic construction."""

    @property
    def model_name(self) -> str: ...

    @model_name.setter
    def model_name(self, value: str) -> None: ...


class BuildDataSettings(Protocol):
    """Data settings read by generic construction."""

    @property
    def error(self) -> str: ...

    @property
    def global_error(self) -> bool: ...

    @property
    def scaled(self) -> bool: ...


class ExperimentConfigContract(Protocol):
    """Common configuration behavior required by generic construction."""

    model: ModelSpec

    @classmethod
    def model_validate(
        cls,
        obj: object,
        *,
        context: object | None = None,
    ) -> Self: ...

    @property
    def conditions(self) -> Conditions: ...

    @property
    def data(self) -> BuildDataSettings: ...

    @property
    def experiment(self) -> BuildExperimentSettings: ...

    @property
    def to_be_fitted(self) -> ToBeFitted: ...


@dataclass(frozen=True)
class ExperimentType[ConfigT: ExperimentConfigContract]:
    """Typed adapter connecting catalog calculations to generic construction."""

    name: str
    config_type: type[ConfigT]
    create_profile_calculation: ProfileCalculationFactory[ConfigT]
    support: ExperimentSupport[ConfigT]


class ExperimentSource:
    """Opaque handle returned by :func:`open` for subsequent construction."""

    __slots__ = ("__weakref__", "_experiment_type_name", "_filename")
    _experiment_type_name: str
    _filename: Path

    def __init__(self) -> None:
        message = "ExperimentSource objects are created by experiment_types.open()"
        raise TypeError(message)

    def __setattr__(self, _name: str, _value: object) -> None:
        message = "ExperimentSource objects are immutable"
        raise AttributeError(message)

    @property
    def filename(self) -> Path:
        """Input filename associated with this source."""
        return self._filename

    @property
    def experiment_type_name(self) -> str:
        """Registered Experiment Type name identified by :func:`open`."""
        return self._experiment_type_name


@dataclass(frozen=True)
class _SourceState:
    raw_config: Mapping[str, object]
    experiment_type: ExperimentType[Any]


@dataclass(frozen=True)
class ExperimentBuildResult:
    """A complete Experiment and nonfatal notices produced while building it."""

    experiment: Experiment
    notices: tuple[NoiseEstimateNotice, ...] = ()


type ExperimentNotice = NoiseEstimateNotice


class ExperimentBuildError(Exception):
    """Base class for expected Experiment construction failures."""


@dataclass(eq=False)
class InvalidExperimentSourceError(ExperimentBuildError):
    filename: Path | None
    explanation: str

    def __post_init__(self) -> None:
        super().__init__(self.explanation)


@dataclass(eq=False)
class ExperimentFileError(ExperimentBuildError):
    filename: Path
    error: FileNotFoundError

    def __post_init__(self) -> None:
        super().__init__(str(self.error))


@dataclass(eq=False)
class ExperimentTomlError(ExperimentBuildError):
    filename: Path
    error: tomllib.TOMLDecodeError | TypeError

    def __post_init__(self) -> None:
        super().__init__(str(self.error))


@dataclass(eq=False)
class ExperimentNameError(ExperimentBuildError):
    filename: Path
    error: ValidationError

    def __post_init__(self) -> None:
        super().__init__(str(self.error))


@dataclass(eq=False)
class UnknownExperimentTypeError(ExperimentBuildError):
    filename: Path
    experiment_type_name: str

    def __post_init__(self) -> None:
        message = f"Unknown Experiment Type {self.experiment_type_name!r}"
        super().__init__(message)


@dataclass(eq=False)
class ExperimentConfigurationError(ExperimentBuildError):
    filename: Path
    error: ValidationError

    def __post_init__(self) -> None:
        super().__init__(str(self.error))


@dataclass(eq=False)
class ExperimentDataError(ExperimentBuildError):
    filename: Path
    error: FileNotFoundError | DatasetLoadError

    def __post_init__(self) -> None:
        super().__init__(str(self.error))


class _NameSettings(BaseSettings):
    name: str


class _NameConfig(BaseSettings):
    experiment: _NameSettings


_REGISTRY: dict[str, ExperimentType[Any]] = {}
_SOURCE_STATES: WeakKeyDictionary[ExperimentSource, _SourceState] = WeakKeyDictionary()


def register_experiment_type(experiment_type: ExperimentType[Any]) -> None:
    """Register an adapter, allowing repeated registration of the same object."""
    registered = _REGISTRY.get(experiment_type.name)
    if registered is experiment_type:
        return
    if registered is not None:
        message = (
            f"Experiment Type {experiment_type.name!r} is already registered "
            "with a different adapter"
        )
        raise ValueError(message)
    _REGISTRY[experiment_type.name] = experiment_type


def registered_experiment_types() -> Mapping[str, ExperimentType[Any]]:
    """Return an immutable snapshot of registered Experiment Types."""
    return MappingProxyType(dict(_REGISTRY))


def _create_source(
    filename: Path,
    experiment_type_name: str,
    raw_config: Mapping[str, object],
    experiment_type_: ExperimentType[Any],
) -> ExperimentSource:
    source = object.__new__(ExperimentSource)
    object.__setattr__(source, "_filename", filename)
    object.__setattr__(source, "_experiment_type_name", experiment_type_name)
    _SOURCE_STATES[source] = _SourceState(deepcopy(raw_config), experiment_type_)
    return source


def _get_source_state(source: ExperimentSource) -> _SourceState:
    if not isinstance(source, ExperimentSource):
        message = "Invalid ExperimentSource; call experiment_types.open() first"
        raise InvalidExperimentSourceError(None, message)

    state = _SOURCE_STATES.get(source)
    if state is None:
        message = (
            "Invalid or stale ExperimentSource; call experiment_types.open() again"
        )
        raise InvalidExperimentSourceError(None, message)

    registered = _REGISTRY.get(state.experiment_type.name)
    if registered is not state.experiment_type:
        message = "Stale ExperimentSource; call experiment_types.open() again"
        raise InvalidExperimentSourceError(source.filename, message)

    return state


def experiment_type[ConfigT: ExperimentConfigContract](
    *,
    name: str,
    config_type: type[ConfigT],
    support: ExperimentSupport[ConfigT],
) -> Callable[[ProfileCalculationFactory[ConfigT]], ExperimentType[ConfigT]]:
    """Create a typed Experiment Type adapter from a profile-calculation function."""

    def decorate(
        create_profile_calculation: ProfileCalculationFactory[ConfigT],
    ) -> ExperimentType[ConfigT]:
        return ExperimentType(
            name=name,
            config_type=config_type,
            create_profile_calculation=create_profile_calculation,
            support=support,
        )

    return decorate


class CpmgBuildExperimentSettings(
    BuildExperimentSettings,
    CpmgExperimentSettings,
    Protocol,
):
    """Generic and scientific settings required by CPMG support."""


class CestBuildExperimentSettings(
    BuildExperimentSettings,
    CestExperimentSettings,
    Protocol,
):
    """Generic and scientific settings required by CEST support."""


class CpmgFamilyConfig(ExperimentConfigContract, Protocol):
    @property
    def data(self) -> RelaxationDataSettings: ...

    @property
    def experiment(self) -> CpmgBuildExperimentSettings: ...


class CestFamilyConfig(ExperimentConfigContract, Protocol):
    @property
    def data(self) -> CestDataSettings: ...

    @property
    def experiment(self) -> CestBuildExperimentSettings: ...


class RelaxationFamilyConfig(ExperimentConfigContract, Protocol):
    @property
    def data(self) -> RelaxationDataSettings: ...

    @property
    def experiment(self) -> BuildExperimentSettings: ...


class ShiftFamilyConfig(ExperimentConfigContract, Protocol):
    @property
    def data(self) -> ShiftDataSettings: ...

    @property
    def experiment(self) -> BuildExperimentSettings: ...


def _cpmg_support[ConfigT: CpmgFamilyConfig]() -> ExperimentSupport[ConfigT]:
    return ExperimentSupport[ConfigT](
        load_dataset=load_relaxation_dataset,
        create_filterer=PlanesFilterer,
        create_printer=CpmgPrinter,
        create_plotter=CpmgPlotter,
    )


def _cest_support[ConfigT: CestFamilyConfig]() -> ExperimentSupport[ConfigT]:
    return ExperimentSupport[ConfigT](
        load_dataset=load_relaxation_dataset,
        create_filterer=CestFilterer,
        create_printer=CestPrinter,
        create_plotter=CestPlotter,
    )


def _relaxation_support[
    ConfigT: RelaxationFamilyConfig,
]() -> ExperimentSupport[ConfigT]:
    return ExperimentSupport[ConfigT](
        load_dataset=load_relaxation_dataset,
        create_filterer=PlanesFilterer,
        create_printer=RelaxationPrinter,
        create_plotter=RelaxationPlotter,
    )


def _shift_support[ConfigT: ShiftFamilyConfig]() -> ExperimentSupport[ConfigT]:
    return ExperimentSupport[ConfigT](
        load_dataset=load_shift_dataset,
        create_filterer=NoFilterer,
        create_printer=ShiftPrinter,
        create_plotter=ShiftPlotter,
    )


def cpmg_type[ConfigT: CpmgFamilyConfig](
    *,
    name: str,
    config_type: type[ConfigT],
) -> Callable[[ProfileCalculationFactory[ConfigT]], ExperimentType[ConfigT]]:
    """Define an Experiment Type using established CPMG support."""
    return experiment_type(
        name=name,
        config_type=config_type,
        support=_cpmg_support(),
    )


def cest_type[ConfigT: CestFamilyConfig](
    *,
    name: str,
    config_type: type[ConfigT],
) -> Callable[[ProfileCalculationFactory[ConfigT]], ExperimentType[ConfigT]]:
    """Define an Experiment Type using established CEST support."""
    return experiment_type(
        name=name,
        config_type=config_type,
        support=_cest_support(),
    )


def relaxation_type[ConfigT: RelaxationFamilyConfig](
    *,
    name: str,
    config_type: type[ConfigT],
) -> Callable[[ProfileCalculationFactory[ConfigT]], ExperimentType[ConfigT]]:
    """Define an Experiment Type using established relaxation support."""
    return experiment_type(
        name=name,
        config_type=config_type,
        support=_relaxation_support(),
    )


def shift_type[ConfigT: ShiftFamilyConfig](
    *,
    name: str,
    config_type: type[ConfigT],
) -> Callable[[ProfileCalculationFactory[ConfigT]], ExperimentType[ConfigT]]:
    """Define an Experiment Type using established shift support."""
    return experiment_type(
        name=name,
        config_type=config_type,
        support=_shift_support(),
    )


def open(filename: Path) -> ExperimentSource:  # noqa: A001
    """Read an Experiment input once and identify its Experiment Type."""
    try:
        raw_config = load_toml(filename)
    except FileNotFoundError as error:
        raise ExperimentFileError(filename, error) from error
    except (tomllib.TOMLDecodeError, TypeError) as error:
        raise ExperimentTomlError(filename, error) from error

    try:
        name_config = _NameConfig.model_validate(raw_config)
    except ValidationError as error:
        raise ExperimentNameError(filename, error) from error

    experiment_type_name = name_config.experiment.name
    experiment_type_ = _REGISTRY.get(experiment_type_name)
    if experiment_type_ is None:
        raise UnknownExperimentTypeError(filename, experiment_type_name)

    return _create_source(
        filename,
        experiment_type_name,
        raw_config,
        experiment_type_,
    )


def _apply_selection(dataset: Dataset, selection: Selection) -> Dataset:
    if (include := selection.include) is not None:
        return [
            (spin_system, data)
            for spin_system, data in dataset
            if spin_system.part_of(include)
        ]

    if (exclude := selection.exclude) is not None:
        return [
            (spin_system, data)
            for spin_system, data in dataset
            if not spin_system.part_of(exclude)
        ]

    return dataset


def _build[ConfigT: GenericConfig](
    source: ExperimentSource,
    raw_config: Mapping[str, object],
    experiment_type_: ExperimentType[ConfigT],
    *,
    selection: Selection,
    model: ModelSpec,
    parameters: ParameterFactory,
) -> ExperimentBuildResult:
    try:
        config = experiment_type_.config_type.model_validate(
            raw_config,
            context={"model": model},
        )
    except ValidationError as error:
        raise ExperimentConfigurationError(source.filename, error) from error

    config.model = model
    config.experiment.model_name = model.name

    printer = experiment_type_.support.create_printer()
    plotter = experiment_type_.support.create_plotter(
        filename=source.filename,
        config=config,
    )

    try:
        dataset = experiment_type_.support.load_dataset(source.filename.parent, config)
    except (FileNotFoundError, DatasetLoadError) as error:
        raise ExperimentDataError(source.filename, error) from error
    selected_dataset = _apply_selection(dataset, selection)

    profiles: list[Profile] = []
    for spin_system, data in selected_dataset:
        calculation = experiment_type_.create_profile_calculation(config, spin_system)
        filterer = experiment_type_.support.create_filterer(
            config=config,
            spectrometer=calculation.spectrometer,
        )
        name_map = parameters.create_parameters(
            config,
            basis=calculation.spectrometer.basis,
            spin_system=calculation.spectrometer.spin_system,
        )
        profiles.append(
            Profile(
                data,
                calculation.spectrometer,
                calculation.pulse_sequence,
                name_map,
                printer,
                filterer,
                config.data.scaled,
            ),
        )

    experiment = Experiment(
        source.filename,
        experiment_type_.name,
        profiles,
        printer,
        plotter,
        parameter_store=parameters.parameter_store,
    )
    experiment.estimate_noise(
        config.data.error,
        global_error=config.data.global_error,
    )
    return ExperimentBuildResult(experiment, experiment.noise_notices)


def build(
    source: ExperimentSource,
    *,
    selection: Selection,
    model: ModelSpec,
    parameters: ParameterFactory,
) -> ExperimentBuildResult:
    """Construct a complete Experiment from an opened input."""
    state = _get_source_state(source)
    return _build(
        source,
        state.raw_config,
        state.experiment_type,
        selection=selection,
        model=model,
        parameters=parameters,
    )

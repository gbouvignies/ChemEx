from __future__ import annotations

from typing import Any, Protocol, TypeVar

from chemex.containers.data import Data
from chemex.nmr.spectrometer import Spectrometer


class NoFilterer:
    """Filterer that does no filtering."""

    def __init__(self, **_extra: Any) -> None:
        return

    def filter(self, data: Data) -> None:
        return


class PlanesDataSettings(Protocol):
    @property
    def filter_planes(self) -> list[int]: ...


class PlanesExperimentConfig(Protocol):
    @property
    def data(self) -> PlanesDataSettings: ...


def _filter_planes(data: Data, planes_to_filter: list[int]) -> None:
    selection = [i for i in planes_to_filter if 0 <= i < data.size]
    data.mask[selection] = False


class PlanesFilterer:
    def __init__(self, config: PlanesExperimentConfig, **_extra: Any) -> None:
        self.config = config

    def filter(self, data: Data) -> None:
        _filter_planes(data, self.config.data.filter_planes)


class CestExperimentSettings(Protocol):
    @property
    def primary_state(self) -> str: ...

    @property
    def sw(self) -> float: ...


class CestDataSettings(Protocol):
    @property
    def filter_planes(self) -> list[int]: ...

    @property
    def filter_offsets(self) -> list[tuple[float, float]]: ...

    @property
    def filter_ref_planes(self) -> bool: ...


class CestExperimentConfig(Protocol):
    @property
    def experiment(self) -> CestExperimentSettings: ...

    @property
    def data(self) -> CestDataSettings: ...


T = TypeVar("T", bound=CestExperimentConfig)


def _filter_ref_planes(data: Data) -> None:
    data.mask[data.refs] = False


def _filter_offsets(
    data: Data,
    config: CestExperimentConfig,
    spectrometer: Spectrometer,
) -> None:
    state = config.experiment.primary_state
    state_ppm = spectrometer.par_values[f"cs_i_{state}"]
    state_offset = spectrometer.ppms_to_offsets(state_ppm)

    offsets = data.metadata
    offset_min = min(offsets[~data.refs])

    for offset_to_filter, filter_bandwidth in config.data.filter_offsets:
        absolute_offset = offset_to_filter + state_offset
        offset_with_aliasing = (
            absolute_offset - offset_min
        ) % config.experiment.sw + offset_min
        shifted_offsets = offsets - offset_with_aliasing
        mask_filter = abs(shifted_offsets) < filter_bandwidth * 0.5
        data.mask[mask_filter] = False


class CestFilterer[T: CestExperimentConfig]:
    def __init__(self, config: T, spectrometer: Spectrometer, **_extra: Any) -> None:
        self.config = config
        self.spectrometer = spectrometer

    def filter(self, data: Data) -> None:
        _filter_planes(data, self.config.data.filter_planes)

        if self.config.data.filter_ref_planes:
            _filter_ref_planes(data)

        _filter_offsets(data, self.config, self.spectrometer)

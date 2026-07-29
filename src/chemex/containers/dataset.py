from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import numpy as np

from chemex.configuration.data import RelaxationDataSettings, ShiftDataSettings
from chemex.containers.data import Data
from chemex.parameters.spin_system import SpinSystem
from chemex.toml import normalize_path
from chemex.typing import Array

# Type aliases
Dataset = list[tuple[SpinSystem, Data]]
ProfilesType = dict[SpinSystem, list[Path]]
StructuredDType = list[tuple[str, str]]


@dataclass(eq=False)
class DatasetLoadError(Exception):
    """Invalid content in a user-supplied dataset file."""

    filename: Path
    explanation: str

    def __post_init__(self) -> None:
        super().__init__(f"Invalid data in '{self.filename}': {self.explanation}")


def _load_text_data(
    filename: Path,
    *,
    dtype: StructuredDType,
    usecols: list[int] | None = None,
) -> Array:
    try:
        return np.loadtxt(filename, dtype=dtype, usecols=usecols)
    except ValueError as error:
        explanation = str(error).splitlines()[0]
        raise DatasetLoadError(filename, explanation) from error


class HasRelaxationData(Protocol):
    @property
    def data(self) -> RelaxationDataSettings: ...


class HasShiftData(Protocol):
    @property
    def data(self) -> ShiftDataSettings: ...


def load_relaxation_dataset(base_path: Path, settings: HasRelaxationData) -> Dataset:
    data_path = normalize_path(base_path, settings.data.path)
    dtype = [("metadata", "f8"), ("exp", "f8"), ("err", "f8")]

    dataset: Dataset = []
    for spin_system, filepaths in settings.data.profiles.items():
        for filepath in filepaths:
            raw_data = _load_text_data(
                data_path / filepath,
                dtype=dtype,
                usecols=[0, 1, 2],
            )
            dataset.append(
                (
                    spin_system,
                    Data(
                        exp=raw_data["exp"],
                        err=raw_data["err"],
                        metadata=raw_data["metadata"],
                    ),
                ),
            )

    return dataset


def load_exsy_dataset(base_path: Path, settings: HasRelaxationData) -> Dataset:
    data_path = normalize_path(base_path, settings.data.path)
    dtype = [
        ("times", "f8"),
        ("states1", "U1"),
        ("states2", "U1"),
        ("exp", "f8"),
        ("err", "f8"),
    ]

    dataset: Dataset = []
    for spin_system, filepaths in settings.data.profiles.items():
        for filepath in filepaths:
            raw_data = _load_text_data(data_path / filepath, dtype=dtype)
            dataset.append(
                (
                    spin_system,
                    Data(
                        exp=raw_data["exp"],
                        err=raw_data["err"],
                        metadata=raw_data[["times", "states1", "states2"]],
                    ),
                ),
            )

    return dataset


def load_shift_dataset(base_path: Path, settings: HasShiftData) -> Dataset:
    data_path = normalize_path(base_path, settings.data.path)

    shifts = _load_text_data(
        data_path,
        dtype=[("spin_system", "U15"), ("exp", "f8"), ("err", "f8")],
    )

    return [
        (SpinSystem(name=spin_system), Data(exp=np.array([exp]), err=np.array([err])))
        for spin_system, exp, err in shifts
    ]

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from chemex.configuration.base import ExperimentConfiguration
from chemex.configuration.data import RelaxationDataSettings, ShiftDataSettings
from chemex.containers.data import Data
from chemex.parameters.spin_system import SpinSystem
from chemex.toml import normalize_path

# Type aliases
Dataset = list[tuple[SpinSystem, Data]]
ProfilesType = dict[SpinSystem, list[Path]]

RelaxationConfig = ExperimentConfiguration[Any, Any, RelaxationDataSettings]
ShiftConfig = ExperimentConfiguration[Any, Any, ShiftDataSettings]


def load_relaxation_dataset(base_path: Path, settings: RelaxationConfig) -> Dataset:
    data_path = normalize_path(base_path, settings.data.path)
    dtype = [("metadata", "f8"), ("exp", "f8"), ("err", "f8")]

    dataset: Dataset = []
    for spin_system, filepaths in settings.data.profiles.items():
        for filepath in filepaths:
            raw_data = np.loadtxt(data_path / filepath, dtype=dtype, usecols=[0, 1, 2])
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


def load_exsy_dataset(base_path: Path, settings: RelaxationConfig) -> Dataset:
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
            raw_data = np.loadtxt(data_path / filepath, dtype=dtype)
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


def load_shift_dataset(base_path: Path, settings: ShiftConfig) -> Dataset:
    data_path = normalize_path(base_path, settings.data.path)

    shifts = np.loadtxt(
        data_path,
        dtype=[("spin_system", "U15"), ("exp", "f8"), ("err", "f8")],
    )

    return [
        (SpinSystem(spin_system), Data(exp=np.array([exp]), err=np.array([err])))
        for spin_system, exp, err in shifts
    ]

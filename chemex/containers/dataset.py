from pathlib import Path

import numpy as np

from chemex.configuration.experiment import ExperimentConfig
from chemex.containers.data import Data
from chemex.parameters.spin_system import SpinSystem
from chemex.toml import normalize_path

# Type aliases
Dataset = list[tuple[SpinSystem, Data]]


def load_relaxation_dataset(base_path: Path, settings: ExperimentConfig) -> Dataset:

    data_path = normalize_path(base_path, settings.data.path)
    dtype = [("metadata", "f8"), ("exp", "f8"), ("err", "f8")]

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
                        metadata=raw_data["metadata"],
                    ),
                )
            )

    return dataset


def load_shift_dataset(base_path: Path, settings: ExperimentConfig) -> Dataset:

    data_path = normalize_path(base_path, settings.data.path)

    shifts = np.loadtxt(
        data_path,
        dtype=[("spin_system", "U15"), ("exp", "f8"), ("err", "f8")],
    )

    return [
        (SpinSystem(spin_system), Data(exp=np.array([exp]), err=np.array([err])))
        for spin_system, exp, err in shifts
    ]

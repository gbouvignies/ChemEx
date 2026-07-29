from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemex.configuration.data import RelaxationDataSettings
from chemex.containers.dataset import DatasetLoadError, load_relaxation_dataset
from chemex.parameters.spin_system import SpinSystem


def _settings(filename: Path) -> SimpleNamespace:
    return SimpleNamespace(
        data=RelaxationDataSettings(
            path=Path(),
            profiles={SpinSystem.from_name("G2N"): [filename]},
        ),
    )


def test_relaxation_dataset_loads_valid_rows(tmp_path: Path) -> None:
    data_file = tmp_path / "profile.dat"
    data_file.write_text("0.0 1.0 0.1\n0.1 0.5 0.1\n", encoding="utf-8")

    dataset = load_relaxation_dataset(tmp_path, _settings(data_file.name))

    assert len(dataset) == 1
    spin_system, data = dataset[0]
    assert spin_system == SpinSystem.from_name("G2N")
    np.testing.assert_allclose(data.metadata, [0.0, 0.1])
    np.testing.assert_allclose(data.exp, [1.0, 0.5])
    np.testing.assert_allclose(data.err, [0.1, 0.1])


@pytest.mark.parametrize(
    ("contents", "explanation"),
    [
        ("not-a-number 1.0 0.1\n", "could not convert"),
        ("0.0 1.0\n", "columns"),
    ],
)
def test_relaxation_dataset_reports_invalid_user_content(
    tmp_path: Path,
    contents: str,
    explanation: str,
) -> None:
    data_file = tmp_path / "profile.dat"
    data_file.write_text(contents, encoding="utf-8")

    with pytest.raises(DatasetLoadError) as error_info:
        load_relaxation_dataset(tmp_path, _settings(data_file.name))

    error = error_info.value
    assert error.filename == data_file
    assert explanation in error.explanation.lower()
    assert isinstance(error.__cause__, ValueError)

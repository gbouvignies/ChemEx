from __future__ import annotations

import pytest
from pydantic import ValidationError

from chemex.experiments.catalog.cpmg_ch3_1h_tq_diff import (
    CpmgCh31HTqDiffSettings,
)


def _settings_data(gradient: float) -> dict[str, object]:
    return {
        "name": "cpmg_ch3_1h_tq_diff",
        "time_t2": 40.0e-3,
        "carrier": 600.25,
        "pw90": 10.5e-6,
        "delta": 1.0e-3,
        "gradient": gradient,
    }


@pytest.mark.parametrize("gradient", [0.0, 30.6e-2])
def test_gradient_accepts_zero_and_positive_values(gradient: float) -> None:
    settings = CpmgCh31HTqDiffSettings.model_validate(_settings_data(gradient))

    assert settings.gradient == gradient


def test_gradient_rejects_negative_values() -> None:
    with pytest.raises(ValidationError):
        CpmgCh31HTqDiffSettings.model_validate(_settings_data(-1.0))

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest

from chemex.containers.data import Data
from chemex.containers.experiment import Experiment, NoDuplicateNoiseNotice
from chemex.containers.profile import Profile
from chemex.uncertainty import estimate_noise_variance


@dataclass
class _NoiseProfile:
    data: Data

    def any_duplicate(self) -> bool:
        return self.data.any_duplicate()

    def set_noise(self, value: float) -> None:
        self.data.err[:] = value
        self.data.mark_dirty()


def _profile(
    experimental: list[float],
    errors: list[float],
    metadata: list[float],
) -> _NoiseProfile:
    return _NoiseProfile(
        Data(
            exp=np.asarray(experimental),
            err=np.asarray(errors),
            metadata=np.asarray(metadata),
        )
    )


def _experiment(*profiles: _NoiseProfile) -> Experiment:
    return Experiment(
        Path("experiment.toml"),
        "test",
        cast("list[Profile]", list(profiles)),
        cast("Any", object()),
        cast("Any", object()),
    )


@pytest.mark.parametrize(
    ("errors", "expected_variance"),
    (([0.04, 0.04], 0.0016), ([0.03, 0.05], 0.0017)),
)
def test_duplicate_variance_fallback_uses_mean_file_variance(
    errors: list[float],
    expected_variance: float,
) -> None:
    profile = _profile([1.0, 2.0], errors, [1.0, 2.0])

    variance = estimate_noise_variance["duplicates"](profile.data)

    assert variance == pytest.approx(expected_variance, rel=1.0e-12, abs=0.0)


@pytest.mark.parametrize("global_error", (False, True))
def test_all_no_duplicate_profiles_preserve_pointwise_file_uncertainties(
    global_error: bool,
) -> None:
    first = _profile([1.0, 2.0], [0.03, 0.05], [1.0, 2.0])
    second = _profile([3.0, 4.0], [0.07, 0.09], [3.0, 4.0])
    experiment = _experiment(first, second)
    before = tuple(np.array(profile.data.err, copy=True) for profile in (first, second))

    experiment.estimate_noise("duplicates", global_error=global_error)

    for profile, expected in zip((first, second), before, strict=True):
        np.testing.assert_array_equal(profile.data.err, expected)
    assert experiment.noise_notices == (NoDuplicateNoiseNotice(experiment.filename),)


@pytest.mark.parametrize(
    ("file_errors", "expected_file_uncertainty"),
    (([0.04, 0.04], 0.04), ([0.03, 0.05], np.sqrt(0.0017))),
)
def test_mixed_duplicate_noise_is_profile_local_in_local_mode(
    file_errors: list[float],
    expected_file_uncertainty: float,
) -> None:
    duplicate = _profile([1.0, 1.2, 3.0], [0.5, 0.5, 0.5], [1.0, 1.0, 2.0])
    no_duplicate = _profile([4.0, 5.0], file_errors, [3.0, 4.0])

    _experiment(duplicate, no_duplicate).estimate_noise(
        "duplicates",
        global_error=False,
    )

    np.testing.assert_allclose(
        duplicate.data.err,
        np.sqrt(0.02),
        rtol=1.0e-12,
        atol=0.0,
    )
    np.testing.assert_allclose(
        no_duplicate.data.err,
        expected_file_uncertainty,
        rtol=1.0e-12,
        atol=0.0,
    )


def test_mixed_duplicate_noise_contributes_file_variance_to_global_pool() -> None:
    duplicate = _profile([1.0, 1.2, 3.0], [0.5, 0.5, 0.5], [1.0, 1.0, 2.0])
    no_duplicate = _profile([4.0, 5.0], [0.04, 0.04], [3.0, 4.0])

    _experiment(duplicate, no_duplicate).estimate_noise(
        "duplicates",
        global_error=True,
    )

    expected = np.sqrt((0.02 + 0.0016) / 2.0)
    np.testing.assert_allclose(
        duplicate.data.err,
        expected,
        rtol=1.0e-12,
        atol=0.0,
    )
    np.testing.assert_allclose(
        no_duplicate.data.err,
        expected,
        rtol=1.0e-12,
        atol=0.0,
    )


@pytest.mark.parametrize(
    ("global_error", "expected"),
    (
        (False, (np.sqrt(0.02), np.sqrt(0.08))),
        (True, (np.sqrt(0.05), np.sqrt(0.05))),
    ),
)
def test_every_profile_with_duplicates_preserves_pooled_variance_behavior(
    global_error: bool,
    expected: tuple[float, float],
) -> None:
    first = _profile([1.0, 1.2], [0.5, 0.5], [1.0, 1.0])
    second = _profile([2.0, 2.4], [0.5, 0.5], [2.0, 2.0])

    _experiment(first, second).estimate_noise(
        "duplicates",
        global_error=global_error,
    )

    np.testing.assert_allclose(
        first.data.err,
        expected[0],
        rtol=1.0e-12,
        atol=0.0,
    )
    np.testing.assert_allclose(
        second.data.err,
        expected[1],
        rtol=1.0e-12,
        atol=0.0,
    )

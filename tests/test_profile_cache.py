from __future__ import annotations

import numpy as np
import pytest
from lmfit import Parameters

import chemex.containers.data as data_module
from chemex.containers.data import Data
from chemex.containers.profile import Profile


class DummyLiouvillian:
    spin_system = "dummy"


class DummySpectrometer:
    def __init__(self) -> None:
        self.liouvillian = DummyLiouvillian()
        self.par_values: dict[str, float] = {}

    def update(self, par_values: dict[str, float]) -> None:
        self.par_values = par_values


class DummyPulseSequence:
    def __init__(self, calc_values: list[float]) -> None:
        self.calc_values = np.array(calc_values, dtype=np.float64)

    def calculate(self, _spectrometer: DummySpectrometer, _data: Data) -> np.ndarray:
        return self.calc_values.copy()

    def is_reference(self, metadata: np.ndarray) -> np.ndarray:
        return metadata < 0


class DummyPrinter:
    header = ""
    simulation = False

    def print(self, _name: str, _data: Data) -> str:
        return ""


class MaskFirstFilterer:
    def filter(self, data: Data) -> None:
        data.mask[0] = False


class StubRng:
    def __init__(self, values: list[float]) -> None:
        self.values = np.array(values, dtype=np.float64)

    def normal(self, _calc: np.ndarray, _err: np.ndarray) -> np.ndarray:
        return self.values.copy()


def make_profile(
    *,
    exp: list[float],
    err: list[float],
    metadata: list[float],
    calc: list[float],
    filterer: MaskFirstFilterer | None = None,
) -> Profile:
    data = Data(
        exp=np.array(exp, dtype=np.float64),
        err=np.array(err, dtype=np.float64),
        metadata=np.array(metadata, dtype=np.float64),
    )
    return Profile(
        data=data,
        spectrometer=DummySpectrometer(),
        pulse_sequence=DummyPulseSequence(calc),
        name_map={"k": "p1"},
        printer=DummyPrinter(),
        filterer=filterer,
        is_scaled=False,
    )


def make_params() -> Parameters:
    params = Parameters()
    params.add("p1", value=1.0)
    return params


def test_set_noise_invalidates_residual_cache() -> None:
    profile = make_profile(
        exp=[2.0, 4.0],
        err=[1.0, 1.0],
        metadata=[0.0, 1.0],
        calc=[1.0, 2.0],
    )
    params = make_params()

    np.testing.assert_allclose(profile.residuals(params), [-1.0, -2.0])
    np.testing.assert_equal(profile.cache.currsize, 1)

    profile.set_noise(0.5)

    np.testing.assert_allclose(profile.residuals(params), [-2.0, -4.0])


def test_filter_invalidates_residual_cache() -> None:
    profile = make_profile(
        exp=[2.0, 4.0],
        err=[1.0, 1.0],
        metadata=[0.0, 1.0],
        calc=[1.0, 2.0],
        filterer=MaskFirstFilterer(),
    )
    params = make_params()

    np.testing.assert_allclose(profile.residuals(params), [-1.0, -2.0])
    np.testing.assert_equal(profile.cache.currsize, 1)

    profile.filter(params)

    np.testing.assert_allclose(profile.residuals(params), [-2.0])


def test_monte_carlo_profile_starts_with_empty_cache_and_fresh_residuals() -> None:
    original_rng = data_module.rng
    data_module.rng = StubRng([10.0, 20.0, 30.0])
    try:
        profile = make_profile(
            exp=[1.0, 2.0, 3.0],
            err=[1.0, 1.0, 1.0],
            metadata=[0.0, 1.0, 2.0],
            calc=[1.0, 2.0, 3.0],
        )
        params = make_params()

        np.testing.assert_allclose(profile.residuals(params), [0.0, 0.0, 0.0])
        np.testing.assert_equal(profile.cache.currsize, 1)

        profile_mc = profile.monte_carlo()

        np.testing.assert_equal(profile_mc.cache.currsize, 0)
        expected = (
            np.array([1.0, 2.0, 3.0]) - profile_mc.data.exp
        ) / profile_mc.data.err
        np.testing.assert_allclose(profile_mc.residuals(params), expected)
    finally:
        data_module.rng = original_rng


def test_bootstrap_profile_starts_with_empty_cache_and_fresh_residuals(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def duplicate_choices(pool: np.ndarray, k: int) -> list[int]:
        return [int(pool[-1])] * k

    monkeypatch.setattr(data_module, "choices", duplicate_choices)

    profile = make_profile(
        exp=[11.0, 12.0, 21.0, 22.0],
        err=[1.0, 2.0, 3.0, 4.0],
        metadata=[-2.0, -1.0, 1.0, 2.0],
        calc=[1.0, 2.0, 3.0, 4.0],
    )
    params = make_params()

    np.testing.assert_allclose(
        profile.residuals(params),
        [-10.0, -5.0, -6.0, -4.5],
    )
    np.testing.assert_equal(profile.cache.currsize, 1)

    profile_bs = profile.bootstrap()

    np.testing.assert_equal(profile_bs.cache.currsize, 0)
    expected = (
        np.array([1.0, 2.0, 3.0, 4.0]) - profile_bs.data.exp
    ) / profile_bs.data.err
    np.testing.assert_allclose(profile_bs.residuals(params), expected)

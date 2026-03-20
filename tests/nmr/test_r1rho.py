from __future__ import annotations

import numpy as np
import pytest

import chemex.nmr.spectrometer_analysis as analysis_module
from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.constants import Distribution
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem


def make_spectrometer() -> Spectrometer:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    spectrometer = Spectrometer.from_spin_system(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    spectrometer.update({"r2_i_a": 4.0})
    spectrometer.b1_i = 2_000.0
    return spectrometer


def test_calculate_r1rho_returns_rate_for_scalar_liouvillian() -> None:
    spectrometer = make_spectrometer()

    r1rho = spectrometer.analysis.calculate_r1rho()

    assert np.isfinite(r1rho)
    assert r1rho >= 0.0


def test_calculate_r1rho_rejects_multi_point_distributions() -> None:
    spectrometer = make_spectrometer()
    spectrometer.set_b1_i_distribution(
        Distribution(
            np.array([1_800.0, 2_000.0, 2_200.0]),
            np.array([0.2, 0.6, 0.2]),
        ),
        nominal=2_000.0,
    )

    with pytest.raises(
        ValueError,
        match="requires single-point B1 and Jeff distributions",
    ):
        spectrometer.analysis.calculate_r1rho()


def test_calculate_r1rho_raises_when_no_nearly_real_eigenvalue_exists(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    spectrometer = make_spectrometer()

    monkeypatch.setattr(
        analysis_module.np.linalg,
        "eigvals",
        lambda _matrix: np.array([1.0 + 1e-3j, -2.0 - 2e-3j]),
    )

    with pytest.raises(ValueError, match="did not find a nearly real eigenvalue"):
        spectrometer.analysis.calculate_r1rho()

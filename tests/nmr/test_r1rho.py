from __future__ import annotations

import numpy as np
import pytest

import chemex.nmr.liouvillian as liouvillian_module
from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.constants import Distribution
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.parameters.spin_system import SpinSystem


def make_liouvillian() -> LiouvillianIS:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    liouvillian = LiouvillianIS(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    liouvillian.update({"r2_i_a": 4.0})
    liouvillian.b1_i = 2_000.0
    return liouvillian


def test_calculate_r1rho_returns_rate_for_scalar_liouvillian() -> None:
    liouvillian = make_liouvillian()

    r1rho = liouvillian.calculate_r1rho()

    assert np.isfinite(r1rho)
    assert r1rho >= 0.0


def test_calculate_r1rho_rejects_multi_point_distributions() -> None:
    liouvillian = make_liouvillian()
    liouvillian.set_b1_i_distribution(
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
        liouvillian.calculate_r1rho()


def test_calculate_r1rho_raises_when_no_nearly_real_eigenvalue_exists(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    liouvillian = make_liouvillian()

    monkeypatch.setattr(
        liouvillian_module.np.linalg,
        "eigvals",
        lambda _matrix: np.array([1.0 + 1e-3j, -2.0 - 2e-3j]),
    )

    with pytest.raises(ValueError, match="did not find a nearly real eigenvalue"):
        liouvillian.calculate_r1rho()

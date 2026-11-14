from __future__ import annotations

import numpy as np
import pytest

import chemex.nmr._engine.analysis as analysis_module
from chemex.configuration.conditions import Conditions
from chemex.models.loader import register_kinetic_settings
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


@pytest.mark.parametrize(
    ("hom_type", "eq_type"),
    [
        ("iz", "iz_eq"),
        ("ixyz", "ixyz_eq"),
        ("ixyzsz", "ixyzsz_eq"),
        ("ixyzsxyz", "ixyzsxyz_eq"),
    ],
)
@pytest.mark.parametrize("model_name", ["1st", "2st"])
def test_calculate_r1rho_on_equilibrium_bases_matches_homogeneous_bases(
    hom_type: str, eq_type: str, model_name: str
) -> None:
    register_kinetic_settings()
    model = ModelSpec.from_name(model_name)
    conditions = Conditions(h_larmor_frq=600.0)
    spin_system = SpinSystem(name="G23N-HN")

    spec_hom = Spectrometer.from_spin_system(
        spin_system,
        Basis(type=hom_type, spin_system="nh", model=model),
        conditions,
    )
    spec_eq = Spectrometer.from_spin_system(
        spin_system,
        Basis(type=eq_type, spin_system="nh", model=model),
        conditions,
    )
    for spec in (spec_hom, spec_eq):
        spec.carrier_i = 120.0
        spec.b1_i = 2_000.0 if "x" in hom_type else 0.0

    params: dict[str, float] = {
        "pa": 0.9 if model_name == "2st" else 1.0,
        "cs_i_a": 122.0,
        "cs_s_a": 8.5,
        "j_is_a": -92.0,
        "r1_i_a": 1.5,
        "r2_i_a": 15.0,
        "r1_s_a": 2.5,
        "r2_s_a": 20.0,
        "r2a_i_a": 16.5,
        "r2a_s_a": 21.5,
        "r1a_is_a": 3.0,
        "r2mq_is_a": 25.0,
        "sigma_is_a": -0.2,
        "etaz_i_a": 0.1,
        "etaz_s_a": 0.1,
        "etaxy_i_a": 0.5,
        "etaxy_s_a": 0.5,
        "mu_is_a": 0.0,
    }
    if model_name == "2st":
        # Choose slow exchange kex = 0.5 s^-1 < r1rho so any unprojected ie_a <-> ie_b
        # bath exchange eigenvalue (-0.5) would also corrupt max(real_eigenvalues).
        params.update(
            {
                "pb": 0.1,
                "kab": 0.05,
                "kba": 0.45,
                "cs_i_b": 125.0,
                "cs_s_b": 8.8,
                "j_is_b": -92.0,
                "r1_i_b": 1.5,
                "r2_i_b": 15.0,
                "r1_s_b": 2.5,
                "r2_s_b": 20.0,
                "r2a_i_b": 16.5,
                "r2a_s_b": 21.5,
                "r1a_is_b": 3.0,
                "r2mq_is_b": 25.0,
                "sigma_is_b": -0.2,
                "etaz_i_b": 0.1,
                "etaz_s_b": 0.1,
                "etaxy_i_b": 0.5,
                "etaxy_s_b": 0.5,
                "mu_is_b": 0.0,
            }
        )

    spec_hom.update(params)
    spec_eq.update(params)

    r1rho_hom = spec_hom.analysis.calculate_r1rho()
    r1rho_eq = spec_eq.analysis.calculate_r1rho()

    assert r1rho_hom > 1.0
    assert r1rho_eq == pytest.approx(r1rho_hom, rel=1e-14, abs=1e-15)


def test_calculate_r1rho_on_ixyz_eq_preserves_zero_relaxation_limit() -> None:
    spectrometer = Spectrometer.from_spin_system(
        SpinSystem(name="G23N-HN"),
        Basis(type="ixyz_eq", spin_system="nh", model=ModelSpec()),
        Conditions(h_larmor_frq=600.0),
    )
    spectrometer.carrier_i = 120.0
    spectrometer.b1_i = 2_000.0
    spectrometer.update({"pa": 1.0, "cs_i_a": 120.0, "r1_i_a": 0.0, "r2_i_a": 0.0})

    assert spectrometer.analysis.calculate_r1rho() == pytest.approx(0.0, abs=1e-14)

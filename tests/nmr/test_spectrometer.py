from __future__ import annotations

import numpy as np

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
    spectrometer.b1_i = 2_500.0
    spectrometer.b1_s = 20_000.0
    return spectrometer


def test_i_pulse_cache_reuses_base_pulses_until_invalidated() -> None:
    spectrometer = make_spectrometer()

    p90_i = spectrometer.p90_i

    assert spectrometer.p90_i is p90_i

    spectrometer.b1_i = 3_000.0

    assert spectrometer.p90_i is not p90_i


def test_b1_i_change_does_not_invalidate_s_pulse_cache() -> None:
    spectrometer = make_spectrometer()

    p180_s = spectrometer.p180_s
    spectrometer.b1_i = 3_000.0

    assert spectrometer.p180_s is p180_s


def test_b1_s_change_does_not_invalidate_i_pulse_cache() -> None:
    spectrometer = make_spectrometer()

    p90_i = spectrometer.p90_i
    spectrometer.b1_s = 25_000.0

    assert spectrometer.p90_i is p90_i


def test_liouvillian_changes_invalidate_both_base_pulse_caches() -> None:
    spectrometer = make_spectrometer()

    p90_i = spectrometer.p90_i
    p180_s = spectrometer.p180_s
    spectrometer.carrier_i = 118.0

    assert spectrometer.p90_i is not p90_i
    assert spectrometer.p180_s is not p180_s


def test_s_pulse_cache_handles_distributed_liouvillian() -> None:
    spectrometer = make_spectrometer()
    spectrometer.jeff_i = Distribution(
        np.array([0.0, 10.0]),
        np.array([0.5, 0.5]),
    )

    p180_s = spectrometer.p180_s

    assert p180_s.shape == (4, 2, 1, 6, 6)
    assert np.isfinite(p180_s).all()


def test_spectrometer_exposes_basis_spin_system_and_ppm_metadata() -> None:
    spectrometer = make_spectrometer()

    assert spectrometer.spin_system.name == "G23N"
    assert spectrometer.basis.spin_system == "nh"
    assert spectrometer.ppm_i != 0.0
    assert spectrometer.ppm_s != 0.0


def test_calculate_r1rho_is_available_on_spectrometer() -> None:
    spectrometer = make_spectrometer()
    spectrometer.update({"r2_i_a": 4.0})

    r1rho = spectrometer.analysis.calculate_r1rho()

    assert np.isfinite(r1rho)
    assert r1rho >= 0.0

from __future__ import annotations

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem


def make_spectrometer() -> Spectrometer:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    liouvillian = LiouvillianIS(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    spectrometer = Spectrometer(liouvillian)
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

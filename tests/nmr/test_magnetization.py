from __future__ import annotations

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr._engine.magnetization import (
    build_equilibrium_magnetization,
    build_start_magnetization,
    collapse_magnetization,
    detect_signal,
    keep_components,
)
from chemex.nmr.basis import Basis
from chemex.nmr.constants import GAMMA, SIGNED_XI_RATIO, XI_RATIO
from chemex.nmr.rates import rate_functions
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spin_system.nucleus import Nucleus


def test_collapse_magnetization_applies_distribution_weights() -> None:
    magnetization = np.array(
        [
            [[1.0], [2.0]],
            [[10.0], [20.0]],
        ]
    )
    weights = np.array(
        [
            [[0.25]],
            [[0.75]],
        ]
    )

    collapsed = collapse_magnetization(magnetization, weights)

    np.testing.assert_allclose(collapsed, [[7.75], [15.5]])


def test_detect_signal_uses_collapsed_magnetization() -> None:
    detection_vector = np.array([[1.0, -1.0]])
    magnetization = np.array(
        [
            [[2.0], [1.0]],
            [[5.0], [2.0]],
        ]
    )
    weights = np.array(
        [
            [[0.25]],
            [[0.75]],
        ]
    )

    detected = detect_signal(detection_vector, magnetization, weights)

    assert detected == 2.5


def test_build_equilibrium_magnetization_uses_populations_and_signed_xi_ratio() -> None:
    basis = Basis(type="ixyzsxyz_eq", spin_system="nh", model=ModelSpec())
    par_values = {"pa": 0.6, "pb": 0.4}

    magnetization = build_equilibrium_magnetization(basis, par_values)
    expected_n_a = 0.6 * SIGNED_XI_RATIO[Nucleus.N15]
    expected_n_b = 0.4 * SIGNED_XI_RATIO[Nucleus.N15]
    expected_h_a = 0.6 * SIGNED_XI_RATIO[Nucleus.H1]
    expected_h_b = 0.4 * SIGNED_XI_RATIO[Nucleus.H1]
    expected = (
        expected_n_a * (basis.vectors["ie_a"] + basis.vectors["iz_a"])
        + expected_n_b * (basis.vectors["ie_b"] + basis.vectors["iz_b"])
        + expected_h_a * (basis.vectors["se_a"] + basis.vectors["sz_a"])
        + expected_h_b * (basis.vectors["se_b"] + basis.vectors["sz_b"])
    )

    assert expected_n_a < 0.0
    assert expected_h_a > 0.0
    np.testing.assert_allclose(magnetization, expected)


def test_build_start_magnetization_keeps_unsigned_xi_ratio_for_negative_gamma() -> None:
    basis = Basis(type="iz", spin_system="nh", model=ModelSpec())
    par_values = {"pa": 0.6, "pb": 0.4}

    magnetization = build_start_magnetization(
        basis,
        par_values,
        terms=["iz"],
        atom=Nucleus.N15,
    )
    expected = XI_RATIO[Nucleus.N15] * (
        0.6 * basis.vectors["iz_a"] + 0.4 * basis.vectors["iz_b"]
    )

    assert float(magnetization.sum()) > 0.0
    np.testing.assert_allclose(magnetization, expected)


@pytest.mark.parametrize(
    ("spin_system", "spin_name", "nucleus_i", "nucleus_s"),
    (
        ("nh", "G23N-HN", Nucleus.N15, Nucleus.H1),
        ("ch", "A10CA-HA", Nucleus.C13, Nucleus.H1),
    ),
)
def test_two_spin_equilibrium_steady_state_noe_matches_solomon_equation(
    spin_system: str,
    spin_name: str,
    nucleus_i: Nucleus,
    nucleus_s: Nucleus,
) -> None:
    model = ModelSpec(name="1st", states="a")
    basis = Basis(type="ixyzsxyz_eq", spin_system=spin_system, model=model)
    conditions = Conditions(h_larmor_frq=600.0)
    spectrometer = Spectrometer.from_spin_system(
        SpinSystem(name=spin_name),
        basis,
        conditions,
    )
    spectrometer.detection = "[iz_a]"

    # 1. Analytical pure-dipolar extreme-narrowing limit (sigma_is = 0.5 * r1_i)
    r1_i = 1.5
    sigma_is = 0.5 * r1_i
    par_values = dict.fromkeys(basis.matrices, 0.0)
    par_values.update(
        {
            "pa": 1.0,
            "r1_i_a": r1_i,
            "r2_i_a": r1_i,
            "r1_s_a": 3.0,
            "r2_s_a": 3.0,
            "r2a_i_a": 2.5,
            "r2a_s_a": 2.5,
            "r2mq_is_a": 3.0,
            "r1a_is_a": 2.5,
            "sigma_is_a": sigma_is,
        }
    )
    spectrometer.update(par_values)
    mag_eq = spectrometer.get_equilibrium()
    i_eq = spectrometer.detect(mag_eq)

    spectrometer.b1_s = 5.0e4
    mag_sat = spectrometer.pulse_s(50.0 / r1_i, phase=0.0) @ mag_eq
    i_sat = spectrometer.detect(mag_sat)

    expected_noe = 1.0 + (SIGNED_XI_RATIO[nucleus_s] / SIGNED_XI_RATIO[nucleus_i]) * (
        sigma_is / r1_i
    )
    assert i_sat / i_eq == pytest.approx(expected_noe, rel=1e-6)
    assert i_sat / i_eq == pytest.approx(
        1.0 + (GAMMA[nucleus_s] / GAMMA[nucleus_i]) * (sigma_is / r1_i),
        rel=1e-3,
    )
    if nucleus_i is Nucleus.N15:
        assert i_sat / i_eq < 0.0
    else:
        assert i_sat / i_eq > 1.0

    # 2. Realistic protein Redfield rates at 600 MHz, tauc = 8.0 ns, S2 = 0.85
    redfield = rate_functions[spin_system](600.0, 8.0, 0.85)
    par_redfield = {f"{key}_a": float(val) for key, val in redfield.items()}
    par_redfield["pa"] = 1.0
    par_redfield["cs_i_a"] = 0.0
    par_redfield["cs_s_a"] = 0.0
    par_redfield["j_is_a"] = -93.0 if spin_system == "nh" else 140.0
    spectrometer.b1_s = 0.0
    spectrometer.update(par_redfield)
    mag_eq_rf = spectrometer.get_equilibrium()
    i_eq_rf = spectrometer.detect(mag_eq_rf)
    # Free relaxation preserves thermal equilibrium strictly
    np.testing.assert_allclose(
        (spectrometer.delays(1.0) @ mag_eq_rf).reshape(mag_eq_rf.shape),
        mag_eq_rf,
        rtol=1e-12,
        atol=1e-14,
    )
    spectrometer.b1_s = 5.0e4
    mag_sat_rf = (
        spectrometer.pulse_s(50.0 / par_redfield["r1_i_a"], phase=0.0) @ mag_eq_rf
    )
    i_sat_rf = spectrometer.detect(mag_sat_rf)
    expected_rf_noe = 1.0 + (
        SIGNED_XI_RATIO[nucleus_s] / SIGNED_XI_RATIO[nucleus_i]
    ) * (par_redfield["sigma_is_a"] / par_redfield["r1_i_a"])
    assert i_sat_rf / i_eq_rf == pytest.approx(expected_rf_noe, rel=1e-5)
    assert i_sat_rf / i_eq_rf == pytest.approx(
        1.0
        + (GAMMA[nucleus_s] / GAMMA[nucleus_i])
        * (par_redfield["sigma_is_a"] / par_redfield["r1_i_a"]),
        rel=1e-3,
    )
    if nucleus_i is Nucleus.N15:
        assert 0.0 < i_sat_rf / i_eq_rf < 1.0


def test_build_start_magnetization_respects_terms_signs_and_populations() -> None:
    basis = Basis(type="iz", spin_system="nh", model=ModelSpec())
    par_values = {"pa": 0.75, "pb": 0.25}

    magnetization = build_start_magnetization(
        basis,
        par_values,
        terms=["iz_a", "-iz_b"],
        atom=Nucleus.H1,
    )
    expected = 0.75 * basis.vectors["iz_a"] - 0.25 * basis.vectors["iz_b"]

    np.testing.assert_allclose(magnetization, expected)


def test_explicit_all_state_start_matches_equilibrium_component_selection() -> None:
    model = ModelSpec(name="3st", states="abc")
    basis = Basis(type="iz", spin_system="nh", model=model)
    par_values = {"pa": 0.6, "pb": 0.3, "pc": 0.1}

    equilibrium_terms = build_start_magnetization(
        basis,
        par_values,
        terms=["iz"],
    )
    explicit_terms = build_start_magnetization(
        basis,
        par_values,
        terms=["iz_a", "iz_b", "iz_c"],
    )

    np.testing.assert_allclose(explicit_terms, equilibrium_terms)


def test_keep_components_masks_unselected_magnetization() -> None:
    basis = Basis(type="izsz", spin_system="nh", model=ModelSpec())
    magnetization = basis.vectors["iz_a"] + 2.0 * basis.vectors["2izsz_b"]

    kept = keep_components(basis, magnetization, components=["iz_a"])

    np.testing.assert_allclose(kept, basis.vectors["iz_a"])

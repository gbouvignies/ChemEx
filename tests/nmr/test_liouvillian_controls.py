from __future__ import annotations

import numpy as np

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr._engine.engine import ISLiouvillianEngine
from chemex.nmr.basis import Basis
from chemex.parameters.spin_system import SpinSystem


def make_liouvillian() -> ISLiouvillianEngine:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    return ISLiouvillianEngine(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )


def test_scalar_control_values_are_stored_as_plain_floats() -> None:
    liouvillian = make_liouvillian()
    liouvillian.carrier_i = 118.0
    liouvillian.carrier_s = -25.0
    liouvillian.offset_i = 150.0
    liouvillian.offset_s = -35.0

    assert isinstance(liouvillian.state.carrier_i, float)
    assert isinstance(liouvillian.state.carrier_s, float)
    assert isinstance(liouvillian.state.offset_i, float)
    assert isinstance(liouvillian.state.offset_s, float)


def test_update_without_parameters_keeps_matrix_shaped_base_liouvillian() -> None:
    liouvillian = make_liouvillian()

    liouvillian.update({})

    assert liouvillian._l_base.shape == (liouvillian.size, liouvillian.size)
    np.testing.assert_allclose(liouvillian._l_base, 0.0)


def test_missing_s_control_terms_use_zero_liouvillian_matrices() -> None:
    liouvillian = make_liouvillian()
    liouvillian.carrier_s = 25.0
    liouvillian.offset_s = -50.0
    liouvillian.b1_s = 20_000.0

    expected_shape = (liouvillian.size, liouvillian.size)
    assert liouvillian._l_carrier_s.shape == expected_shape
    assert liouvillian._l_offset_s.shape == expected_shape
    assert liouvillian.l_b1x_s.shape == expected_shape
    assert liouvillian.l_b1y_s.shape == expected_shape
    np.testing.assert_allclose(liouvillian._l_carrier_s, 0.0)
    np.testing.assert_allclose(liouvillian._l_offset_s, 0.0)
    np.testing.assert_allclose(liouvillian.l_b1x_s, 0.0)
    np.testing.assert_allclose(liouvillian.l_b1y_s, 0.0)


def test_offset_sign_matches_gyromagnetic_convention_without_h_larmor_frq() -> None:
    """offset_i and offset_s remain active with the physical sign when ppm_i == ppm_s == 0."""
    basis_nh = Basis(type="ixyzsxyz", spin_system="nh", model=ModelSpec())
    engine_no_b0 = ISLiouvillianEngine(
        SpinSystem(name="G23N-HN"),
        basis_nh,
        Conditions(),
    )
    engine_with_b0 = ISLiouvillianEngine(
        SpinSystem(name="G23N-HN"),
        basis_nh,
        Conditions(h_larmor_frq=600.0),
    )

    assert engine_no_b0.ppm_i == 0.0
    assert engine_no_b0.ppm_s == 0.0

    engine_no_b0.offset_i = 125.0
    engine_no_b0.offset_s = -80.0
    engine_with_b0.offset_i = 125.0
    engine_with_b0.offset_s = -80.0

    assert np.any(engine_no_b0._l_offset_i != 0.0)
    assert np.any(engine_no_b0._l_offset_s != 0.0)
    np.testing.assert_allclose(
        engine_no_b0._l_offset_i,
        engine_with_b0._l_offset_i,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        engine_no_b0._l_offset_s,
        engine_with_b0._l_offset_s,
        rtol=0.0,
        atol=0.0,
    )

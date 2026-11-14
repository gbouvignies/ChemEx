from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import expm

from chemex.configuration.conditions import Conditions
from chemex.containers.data import Data
from chemex.experiments.catalog.cpmg_ch3_1h_tq_diff import (
    CpmgCh31HTqDiffConfig,
    CpmgCh31HTqDiffSequence,
)
from chemex.experiments.catalog.cpmg_ch3_1h_tq_diff import (
    build_spectrometer as build_tq_diff_spectrometer,
)
from chemex.experiments.catalog.cpmg_ch3_mq import (
    CpmgCh3MqConfig,
    CpmgCh3MqSequence,
)
from chemex.experiments.catalog.cpmg_ch3_mq import (
    build_spectrometer as build_ch3_mq_spectrometer,
)
from chemex.experiments.catalog.relaxation_hznz import (
    RelaxationHzNzConfig,
    RelaxationHzNzSequence,
)
from chemex.experiments.catalog.relaxation_hznz import (
    build_spectrometer as build_hznz_spectrometer,
)
from chemex.experiments.catalog.relaxation_nz import (
    RelaxationNzConfig,
    RelaxationNzSequence,
)
from chemex.experiments.catalog.relaxation_nz import (
    build_spectrometer as build_nz_spectrometer,
)
from chemex.models.loader import register_kinetic_settings
from chemex.models.model import ModelSpec
from chemex.nmr._pulses.propagators import calculate_propagators
from chemex.nmr.basis import Basis
from chemex.nmr.constants import Distribution
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem


def test_calculate_propagators_single_delay_matches_expm() -> None:
    liouv = np.array([[0.0, 1.0], [-2.0, -3.0]])
    delay = 0.25

    propagator = calculate_propagators(liouv, delay)
    expected = expm(liouv * delay)

    np.testing.assert_allclose(propagator, expected)


def test_calculate_propagators_multiple_delays_match_expm() -> None:
    liouv = np.array([[0.0, 1.0], [-2.0, -3.0]])
    delays = np.array([0.1, 0.25, 0.5])

    propagators = calculate_propagators(liouv, delays)
    expected = np.array([expm(liouv * delay) for delay in delays])

    np.testing.assert_allclose(propagators, expected)


def test_calculate_propagators_handles_stacked_liouvillians() -> None:
    liouvillians = np.array(
        [
            [[0.0, 1.0], [-2.0, -3.0]],
            [[-1.0, 0.5], [-0.25, -2.0]],
        ]
    )
    delays = np.array([0.1, 0.25])

    propagators = calculate_propagators(liouvillians, delays)
    expected = np.array(
        [[expm(liouv * delay) for liouv in liouvillians] for delay in delays]
    )

    np.testing.assert_allclose(propagators, expected)


def test_calculate_propagators_handles_stacked_liouvillians_for_scalar_delay() -> None:
    liouvillians = np.array(
        [
            [[0.0, 1.0], [-2.0, -3.0]],
            [[-1.0, 0.5], [-0.25, -2.0]],
        ]
    )
    delay = 0.25

    propagators = calculate_propagators(liouvillians, delay)
    expected = np.array([expm(liouv * delay) for liouv in liouvillians])

    np.testing.assert_allclose(propagators, expected)


def test_complete_dephasing_removes_undamped_oscillatory_modes() -> None:
    liouvillians = np.array(
        [
            [[0.0, 1.0], [-1.0, 0.0]],
            [[0.0, 2.0], [-2.0, 0.0]],
        ]
    )
    delays = [0.25, 0.5]

    propagators = calculate_propagators(liouvillians, delays, dephasing=True)

    assert propagators.shape == (2, 2, 2, 2)
    assert propagators.dtype == np.float64
    np.testing.assert_array_equal(propagators, np.zeros_like(propagators))


@pytest.mark.parametrize("damping", [0.0, 1.0e-6, 3.0])
def test_complete_dephasing_agrees_with_broad_uniform_b1_limit(
    damping: float,
) -> None:
    delay = 0.25
    liouv = np.array([[-damping, 1.0], [-1.0, -damping]])
    # The finite uniform ensemble tends to zero as its B1 width increases.
    frequencies = np.linspace(-1000.0, 1000.0, 20001)
    phases = frequencies * delay
    ensemble = np.exp(-damping * delay) * np.array(
        [
            [np.cos(phases).mean(), np.sin(phases).mean()],
            [-np.sin(phases).mean(), np.cos(phases).mean()],
        ]
    )

    np.testing.assert_array_equal(
        calculate_propagators(liouv, delay, dephasing=True), np.zeros((2, 2))
    )
    assert np.linalg.norm(ensemble) < 0.01


def test_complete_dephasing_matches_damped_legacy_limit() -> None:
    liouv = np.array([[-3.0, 1.0, 0.0], [-1.0, -3.0, 0.0], [0.0, 0.0, -2.0]])
    expected = np.diag([0.0, 0.0, np.exp(-0.5)])

    np.testing.assert_allclose(
        calculate_propagators(liouv, 0.25, dephasing=True),
        expected,
        rtol=0.0,
        atol=2.0e-15,
    )


@pytest.mark.parametrize(("frequency", "retained"), [(0.5e-6, True), (2.0e-6, False)])
def test_complete_dephasing_uses_historical_strict_imaginary_cutoff(
    frequency: float, retained: bool
) -> None:
    liouv = np.array([[0.0, frequency], [-frequency, 0.0]])
    expected = expm(liouv * 0.25) if retained else np.zeros((2, 2))

    np.testing.assert_allclose(
        calculate_propagators(liouv, 0.25, dephasing=True),
        expected,
        rtol=0.0,
        atol=2.0e-15,
    )


def test_complete_dephasing_at_zero_delay_is_identity() -> None:
    liouv = np.array([[0.0, 1.0], [-1.0, 0.0]])

    np.testing.assert_allclose(
        calculate_propagators(liouv, [0.0, 0.25], dephasing=True),
        [np.eye(2), np.zeros((2, 2))],
        rtol=0.0,
        atol=2.0e-15,
    )


def test_complete_dephasing_handles_repeated_modes_and_non_normal_basis() -> None:
    base = np.array(
        [
            [0.0, 2.0, 0.0, 0.0, 0.0],
            [-2.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 2.0, 0.0],
            [0.0, 0.0, -2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, -1.0],
        ]
    )
    similarity = np.eye(5)
    similarity[0, 4] = 3.0
    liouv = similarity @ base @ np.linalg.inv(similarity)
    expected = (
        similarity
        @ np.diag([0.0, 0.0, 0.0, 0.0, np.exp(-0.25)])
        @ np.linalg.inv(similarity)
    )

    np.testing.assert_allclose(
        calculate_propagators(liouv, 0.25, dephasing=True),
        expected,
        rtol=0.0,
        atol=1.0e-14,
    )


def test_complete_dephasing_handles_nearly_defective_oscillatory_block() -> None:
    frequency = 2.0e-6
    liouv = np.array([[0.0, 1.0], [-(frequency**2), 0.0]])
    assert np.linalg.cond(np.linalg.eig(liouv).eigenvectors) > 1.0e5

    np.testing.assert_allclose(
        calculate_propagators(liouv, 0.25, dephasing=True),
        np.zeros((2, 2)),
        rtol=0.0,
        atol=1.0e-12,
    )


def test_calculate_propagators_preserves_secular_transfer_at_jordan_blocks() -> None:
    # 2x2 defective Jordan block in 4D Spectrometer shape (n_jeff=1, n_b1=1, 2, 2)
    liouv_2x2 = np.array([[-110.0, 100.0], [0.0, -110.0]])[np.newaxis, np.newaxis, :, :]
    delay = 0.02
    expected_single = np.exp(-110.0 * delay) * np.array(
        [[1.0, 100.0 * delay], [0.0, 1.0]]
    )

    got_single = calculate_propagators(liouv_2x2, delay)
    assert got_single.shape == (1, 1, 2, 2)
    np.testing.assert_allclose(
        got_single[0, 0], expected_single, rtol=1e-14, atol=1e-15
    )

    delays = [0.0, 0.02, 0.05]
    got_multi = calculate_propagators(liouv_2x2, delays)
    assert got_multi.shape == (3, 1, 1, 2, 2)
    for idx, t in enumerate(delays):
        expected = np.exp(-110.0 * t) * np.array([[1.0, 100.0 * t], [0.0, 1.0]])
        np.testing.assert_allclose(
            got_multi[idx, 0, 0], expected, rtol=1e-14, atol=1e-15
        )

    # 3x3 Jordan block verifying quadratic secular transfer (k * t)^2 / 2
    rate = -25.0
    coupling = 40.0
    liouv_3x3 = np.array(
        [[rate, coupling, 0.0], [0.0, rate, coupling], [0.0, 0.0, rate]]
    )[np.newaxis, np.newaxis, :, :]
    expected_3x3 = np.exp(rate * delay) * np.array(
        [
            [1.0, coupling * delay, 0.5 * (coupling * delay) ** 2],
            [0.0, 1.0, coupling * delay],
            [0.0, 0.0, 1.0],
        ]
    )
    np.testing.assert_allclose(
        calculate_propagators(liouv_3x3, delay)[0, 0],
        expected_3x3,
        rtol=1e-14,
        atol=1e-15,
    )


def test_spectrometer_delays_accurate_at_coalescence_and_pb_zero_boundary() -> None:
    register_kinetic_settings()
    model_2st = ModelSpec.from_name("2st")
    basis = Basis(type="ixy", spin_system="nh", model=model_2st)
    spectrometer = Spectrometer.from_spin_system(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    spectrometer.carrier_i = 120.0

    # 1. Exact two-state coalescence: pa = pb = 0.5, r2_a = r2_b, kex = |dw_rad_s|
    k = 100.0
    dw_ppm = (2.0 * k) / abs(spectrometer.ppm_i)
    for eps in (0.0, 1.0e-12, -1.0e-12, 1.0e-8):
        spectrometer.update(
            {
                "pa": 0.5,
                "pb": 0.5,
                "kab": k * (1.0 + eps),
                "kba": k * (1.0 + eps),
                "cs_i_a": 120.0 - 0.5 * dw_ppm,
                "cs_i_b": 120.0 + 0.5 * dw_ppm,
                "r2_i_a": 15.0,
                "r2_i_b": 15.0,
            }
        )
        expected = expm(spectrometer._engine.l_free[0, 0] * 0.02)
        np.testing.assert_allclose(
            spectrometer.delays(0.02)[0, 0],
            expected,
            rtol=1e-13,
            atol=1e-14,
        )

    # 2. Defective pb = 0 boundary: kab = 0, dw = 0, r2_a - r2_b = kba
    spectrometer.update(
        {
            "pa": 1.0,
            "pb": 0.0,
            "kab": 0.0,
            "kba": 100.0,
            "cs_i_a": 120.0,
            "cs_i_b": 120.0,
            "r2_i_a": 115.0,
            "r2_i_b": 15.0,
        }
    )
    expected_pb0 = expm(spectrometer._engine.l_free[0, 0] * 0.02)
    np.testing.assert_allclose(
        spectrometer.delays(0.02)[0, 0],
        expected_pb0,
        rtol=1e-14,
        atol=1e-15,
    )


def test_calculate_propagators_handles_mixed_norm_batches_without_state_bleed() -> None:
    # Guards against SciPy 1.18.0 batched expm squaring-count bleed when a large-norm
    # slice (requiring squarings) precedes a small-norm or nilpotent slice.
    liouv_big = np.array([[0.0, -20.0], [20.0, 0.0]])
    liouv_small = np.array([[0.0, -0.5], [0.5, 0.0]])
    liouv_nilpotent = np.array([[0.0, 100.0], [0.0, 0.0]])
    batch = np.stack([liouv_big, liouv_small, liouv_nilpotent])[:, np.newaxis, :, :]

    got = calculate_propagators(batch, 1.0)
    expected = np.stack(
        [
            np.array([[np.cos(20.0), -np.sin(20.0)], [np.sin(20.0), np.cos(20.0)]]),
            np.array([[np.cos(0.5), -np.sin(0.5)], [np.sin(0.5), np.cos(0.5)]]),
            np.array([[1.0, 100.0], [0.0, 1.0]]),
        ]
    )[:, np.newaxis, :, :]

    assert got.shape == expected.shape
    np.testing.assert_allclose(got, expected, rtol=1e-13, atol=1e-14)


def test_complete_dephasing_preserves_secular_transfer_on_real_jordan_blocks() -> None:
    # Pure real defective Jordan block must retain secular transfer under dephasing=True
    liouv_nilpotent = np.array([[0.0, 1.0], [0.0, 0.0]])
    np.testing.assert_allclose(
        calculate_propagators(liouv_nilpotent, 0.25, dephasing=True),
        np.array([[1.0, 0.25], [0.0, 1.0]]),
        rtol=0.0,
        atol=1.0e-14,
    )

    # Realistic 2-state CEST spin-lock Liouvillian at pb=0, dw=0, r2_a - r2_b = kba:
    # ix_a / ix_b form a real non-oscillatory Jordan block along B1x while iy/iz oscillate.
    register_kinetic_settings()
    model_2st = ModelSpec.from_name("2st")
    basis = Basis(type="ixyz_eq", spin_system="nh", model=model_2st)
    spectrometer = Spectrometer.from_spin_system(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    spectrometer.carrier_i = 120.0
    spectrometer.b1_i = 25.0
    spectrometer.update(
        {
            "pa": 1.0,
            "pb": 0.0,
            "kab": 0.0,
            "kba": 100.0,
            "cs_i_a": 120.0,
            "cs_i_b": 120.0,
            "r2_i_a": 115.0,
            "r2_i_b": 15.0,
            "r1_i_a": 101.5,
            "r1_i_b": 1.5,
        }
    )
    delay = 0.02
    liouv = spectrometer._engine.l_free + spectrometer._engine.l_b1x_i
    propagator = calculate_propagators(liouv, delay, dephasing=True)
    # In ixyz_eq (components: ie, ix, iy, iz), ix_a is index 1 and ix_b is index 5
    expected_ix_transfer = delay * 100.0 * np.exp(-115.0 * delay)
    assert propagator.shape == (1, 1, 8, 8)
    np.testing.assert_allclose(
        propagator[0, 0, 1, 5],
        expected_ix_transfer,
        rtol=1e-13,
        atol=1e-14,
    )


@pytest.mark.parametrize("dephasing", [False, True])
def test_calculate_propagators_squeezes_only_scalar_delays(dephasing: bool) -> None:
    liouv_2d = np.array([[-2.0, 0.0], [0.0, -3.0]])
    liouv_4d = np.broadcast_to(liouv_2d, (2, 3, 2, 2)).copy()
    liouv_5d = np.broadcast_to(liouv_2d, (5, 2, 3, 2, 2)).copy()

    # Scalar inputs drop the leading delay axis
    assert calculate_propagators(liouv_2d, 0.25, dephasing=dephasing).shape == (2, 2)
    assert calculate_propagators(
        liouv_4d, np.float64(0.25), dephasing=dephasing
    ).shape == (2, 3, 2, 2)
    assert calculate_propagators(
        liouv_4d, np.array(0.25), dephasing=dephasing
    ).shape == (2, 3, 2, 2)
    assert calculate_propagators(liouv_5d, 0.25, dephasing=dephasing).shape == (
        5,
        2,
        3,
        2,
        2,
    )

    # 1-element 1D sequences / iterables preserve the leading delay axis of length 1
    assert calculate_propagators(liouv_2d, [0.25], dephasing=dephasing).shape == (
        1,
        2,
        2,
    )
    assert calculate_propagators(liouv_4d, [0.25], dephasing=dephasing).shape == (
        1,
        2,
        3,
        2,
        2,
    )
    assert calculate_propagators(liouv_4d, (0.25,), dephasing=dephasing).shape == (
        1,
        2,
        3,
        2,
        2,
    )
    assert calculate_propagators(
        liouv_4d, np.array([0.25]), dephasing=dephasing
    ).shape == (1, 2, 3, 2, 2)
    assert calculate_propagators(
        liouv_4d, (t for t in [0.25]), dephasing=dephasing
    ).shape == (1, 2, 3, 2, 2)
    assert calculate_propagators(
        liouv_4d, {0: 0.25}.values(), dephasing=dephasing
    ).shape == (1, 2, 3, 2, 2)
    assert calculate_propagators(liouv_4d, [], dephasing=dephasing).shape == (
        0,
        2,
        3,
        2,
        2,
    )


def test_single_delay_sequences_work_under_multi_point_distributions() -> None:
    register_kinetic_settings()
    jeff_dist = Distribution(np.array([-17.5, 17.5]), np.array([0.5, 0.5]))

    # 1. cpmg_ch3_mq with reference-only metadata (all_delays has length 1)
    mq_cfg = CpmgCh3MqConfig.model_validate(
        {
            "experiment": {
                "name": "cpmg_ch3_mq",
                "time_t2": 0.04,
                "t_zeta": 1.0 / (8.0 * 125.0),
            },
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    mq_cfg.model = ModelSpec.from_name("1st")
    mq_spec = build_ch3_mq_spectrometer(mq_cfg, SpinSystem.from_name("L18CD1-HD1"))
    mq_spec.jeff_i = jeff_dist
    mq_spec.update(
        {"pa": 1.0, "r2mq_is_a": 10.0, "mu_is_a": 0.0, "cs_i_a": 22.0, "cs_s_a": 0.8}
    )
    mq_out = CpmgCh3MqSequence(mq_cfg.experiment).calculate(
        mq_spec,
        Data(exp=np.array([1.0]), err=np.array([0.01]), metadata=np.array([0.0])),
    )
    assert mq_out.shape == (1,)
    assert np.isfinite(mq_out).all()

    # 2. relaxation_nz & relaxation_hznz with a 1-point time array under distributed jeff_i
    nz_cfg = RelaxationNzConfig.model_validate(
        {
            "experiment": {"name": "relaxation_nz"},
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    nz_cfg.model = ModelSpec.from_name("1st")
    nz_spec = build_nz_spectrometer(nz_cfg, SpinSystem.from_name("G23N-HN"))
    nz_spec.jeff_i = jeff_dist
    nz_spec.update({"pa": 1.0, "r1_i_a": 1.5})
    nz_out = RelaxationNzSequence(nz_cfg.experiment).calculate(
        nz_spec,
        Data(exp=np.array([0.8]), err=np.array([0.01]), metadata=np.array([0.1])),
    )
    assert nz_out.shape == (1,)
    expected_nz = nz_spec.get_start_magnetization(nz_cfg.experiment.start_terms)[0, 0]
    np.testing.assert_allclose(nz_out, [expected_nz * np.exp(-0.15)])

    hznz_cfg = RelaxationHzNzConfig.model_validate(
        {
            "experiment": {"name": "relaxation_hznz"},
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    hznz_cfg.model = ModelSpec.from_name("1st")
    hznz_spec = build_hznz_spectrometer(hznz_cfg, SpinSystem.from_name("G23N-HN"))
    hznz_spec.jeff_i = jeff_dist
    hznz_spec.update({"pa": 1.0, "r1_i_a": 1.5, "r1a_is_a": 2.5, "etaz_i_a": 0.0})
    hznz_out = RelaxationHzNzSequence(hznz_cfg.experiment).calculate(
        hznz_spec,
        Data(exp=np.array([0.8]), err=np.array([0.01]), metadata=np.array([0.1])),
    )
    assert hznz_out.shape == (1,)
    assert np.isfinite(hznz_out).all()

    # 3. cpmg_ch3_1h_tq_diff under distributed jeff_i (scalar delays must not unpack axis 0)
    tq_cfg = CpmgCh31HTqDiffConfig.model_validate(
        {
            "experiment": {
                "name": "cpmg_ch3_1h_tq_diff",
                "time_t2": 0.04,
                "carrier": 0.8,
                "pw90": 10.0e-6,
                "delta": 1.5e-3,
                "gradient": 20.0,
                "tauc": 0.001,
            },
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    tq_cfg.model = ModelSpec.from_name("1st")
    tq_spec = build_tq_diff_spectrometer(tq_cfg, SpinSystem.from_name("L18HD1-CD1"))
    tq_spec.jeff_i = jeff_dist
    tq_spec.update(
        {
            "pa": 1.0,
            "r1_i_a": 1.5,
            "r2_i_a": 15.0,
            "r2a_i_a": 16.0,
            "r1a_is_a": 2.0,
            "etaxy_i_a": 0.0,
            "etaz_i_a": 0.0,
            "cs_i_a": 0.8,
            "j_is_a": 125.0,
            "d_a": 0.0,
        }
    )
    tq_out = CpmgCh31HTqDiffSequence(tq_cfg.experiment).calculate(
        tq_spec,
        Data(exp=np.array([1.0]), err=np.array([0.01]), metadata=np.array([2.0])),
    )
    assert tq_out.shape == (1,)
    assert np.isfinite(tq_out).all()


def test_calculate_propagators_rejects_invalid_shapes() -> None:
    with pytest.raises(ValueError, match="square matrix dimensions"):
        calculate_propagators(np.zeros((2, 3)), 0.1)

    with pytest.raises(ValueError, match="scalar or a 1D iterable"):
        calculate_propagators(np.zeros((2, 2)), np.zeros((1, 2)))

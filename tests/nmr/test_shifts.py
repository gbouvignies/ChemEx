from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.containers.data import Data
from chemex.experiments.catalog.shift_13c_sq import (
    Shift13CSqConfig,
    Shift13CSqSequence,
)
from chemex.experiments.catalog.shift_13c_sq import (
    build_spectrometer as build_13c_sq_spectrometer,
)
from chemex.experiments.catalog.shift_15n_sq import (
    Shift15NSqConfig,
    Shift15NSqSequence,
)
from chemex.experiments.catalog.shift_15n_sq import (
    build_spectrometer as build_15n_sq_spectrometer,
)
from chemex.experiments.catalog.shift_15n_sqmq import (
    Shift15NSqMqConfig,
    Shift15NSqMqSequence,
)
from chemex.experiments.catalog.shift_15n_sqmq import (
    build_spectrometer as build_15n_sqmq_spectrometer,
)
from chemex.models.loader import register_kinetic_settings
from chemex.models.model import ModelSpec
from chemex.nmr._engine.analysis import _select_observed_frequency_offset
from chemex.nmr._engine.engine import ISLiouvillianEngine
from chemex.nmr.basis import Basis
from chemex.nmr.constants import XI_RATIO, Distribution
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spin_system.nucleus import Nucleus

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module", autouse=True)
def _register_models() -> None:
    register_kinetic_settings()


def _analytical_two_state_frequency(
    *,
    omega_a: float,
    omega_b: float,
    r2_a: float,
    r2_b: float,
    kab: float,
    kba: float,
    observed_state: str = "a",
) -> float:
    """Closed-form scalar two-state Bloch-McConnell precession frequency (rad/s)."""
    alpha = complex(-r2_a - kab, omega_a)
    beta = complex(-r2_b - kba, omega_b)
    diff = alpha - beta
    disc = np.lib.scimath.sqrt(diff * diff + 4.0 * kab * kba)
    eig_plus = 0.5 * (alpha + beta + disc)
    eig_minus = 0.5 * (alpha + beta - disc)

    if abs(eig_plus - eig_minus) <= 1e-12 * max(1.0, abs(alpha), abs(beta)):
        return float(0.5 * (omega_a + omega_b))

    liouv = np.array([[alpha, kba], [kab, beta]], dtype=np.complex128)
    eye = np.eye(2, dtype=np.complex128)
    obs_idx = 0 if observed_state == "a" else 1

    # Exact Cayley-Hamilton Riesz projectors P_+ and P_-
    proj_plus_col = (liouv - eig_minus * eye)[:, obs_idx] / (eig_plus - eig_minus)
    proj_minus_col = (liouv - eig_plus * eye)[:, obs_idx] / (eig_minus - eig_plus)
    amp_plus = proj_plus_col.sum()
    amp_minus = proj_minus_col.sum()

    return float(eig_plus.imag if abs(amp_plus) >= abs(amp_minus) else eig_minus.imag)


def make_engine() -> ISLiouvillianEngine:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    engine = ISLiouvillianEngine(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    engine.update({"r2_i_a": 4.0, "cs_i_a": 1.0})
    return engine


def make_spectrometer() -> Spectrometer:
    return Spectrometer(make_engine())


def test_calculate_shifts_returns_signed_coherence_subspace_frequencies() -> None:
    engine = make_engine()
    spectrometer = Spectrometer(engine)
    expected = np.array([0.0, 1.0 * engine.ppm_i])

    shifts = spectrometer.analysis.calculate_shifts()

    assert shifts.shape == (2,)
    np.testing.assert_allclose(
        np.sort(shifts), np.sort(expected), rtol=1e-13, atol=1e-12
    )


def test_calculate_shifts_rejects_multi_point_distributions() -> None:
    spectrometer = make_spectrometer()
    spectrometer.jeff_i = Distribution(
        np.array([-15.0, 0.0, 15.0]),
        np.array([0.2, 0.6, 0.2]),
    )

    with pytest.raises(
        ValueError,
        match="requires single-point B1 and Jeff distributions",
    ):
        spectrometer.analysis.calculate_shifts()

    with pytest.raises(
        ValueError,
        match="requires single-point B1 and Jeff distributions",
    ):
        spectrometer.analysis.calculate_shift("a")


def test_calculate_shift_validates_coherence_subspace_and_state() -> None:
    model = ModelSpec()
    iz_spec = Spectrometer.from_spin_system(
        SpinSystem(name="G23N-HN"),
        Basis(type="iz", spin_system="nh", model=model),
        Conditions(h_larmor_frq=600.0),
    )
    with pytest.raises(ValueError, match="lacks components"):
        iz_spec.analysis.calculate_shift("a", coherence="sq_i")

    ixy_spec = Spectrometer.from_spin_system(
        SpinSystem(name="G23N-HN"),
        Basis(type="ixy", spin_system="nh", model=model),
        Conditions(h_larmor_frq=600.0),
    )
    with pytest.raises(ValueError, match="Unknown observed_state"):
        ixy_spec.analysis.calculate_shift("z", coherence="sq_i")

    coupled_spec = Spectrometer.from_spin_system(
        SpinSystem(name="G23N-HN"),
        Basis(type="ixyzsz", spin_system="nh", model=model),
        Conditions(h_larmor_frq=600.0),
    )
    coupled_spec.update({"j_is_a": -92.0, "j_is_b": -92.0})
    with pytest.raises(ValueError, match="coupled outside its scalar subspace"):
        coupled_spec.analysis.calculate_shift("a", coherence="sq_i")


@pytest.mark.parametrize(
    ("kex", "pb", "dw_ppm", "r2_a", "r2_b"),
    [
        (20.0, 0.10, 3.0, 8.0, 12.0),  # slow exchange
        (800.0, 0.15, -2.5, 10.0, 18.0),  # intermediate exchange
        (20000.0, 0.10, 4.0, 9.0, 9.0),  # fast exchange (major a)
        (20000.0, 0.85, -3.5, 9.0, 11.0),  # fast exchange past 50% titration (pa < 0.5)
    ],
)
def test_sq_shift_matches_analytical_two_state_bloch_mcconnell_for_15n_and_13c(
    kex: float,
    pb: float,
    dw_ppm: float,
    r2_a: float,
    r2_b: float,
) -> None:
    model = ModelSpec.from_name("2st")
    pa = 1.0 - pb
    kab = pb * kex
    kba = pa * kex
    dummy_data = Data(exp=np.zeros(1), err=np.ones(1))

    # 15N (gamma < 0, ppm_i > 0)
    cfg_n = Shift15NSqConfig.model_validate(
        {
            "experiment": {"name": "shift_15n_sq"},
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    cfg_n.model = model
    spec_n = build_15n_sq_spectrometer(cfg_n, SpinSystem.from_name("G23N-HN"))
    seq_n = Shift15NSqSequence(cfg_n.experiment)
    cs_n_a = 118.5
    spec_n.update(
        {
            "pa": pa,
            "pb": pb,
            "kab": kab,
            "kba": kba,
            "r2_i_a": r2_a,
            "r2_i_b": r2_b,
            "cs_i_a": cs_n_a,
            "cs_i_b": cs_n_a + dw_ppm,
        }
    )
    expected_n = (
        _analytical_two_state_frequency(
            omega_a=cs_n_a * spec_n.ppm_i,
            omega_b=(cs_n_a + dw_ppm) * spec_n.ppm_i,
            r2_a=r2_a,
            r2_b=r2_b,
            kab=kab,
            kba=kba,
        )
        / spec_n.ppm_i
    )
    np.testing.assert_allclose(
        seq_n.calculate(spec_n, dummy_data)[0],
        expected_n,
        rtol=0.0,
        atol=1e-12,
    )

    # 13C (gamma > 0, ppm_i < 0)
    cfg_c = Shift13CSqConfig.model_validate(
        {
            "experiment": {"name": "shift_13c_sq"},
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    cfg_c.model = model
    spec_c = build_13c_sq_spectrometer(cfg_c, SpinSystem.from_name("A5CB-HB"))
    seq_c = Shift13CSqSequence(cfg_c.experiment)
    cs_c_a = 18.2
    spec_c.update(
        {
            "pa": pa,
            "pb": pb,
            "kab": kab,
            "kba": kba,
            "r2_i_a": r2_a,
            "r2_i_b": r2_b,
            "cs_i_a": cs_c_a,
            "cs_i_b": cs_c_a + dw_ppm,
        }
    )
    expected_c = (
        _analytical_two_state_frequency(
            omega_a=cs_c_a * spec_c.ppm_i,
            omega_b=(cs_c_a + dw_ppm) * spec_c.ppm_i,
            r2_a=r2_a,
            r2_b=r2_b,
            kab=kab,
            kba=kba,
        )
        / spec_c.ppm_i
    )
    np.testing.assert_allclose(
        seq_c.calculate(spec_c, dummy_data)[0],
        expected_c,
        rtol=0.0,
        atol=1e-12,
    )


def test_13c_and_15n_sq_shifts_do_not_bounce_across_zero_ppm() -> None:
    model = ModelSpec.from_name("2st")
    dummy_data = Data(exp=np.zeros(1), err=np.ones(1))
    cfg_c = Shift13CSqConfig.model_validate(
        {
            "experiment": {"name": "shift_13c_sq"},
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    cfg_c.model = model
    spec_c = build_13c_sq_spectrometer(cfg_c, SpinSystem.from_name("A5CB"))
    seq_c = Shift13CSqSequence(cfg_c.experiment)
    base = {
        "pa": 0.90,
        "pb": 0.10,
        "kab": 500.0,
        "kba": 4500.0,
        "r2_i_a": 10.0,
        "r2_i_b": 10.0,
    }
    dw = -2.0
    offsets = []
    for cs_a in (15.0, 0.05, 0.0, -0.05, -3.0):
        spec_c.update({**base, "cs_i_a": cs_a, "cs_i_b": cs_a + dw})
        offsets.append(seq_c.calculate(spec_c, dummy_data)[0] - cs_a)

    for offset in offsets[1:]:
        assert offset < 0.0
        np.testing.assert_allclose(offset, offsets[0], rtol=0.0, atol=1e-13)


def test_shift_15n_sqmq_invariance_to_cs_a_glycine_dq_zero_and_skrynnikov_formula() -> (
    None
):
    model = ModelSpec.from_name("2st")
    cfg = Shift15NSqMqConfig.model_validate(
        {
            "experiment": {"name": "shift_15n_sqmq"},
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    cfg.model = model
    spec = build_15n_sqmq_spectrometer(cfg, SpinSystem.from_name("G10N-HN"))
    seq = Shift15NSqMqSequence(cfg.experiment)
    dummy_data = Data(exp=np.zeros(1), err=np.ones(1))

    pa, pb = 0.90, 0.10
    kex = 5000.0
    kab, kba = pb * kex, pa * kex
    dw_n, dw_h = -2.0, 0.20
    r2_i_a, r2_i_b = 6.0, 8.0
    r2mq_a, r2mq_b = 11.0, 14.0
    mu_a, mu_b = 1.5, -0.5

    ratio = XI_RATIO[Nucleus.N15] / XI_RATIO[Nucleus.H1]
    gly_cs_n = 102.0
    gly_cs_h = gly_cs_n * ratio - 0.005  # places omega_DQ right near 0 rad/s

    shifts_ppb = []
    for cs_n_a, cs_h_a in ((0.0, 0.0), (120.0, 8.5), (gly_cs_n, gly_cs_h)):
        spec.update(
            {
                "pa": pa,
                "pb": pb,
                "kab": kab,
                "kba": kba,
                "cs_i_a": cs_n_a,
                "cs_i_b": cs_n_a + dw_n,
                "cs_s_a": cs_h_a,
                "cs_s_b": cs_h_a + dw_h,
                "r2_i_a": r2_i_a,
                "r2_i_b": r2_i_b,
                "r2mq_is_a": r2mq_a,
                "r2mq_is_b": r2mq_b,
                "mu_is_a": mu_a,
                "mu_is_b": mu_b,
            }
        )
        shifts_ppb.append(seq.calculate(spec, dummy_data)[0])

    # Exact invariance to CS_A (including omitted CS_A = 0.0 and Glycine omega_DQ ~= 0)
    # Within binary64 ULP of `(cs_s_a + dw_h) * ppm_s - cs_s_a * ppm_s` (~1e-12 ppb)
    np.testing.assert_allclose(shifts_ppb[0], shifts_ppb[1], rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(shifts_ppb[2], shifts_ppb[1], rtol=0.0, atol=1e-10)

    # Exact match to analytical Bloch-McConnell with DQ/ZQ cross-correlated mu_is
    dw_i_rad = dw_n * spec.ppm_i
    dw_s_rad = dw_h * spec.ppm_s
    w_sq = _analytical_two_state_frequency(
        omega_a=0.0,
        omega_b=dw_i_rad,
        r2_a=r2_i_a,
        r2_b=r2_i_b,
        kab=kab,
        kba=kba,
    )
    w_dq = _analytical_two_state_frequency(
        omega_a=0.0,
        omega_b=dw_i_rad + dw_s_rad,
        r2_a=r2mq_a + mu_a,
        r2_b=r2mq_b + mu_b,
        kab=kab,
        kba=kba,
    )
    w_zq = _analytical_two_state_frequency(
        omega_a=0.0,
        omega_b=dw_i_rad - dw_s_rad,
        r2_a=r2mq_a - mu_a,
        r2_b=r2mq_b - mu_b,
        kab=kab,
        kba=kba,
    )
    expected_ppb = 1e3 * (w_sq - 0.5 * (w_dq + w_zq)) / spec.ppm_i
    np.testing.assert_allclose(shifts_ppb[0], expected_ppb, rtol=0.0, atol=1e-10)

    # Fast-exchange asymptotic verification against Skrynnikov, Dahlquist & Kay (2002):
    # Delta_delta_SQMQ (ppb) -> 3000 * pa * pb * (pa - pb) * dw_N * (dw_H * ppm_s)^2 / kex^2
    kex_fast = 80000.0
    spec.update(
        {
            "pa": pa,
            "pb": pb,
            "kab": pb * kex_fast,
            "kba": pa * kex_fast,
            "cs_i_a": 0.0,
            "cs_i_b": dw_n,
            "cs_s_a": 0.0,
            "cs_s_b": dw_h,
            "r2_i_a": 5.0,
            "r2_i_b": 5.0,
            "r2mq_is_a": 10.0,
            "r2mq_is_b": 10.0,
            "mu_is_a": 0.0,
            "mu_is_b": 0.0,
        }
    )
    calc_fast_ppb = seq.calculate(spec, dummy_data)[0]
    sdk_asymptotic_ppb = (
        3000.0 * pa * pb * (pa - pb) * dw_n * (dw_s_rad**2) / (kex_fast**2)
    )
    np.testing.assert_allclose(calc_fast_ppb, sdk_asymptotic_ppb, rtol=2e-3)


def test_coalescence_exceptional_points_and_jordan_blocks_are_stable() -> None:
    # 1. Symmetric coalescence exceptional point (pa = pb = 0.5, dR2 = 0, kex = |dw|)
    model_2st = ModelSpec.from_name("2st")
    cfg = Shift15NSqConfig.model_validate(
        {
            "experiment": {"name": "shift_15n_sq"},
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    cfg.model = model_2st
    spec = build_15n_sq_spectrometer(cfg, SpinSystem.from_name("G23N-HN"))
    seq = Shift15NSqSequence(cfg.experiment)
    dummy_data = Data(exp=np.zeros(1), err=np.ones(1))

    dw_ppm = 2.0
    kex_ep = abs(dw_ppm * spec.ppm_i)
    for rel_dk in (-1e-3, -1e-7, 0.0, 1e-7, 1e-3):
        kex = kex_ep * (1.0 + rel_dk)
        spec.update(
            {
                "pa": 0.5,
                "pb": 0.5,
                "kab": 0.5 * kex,
                "kba": 0.5 * kex,
                "r2_i_a": 10.0,
                "r2_i_b": 10.0,
                "cs_i_a": 115.0,
                "cs_i_b": 115.0 + dw_ppm,
            }
        )
        calc = seq.calculate(spec, dummy_data)[0]
        if rel_dk == 0.0:
            expected = 115.0 + 0.5 * dw_ppm
        else:
            expected = 115.0 + (
                _analytical_two_state_frequency(
                    omega_a=0.0,
                    omega_b=dw_ppm * spec.ppm_i,
                    r2_a=10.0,
                    r2_b=10.0,
                    kab=0.5 * kex,
                    kba=0.5 * kex,
                )
                / spec.ppm_i
            )
        np.testing.assert_allclose(calc, expected, rtol=0.0, atol=1e-9)

    # 2. Asymmetric coalescence exceptional point (pa = 0.8, pb = 0.2, dR2 != 0)
    pa, pb = 0.8, 0.2
    kex_asym = 1000.0
    r2_a = 700.0
    r2_b = r2_a + (pb - pa) * kex_asym  # 100.0 -> equal diagonal real parts
    dw_asym_ppm = (2.0 * np.sqrt(pa * pb) * kex_asym) / spec.ppm_i
    spec.update(
        {
            "pa": pa,
            "pb": pb,
            "kab": pb * kex_asym,
            "kba": pa * kex_asym,
            "r2_i_a": r2_a,
            "r2_i_b": r2_b,
            "cs_i_a": 120.0,
            "cs_i_b": 120.0 + dw_asym_ppm,
        }
    )
    calc_asym_ep = seq.calculate(spec, dummy_data)[0]
    np.testing.assert_allclose(
        calc_asym_ep,
        120.0 + 0.5 * dw_asym_ppm,
        rtol=0.0,
        atol=1e-12,
    )

    # 3. 3x3 matrix with exact 2x2 Jordan block in (b, c) subspace
    liouv_jordan = np.array(
        [
            [-10.0 + 1000.0j, 0.0, 0.0],
            [0.0, -110.0 + 5000.0j, 100.0],
            [0.0, 0.0, -110.0 + 5000.0j],
        ],
        dtype=np.complex128,
    )
    assert _select_observed_frequency_offset(liouv_jordan, 0) == pytest.approx(0.0)
    assert _select_observed_frequency_offset(liouv_jordan, 1) == pytest.approx(0.0)
    assert _select_observed_frequency_offset(liouv_jordan, 2) == pytest.approx(0.0)


def test_three_state_dark_mode_zero_population_states_and_fd_derivative() -> None:
    model_3st = ModelSpec.from_name("3st")
    dummy_data = Data(exp=np.zeros(1), err=np.ones(1))

    # 1. Symmetric fork a <-> c <-> b with dark antisymmetric mode (1^T v_dark == 0)
    cfg_a = Shift15NSqConfig.model_validate(
        {
            "experiment": {"name": "shift_15n_sq", "observed_state": "a"},
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    cfg_a.model = model_3st
    spec_3 = build_15n_sq_spectrometer(cfg_a, SpinSystem.from_name("G23N-HN"))
    seq_a = Shift15NSqSequence(cfg_a.experiment)
    spec_3.update(
        {
            "pa": 0.4,
            "pb": 0.4,
            "pc": 0.2,
            "kab": 0.0,
            "kba": 0.0,
            "kac": 50.0,
            "kca": 100.0,
            "kbc": 50.0,
            "kcb": 100.0,
            "r2_i_a": 10.0,
            "r2_i_b": 10.0,
            "r2_i_c": 10.0,
            "cs_i_a": 120.0,
            "cs_i_b": 120.0,
            "cs_i_c": 123.0,
        }
    )
    expected_fork = (
        _analytical_two_state_frequency(
            omega_a=120.0 * spec_3.ppm_i,
            omega_b=123.0 * spec_3.ppm_i,
            r2_a=10.0,
            r2_b=10.0,
            kab=50.0,
            kba=200.0,
        )
        / spec_3.ppm_i
    )
    np.testing.assert_allclose(
        seq_a.calculate(spec_3, dummy_data)[0],
        expected_fork,
        rtol=0.0,
        atol=1e-12,
    )

    # 2. Zero-population states (pb = pc = 0, kex = 0) with observed_state a, b, c
    spec_3.update(
        {
            "pa": 1.0,
            "pb": 0.0,
            "pc": 0.0,
            "kab": 0.0,
            "kba": 0.0,
            "kac": 0.0,
            "kca": 0.0,
            "kbc": 0.0,
            "kcb": 0.0,
            "r2_i_a": 8.0,
            "r2_i_b": 12.0,
            "r2_i_c": 16.0,
            "cs_i_a": 115.0,
            "cs_i_b": 118.5,
            "cs_i_c": 122.25,
        }
    )
    for state, expected_cs in (("a", 115.0), ("b", 118.5), ("c", 122.25)):
        cfg_s = Shift15NSqConfig.model_validate(
            {
                "experiment": {"name": "shift_15n_sq", "observed_state": state},
                "conditions": {"h_larmor_frq": 600.0},
                "data": {},
            }
        )
        cfg_s.model = model_3st
        seq_s = Shift15NSqSequence(cfg_s.experiment)
        assert seq_s.calculate(spec_3, dummy_data)[0] == pytest.approx(
            expected_cs, abs=1e-12
        )

    # 3. Forward finite-difference derivative at DW_AB = 0 matches pb in fast exchange
    cfg_2 = Shift15NSqConfig.model_validate(
        {
            "experiment": {"name": "shift_15n_sq"},
            "conditions": {"h_larmor_frq": 600.0},
            "data": {},
        }
    )
    cfg_2.model = ModelSpec.from_name("2st")
    spec_2 = build_15n_sq_spectrometer(cfg_2, SpinSystem.from_name("G23N-HN"))
    seq_2 = Shift15NSqSequence(cfg_2.experiment)
    step = float(np.sqrt(np.finfo(np.float64).eps))
    shifts_fd = []
    for dw in (0.0, step):
        spec_2.update(
            {
                "pa": 0.9,
                "pb": 0.1,
                "kab": 50.0,
                "kba": 450.0,
                "r2_i_a": 10.0,
                "r2_i_b": 10.0,
                "cs_i_a": 120.0,
                "cs_i_b": 120.0 + dw,
            }
        )
        shifts_fd.append(seq_2.calculate(spec_2, dummy_data)[0])
    fd_slope = (shifts_fd[1] - shifts_fd[0]) / step
    np.testing.assert_allclose(fd_slope, 0.1, rtol=1e-5)

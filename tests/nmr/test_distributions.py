"""Scientific and numerical tests for the B1 distribution kernel and quadrature."""

from __future__ import annotations

import math
from typing import Any, Literal

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy import integrate, stats
from scipy.special import log_ndtr, ndtr, ndtri

from chemex.configuration.b1_config import parse_distribution_config
from chemex.containers.data import Data
from chemex.experiments.catalog import (
    cest_1hn_ip_ap,
    cest_15n,
    coscest_1hn_ip_ap,
    dcest_15n,
)
from chemex.models.model import ModelSpec
from chemex.nmr.constants import Distribution
from chemex.nmr.distributions import (
    custom,
    dephasing,
    gaussian,
    homogeneous,
    skewed,
)
from chemex.nmr.distributions._quadrature import (
    MAX_SKEWNESS,
    compute_skewnorm_mode_z,
    skewness_to_alpha_and_sigma_z,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import ensure_plugins_registered

_HALF_LOG_2PI = 0.5 * math.log(2.0 * math.pi)
_SQRT_2_OVER_PI = math.sqrt(2.0 / math.pi)


def _reference_cest_distribution(
    nominal: float,
    scale: float,
    skewness: float,
    center: str,
    *,
    panels: int = 200,
    points_per_panel: int = 10,
) -> Distribution:
    """Build a 2,000-node composite Gauss-Legendre oracle on p_+(B1)."""
    alpha, sigma_z = skewness_to_alpha_and_sigma_z(skewness)
    if alpha == 0.0:
        m_anchor = 0.0
        z0 = -1.0 / scale
        f0 = float(ndtr(np.float64(z0)))
        s0 = float(ndtr(np.float64(-z0)))
        z_lo = float(ndtri(np.float64(f0 + 1.0e-9 * s0)))
        z_hi = -float(ndtri(np.float64(1.0e-9 * s0)))
    else:
        delta = alpha / math.sqrt(1.0 + alpha * alpha)
        m_anchor = (
            delta * _SQRT_2_OVER_PI
            if center == "mean"
            else compute_skewnorm_mode_z(alpha)
        )
        z0 = m_anchor - sigma_z / scale
        f0 = float(stats.skewnorm.cdf(z0, alpha))
        s0 = float(stats.skewnorm.sf(z0, alpha))
        z_lo = float(stats.skewnorm.ppf(f0 + 1.0e-9 * s0, alpha))
        z_hi = float(stats.skewnorm.isf(1.0e-9 * s0, alpha))

    gl_x, gl_w = leggauss(points_per_panel)
    edges = np.linspace(z_lo, z_hi, panels + 1)
    z_parts: list[np.ndarray] = []
    lw_parts: list[np.ndarray] = []
    for idx in range(panels):
        lo, hi = edges[idx], edges[idx + 1]
        z = 0.5 * (hi - lo) * gl_x + 0.5 * (hi + lo)
        lw = np.log(gl_w * 0.5 * (hi - lo)) - 0.5 * z * z
        if alpha != 0.0:
            lw = lw + log_ndtr(alpha * z)
        z_parts.append(z)
        lw_parts.append(lw)

    z_all = np.concatenate(z_parts)
    lw_all = np.concatenate(lw_parts)
    weights = np.exp(lw_all - np.max(lw_all))
    weights /= np.sum(weights)
    scales = np.maximum(0.0, 1.0 + scale * ((z_all - m_anchor) / sigma_z))
    return Distribution(float(nominal) * scales, weights)


class TestFiveModelInvariants:
    """Validate physical support, normalization, flags, and edge inputs across all 5 models."""

    def test_homogeneous_and_dephasing_single_point_contracts(self) -> None:
        dist_h = homogeneous.generate(25.0)
        np.testing.assert_allclose(dist_h.values, [25.0])
        np.testing.assert_allclose(dist_h.weights, [1.0])
        assert dist_h.dephasing is False

        dist_h_zero = homogeneous.generate(0.0)
        np.testing.assert_allclose(dist_h_zero.values, [0.0])
        np.testing.assert_allclose(dist_h_zero.weights, [1.0])

        dist_d = dephasing.generate(25.0)
        np.testing.assert_allclose(dist_d.values, [25.0])
        np.testing.assert_allclose(dist_d.weights, [1.0])
        assert dist_d.dephasing is True

        for fn in (homogeneous.generate, dephasing.generate):
            for bad in (-1.0, float("nan"), float("inf")):
                with pytest.raises(ValueError, match="non-negative"):
                    fn(bad)

    @pytest.mark.parametrize(
        ("dist_type", "extra"),
        [
            ("gaussian", {}),
            ("skewed", {"skewness": -0.75, "center": "mode"}),
            ("skewed", {"skewness": -0.75, "center": "mean"}),
            ("skewed", {"skewness": 0.75, "center": "mode"}),
            ("skewed", {"skewness": 0.75, "center": "mean"}),
        ],
    )
    @pytest.mark.parametrize("scale", [0.05, 0.25, 0.5, 1.0])
    @pytest.mark.parametrize("res", [5, 13, 25])
    def test_quadrature_nodes_non_negative_ordered_and_normalized(
        self,
        dist_type: str,
        extra: dict[str, object],
        scale: float,
        res: int,
    ) -> None:
        config = parse_distribution_config(
            {"type": dist_type, "scale": scale, "res": res, **extra}
        )
        dist = config.get_distribution(25.0)
        assert dist.dephasing is False
        assert dist.values.size == res
        assert dist.weights.size == res
        assert np.all(dist.values >= 0.0)
        assert np.all(np.diff(dist.values) > 0.0)
        assert np.all(dist.weights > 0.0)
        np.testing.assert_allclose(dist.weights.sum(), 1.0, rtol=0.0, atol=1.0e-13)

    @pytest.mark.parametrize("tiny_scale", [0.0, 1.0e-12, 1.0e-100])
    def test_zero_and_subnormal_scales_remain_finite_and_centered(
        self,
        tiny_scale: float,
    ) -> None:
        for fn in (
            lambda s: gaussian.generate(25.0, s, res=15),
            lambda s: skewed.generate(25.0, s, -0.6, res=15, center="mode"),
            lambda s: skewed.generate(25.0, s, -0.6, res=15, center="mean"),
        ):
            dist = fn(tiny_scale)
            assert np.all(np.isfinite(dist.values))
            assert np.all(np.isfinite(dist.weights))
            assert np.all(dist.weights > 0.0)
            np.testing.assert_allclose(dist.weights.sum(), 1.0, rtol=0.0, atol=1.0e-13)
            np.testing.assert_allclose(
                np.sum(dist.values * dist.weights), 25.0, rtol=1.0e-8, atol=0.0
            )

    def test_zero_nominal_or_res_one_collapses_to_single_point(self) -> None:
        for dist in (
            gaussian.generate(0.0, 0.1, res=15),
            gaussian.generate(20.0, 0.1, res=1),
            skewed.generate(0.0, 0.1, -0.5, res=15, center="mode"),
            skewed.generate(20.0, 0.1, -0.5, res=1, center="mean"),
        ):
            assert dist.values.size == 1
            np.testing.assert_allclose(dist.weights, [1.0])

    @pytest.mark.parametrize(
        "gen_fn",
        [
            lambda v, s: gaussian.generate(v, s, 11),
            lambda v, s: skewed.generate(v, s, 0.2, 11, center="mode"),
        ],
    )
    def test_negative_or_non_finite_nominal_or_scale_rejected(self, gen_fn) -> None:
        for v, s in [
            (-1.0, 0.1),
            (10.0, -0.1),
            (float("nan"), 0.1),
            (float("inf"), 0.1),
            (10.0, float("nan")),
        ]:
            with pytest.raises(ValueError, match="non-negative"):
                gen_fn(v, s)

    @pytest.mark.parametrize(
        "invalid_skewness",
        [-2.0, -1.0, -MAX_SKEWNESS, MAX_SKEWNESS, 1.0, 2.0, float("nan")],
    )
    def test_skewed_rejects_out_of_domain_skewness(
        self, invalid_skewness: float
    ) -> None:
        with pytest.raises(ValueError, match="skew-normal range"):
            skewed.generate(20.0, 0.1, invalid_skewness, 11, center="mode")

    def test_custom_normalizes_and_validates_standalone_calls(self) -> None:
        dist = custom.generate(15.0, [1.1, 0.9, 1.0], [1.0, 2.0, 1.0])
        assert dist.dephasing is False
        np.testing.assert_allclose(dist.values, [16.5, 13.5, 15.0])
        np.testing.assert_allclose(dist.weights, [0.25, 0.5, 0.25])

        for args in [
            (-1.0, [1.0], [1.0]),
            (15.0, [], []),
            (15.0, [0.9, 1.0], [1.0]),
            (15.0, [0.0, 1.0], [0.5, 0.5]),
            (15.0, [-0.1, 1.0], [0.5, 0.5]),
            (15.0, [float("nan")], [1.0]),
            (15.0, [1.0], [0.0]),
            (15.0, [1.0], [float("inf")]),
        ]:
            with pytest.raises(ValueError):
                custom.generate(*args)


class TestPositiveConditionedMomentsAndModeSolver:
    """Compare quadrature moments to scipy.integrate.quad and validate the mode solver."""

    @pytest.mark.parametrize("scale", [0.05, 0.25, 0.6])
    def test_gaussian_moments_match_positive_conditioned_quad(
        self,
        scale: float,
    ) -> None:
        nominal = 20.0
        dist = gaussian.generate(nominal, scale=scale, res=35)
        q_mean = float(np.sum(dist.values * dist.weights))
        q_var = float(np.sum(((dist.values - q_mean) ** 2) * dist.weights))

        sigma = scale * nominal
        surv_zero = float(stats.norm.sf(0.0, loc=nominal, scale=sigma))

        def pdf_pos(b: float) -> float:
            return float(stats.norm.pdf(b, loc=nominal, scale=sigma)) / surv_zero

        b_upper = float(stats.norm.isf(1.0e-12 * surv_zero, loc=nominal, scale=sigma))
        b_lo = float(
            stats.norm.ppf(
                1.0 - surv_zero + 1.0e-5 * surv_zero, loc=nominal, scale=sigma
            )
        )
        b_hi = float(stats.norm.isf(1.0e-5 * surv_zero, loc=nominal, scale=sigma))
        norm_bracket = 1.0 - 2.0e-5

        b_mean, _ = integrate.quad(
            lambda b: b * pdf_pos(b) / norm_bracket,
            b_lo,
            b_hi,
            epsabs=1.0e-12,
            epsrel=1.0e-12,
        )
        b_var, _ = integrate.quad(
            lambda b: ((b - b_mean) ** 2) * pdf_pos(b) / norm_bracket,
            b_lo,
            b_hi,
            epsabs=1.0e-12,
            epsrel=1.0e-12,
        )
        np.testing.assert_allclose(q_mean, b_mean, rtol=1.0e-9, atol=0.0)
        np.testing.assert_allclose(q_var, b_var, rtol=1.0e-8, atol=0.0)

        full_mean, _ = integrate.quad(
            lambda b: b * pdf_pos(b), 0.0, b_upper, epsabs=1.0e-12, epsrel=1.0e-12
        )
        full_var, _ = integrate.quad(
            lambda b: ((b - full_mean) ** 2) * pdf_pos(b),
            0.0,
            b_upper,
            epsabs=1.0e-12,
            epsrel=1.0e-12,
        )
        np.testing.assert_allclose(q_mean, full_mean, rtol=1.0e-4, atol=0.0)
        np.testing.assert_allclose(q_var, full_var, rtol=5.0e-4, atol=0.0)

    @pytest.mark.parametrize("center", ["mode", "mean"])
    @pytest.mark.parametrize("skewness", [-0.85, -0.4, 0.4, 0.85])
    @pytest.mark.parametrize("scale", [0.05, 0.25])
    def test_skewed_moments_match_positive_conditioned_quad(
        self,
        center: Literal["mode", "mean"],
        skewness: float,
        scale: float,
    ) -> None:
        nominal = 25.0
        dist = skewed.generate(
            nominal,
            scale=scale,
            skewness=skewness,
            res=35,
            center=center,
        )
        q_mean = float(np.sum(dist.values * dist.weights))
        q_var = float(np.sum(((dist.values - q_mean) ** 2) * dist.weights))

        alpha, sigma_z = skewness_to_alpha_and_sigma_z(skewness)
        omega = scale * nominal / sigma_z
        delta = alpha / math.sqrt(1.0 + alpha * alpha)
        anchor_z = (
            delta * _SQRT_2_OVER_PI
            if center == "mean"
            else compute_skewnorm_mode_z(alpha)
        )
        xi = nominal - omega * anchor_z

        surv_zero = float(stats.skewnorm.sf(0.0, alpha, loc=xi, scale=omega))

        def pdf_pos(b: float) -> float:
            return float(stats.skewnorm.pdf(b, alpha, loc=xi, scale=omega)) / surv_zero

        b_upper = float(
            stats.skewnorm.isf(1.0e-12 * surv_zero, alpha, loc=xi, scale=omega)
        )
        b_lo = float(
            stats.skewnorm.ppf(
                1.0 - surv_zero + 1.0e-5 * surv_zero, alpha, loc=xi, scale=omega
            )
        )
        b_hi = float(stats.skewnorm.isf(1.0e-5 * surv_zero, alpha, loc=xi, scale=omega))
        norm_bracket = 1.0 - 2.0e-5

        b_mean, _ = integrate.quad(
            lambda b: b * pdf_pos(b) / norm_bracket,
            b_lo,
            b_hi,
            epsabs=1.0e-11,
            epsrel=1.0e-11,
        )
        b_var, _ = integrate.quad(
            lambda b: ((b - b_mean) ** 2) * pdf_pos(b) / norm_bracket,
            b_lo,
            b_hi,
            epsabs=1.0e-11,
            epsrel=1.0e-11,
        )
        np.testing.assert_allclose(q_mean, b_mean, rtol=1.0e-8, atol=0.0)
        np.testing.assert_allclose(q_var, b_var, rtol=1.0e-7, atol=0.0)

        full_mean, _ = integrate.quad(
            lambda b: b * pdf_pos(b), 0.0, b_upper, epsabs=1.0e-11, epsrel=1.0e-11
        )
        full_var, _ = integrate.quad(
            lambda b: ((b - full_mean) ** 2) * pdf_pos(b),
            0.0,
            b_upper,
            epsabs=1.0e-11,
            epsrel=1.0e-11,
        )
        np.testing.assert_allclose(q_mean, full_mean, rtol=1.0e-4, atol=0.0)
        np.testing.assert_allclose(q_var, full_var, rtol=5.0e-4, atol=0.0)

    @pytest.mark.parametrize(
        "alpha",
        [
            1.0e-10,
            1.0e-6,
            1.0e-3,
            0.05,
            0.166992,
            0.5,
            1.814141,
            2.0,
            10.0,
            50.0,
            200.0,
        ],
    )
    def test_skewnorm_mode_solver_machine_precision_and_antisymmetry(
        self, alpha: float
    ) -> None:
        z_mode = compute_skewnorm_mode_z(alpha)
        az = alpha * z_mode
        ratio = math.exp(
            -0.5 * az * az - _HALF_LOG_2PI - float(log_ndtr(np.float64(az)))
        )
        residual = -z_mode + alpha * ratio
        assert abs(residual) < 1.0e-14
        assert compute_skewnorm_mode_z(-alpha) == -z_mode

    @pytest.mark.parametrize("center", ["mode", "mean"])
    @pytest.mark.parametrize("scale", [0.05, 0.15, 0.5])
    def test_skewed_zero_skewness_reduces_identically_to_gaussian(
        self,
        center: Literal["mode", "mean"],
        scale: float,
    ) -> None:
        dist_s = skewed.generate(
            18.5,
            scale,
            0.0,
            res=15,
            center=center,
        )
        dist_g = gaussian.generate(18.5, scale, res=15)
        np.testing.assert_allclose(dist_s.values, dist_g.values, rtol=0.0, atol=1.0e-15)
        np.testing.assert_allclose(
            dist_s.weights, dist_g.weights, rtol=0.0, atol=1.0e-15
        )


class TestFiniteDifferenceSmoothnessForFitting:
    """Verify TRF finite-difference smoothness for future joint B1 fitting."""

    @staticmethod
    def _eval_cest_normalized(skewness: float, center: str) -> np.ndarray:
        ensure_plugins_registered()
        model = ModelSpec.from_name("2st")
        cfg = cest_15n.Cest15NConfig.model_validate(
            {
                "experiment": {
                    "name": "cest_15n",
                    "time_t1": 0.5,
                    "carrier": 118.0,
                    "b1_frq": 25.0,
                    "b1_distribution": {
                        "type": "skewed",
                        "scale": 0.1,
                        "skewness": skewness,
                        "center": center,
                        "res": 25,
                    },
                },
                "conditions": {"h_larmor_frq": 600.0},
                "data": {},
            },
            context={"model": model},
        )
        cfg.model = model
        cfg.experiment.model_name = model.name
        calc = cest_15n.create_profile_calculation(cfg, SpinSystem(name="G23N-HN"))
        calc.spectrometer.update(
            {
                "pa": 0.9,
                "pb": 0.1,
                "kab": 20.0,
                "kba": 180.0,
                "cs_i_a": 118.0,
                "cs_i_b": 114.0,
                "r1_i_a": 1.5,
                "r1_i_b": 1.5,
                "r2_i_a": 2.0,
                "r2_i_b": 2.0,
            }
        )
        offsets = np.array([-2.0e4, -300.0, -240.0, -100.0, -15.0, 0.0, 15.0])
        data = Data(
            exp=np.ones_like(offsets), err=np.ones_like(offsets), metadata=offsets
        )
        intensities = calc.pulse_sequence.calculate(calc.spectrometer, data)
        return intensities[1:] / intensities[0]

    @pytest.mark.parametrize("center", ["mode", "mean"])
    @pytest.mark.parametrize("skewness", [-0.85, -0.4, -0.02, 0.0, 0.02, 0.4, 0.85])
    def test_trf_forward_derivative_matches_central_difference(
        self,
        center: str,
        skewness: float,
    ) -> None:
        h_trf = 1.4901161193847656e-8
        h_central = 1.0e-5
        y0 = self._eval_cest_normalized(skewness, center)
        df_trf = (self._eval_cest_normalized(skewness + h_trf, center) - y0) / h_trf
        df_central = (
            self._eval_cest_normalized(skewness + h_central, center)
            - self._eval_cest_normalized(skewness - h_central, center)
        ) / (2.0 * h_central)
        np.testing.assert_allclose(df_trf, df_central, rtol=2.0e-3, atol=2.0e-6)

    @pytest.mark.parametrize(
        ("center", "min_slope"), [("mode", 1.0e-5), ("mean", 1.0e-7)]
    )
    def test_skewness_parameterization_is_linearly_regular_at_zero(
        self,
        center: str,
        min_slope: float,
    ) -> None:
        d1 = self._eval_cest_normalized(1.0e-3, center) - self._eval_cest_normalized(
            -1.0e-3, center
        )
        d2 = self._eval_cest_normalized(2.0e-3, center) - self._eval_cest_normalized(
            -2.0e-3, center
        )
        max_d1 = float(np.max(np.abs(d1)))
        assert max_d1 > min_slope
        assert float(np.max(np.abs(d2))) / max_d1 == pytest.approx(2.0, rel=2.0e-2)


class TestIndependentReferenceQuadratureAccuracy:
    """Verify automatic resolution against 2,000-node reference across CEST families."""

    _N_PARAMS = {
        "pa": 0.9,
        "pb": 0.1,
        "kab": 20.0,
        "kba": 180.0,
        "cs_i_a": 118.0,
        "cs_i_b": 114.0,
        "r1_i_a": 1.5,
        "r1_i_b": 1.5,
        "r2_i_a": 2.0,
        "r2_i_b": 2.0,
    }
    _HN_PARAMS = {
        "pa": 0.9,
        "pb": 0.1,
        "kab": 20.0,
        "kba": 180.0,
        "cs_i_a": 8.0,
        "cs_i_b": 7.5,
        "cs_s_a": 118.0,
        "cs_s_b": 115.0,
        "j_is_a": -92.0,
        "j_is_b": -92.0,
        "r1_i_a": 1.5,
        "r1_i_b": 1.5,
        "r1_s_a": 1.2,
        "r1_s_b": 1.2,
        "r1a_is_a": 2.0,
        "r1a_is_b": 2.0,
        "r2_i_a": 2.0,
        "r2_i_b": 2.0,
        "r2a_i_a": 2.5,
        "r2a_i_b": 2.5,
        "etaz_i_a": 0.2,
        "etaz_i_b": 0.2,
        "etaxy_i_a": 0.5,
        "etaxy_i_b": 0.5,
    }
    _N_OFFSETS = np.array([-2.0e4, -300.0, -240.0, -100.0, -15.0, 0.0, 15.0])
    _H_OFFSETS = np.array([-2.0e4, -400.0, -300.0, -120.0, -20.0, 0.0, 20.0])

    @pytest.mark.parametrize(
        ("module", "spin_name", "exp_cfg", "params", "offsets", "b1_cfg"),
        [
            (
                cest_15n,
                "G23N-HN",
                {"name": "cest_15n", "time_t1": 0.5, "carrier": 118.0, "b1_frq": 15.0},
                _N_PARAMS,
                _N_OFFSETS,
                {"type": "gaussian", "scale": 0.1},
            ),
            (
                cest_15n,
                "G23N-HN",
                {"name": "cest_15n", "time_t1": 0.5, "carrier": 118.0, "b1_frq": 60.0},
                _N_PARAMS,
                _N_OFFSETS,
                {"type": "gaussian", "scale": 0.1},
            ),
            (
                cest_15n,
                "G23N-HN",
                {"name": "cest_15n", "time_t1": 0.5, "carrier": 118.0, "b1_frq": 150.0},
                _N_PARAMS,
                _N_OFFSETS,
                {"type": "gaussian", "scale": 0.1},
            ),
            (
                cest_15n,
                "G23N-HN",
                {"name": "cest_15n", "time_t1": 0.5, "carrier": 118.0, "b1_frq": 25.0},
                _N_PARAMS,
                _N_OFFSETS,
                {"type": "skewed", "scale": 0.1, "skewness": -0.90, "center": "mode"},
            ),
            (
                cest_15n,
                "G23N-HN",
                {"name": "cest_15n", "time_t1": 0.5, "carrier": 118.0, "b1_frq": 25.0},
                _N_PARAMS,
                _N_OFFSETS,
                {"type": "skewed", "scale": 0.1, "skewness": -0.99, "center": "mean"},
            ),
            (
                cest_1hn_ip_ap,
                "G23HN-N",
                {
                    "name": "cest_1hn_ip_ap",
                    "time_t1": 0.4,
                    "carrier": 8.0,
                    "b1_frq": 25.0,
                    "d1": 0.6,
                    "eta_block": 4,
                },
                _HN_PARAMS,
                _H_OFFSETS,
                {"type": "skewed", "scale": 0.1, "skewness": -0.90, "center": "mode"},
            ),
            (
                coscest_1hn_ip_ap,
                "G23HN-N",
                {
                    "name": "coscest_1hn_ip_ap",
                    "time_t1": 0.4,
                    "carrier": 8.0,
                    "b1_frq": 25.0,
                    "d1": 0.6,
                    "cos_n": 3,
                    "sw": 800.0,
                },
                _HN_PARAMS,
                _H_OFFSETS,
                {"type": "skewed", "scale": 0.1, "skewness": -0.90, "center": "mode"},
            ),
            (
                dcest_15n,
                "G23N-HN",
                {
                    "name": "dcest_15n",
                    "time_t1": 0.3,
                    "carrier": 118.0,
                    "b1_eff": 20.0,
                    "pw90": 40.0e-6,
                    "sw": 700.0,
                },
                _N_PARAMS,
                _N_OFFSETS,
                {"type": "skewed", "scale": 0.1, "skewness": -0.90, "center": "mode"},
            ),
        ],
    )
    def test_auto_res_matches_2000_node_reference_within_1e_4(
        self,
        module: Any,
        spin_name: str,
        exp_cfg: dict[str, object],
        params: dict[str, float],
        offsets: np.ndarray,
        b1_cfg: dict[str, Any],
    ) -> None:
        ensure_plugins_registered()
        model = ModelSpec.from_name("2st")
        config = module.EXPERIMENT_TYPE.config_type.model_validate(
            {
                "experiment": {**exp_cfg, "b1_distribution": b1_cfg},
                "conditions": {"h_larmor_frq": 600.0},
                "data": {},
            },
            context={"model": model},
        )
        config.model = model
        config.experiment.model_name = model.name

        calc_auto = module.EXPERIMENT_TYPE.create_profile_calculation(
            config, SpinSystem(name=spin_name)
        )
        calc_auto.spectrometer.update(params)
        data = Data(
            exp=np.ones_like(offsets), err=np.ones_like(offsets), metadata=offsets
        )
        y_auto = calc_auto.pulse_sequence.calculate(calc_auto.spectrometer, data)
        ratio_auto = y_auto[1:] / y_auto[0]

        nominal = config.experiment.get_b1_nominal()
        dist_ref = _reference_cest_distribution(
            nominal,
            float(b1_cfg["scale"]),
            float(b1_cfg.get("skewness", 0.0)),
            str(b1_cfg.get("center", "mean")),
        )
        calc_ref = module.EXPERIMENT_TYPE.create_profile_calculation(
            config, SpinSystem(name=spin_name)
        )
        calc_ref.spectrometer.set_b1_i_distribution(dist_ref, nominal=nominal)
        calc_ref.spectrometer.update(params)
        y_ref = calc_ref.pulse_sequence.calculate(calc_ref.spectrometer, data)
        ratio_ref = y_ref[1:] / y_ref[0]

        assert float(np.max(np.abs(ratio_auto - ratio_ref))) <= 1.0e-4

"""Tests for B1 field inhomogeneity distributions."""

from __future__ import annotations

import numpy as np
import pytest
from pydantic import ValidationError
from scipy import stats

from chemex.configuration.b1_config import parse_distribution_config
from chemex.nmr.distributions import (
    gaussian,
    hermite,
    registry,
    skewed,
    truncated_skewed,
)


class TestDistributionRegistry:
    """Test the distribution registry and plugin system."""

    def test_registry_has_distributions(self):
        """Test that distributions are auto-registered."""
        generators = registry.get_all_generators()
        assert len(generators) > 0
        assert "gaussian" in generators
        assert "custom" in generators

    def test_registry_has_config_classes(self):
        """Test that config classes are registered."""
        configs = registry.get_all_config_classes()
        assert len(configs) > 0

    def test_get_gaussian_config_class(self):
        """Test retrieving a specific config class."""
        gaussian_config = registry.get_config_class("gaussian")
        assert gaussian_config is not None
        instance = gaussian_config()
        assert instance.type == "gaussian"

    def test_get_nonexistent_generator_fails(self):
        """Test that requesting an unknown distribution raises an error."""
        with pytest.raises(ValueError, match="Unknown B1 distribution type"):
            registry.get_generator("nonexistent")


class TestDistributionConfigParsing:
    """Test parsing distribution configs via the registered plugin union."""

    def test_flat_schema_gaussian(self):
        """Test flat TOML schema for Gaussian distribution."""
        config = parse_distribution_config(
            {
                "type": "gaussian",
                "scale": 0.1,
                "res": 11,
            }
        )
        assert config.type == "gaussian"
        assert config.scale == 0.1
        assert config.res == 11

    def test_flat_schema_skewed(self):
        """Test flat TOML schema for skewed distribution."""
        config = parse_distribution_config(
            {
                "type": "skewed",
                "scale": 0.15,
                "skewness": -0.5,
                "res": 11,
            }
        )
        assert config.type == "skewed"
        assert config.skewness == -0.5

    def test_flat_schema_custom(self):
        """Test flat TOML schema for custom distribution."""
        config = parse_distribution_config(
            {
                "type": "custom",
                "scales": [0.95, 1.0, 1.03],
                "weights": [0.25, 0.5, 0.25],
            }
        )
        assert config.type == "custom"
        assert config.scales == [0.95, 1.0, 1.03]
        assert config.weights == [0.25, 0.5, 0.25]

    def test_unknown_distribution_option_rejected(self):
        """Test that mistyped distribution options are rejected."""
        with pytest.raises(ValidationError, match="sclae"):
            parse_distribution_config(
                {
                    "type": "gaussian",
                    "sclae": 0.2,
                }
            )


class TestGaussianDistribution:
    """Test Gaussian distribution generation."""

    def test_gaussian_generates_correct_shape(self):
        """Test that Gaussian distribution has correct number of points."""
        config = parse_distribution_config(
            {
                "type": "gaussian",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution(15.0)
        assert len(dist.values) == 11
        assert len(dist.weights) == 11

    def test_gaussian_weights_normalized(self):
        """Test that Gaussian weights sum to 1.0."""
        config = parse_distribution_config(
            {
                "type": "gaussian",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution(15.0)
        assert np.isclose(dist.weights.sum(), 1.0)

    def test_gaussian_centered_at_nominal(self):
        """Test that Gaussian is approximately centered at nominal value."""
        config = parse_distribution_config(
            {
                "type": "gaussian",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution(15.0)
        # Mean should be close to nominal value
        weighted_mean = np.sum(dist.values * dist.weights)
        assert np.isclose(weighted_mean, 15.0, rtol=0.01)

    def test_gaussian_scale_zero_gives_single_point(self):
        """Test that scale=0 gives a single point distribution."""
        config = parse_distribution_config(
            {
                "type": "gaussian",
                "scale": 0.0,
                "res": 11,
            }
        )
        dist = config.get_distribution(15.0)
        assert len(dist.values) == 1
        assert dist.values[0] == 15.0
        assert dist.weights[0] == 1.0


class TestHermiteDistribution:
    """Test Hermite quadrature distribution."""

    def test_hermite_generates_correct_shape(self):
        """Test that Hermite distribution has correct number of points."""
        config = parse_distribution_config(
            {
                "type": "hermite",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution(15.0)
        assert len(dist.values) == 11
        assert len(dist.weights) == 11

    def test_hermite_weights_normalized(self):
        """Test that Hermite weights sum to 1.0."""
        config = parse_distribution_config(
            {
                "type": "hermite",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution(15.0)
        assert np.isclose(dist.weights.sum(), 1.0)


class TestSkewedDistribution:
    """Test skewed distribution generation."""

    def test_skewed_generates_correct_shape(self):
        """Test that skewed distribution has correct number of points."""
        config = parse_distribution_config(
            {
                "type": "skewed",
                "scale": 0.15,
                "skewness": -0.5,
                "res": 11,
            }
        )
        dist = config.get_distribution(15.95)
        assert len(dist.values) == 11
        assert len(dist.weights) == 11

    def test_skewed_weights_normalized(self):
        """Test that skewed weights sum to 1.0."""
        config = parse_distribution_config(
            {
                "type": "skewed",
                "scale": 0.15,
                "skewness": -0.5,
                "res": 11,
            }
        )
        dist = config.get_distribution(15.95)
        assert np.isclose(dist.weights.sum(), 1.0)

    def test_skewed_zero_skewness_uses_hermite(self):
        """Test that zero skewness falls back to Hermite."""
        config_skewed = parse_distribution_config(
            {
                "type": "skewed",
                "scale": 0.1,
                "skewness": 0.0,
                "res": 11,
            }
        )
        config_hermite = parse_distribution_config(
            {
                "type": "hermite",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist_skewed = config_skewed.get_distribution(15.0)
        dist_hermite = config_hermite.get_distribution(15.0)
        # Should be identical
        assert np.allclose(dist_skewed.values, dist_hermite.values)
        assert np.allclose(dist_skewed.weights, dist_hermite.weights)


class TestCustomDistribution:
    """Test user-defined custom distribution."""

    def test_custom_basic(self):
        """Test basic custom distribution creation."""
        config = parse_distribution_config(
            {
                "type": "custom",
                "scales": [0.95, 1.0, 1.03],
                "weights": [0.25, 0.5, 0.25],
            }
        )
        dist = config.get_distribution(15.95)
        assert len(dist.values) == 3
        assert len(dist.weights) == 3

    def test_custom_correct_values(self):
        """Test that custom distribution computes correct B1 values."""
        config = parse_distribution_config(
            {
                "type": "custom",
                "scales": [0.9, 1.0, 1.1],
                "weights": [1.0, 1.0, 1.0],
            }
        )
        dist = config.get_distribution(10.0)
        expected_values = np.array([9.0, 10.0, 11.0])
        assert np.allclose(dist.values, expected_values)

    def test_custom_weights_normalized(self):
        """Test that custom weights are auto-normalized."""
        config = parse_distribution_config(
            {
                "type": "custom",
                "scales": [0.9, 1.0, 1.1],
                "weights": [1.0, 2.0, 1.0],  # Sum = 4.0
            }
        )
        dist = config.get_distribution(15.0)
        assert np.isclose(dist.weights.sum(), 1.0)
        expected_weights = np.array([0.25, 0.5, 0.25])
        assert np.allclose(dist.weights, expected_weights)

    def test_custom_single_point(self):
        """Test custom distribution with single point."""
        config = parse_distribution_config(
            {
                "type": "custom",
                "scales": [1.0],
                "weights": [1.0],
            }
        )
        dist = config.get_distribution(15.0)
        assert len(dist.values) == 1
        assert dist.values[0] == 15.0
        assert dist.weights[0] == 1.0

    def test_custom_mismatched_lengths_rejected(self):
        """Test that mismatched scales and weights are rejected."""
        with pytest.raises(ValidationError, match="same length"):
            parse_distribution_config(
                {
                    "type": "custom",
                    "scales": [0.95, 1.0, 1.05],
                    "weights": [0.5, 0.5],  # Wrong length!
                }
            )

    def test_custom_negative_scale_rejected(self):
        """Test that negative scales are rejected."""
        with pytest.raises(ValidationError, match="positive"):
            parse_distribution_config(
                {
                    "type": "custom",
                    "scales": [-0.95, 1.0],
                    "weights": [0.5, 0.5],
                }
            )

    def test_custom_negative_weight_rejected(self):
        """Test that negative weights are rejected."""
        with pytest.raises(ValidationError, match="positive"):
            parse_distribution_config(
                {
                    "type": "custom",
                    "scales": [0.95, 1.0],
                    "weights": [-0.5, 1.5],
                }
            )


class TestOtherDistributions:
    """Test other distribution types."""

    @pytest.mark.parametrize("dist_type", ["beta"])
    def test_distribution_generates(self, dist_type):
        """Test that distribution generates and has normalized weights."""
        config = parse_distribution_config(
            {
                "type": dist_type,
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution(15.0)
        assert len(dist.values) > 0
        assert len(dist.weights) > 0
        assert np.isclose(dist.weights.sum(), 1.0)


class TestSkewedAnalyticalInvariants:
    """Scientific reference tests for Azzalini skew-normal B1 distributions."""

    @pytest.mark.parametrize(
        "target_skewness",
        [-0.98, -0.8, -0.5, -0.1, -0.005, 0.005, 0.1, 0.5, 0.8, 0.98],
    )
    def test_skewed_analytical_inversion_reproduces_azzalini_moments(
        self,
        target_skewness: float,
    ) -> None:
        """Analytical (alpha, omega) inversion matches scipy.stats.skewnorm and quadrature moments."""
        nominal = 25.0
        scale = 0.04
        alpha, omega = skewed._skewness_to_alpha_and_omega(target_skewness, scale)
        ref_var, ref_skew = stats.skewnorm.stats(alpha, scale=omega, moments="vs")
        np.testing.assert_allclose(float(ref_var), scale**2, rtol=1e-12, atol=0.0)
        np.testing.assert_allclose(
            float(ref_skew),
            target_skewness,
            rtol=1e-12,
            atol=1e-14,
        )

        dist = skewed.generate(nominal, scale, target_skewness, res=31)
        rel_values = dist.values / nominal
        mean = float(np.sum(dist.weights * rel_values))
        centered = rel_values - mean
        std = float(np.sqrt(np.sum(dist.weights * centered**2)))
        empirical_skew = float(np.sum(dist.weights * centered**3) / std**3)

        np.testing.assert_allclose(std, scale, rtol=5e-3)
        np.testing.assert_allclose(empirical_skew, target_skewness, atol=1.5e-2)

    def test_skewed_smooth_across_small_nonzero_skewness(self) -> None:
        """Quadrature nodes and weights vary smoothly across |skewness| = 0.01."""
        d_below = skewed.generate(20.0, 0.1, 0.008, res=15)
        d_at = skewed.generate(20.0, 0.1, 0.010, res=15)
        d_above = skewed.generate(20.0, 0.1, 0.012, res=15)

        assert np.max(np.abs(d_at.values - d_below.values)) < 2e-2
        assert np.max(np.abs(d_above.values - d_at.values)) < 2e-2
        assert (
            np.max(np.abs(d_above.values - 2.0 * d_at.values + d_below.values)) < 5e-3
        )
        assert np.max(np.abs(d_at.weights - d_below.weights)) < 1e-3
        assert np.max(np.abs(d_above.weights - d_at.weights)) < 1e-3

    @pytest.mark.parametrize(
        "invalid_skewness",
        [
            -2.0,
            -1.0,
            -skewed.MAX_SKEWNESS,
            skewed.MAX_SKEWNESS,
            1.0,
            2.0,
        ],
    )
    def test_skewed_rejects_out_of_domain_skewness(
        self,
        invalid_skewness: float,
    ) -> None:
        """Skewness outside (-MAX_SKEWNESS, MAX_SKEWNESS) is rejected."""
        with pytest.raises(ValidationError, match="skew-normal range"):
            parse_distribution_config(
                {
                    "type": "skewed",
                    "scale": 0.1,
                    "skewness": invalid_skewness,
                    "res": 11,
                }
            )
        with pytest.raises(ValueError, match="skew-normal range"):
            skewed.generate(20.0, 0.1, invalid_skewness, 11)


class TestTruncatedSkewedDistribution:
    """Tests for truncated skew-normal B1 distribution."""

    def test_dead_broken_cdf_helpers_are_absent(self) -> None:
        assert not hasattr(truncated_skewed, "_skewnorm_cdf")
        assert not hasattr(truncated_skewed, "_owens_t_approx")

    @pytest.mark.parametrize("alpha", [-8.0, -0.8, 0.0, 5.0])
    @pytest.mark.parametrize("scale", [0.1, 0.5, 1.0])
    def test_truncated_skewed_respects_bounds_and_slant_parameter(
        self,
        alpha: float,
        scale: float,
    ) -> None:
        nominal = 20.0
        config = parse_distribution_config(
            {
                "type": "truncated_skewed",
                "scale": scale,
                "skewness": alpha,
                "res": 15,
            }
        )
        dist = config.get_distribution(nominal)
        assert len(dist.values) == 15
        assert np.all(dist.values >= 0.0)
        assert np.all(dist.values <= nominal)
        assert np.all(dist.weights > 0.0)
        np.testing.assert_allclose(dist.weights.sum(), 1.0, rtol=0.0, atol=1e-13)


class TestNonNegativeB1Invariants:
    """Ensure all B1 distributions produce non-negative nodes across allowed scales."""

    @pytest.mark.parametrize(
        ("dist_type", "extra"),
        [
            ("gaussian", {}),
            ("hermite", {}),
            ("skewed", {"skewness": -0.7}),
            ("skewed", {"skewness": 0.7}),
            ("truncated_skewed", {"skewness": -2.0}),
        ],
    )
    @pytest.mark.parametrize("scale", [0.05, 0.25, 0.6, 1.0])
    @pytest.mark.parametrize("res", [5, 11, 25])
    def test_distribution_nodes_are_non_negative_and_ordered(
        self,
        dist_type: str,
        extra: dict[str, float],
        scale: float,
        res: int,
    ) -> None:
        config = parse_distribution_config(
            {
                "type": dist_type,
                "scale": scale,
                "res": res,
                **extra,
            }
        )
        dist = config.get_distribution(25.0)
        assert len(dist.values) == res
        assert np.all(dist.values >= 0.0)
        assert np.all(np.diff(dist.values) > 0.0)
        assert np.all(dist.weights >= 0.0)
        np.testing.assert_allclose(dist.weights.sum(), 1.0, rtol=0.0, atol=1e-13)

    @pytest.mark.parametrize(
        "gen_fn",
        [
            lambda v, s: gaussian.generate(v, s, 11),
            lambda v, s: hermite.generate(v, s, 11),
            lambda v, s: skewed.generate(v, s, 0.2, 11),
            lambda v, s: truncated_skewed.generate(v, s, -0.8, 11),
        ],
    )
    def test_negative_nominal_or_scale_rejected(self, gen_fn) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            gen_fn(-1.0, 0.1)
        with pytest.raises(ValueError, match="non-negative"):
            gen_fn(10.0, -0.1)

"""Tests for B1 field inhomogeneity distributions."""

from __future__ import annotations

import numpy as np
import pytest
from pydantic import ValidationError

from chemex.configuration.b1_config import parse_distribution_config
from chemex.nmr.distributions import registry


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

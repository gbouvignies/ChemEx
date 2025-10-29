"""Tests for B1 field inhomogeneity distributions."""

from __future__ import annotations

import numpy as np
import pytest
from pydantic import ValidationError

from chemex.configuration.b1_config import B1FieldConfig
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


class TestB1FieldConfig:
    """Test the B1FieldConfig model and schema validation."""

    def test_default_uses_gaussian(self):
        """Test that default distribution is Gaussian."""
        config = B1FieldConfig.model_validate({"value": 15.0})
        assert config.distribution.type == "gaussian"

    def test_flat_schema_gaussian(self):
        """Test flat TOML schema for Gaussian distribution."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "gaussian",
                "scale": 0.1,
                "res": 11,
            }
        )
        assert config.value == 15.0
        assert config.distribution.type == "gaussian"
        assert config.distribution.scale == 0.1
        assert config.distribution.res == 11

    def test_flat_schema_skewed(self):
        """Test flat TOML schema for skewed distribution."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.95,
                "type": "skewed",
                "scale": 0.15,
                "skewness": -0.5,
                "res": 11,
            }
        )
        assert config.distribution.type == "skewed"
        assert config.distribution.skewness == -0.5

    def test_flat_schema_custom(self):
        """Test flat TOML schema for custom distribution."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.95,
                "type": "custom",
                "scales": [0.95, 1.0, 1.03],
                "weights": [0.25, 0.5, 0.25],
            }
        )
        assert config.distribution.type == "custom"
        assert config.distribution.scales == [0.95, 1.0, 1.03]
        assert config.distribution.weights == [0.25, 0.5, 0.25]

    def test_nested_schema_rejected(self):
        """Test that nested schema is rejected."""
        with pytest.raises(ValueError, match="Nested distribution schema"):
            B1FieldConfig.model_validate(
                {
                    "value": 15.0,
                    "distribution": {
                        "type": "gaussian",
                        "scale": 0.1,
                        "res": 11,
                    },
                }
            )

    def test_invalid_value_rejected(self):
        """Test that invalid nominal values are rejected."""
        with pytest.raises(ValidationError):
            B1FieldConfig.model_validate({"value": -5.0})
        with pytest.raises(ValidationError):
            B1FieldConfig.model_validate({"value": 0.0})


class TestGaussianDistribution:
    """Test Gaussian distribution generation."""

    def test_gaussian_generates_correct_shape(self):
        """Test that Gaussian distribution has correct number of points."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "gaussian",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution()
        assert len(dist.values) == 11
        assert len(dist.weights) == 11

    def test_gaussian_weights_normalized(self):
        """Test that Gaussian weights sum to 1.0."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "gaussian",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution()
        assert np.isclose(dist.weights.sum(), 1.0)

    def test_gaussian_centered_at_nominal(self):
        """Test that Gaussian is approximately centered at nominal value."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "gaussian",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution()
        # Mean should be close to nominal value
        weighted_mean = np.sum(dist.values * dist.weights)
        assert np.isclose(weighted_mean, 15.0, rtol=0.01)

    def test_gaussian_scale_zero_gives_single_point(self):
        """Test that scale=0 gives a single point distribution."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "gaussian",
                "scale": 0.0,
                "res": 11,
            }
        )
        dist = config.get_distribution()
        assert len(dist.values) == 1
        assert dist.values[0] == 15.0
        assert dist.weights[0] == 1.0


class TestHermiteDistribution:
    """Test Hermite quadrature distribution."""

    def test_hermite_generates_correct_shape(self):
        """Test that Hermite distribution has correct number of points."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "hermite",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution()
        assert len(dist.values) == 11
        assert len(dist.weights) == 11

    def test_hermite_weights_normalized(self):
        """Test that Hermite weights sum to 1.0."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "hermite",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution()
        assert np.isclose(dist.weights.sum(), 1.0)


class TestSkewedDistribution:
    """Test skewed distribution generation."""

    def test_skewed_generates_correct_shape(self):
        """Test that skewed distribution has correct number of points."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.95,
                "type": "skewed",
                "scale": 0.15,
                "skewness": -0.5,
                "res": 11,
            }
        )
        dist = config.get_distribution()
        assert len(dist.values) == 11
        assert len(dist.weights) == 11

    def test_skewed_weights_normalized(self):
        """Test that skewed weights sum to 1.0."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.95,
                "type": "skewed",
                "scale": 0.15,
                "skewness": -0.5,
                "res": 11,
            }
        )
        dist = config.get_distribution()
        assert np.isclose(dist.weights.sum(), 1.0)

    def test_skewed_zero_skewness_uses_hermite(self):
        """Test that zero skewness falls back to Hermite."""
        config_skewed = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "skewed",
                "scale": 0.1,
                "skewness": 0.0,
                "res": 11,
            }
        )
        config_hermite = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "hermite",
                "scale": 0.1,
                "res": 11,
            }
        )
        dist_skewed = config_skewed.get_distribution()
        dist_hermite = config_hermite.get_distribution()
        # Should be identical
        assert np.allclose(dist_skewed.values, dist_hermite.values)
        assert np.allclose(dist_skewed.weights, dist_hermite.weights)


class TestCustomDistribution:
    """Test user-defined custom distribution."""

    def test_custom_basic(self):
        """Test basic custom distribution creation."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.95,
                "type": "custom",
                "scales": [0.95, 1.0, 1.03],
                "weights": [0.25, 0.5, 0.25],
            }
        )
        dist = config.get_distribution()
        assert len(dist.values) == 3
        assert len(dist.weights) == 3

    def test_custom_correct_values(self):
        """Test that custom distribution computes correct B1 values."""
        config = B1FieldConfig.model_validate(
            {
                "value": 10.0,
                "type": "custom",
                "scales": [0.9, 1.0, 1.1],
                "weights": [1.0, 1.0, 1.0],
            }
        )
        dist = config.get_distribution()
        expected_values = np.array([9.0, 10.0, 11.0])
        assert np.allclose(dist.values, expected_values)

    def test_custom_weights_normalized(self):
        """Test that custom weights are auto-normalized."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "custom",
                "scales": [0.9, 1.0, 1.1],
                "weights": [1.0, 2.0, 1.0],  # Sum = 4.0
            }
        )
        dist = config.get_distribution()
        assert np.isclose(dist.weights.sum(), 1.0)
        expected_weights = np.array([0.25, 0.5, 0.25])
        assert np.allclose(dist.weights, expected_weights)

    def test_custom_single_point(self):
        """Test custom distribution with single point."""
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": "custom",
                "scales": [1.0],
                "weights": [1.0],
            }
        )
        dist = config.get_distribution()
        assert len(dist.values) == 1
        assert dist.values[0] == 15.0
        assert dist.weights[0] == 1.0

    def test_custom_mismatched_lengths_rejected(self):
        """Test that mismatched scales and weights are rejected."""
        with pytest.raises(ValidationError, match="same length"):
            B1FieldConfig.model_validate(
                {
                    "value": 15.0,
                    "type": "custom",
                    "scales": [0.95, 1.0, 1.05],
                    "weights": [0.5, 0.5],  # Wrong length!
                }
            )

    def test_custom_negative_scale_rejected(self):
        """Test that negative scales are rejected."""
        with pytest.raises(ValidationError, match="positive"):
            B1FieldConfig.model_validate(
                {
                    "value": 15.0,
                    "type": "custom",
                    "scales": [-0.95, 1.0],
                    "weights": [0.5, 0.5],
                }
            )

    def test_custom_negative_weight_rejected(self):
        """Test that negative weights are rejected."""
        with pytest.raises(ValidationError, match="positive"):
            B1FieldConfig.model_validate(
                {
                    "value": 15.0,
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
        config = B1FieldConfig.model_validate(
            {
                "value": 15.0,
                "type": dist_type,
                "scale": 0.1,
                "res": 11,
            }
        )
        dist = config.get_distribution()
        assert len(dist.values) > 0
        assert len(dist.weights) > 0
        assert np.isclose(dist.weights.sum(), 1.0)


class TestBackwardCompatibility:
    """Test backward compatibility with legacy float format."""

    def test_float_value_creates_gaussian_default(self):
        """Test that a simple float value defaults to Gaussian."""
        config = B1FieldConfig.model_validate({"value": 13.0})
        assert config.value == 13.0
        assert config.distribution.type == "gaussian"
        # Should use default scale and res
        assert config.distribution.scale == 0.1
        assert config.distribution.res == 11

    def test_get_distribution_from_default(self):
        """Test getting distribution from default config."""
        config = B1FieldConfig.model_validate({"value": 13.0})
        dist = config.get_distribution()
        assert len(dist.values) == 11
        assert np.isclose(dist.weights.sum(), 1.0)

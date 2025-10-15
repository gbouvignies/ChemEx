"""Tests for 3-state Eyring kinetic model."""

import numpy as np
import pytest
from scipy import constants

from chemex.configuration.conditions import Conditions
from chemex.models.kinetic.settings_3st_eyring import (
    MAX_RATE_CONSTANT,
    calculate_kij_3st_eyring,
    make_settings_3st_eyring,
)


class TestCalculateKij3stEyring:
    """Test suite for calculate_kij_3st_eyring function."""

    @pytest.fixture
    def default_params(self):
        """Default thermodynamic parameters for testing."""
        return {
            "dh_b": 8000.0,
            "ds_b": 0.0,
            "dh_c": 12000.0,
            "ds_c": 0.0,
            "dh_ab": 75000.0,
            "ds_ab": 0.0,
            "dh_ac": 80000.0,
            "ds_ac": 0.0,
            "dh_bc": 70000.0,
            "ds_bc": 0.0,
            "temperature": 25.0,
        }

    def test_basic_functionality(self, default_params):
        """Test basic functionality of rate constant calculation."""
        rates = calculate_kij_3st_eyring(**default_params)

        # Check that all expected rate constants are present
        expected_keys = [f"k{i}{j}" for i in "abc" for j in "abc" if i != j]
        assert set(rates.keys()) == set(expected_keys)
        assert len(rates) == 6  # 3 states with 2 transitions each

        # Check that all rates are positive and finite
        for key, value in rates.items():
            assert value > 0, f"Rate {key} should be positive, got {value}"
            assert np.isfinite(value), f"Rate {key} should be finite, got {value}"
            assert value <= MAX_RATE_CONSTANT, (
                f"Rate {key} exceeds maximum, got {value}"
            )

    def test_temperature_dependence(self, default_params):
        """Test that rate constants change appropriately with temperature."""
        params_25 = default_params.copy()
        params_50 = default_params.copy()
        params_50["temperature"] = 50.0

        rates_25 = calculate_kij_3st_eyring(**params_25)
        rates_50 = calculate_kij_3st_eyring(**params_50)

        # Higher temperature should generally increase rates
        for key in ["kab", "kba", "kbc", "kcb"]:
            assert rates_50[key] > rates_25[key], (
                f"Rate {key} should increase with temperature: "
                f"{rates_25[key]:.2e} at 25°C vs {rates_50[key]:.2e} at 50°C"
            )

    def test_eyring_equation_implementation(self, default_params):
        """Test that the Eyring equation is correctly implemented."""
        rates = calculate_kij_3st_eyring(**default_params)

        # Manually calculate kab to verify implementation
        T = constants.convert_temperature(default_params["temperature"], "C", "K")
        kbt_h = constants.k * T / constants.h
        RT = constants.R * T

        dh_a, ds_a = 0.0, 0.0  # Reference state
        dh_ab = default_params["dh_ab"]
        ds_ab = default_params["ds_ab"]

        ddg_ab = dh_ab - dh_a - T * (ds_ab - ds_a)
        expected_kab = kbt_h * np.exp(-ddg_ab / RT)
        expected_kab = np.clip(expected_kab, 0.0, MAX_RATE_CONSTANT)

        assert np.isclose(rates["kab"], expected_kab, rtol=1e-10), (
            f"Calculated kab {rates['kab']:.2e} doesn't match expected {expected_kab:.2e}"
        )

    def test_detailed_balance(self, default_params):
        """Test that detailed balance is maintained for all pairs."""
        rates = calculate_kij_3st_eyring(**default_params)

        T = constants.convert_temperature(default_params["temperature"], "C", "K")
        RT = constants.R * T

        # Test A <-> B
        dh_b = default_params["dh_b"]
        ds_b = default_params["ds_b"]
        dg_ab = dh_b - T * ds_b
        expected_ratio_ab = np.exp(-dg_ab / RT)
        actual_ratio_ab = rates["kab"] / rates["kba"]
        assert np.isclose(actual_ratio_ab, expected_ratio_ab, rtol=1e-10)

        # Test A <-> C
        dh_c = default_params["dh_c"]
        ds_c = default_params["ds_c"]
        dg_ac = dh_c - T * ds_c
        expected_ratio_ac = np.exp(-dg_ac / RT)
        actual_ratio_ac = rates["kac"] / rates["kca"]
        assert np.isclose(actual_ratio_ac, expected_ratio_ac, rtol=1e-10)

        # Test B <-> C
        dg_bc = (dh_c - dh_b) - T * (ds_c - ds_b)
        expected_ratio_bc = np.exp(-dg_bc / RT)
        actual_ratio_bc = rates["kbc"] / rates["kcb"]
        assert np.isclose(actual_ratio_bc, expected_ratio_bc, rtol=1e-10)

    def test_entropy_effects(self, default_params):
        """Test that entropy parameters affect rate constants correctly."""
        rates_zero_entropy = calculate_kij_3st_eyring(**default_params)

        # Test with positive entropy for state B
        params_with_entropy = default_params.copy()
        params_with_entropy["ds_b"] = 50.0  # J/mol/K
        rates_with_entropy = calculate_kij_3st_eyring(**params_with_entropy)

        # Entropy should affect equilibrium
        ratio_zero = rates_zero_entropy["kab"] / rates_zero_entropy["kba"]
        ratio_entropy = rates_with_entropy["kab"] / rates_with_entropy["kba"]

        assert ratio_zero != ratio_entropy, "Entropy should affect rate ratios"

    def test_numerical_stability(self, default_params):
        """Test numerical stability with extreme values."""
        # Test with very high barriers
        params_high_barrier = default_params.copy()
        params_high_barrier.update(
            {
                "dh_ab": 200000.0,
                "dh_ac": 200000.0,
                "dh_bc": 200000.0,
            }
        )

        rates_high = calculate_kij_3st_eyring(**params_high_barrier)

        for key, value in rates_high.items():
            assert value > 0, f"Rate {key} should be positive with high barriers"
            assert np.isfinite(value), f"Rate {key} should be finite"
            assert value < 1e-10, f"Rate {key} should be very small with high barriers"

        # Test with very low barriers
        params_low_barrier = default_params.copy()
        params_low_barrier.update(
            {
                "dh_ab": 1000.0,
                "dh_ac": 1000.0,
                "dh_bc": 1000.0,
            }
        )

        rates_low = calculate_kij_3st_eyring(**params_low_barrier)

        for key, value in rates_low.items():
            assert value > 0, f"Rate {key} should be positive"
            assert np.isfinite(value), f"Rate {key} should be finite"
            assert value <= MAX_RATE_CONSTANT, (
                f"Rate {key} should be clipped at maximum"
            )

    def test_caching(self, default_params):
        """Test that function caching works correctly."""
        rates1 = calculate_kij_3st_eyring(**default_params)
        rates2 = calculate_kij_3st_eyring(**default_params)

        # Results should be identical (same object due to caching)
        assert rates1 is rates2, "Caching should return the same object"

        # Different parameters should give different objects
        params_diff = default_params.copy()
        params_diff["temperature"] = 30.0
        rates3 = calculate_kij_3st_eyring(**params_diff)

        assert rates1 is not rates3, (
            "Different parameters should give different objects"
        )

    def test_parameter_types(self, default_params):
        """Test that function handles different numeric types correctly."""
        # Test with integer values
        params_int = {
            k: int(v) if k != "temperature" else v for k, v in default_params.items()
        }
        rates_int = calculate_kij_3st_eyring(**params_int)
        assert len(rates_int) == 6

        # Test with numpy types
        params_numpy = {k: np.float64(v) for k, v in default_params.items()}
        rates_numpy = calculate_kij_3st_eyring(**params_numpy)
        assert len(rates_numpy) == 6


class TestMakeSettings3stEyring:
    """Test suite for make_settings_3st_eyring function."""

    def test_basic_settings_creation(self):
        """Test basic settings creation with valid conditions."""
        conditions = Conditions(temperature=25.0)
        settings = make_settings_3st_eyring(conditions)

        # Check that all expected parameters are present
        expected_params = [
            "dh_b",
            "dh_c",
            "ds_b",
            "ds_c",
            "dh_ab",
            "dh_ac",
            "dh_bc",
            "ds_ab",
            "ds_ac",
            "ds_bc",
        ]
        for param in expected_params:
            assert param in settings, f"Parameter {param} missing from settings"

        # Check that rate constants are created
        rate_constants = [f"k{i}{j}" for i in "abc" for j in "abc" if i != j]
        for k in rate_constants:
            assert k in settings, f"Rate constant {k} missing from settings"

        # Check that populations are created
        for pop in ["pa", "pb", "pc"]:
            assert pop in settings, f"Population {pop} missing from settings"

    def test_none_temperature_error(self):
        """Test that None temperature raises appropriate error."""
        conditions = Conditions(temperature=None)

        with pytest.raises(ValueError, match="The 'temperature' is None"):
            make_settings_3st_eyring(conditions)

    def test_default_parameter_values(self):
        """Test that default parameter values are reasonable."""
        conditions = Conditions(temperature=25.0)
        settings = make_settings_3st_eyring(conditions)

        assert settings["dh_b"].value == 8e3, "Default DH_B should be 8000 J/mol"
        assert settings["dh_c"].value == 8e3, "Default DH_C should be 8000 J/mol"
        assert settings["dh_ab"].value == 6.5e4, "Default DH_AB should be 65000 J/mol"
        assert settings["ds_b"].value == 0.0, "Default DS_B should be 0"

        # Check vary flags
        assert settings["dh_b"].vary is True, "DH_B should be variable"
        assert settings["dh_ab"].vary is True, "DH_AB should be variable"
        assert settings["ds_b"].vary is False, "DS_B should be fixed by default"

    def test_temperature_incorporation(self):
        """Test that temperature is correctly incorporated into expressions."""
        temperature = 37.0
        conditions = Conditions(temperature=temperature)
        settings = make_settings_3st_eyring(conditions)

        kab_expr = settings["kab"].expr
        assert str(temperature) in kab_expr, (
            f"Temperature {temperature} should appear in kab expression"
        )

    def test_parameter_bounds(self):
        """Test that parameters have appropriate bounds."""
        conditions = Conditions(temperature=25.0)
        settings = make_settings_3st_eyring(conditions)

        # Rate constants should have minimum of 0
        for k in ["kab", "kba", "kac", "kca", "kbc", "kcb"]:
            assert settings[k].min == 0.0, f"{k} should have min=0.0"

        # Populations should have bounds [0, 1]
        for pop in ["pa", "pb", "pc"]:
            assert settings[pop].min == 0.0, f"{pop} should have min=0.0"
            assert settings[pop].max == 1.0, f"{pop} should have max=1.0"


class TestPhysicalRealism:
    """Test suite for physical realism of the model."""

    def test_arrhenius_behavior(self):
        """Test that the model exhibits Arrhenius behavior."""
        params = {
            "dh_b": 5000.0,
            "ds_b": 0.0,
            "dh_c": 10000.0,
            "ds_c": 0.0,
            "dh_ab": 60000.0,
            "ds_ab": 0.0,
            "dh_ac": 65000.0,
            "ds_ac": 0.0,
            "dh_bc": 55000.0,
            "ds_bc": 0.0,
        }

        temperatures = [10.0, 25.0, 40.0, 55.0]
        rates_kab = []

        for temp in temperatures:
            params["temperature"] = temp
            rates = calculate_kij_3st_eyring(**params)
            rates_kab.append(rates["kab"])

        # Check that rates increase monotonically with temperature
        for i in range(1, len(rates_kab)):
            assert rates_kab[i] > rates_kab[i - 1], (
                f"Rate at {temperatures[i]}°C should be higher than "
                f"at {temperatures[i - 1]}°C"
            )

        # Check Arrhenius plot linearity (ln(k) vs 1/T)
        ln_rates = [np.log(r) for r in rates_kab]
        inv_temps = [
            1 / constants.convert_temperature(t, "C", "K") for t in temperatures
        ]

        # Simple linear regression
        n = len(ln_rates)
        sum_x = sum(inv_temps)
        sum_y = sum(ln_rates)
        sum_xy = sum(x * y for x, y in zip(inv_temps, ln_rates, strict=False))
        sum_x2 = sum(x * x for x in inv_temps)

        # Calculate correlation coefficient
        r_squared = ((n * sum_xy - sum_x * sum_y) ** 2) / (
            (n * sum_x2 - sum_x**2) * (n * sum(y * y for y in ln_rates) - sum_y**2)
        )

        # Should have good linear correlation (R² > 0.99 for Arrhenius behavior)
        assert r_squared > 0.99, f"Arrhenius plot should be linear, R² = {r_squared}"

    def test_cyclic_consistency(self):
        """Test that cyclic paths give consistent thermodynamics.

        For the cycle A -> B -> C -> A, the free energy change should be zero.
        This means: kab/kba * kbc/kcb * kca/kac = 1
        """
        params = {
            "dh_b": 8000.0,
            "ds_b": 10.0,
            "dh_c": 12000.0,
            "ds_c": -5.0,
            "dh_ab": 75000.0,
            "ds_ab": 50.0,
            "dh_ac": 80000.0,
            "ds_ac": 30.0,
            "dh_bc": 70000.0,
            "ds_bc": 20.0,
            "temperature": 25.0,
        }

        rates = calculate_kij_3st_eyring(**params)

        # Calculate the product of equilibrium constants around the cycle
        cycle_product = (
            (rates["kab"] / rates["kba"])
            * (rates["kbc"] / rates["kcb"])
            * (rates["kca"] / rates["kac"])
        )

        # Should be 1.0 (within numerical precision)
        assert np.isclose(cycle_product, 1.0, rtol=1e-10), (
            f"Cyclic thermodynamic consistency violated: "
            f"cycle product = {cycle_product}, expected 1.0"
        )


if __name__ == "__main__":
    pytest.main([__file__])

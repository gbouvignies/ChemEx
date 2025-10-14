"""Tests for 4-state Eyring kinetic model."""

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.kinetic.settings_4st_eyring import (
    CELSIUS_TO_KELVIN,
    MAX_RATE_CONSTANT,
    MAX_TEMPERATURE,
    MIN_TEMPERATURE,
    calculate_kij_4st_eyring,
    make_settings_4st_eyring,
)


class TestCalculateKij4stEyring:
    """Test suite for calculate_kij_4st_eyring function."""

    @pytest.fixture
    def default_params(self):
        """Default thermodynamic parameters for testing."""
        return {
            "dh_b": 8000.0,
            "ds_b": 0.0,
            "dh_c": 12000.0,
            "ds_c": 0.0,
            "dh_d": 15000.0,
            "ds_d": 0.0,
            "dh_ab": 75000.0,
            "ds_ab": 0.0,
            "dh_ac": 80000.0,
            "ds_ac": 0.0,
            "dh_ad": 85000.0,
            "ds_ad": 0.0,
            "dh_bc": 70000.0,
            "ds_bc": 0.0,
            "dh_bd": 77000.0,
            "ds_bd": 0.0,
            "dh_cd": 72000.0,
            "ds_cd": 0.0,
            "temperature": 25.0,
        }

    def test_basic_functionality(self, default_params):
        """Test basic functionality of rate constant calculation."""
        rates = calculate_kij_4st_eyring(**default_params)

        # Check that all expected rate constants are present
        expected_keys = [f"k{i}{j}" for i in "abcd" for j in "abcd" if i != j]
        assert set(rates.keys()) == set(expected_keys)
        assert len(rates) == 12  # 4 states with 3 transitions each

        # Check that all rates are positive and finite
        for key, value in rates.items():
            assert value > 0, f"Rate {key} should be positive, got {value}"
            assert np.isfinite(value), f"Rate {key} should be finite, got {value}"
            assert value <= MAX_RATE_CONSTANT, (
                f"Rate {key} exceeds maximum, got {value}"
            )

    def test_temperature_validation(self, default_params):
        """Test temperature range validation."""
        # Test too low temperature
        with pytest.raises(ValueError, match="Temperature.*outside reasonable range"):
            params = default_params.copy()
            params["temperature"] = MIN_TEMPERATURE - 10
            calculate_kij_4st_eyring(**params)

        # Test too high temperature
        with pytest.raises(ValueError, match="Temperature.*outside reasonable range"):
            params = default_params.copy()
            params["temperature"] = MAX_TEMPERATURE + 10
            calculate_kij_4st_eyring(**params)

        # Test boundary values should work
        params_min = default_params.copy()
        params_min["temperature"] = MIN_TEMPERATURE
        rates_min = calculate_kij_4st_eyring(**params_min)
        assert len(rates_min) == 12

        params_max = default_params.copy()
        params_max["temperature"] = MAX_TEMPERATURE
        rates_max = calculate_kij_4st_eyring(**params_max)
        assert len(rates_max) == 12

    def test_temperature_dependence(self, default_params):
        """Test that rate constants change appropriately with temperature."""
        # Calculate rates at two different temperatures
        params_25 = default_params.copy()
        params_50 = default_params.copy()
        params_50["temperature"] = 50.0

        rates_25 = calculate_kij_4st_eyring(**params_25)
        rates_50 = calculate_kij_4st_eyring(**params_50)

        # Higher temperature should generally increase rates (Arrhenius behavior)
        # Test a few representative transitions
        for key in ["kab", "kba", "kbc", "kcb"]:
            assert rates_50[key] > rates_25[key], (
                f"Rate {key} should increase with temperature: "
                f"{rates_25[key]:.2e} at 25°C vs {rates_50[key]:.2e} at 50°C"
            )

    def test_eyring_equation_implementation(self, default_params):
        """Test that the Eyring equation is correctly implemented."""
        from scipy.constants import R, h, k

        rates = calculate_kij_4st_eyring(**default_params)

        # Manually calculate one rate constant to verify implementation
        T = default_params["temperature"] + CELSIUS_TO_KELVIN
        kbt_h = k * T / h
        RT = R * T

        # Calculate kab: A -> B transition
        dh_a, ds_a = 0.0, 0.0  # Reference state
        dh_ab = default_params["dh_ab"]
        ds_ab = default_params["ds_ab"]

        ddg_ab = dh_ab - dh_a - T * (ds_ab - ds_a)
        expected_kab = kbt_h * np.exp(-ddg_ab / RT)

        assert np.isclose(rates["kab"], expected_kab, rtol=1e-10), (
            f"Calculated kab {rates['kab']:.2e} doesn't match expected {expected_kab:.2e}"
        )

    def test_reverse_transition_consistency(self, default_params):
        """Test that forward and reverse transitions use correct barriers."""
        from scipy.constants import R

        rates = calculate_kij_4st_eyring(**default_params)

        # The ratio of forward to reverse rates should depend on state energy differences
        T = default_params["temperature"] + CELSIUS_TO_KELVIN
        RT = R * T

        # For A <-> B transition
        dh_b = default_params["dh_b"]
        ds_b = default_params["ds_b"]
        dg_ab = dh_b - T * ds_b  # Free energy difference A -> B

        # Detailed balance: kab/kba = exp(-ΔG_AB/RT)
        expected_ratio = np.exp(-dg_ab / RT)
        actual_ratio = rates["kab"] / rates["kba"]

        assert np.isclose(actual_ratio, expected_ratio, rtol=1e-10), (
            f"Detailed balance violated for A<->B: "
            f"expected ratio {expected_ratio:.2e}, got {actual_ratio:.2e}"
        )

    def test_entropy_effects(self, default_params):
        """Test that entropy parameters affect rate constants correctly."""
        # Test with zero entropy (default case)
        rates_zero_entropy = calculate_kij_4st_eyring(**default_params)

        # Test with positive entropy for state B
        params_with_entropy = default_params.copy()
        params_with_entropy["ds_b"] = 50.0  # J/mol/K
        rates_with_entropy = calculate_kij_4st_eyring(**params_with_entropy)

        # Positive entropy should affect the equilibrium
        # The ratio kab/kba should change
        ratio_zero = rates_zero_entropy["kab"] / rates_zero_entropy["kba"]
        ratio_entropy = rates_with_entropy["kab"] / rates_with_entropy["kba"]

        assert ratio_zero != ratio_entropy, "Entropy should affect rate ratios"

    def test_numerical_stability(self, default_params):
        """Test numerical stability with extreme values."""
        # Test with very high barriers (should give very small rates)
        params_high_barrier = default_params.copy()
        params_high_barrier.update(
            {
                "dh_ab": 200000.0,  # Very high barrier
                "dh_ac": 200000.0,
                "dh_ad": 200000.0,
                "dh_bc": 200000.0,
                "dh_bd": 200000.0,
                "dh_cd": 200000.0,
            }
        )

        rates_high = calculate_kij_4st_eyring(**params_high_barrier)

        # All rates should be very small but still positive and finite
        for key, value in rates_high.items():
            assert value > 0, f"Rate {key} should be positive even with high barriers"
            assert np.isfinite(value), f"Rate {key} should be finite"
            assert value < 1e-10, f"Rate {key} should be very small with high barriers"

    def test_caching(self, default_params):
        """Test that function caching works correctly."""
        # Call function twice with same parameters
        rates1 = calculate_kij_4st_eyring(**default_params)
        rates2 = calculate_kij_4st_eyring(**default_params)

        # Results should be identical (same object due to caching)
        assert rates1 is rates2, "Caching should return the same object"

        # Call with slightly different parameters
        params_diff = default_params.copy()
        params_diff["temperature"] = 26.0
        rates3 = calculate_kij_4st_eyring(**params_diff)

        # Should be different object
        assert rates1 is not rates3, (
            "Different parameters should give different objects"
        )

    def test_parameter_types(self, default_params):
        """Test that function handles different numeric types correctly."""
        # Test with integer values
        params_int = {
            k: int(v) if k != "temperature" else v for k, v in default_params.items()
        }
        rates_int = calculate_kij_4st_eyring(**params_int)
        assert len(rates_int) == 12

        # Test with numpy types
        params_numpy = {k: np.float64(v) for k, v in default_params.items()}
        rates_numpy = calculate_kij_4st_eyring(**params_numpy)
        assert len(rates_numpy) == 12


class TestMakeSettings4stEyring:
    """Test suite for make_settings_4st_eyring function."""

    def test_basic_settings_creation(self):
        """Test basic settings creation with valid conditions."""
        conditions = Conditions(temperature=25.0)
        settings = make_settings_4st_eyring(conditions)

        # Check that all expected parameters are present
        expected_params = [
            # State enthalpies and entropies
            "dh_b",
            "dh_c",
            "dh_d",
            "ds_b",
            "ds_c",
            "ds_d",
            # Transition barriers
            "dh_ab",
            "dh_ac",
            "dh_ad",
            "dh_bc",
            "dh_bd",
            "dh_cd",
            "ds_ab",
            "ds_ac",
            "ds_ad",
            "ds_bc",
            "ds_bd",
            "ds_cd",
        ]

        for param in expected_params:
            assert param in settings, f"Parameter {param} missing from settings"

        # Check that rate constants are also created
        rate_constants = [f"k{i}{j}" for i in "abcd" for j in "abcd" if i != j]
        for k in rate_constants:
            assert k in settings, f"Rate constant {k} missing from settings"

        # Check that populations are created
        populations = [f"p{state}" for state in "abcd"]
        for pop in populations:
            assert pop in settings, f"Population {pop} missing from settings"

    def test_none_temperature_error(self):
        """Test that None temperature raises appropriate error."""
        conditions = Conditions(temperature=None)

        with pytest.raises(ValueError, match="The 'temperature' is None"):
            make_settings_4st_eyring(conditions)

    def test_default_parameter_values(self):
        """Test that default parameter values are reasonable."""
        conditions = Conditions(temperature=25.0)
        settings = make_settings_4st_eyring(conditions)

        # Check some default values
        assert settings["dh_b"].value == 8e3, "Default DH_B should be 8000 J/mol"
        assert settings["dh_ab"].value == 7.5e4, "Default DH_AB should be 75000 J/mol"
        assert settings["ds_b"].value == 0.0, "Default DS_B should be 0"

        # Check vary flags
        assert settings["dh_b"].vary is True, "DH_B should be variable by default"
        assert settings["ds_b"].vary is False, "DS_B should be fixed by default"
        assert settings["dh_ab"].vary is True, "DH_AB should be variable by default"

    def test_temperature_incorporation(self):
        """Test that temperature is correctly incorporated into expressions."""
        temperature = 37.0  # Different from default
        conditions = Conditions(temperature=temperature)
        settings = make_settings_4st_eyring(conditions)

        # Check that temperature appears in rate constant expressions
        kab_expr = settings["kab"].expr
        assert str(temperature) in kab_expr, (
            f"Temperature {temperature} should appear in kab expression"
        )


class TestPhysicalRealism:
    """Test suite for physical realism of the model."""

    def test_arrhenius_behavior(self):
        """Test that the model exhibits Arrhenius behavior."""
        params = {
            "dh_b": 5000.0,
            "ds_b": 0.0,
            "dh_c": 10000.0,
            "ds_c": 0.0,
            "dh_d": 8000.0,
            "ds_d": 0.0,
            "dh_ab": 60000.0,
            "ds_ab": 0.0,
            "dh_ac": 65000.0,
            "ds_ac": 0.0,
            "dh_ad": 70000.0,
            "ds_ad": 0.0,
            "dh_bc": 55000.0,
            "ds_bc": 0.0,
            "dh_bd": 62000.0,
            "ds_bd": 0.0,
            "dh_cd": 58000.0,
            "ds_cd": 0.0,
        }

        temperatures = [10.0, 25.0, 40.0, 55.0]
        rates_kab = []

        for temp in temperatures:
            params["temperature"] = temp
            rates = calculate_kij_4st_eyring(**params)
            rates_kab.append(rates["kab"])

        # Check that rates increase monotonically with temperature
        for i in range(1, len(rates_kab)):
            assert rates_kab[i] > rates_kab[i - 1], (
                f"Rate at {temperatures[i]}°C should be higher than at {temperatures[i - 1]}°C"
            )

        # Check Arrhenius plot linearity (ln(k) vs 1/T should be approximately linear)
        ln_rates = [np.log(r) for r in rates_kab]
        inv_temps = [1 / (t + CELSIUS_TO_KELVIN) for t in temperatures]

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

    def test_detailed_balance(self):
        """Test that detailed balance is maintained."""
        params = {
            "dh_b": 8000.0,
            "ds_b": 10.0,
            "dh_c": 12000.0,
            "ds_c": -5.0,
            "dh_d": 15000.0,
            "ds_d": 15.0,
            "dh_ab": 75000.0,
            "ds_ab": 50.0,
            "dh_ac": 80000.0,
            "ds_ac": 30.0,
            "dh_ad": 85000.0,
            "ds_ad": 70.0,
            "dh_bc": 70000.0,
            "ds_bc": 20.0,
            "dh_bd": 77000.0,
            "ds_bd": 40.0,
            "dh_cd": 72000.0,
            "ds_cd": 25.0,
            "temperature": 25.0,
        }

        rates = calculate_kij_4st_eyring(**params)

        from scipy.constants import R

        T = params["temperature"] + CELSIUS_TO_KELVIN
        RT = R * T

        # Test detailed balance for all pairs
        pairs = [("a", "b"), ("a", "c"), ("a", "d"), ("b", "c"), ("b", "d"), ("c", "d")]

        for i, j in pairs:
            # Calculate free energy difference
            if i == "a":
                dh_i, ds_i = 0.0, 0.0
            else:
                dh_i = params[f"dh_{i}"]
                ds_i = params[f"ds_{i}"]

            if j == "a":
                dh_j, ds_j = 0.0, 0.0
            else:
                dh_j = params[f"dh_{j}"]
                ds_j = params[f"ds_{j}"]

            dg_ij = (dh_j - dh_i) - T * (ds_j - ds_i)

            # Detailed balance: kij/kji = exp(-ΔG_ij/RT)
            expected_ratio = np.exp(-dg_ij / RT)
            actual_ratio = rates[f"k{i}{j}"] / rates[f"k{j}{i}"]

            assert np.isclose(actual_ratio, expected_ratio, rtol=1e-10), (
                f"Detailed balance violated for {i}<->{j}: "
                f"expected ratio {expected_ratio:.2e}, got {actual_ratio:.2e}"
            )


if __name__ == "__main__":
    pytest.main([__file__])

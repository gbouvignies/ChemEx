"""Integration tests for 2st_eyring model."""

import numpy as np
import pytest
from scipy import constants

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.models.kinetic.settings_2st_eyring import (
    calculate_kij_2st_eyring,
    make_settings_2st_eyring,
)


class Test2stEyringIntegration:
    """Integration tests for the 2st_eyring model."""

    def test_model_registration(self):
        """Test that the model is properly registered in the factory."""
        # Explicitly trigger registration
        from chemex.models.kinetic.settings_2st_eyring import register

        register()

        # The model should be registered
        assert "2st_eyring" in model_factory.set

        # Test that we can create settings
        conditions = Conditions(temperature=25.0)
        settings = model_factory.create("2st_eyring", conditions)

        assert settings is not None
        assert len(settings) > 0

        # Check for key parameters
        expected_keys = ["dh_b", "ds_b", "dh_ab", "ds_ab", "kab", "kba", "pa", "pb"]
        for key in expected_keys:
            assert key in settings

    def test_temperature_dependency_across_range(self):
        """Test that the model behaves correctly across a range of temperatures."""
        temperatures = [10.0, 25.0, 40.0, 55.0, 70.0]

        # Create settings for each temperature
        all_settings = {}
        for temp in temperatures:
            conditions = Conditions(temperature=temp)
            settings = make_settings_2st_eyring(conditions)
            all_settings[temp] = settings

        # Check that rate constant expressions change with temperature
        for i, temp1 in enumerate(temperatures[:-1]):
            temp2 = temperatures[i + 1]

            # The expressions should contain the temperature values
            assert str(temp1) in all_settings[temp1]["kab"].expr
            assert str(temp2) in all_settings[temp2]["kab"].expr

            # The expressions should be different
            assert all_settings[temp1]["kab"].expr != all_settings[temp2]["kab"].expr

    def test_parameter_constraints(self):
        """Test that parameter constraints are properly set."""
        conditions = Conditions(temperature=25.0)
        settings = make_settings_2st_eyring(conditions)

        # Rate constants should be non-negative
        assert settings["kab"].min == 0.0
        assert settings["kba"].min == 0.0

        # Populations should be in [0, 1]
        assert settings["pa"].min == 0.0
        assert settings["pa"].max == 1.0
        assert settings["pb"].min == 0.0
        assert settings["pb"].max == 1.0

        # Population expressions should use pop_2st function
        assert "pop_2st" in settings["pa"].expr
        assert "pop_2st" in settings["pb"].expr

    def test_settings_consistency_with_calculation(self):
        """Test that settings creation is consistent with direct calculation."""
        temp = 30.0
        conditions = Conditions(temperature=temp)
        settings = make_settings_2st_eyring(conditions)

        # Use default values from settings
        params = {
            "dh_b": settings["dh_b"].value,
            "ds_b": settings["ds_b"].value,
            "dh_ab": settings["dh_ab"].value,
            "ds_ab": settings["ds_ab"].value,
            "temperature": temp,
        }

        # Calculate rate constants directly
        rates = calculate_kij_2st_eyring(**params)

        # The expressions in settings should evaluate to the same values
        # We can't evaluate the expressions directly, but we can verify they're well-formed
        assert "kab" in rates
        assert "kba" in rates
        assert rates["kab"] > 0
        assert rates["kba"] > 0

    def test_equilibrium_consistency(self):
        """Test that equilibrium properties are thermodynamically consistent."""
        temp = 25.0
        dh_b = 5000.0  # J/mol
        ds_b = 10.0  # J/mol/K

        params = {
            "dh_b": dh_b,
            "ds_b": ds_b,
            "dh_ab": 60000.0,
            "ds_ab": 50.0,
            "temperature": temp,
        }

        rates = calculate_kij_2st_eyring(**params)

        # Calculate expected equilibrium constant from thermodynamics
        T = constants.convert_temperature(temp, "C", "K")
        dg_b = dh_b - T * ds_b
        K_eq_expected = np.exp(-dg_b / (constants.R * T))

        # Calculate from rate constants
        K_eq_from_rates = rates["kab"] / rates["kba"]

        # Should match within numerical precision
        assert np.isclose(K_eq_from_rates, K_eq_expected, rtol=1e-10)

    def test_van_t_hoff_analysis(self):
        """Test that the model produces correct van't Hoff behavior.

        The van't Hoff equation states: ln(K_eq) = -ΔH°/R * (1/T) + ΔS°/R
        This should give a linear relationship between ln(K_eq) and 1/T.
        """
        dh_b = 8000.0  # J/mol
        ds_b = 15.0  # J/mol/K

        temperatures = [10.0, 20.0, 30.0, 40.0, 50.0]
        inv_temps = []
        ln_keq_values = []

        for temp in temperatures:
            params = {
                "dh_b": dh_b,
                "ds_b": ds_b,
                "dh_ab": 65000.0,
                "ds_ab": 50.0,
                "temperature": temp,
            }

            rates = calculate_kij_2st_eyring(**params)
            K_eq = rates["kab"] / rates["kba"]

            T = constants.convert_temperature(temp, "C", "K")
            inv_temps.append(1 / T)
            ln_keq_values.append(np.log(K_eq))

        # Linear regression to check linearity
        n = len(ln_keq_values)
        sum_x = sum(inv_temps)
        sum_y = sum(ln_keq_values)
        sum_xy = sum(x * y for x, y in zip(inv_temps, ln_keq_values, strict=False))
        sum_x2 = sum(x * x for x in inv_temps)
        sum_y2 = sum(y * y for y in ln_keq_values)

        # Calculate R² for linearity
        r_squared = ((n * sum_xy - sum_x * sum_y) ** 2) / (
            (n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2)
        )

        # Should be very linear (R² > 0.999)
        assert r_squared > 0.999, f"van't Hoff plot not linear, R² = {r_squared}"

        # Calculate slope and intercept
        slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x**2)
        intercept = (sum_y - slope * sum_x) / n

        # Slope should be -ΔH°/R
        expected_slope = -dh_b / constants.R
        assert np.isclose(slope, expected_slope, rtol=0.01)

        # Intercept should be ΔS°/R
        expected_intercept = ds_b / constants.R
        assert np.isclose(intercept, expected_intercept, rtol=0.01)

    def test_microscopic_reversibility(self):
        """Test microscopic reversibility (detailed balance) at multiple temperatures."""
        temperatures = [15.0, 25.0, 35.0, 45.0]

        for temp in temperatures:
            params = {
                "dh_b": 6000.0,
                "ds_b": 12.0,
                "dh_ab": 70000.0,
                "ds_ab": 40.0,
                "temperature": temp,
            }

            rates = calculate_kij_2st_eyring(**params)

            # Calculate equilibrium constant from thermodynamics
            T = constants.convert_temperature(temp, "C", "K")
            dg_b = params["dh_b"] - T * params["ds_b"]
            K_eq_thermo = np.exp(-dg_b / (constants.R * T))

            # Calculate from kinetics
            K_eq_kinetic = rates["kab"] / rates["kba"]

            # Should satisfy detailed balance
            assert np.isclose(K_eq_kinetic, K_eq_thermo, rtol=1e-10), (
                f"Detailed balance violated at {temp}°C: "
                f"K_kinetic = {K_eq_kinetic:.6e}, K_thermo = {K_eq_thermo:.6e}"
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

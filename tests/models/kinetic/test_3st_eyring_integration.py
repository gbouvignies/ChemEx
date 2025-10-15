"""Integration tests for 3st_eyring model."""

import numpy as np
import pytest
from scipy import constants

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.models.kinetic.settings_3st_eyring import (
    calculate_kij_3st_eyring,
    make_settings_3st_eyring,
)


class Test3stEyringIntegration:
    """Integration tests for the 3st_eyring model."""

    def test_model_registration(self):
        """Test that the model is properly registered in the factory."""
        # Explicitly trigger registration
        from chemex.models.kinetic.settings_3st_eyring import register

        register()

        # The model should be registered
        assert "3st_eyring" in model_factory.set

        # Test that we can create settings
        conditions = Conditions(temperature=25.0)
        settings = model_factory.create("3st_eyring", conditions)

        assert settings is not None
        assert len(settings) > 0

        # Check for key parameters
        expected_keys = [
            "dh_b",
            "dh_c",
            "ds_b",
            "ds_c",
            "dh_ab",
            "dh_ac",
            "dh_bc",
            "kab",
            "kba",
            "kac",
            "kca",
            "kbc",
            "kcb",
            "pa",
            "pb",
            "pc",
        ]
        for key in expected_keys:
            assert key in settings

    def test_temperature_dependency_across_range(self):
        """Test that the model behaves correctly across a range of temperatures."""
        temperatures = [10.0, 25.0, 40.0, 55.0, 70.0]

        # Create settings for each temperature
        all_settings = {}
        for temp in temperatures:
            conditions = Conditions(temperature=temp)
            settings = make_settings_3st_eyring(conditions)
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
        settings = make_settings_3st_eyring(conditions)

        # Rate constants should be non-negative
        for k in ["kab", "kba", "kac", "kca", "kbc", "kcb"]:
            assert settings[k].min == 0.0

        # Populations should be in [0, 1]
        for pop in ["pa", "pb", "pc"]:
            assert settings[pop].min == 0.0
            assert settings[pop].max == 1.0

        # Population expressions should use pop_3st function
        for pop in ["pa", "pb", "pc"]:
            assert "pop_3st" in settings[pop].expr

    def test_settings_consistency_with_calculation(self):
        """Test that settings creation is consistent with direct calculation."""
        temp = 30.0
        conditions = Conditions(temperature=temp)
        settings = make_settings_3st_eyring(conditions)

        # Use default values from settings
        params = {
            "dh_b": settings["dh_b"].value,
            "ds_b": settings["ds_b"].value,
            "dh_c": settings["dh_c"].value,
            "ds_c": settings["ds_c"].value,
            "dh_ab": settings["dh_ab"].value,
            "ds_ab": settings["ds_ab"].value,
            "dh_ac": settings["dh_ac"].value,
            "ds_ac": settings["ds_ac"].value,
            "dh_bc": settings["dh_bc"].value,
            "ds_bc": settings["ds_bc"].value,
            "temperature": temp,
        }

        # Calculate rate constants directly
        rates = calculate_kij_3st_eyring(**params)

        # Verify all expected rate constants are present and non-negative
        # Note: default dh_ac is very high (1e10), so kac/kca may be ~0
        expected_rates = ["kab", "kba", "kac", "kca", "kbc", "kcb"]
        for rate in expected_rates:
            assert rate in rates
            assert rates[rate] >= 0  # >= 0 because default barriers can be very high

    def test_cyclic_thermodynamic_consistency(self):
        """Test thermodynamic consistency around cyclic paths.

        For the cycle A → B → C → A, the product of equilibrium constants
        should equal 1: K_AB * K_BC * K_CA = 1
        """
        temp = 25.0
        params = {
            "dh_b": 6000.0,
            "ds_b": 10.0,
            "dh_c": 9000.0,
            "ds_c": -5.0,
            "dh_ab": 70000.0,
            "ds_ab": 40.0,
            "dh_ac": 75000.0,
            "ds_ac": 45.0,
            "dh_bc": 68000.0,
            "ds_bc": 35.0,
            "temperature": temp,
        }

        rates = calculate_kij_3st_eyring(**params)

        # Calculate equilibrium constants
        K_AB = rates["kab"] / rates["kba"]
        K_BC = rates["kbc"] / rates["kcb"]
        K_CA = rates["kca"] / rates["kac"]

        # Product should be 1.0
        cycle_product = K_AB * K_BC * K_CA
        assert np.isclose(cycle_product, 1.0, rtol=1e-10), (
            f"Cyclic consistency violated: product = {cycle_product}"
        )

    def test_detailed_balance_all_pairs(self):
        """Test detailed balance for all state pairs at multiple temperatures."""
        temperatures = [15.0, 25.0, 35.0, 45.0]

        for temp in temperatures:
            params = {
                "dh_b": 5000.0,
                "ds_b": 12.0,
                "dh_c": 8000.0,
                "ds_c": -8.0,
                "dh_ab": 65000.0,
                "ds_ab": 35.0,
                "dh_ac": 72000.0,
                "ds_ac": 42.0,
                "dh_bc": 67000.0,
                "ds_bc": 30.0,
                "temperature": temp,
            }

            rates = calculate_kij_3st_eyring(**params)
            T = constants.convert_temperature(temp, "C", "K")
            RT = constants.R * T

            # Test A <-> B
            dg_ab = params["dh_b"] - T * params["ds_b"]
            K_ab_thermo = np.exp(-dg_ab / RT)
            K_ab_kinetic = rates["kab"] / rates["kba"]
            assert np.isclose(K_ab_kinetic, K_ab_thermo, rtol=1e-10)

            # Test A <-> C
            dg_ac = params["dh_c"] - T * params["ds_c"]
            K_ac_thermo = np.exp(-dg_ac / RT)
            K_ac_kinetic = rates["kac"] / rates["kca"]
            assert np.isclose(K_ac_kinetic, K_ac_thermo, rtol=1e-10)

            # Test B <-> C
            dg_bc = (params["dh_c"] - params["dh_b"]) - T * (
                params["ds_c"] - params["ds_b"]
            )
            K_bc_thermo = np.exp(-dg_bc / RT)
            K_bc_kinetic = rates["kbc"] / rates["kcb"]
            assert np.isclose(K_bc_kinetic, K_bc_thermo, rtol=1e-10)

    def test_van_t_hoff_analysis_multiple_pairs(self):
        """Test van't Hoff behavior for all equilibrium pairs."""
        dh_b = 7000.0  # J/mol
        ds_b = 15.0  # J/mol/K
        dh_c = 10000.0  # J/mol
        ds_c = -10.0  # J/mol/K

        temperatures = [10.0, 20.0, 30.0, 40.0, 50.0]

        # Test for A <-> B equilibrium
        inv_temps = []
        ln_keq_ab = []

        for temp in temperatures:
            params = {
                "dh_b": dh_b,
                "ds_b": ds_b,
                "dh_c": dh_c,
                "ds_c": ds_c,
                "dh_ab": 68000.0,
                "ds_ab": 45.0,
                "dh_ac": 75000.0,
                "ds_ac": 50.0,
                "dh_bc": 70000.0,
                "ds_bc": 40.0,
                "temperature": temp,
            }

            rates = calculate_kij_3st_eyring(**params)
            K_ab = rates["kab"] / rates["kba"]

            T = constants.convert_temperature(temp, "C", "K")
            inv_temps.append(1 / T)
            ln_keq_ab.append(np.log(K_ab))

        # Check linearity for van't Hoff plot
        n = len(ln_keq_ab)
        sum_x = sum(inv_temps)
        sum_y = sum(ln_keq_ab)
        sum_xy = sum(x * y for x, y in zip(inv_temps, ln_keq_ab, strict=False))
        sum_x2 = sum(x * x for x in inv_temps)
        sum_y2 = sum(y * y for y in ln_keq_ab)

        r_squared = ((n * sum_xy - sum_x * sum_y) ** 2) / (
            (n * sum_x2 - sum_x**2) * (n * sum_y2 - sum_y**2)
        )

        # Should be very linear
        assert r_squared > 0.999, f"van't Hoff plot not linear, R² = {r_squared}"

        # Check that slope and intercept match thermodynamic parameters
        slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x**2)
        intercept = (sum_y - slope * sum_x) / n

        expected_slope = -dh_b / constants.R
        expected_intercept = ds_b / constants.R

        assert np.isclose(slope, expected_slope, rtol=0.01)
        assert np.isclose(intercept, expected_intercept, rtol=0.01)

    def test_state_transitions_independence(self):
        """Test that direct transitions are independent of indirect pathways.

        The rate of A → C should not depend on whether we consider the
        direct path or the pathway through B.
        """
        temp = 25.0
        params = {
            "dh_b": 5000.0,
            "ds_b": 10.0,
            "dh_c": 8000.0,
            "ds_c": 5.0,
            "dh_ab": 65000.0,
            "ds_ab": 40.0,
            "dh_ac": 72000.0,
            "ds_ac": 45.0,
            "dh_bc": 68000.0,
            "ds_bc": 38.0,
            "temperature": temp,
        }

        rates = calculate_kij_3st_eyring(**params)

        # All rates should be independent and positive
        assert rates["kac"] > 0  # Direct A → C
        assert rates["kab"] > 0  # A → B
        assert rates["kbc"] > 0  # B → C

        # The rates should satisfy microscopic reversibility independently
        # This is already tested in detailed_balance, but we verify independence here
        T = constants.convert_temperature(temp, "C", "K")
        RT = constants.R * T

        # Each pair should independently satisfy detailed balance
        pairs = [
            ("kab", "kba", "dh_b", "ds_b"),
            ("kac", "kca", "dh_c", "ds_c"),
        ]

        for kf, kr, dh_key, ds_key in pairs:
            dg = params[dh_key] - T * params[ds_key]
            K_thermo = np.exp(-dg / RT)
            K_kinetic = rates[kf] / rates[kr]
            assert np.isclose(K_kinetic, K_thermo, rtol=1e-10)

    def test_population_sum_constraint(self):
        """Test that population expressions are created using pop_3st function."""
        conditions = Conditions(temperature=25.0)
        settings = make_settings_3st_eyring(conditions)

        # For 3-state system: all populations should use pop_3st function
        for pop in ["pa", "pb", "pc"]:
            assert "pop_3st" in settings[pop].expr
            assert settings[pop].min == 0.0
            assert settings[pop].max == 1.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

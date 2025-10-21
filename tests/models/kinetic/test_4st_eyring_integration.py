"""Integration tests for 4st_eyring model."""

import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.models.kinetic.settings_4st_eyring import make_settings_4st_eyring


class TestEyring4stIntegration:
    """Integration tests for the 4st_eyring model."""

    def test_model_registration(self):
        """Test that the model is properly registered in the factory."""
        # Explicitly trigger registration
        from chemex.models.kinetic.settings_4st_eyring import register

        register()

        # The model should be registered
        assert "4st_eyring" in model_factory.set

        # Test that we can create settings
        conditions = Conditions(temperature=25.0)
        settings = model_factory.create("4st_eyring", conditions)

        assert settings is not None
        assert len(settings) > 0

        # Check for key parameters
        expected_keys = ["dh_b", "dh_c", "dh_d", "kab", "kba", "pa", "pb"]
        for key in expected_keys:
            assert key in settings

    def test_temperature_dependency_across_range(self):
        """Test that the model behaves correctly across a range of temperatures."""
        temperatures = [10.0, 25.0, 40.0, 55.0, 70.0]

        # Create settings for each temperature
        all_settings = {}
        for temp in temperatures:
            conditions = Conditions(temperature=temp)
            settings = make_settings_4st_eyring(conditions)
            all_settings[temp] = settings

        # Check that rate constant expressions change with temperature
        for i, temp1 in enumerate(temperatures[:-1]):
            temp2 = temperatures[i + 1]

            # The expressions should contain the temperature values
            expr1 = all_settings[temp1]["kab"].expr
            expr2 = all_settings[temp2]["kab"].expr

            assert str(temp1) in expr1
            assert str(temp2) in expr2
            assert expr1 != expr2

    def test_parameter_bounds_and_defaults(self):
        """Test that default parameter values and bounds are reasonable."""
        conditions = Conditions(temperature=25.0)
        settings = make_settings_4st_eyring(conditions)

        # Check enthalpy parameters
        dh_params = ["dh_b", "dh_c", "dh_d"]
        for param in dh_params:
            value = settings[param].value
            assert value is not None
            assert value > 0  # Should be positive
            assert value < 100000  # Should be reasonable (< 100 kJ/mol)
            assert settings[param].vary is True  # Should be variable

        # Check activation enthalpies
        activation_params = ["dh_ab", "dh_ac", "dh_ad", "dh_bc", "dh_bd", "dh_cd"]
        for param in activation_params:
            value = settings[param].value
            assert value is not None
            assert value > 50000  # Should be reasonable activation energy
            assert value < 200000  # Should not be too high

        # Check entropy parameters (typically zero by default)
        entropy_params = [
            "ds_b",
            "ds_c",
            "ds_d",
            "ds_ab",
            "ds_ac",
            "ds_ad",
            "ds_bc",
            "ds_bd",
            "ds_cd",
        ]
        for param in entropy_params:
            assert settings[param].value == 0.0  # Default to zero
            assert settings[param].vary is False  # Fixed by default

    def test_rate_constant_generation(self):
        """Test that all required rate constants are generated."""
        conditions = Conditions(temperature=30.0)
        settings = make_settings_4st_eyring(conditions)

        # Check that all 12 rate constants are present
        states = "abcd"
        for i in states:
            for j in states:
                if i != j:
                    rate_key = f"k{i}{j}"
                    assert rate_key in settings

                    # Check that the expression is properly formatted
                    expr = settings[rate_key].expr
                    assert "kij_4st_eyring" in expr
                    assert f"['k{i}{j}']" in expr
                    assert "30.0" in expr  # Temperature should be embedded

    def test_population_generation(self):
        """Test that population parameters are correctly generated."""
        conditions = Conditions(temperature=25.0)
        settings = make_settings_4st_eyring(conditions)

        # Check that all 4 populations are present
        for state in "abcd":
            pop_key = f"p{state}"
            assert pop_key in settings

            # Check that the expression uses the pop_4st function
            expr = settings[pop_key].expr
            assert "pop_4st(" in expr
            assert f"['p{state}']" in expr

    def test_consistency_with_2st_eyring(self):
        """Test consistency with simpler 2-state Eyring model for A<->B only."""
        # This test would require importing the 2st_eyring model
        # For now, we just check that A<->B transitions are properly handled
        conditions = Conditions(temperature=25.0)
        settings = make_settings_4st_eyring(conditions)

        # Check that A<->B parameters are present and reasonable
        assert "kab" in settings
        assert "kba" in settings
        assert "pa" in settings
        assert "pb" in settings

        # Check that the rate expressions contain the relevant thermodynamic parameters
        kab_expr = settings["kab"].expr
        assert "dh_b" in kab_expr or "DH_B" in kab_expr
        assert "dh_ab" in kab_expr or "DH_AB" in kab_expr


if __name__ == "__main__":
    pytest.main([__file__])

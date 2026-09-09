"""Cross-model thermodynamic invariants for the public Eyring models."""

from __future__ import annotations

import math
import sys
from collections.abc import Callable
from decimal import Decimal, localcontext
from types import SimpleNamespace

import pytest
from scipy import constants

from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Method
from chemex.configuration.parameters import DefaultSetting
from chemex.models.factory import model_factory
from chemex.models.kinetic._eyring import calculate_rate_from_coordinates
from chemex.models.kinetic.settings_2st_eyring import (
    calculate_kij_2st_eyring,
    calculate_populations_2st_eyring,
    make_settings_2st_eyring,
)
from chemex.models.kinetic.settings_2st_eyring import (
    register as register_2st_eyring,
)
from chemex.models.kinetic.settings_3st_eyring import (
    calculate_kij_3st_eyring_fork,
    calculate_kij_3st_eyring_linear,
    calculate_populations_3st_eyring,
    make_settings_3st_eyring_fork,
    make_settings_3st_eyring_linear,
)
from chemex.models.kinetic.settings_4st_eyring import (
    calculate_kij_4st_eyring,
    calculate_populations_4st_eyring,
    make_settings_4st_eyring,
)
from chemex.nmr.basis import Basis
from chemex.optimize.deterministic_uncertainty import (
    compile_model_constraint_linearization_capabilities,
)
from chemex.optimize.uncertainty import (
    FunctionAnalyticPartialCapability,
)
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import ConstraintDomainError
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.userfunctions import (
    AnalyticFunctionLinearization,
    function_linearization_registry,
    user_function_registry,
)
from chemex.runtime import AnalysisSession

RateCalculator = Callable[..., dict[str, float]]


def _decimal_directional_rate(
    celsius: str,
    initial: tuple[str, str],
    transition_state: tuple[str, str],
) -> float:
    """Evaluate one directional rate independently from exact SI literals."""
    with localcontext() as context:
        context.prec = 100
        gas_constant = Decimal("8.31446261815324")
        frequency_factor = Decimal("1.380649e-23") / Decimal("6.62607015e-34")
        kelvin = Decimal(celsius) + Decimal("273.15")
        activation_enthalpy = Decimal(transition_state[0]) - Decimal(initial[0])
        activation_entropy = Decimal(transition_state[1]) - Decimal(initial[1])
        exponent = activation_entropy / gas_constant - activation_enthalpy / (
            gas_constant * kelvin
        )
        return float(frequency_factor * kelvin * exponent.exp())


def _prepare_eyring_model(
    model_name: str,
    *,
    temperature: float = 25.0,
    defaults: dict[str, float] | None = None,
) -> tuple[AnalysisSession, dict[str, str]]:
    conditions = Conditions(
        h_larmor_frq=600.0,
        temperature=temperature,
        p_total=1.0e-3,
        l_total=2.0e-3,
    )
    session = AnalysisSession.create()
    session.set_model(model_name)
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    spin_system = SpinSystem.from_name("G23N-HN")
    settings = model_factory.create(model_name, conditions)
    local_ids = {
        name: setting.name_setting.get_param_name(spin_system, conditions).id_
        for name, setting in settings.items()
    }
    config = SimpleNamespace(
        conditions=conditions,
        to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
    )
    session.parameter_factory.create_parameters(
        config,  # ty: ignore[invalid-argument-type]
        basis=basis,
        spin_system=spin_system,
    )
    session.parameters.set_defaults(
        [
            (ParamName.from_section(name), DefaultSetting(value))
            for name, value in (defaults or {}).items()
        ]
    )
    assert session.parameter_factory.try_seal_definitions(), repr(
        session.parameter_factory.native_construction_error
    )
    return session, local_ids


def _assert_log_detailed_balance(
    rates: dict[str, float],
    states: dict[str, tuple[float, float]],
    edges: tuple[str, ...],
    celsius: float,
) -> None:
    kelvin = celsius + constants.zero_Celsius
    for edge in edges:
        initial, final = edge
        dh_i, ds_i = states[initial]
        dh_j, ds_j = states[final]
        actual = math.log(rates[f"k{edge}"]) - math.log(rates[f"k{edge[::-1]}"])
        expected = -((dh_j - dh_i) - kelvin * (ds_j - ds_i)) / (constants.R * kelvin)
        assert actual == pytest.approx(expected, rel=2.0e-14, abs=2.0e-14)


@pytest.mark.parametrize("celsius", (-241.15, 25.0, 50.0, 1.0e6))
def test_two_state_log_detailed_balance_across_temperature(celsius: float) -> None:
    rates = calculate_kij_2st_eyring(8_000.0, 10.0, 65_000.0, 20.0, celsius)

    _assert_log_detailed_balance(
        rates,
        {"a": (0.0, 0.0), "b": (8_000.0, 10.0)},
        ("ab",),
        celsius,
    )


def test_subnormal_representable_edge_preserves_log_detailed_balance() -> None:
    kelvin = 32.0
    state_enthalpy = 10.0 * constants.R * kelvin
    transition_enthalpy = 750.0 * constants.R * kelvin
    rates = calculate_kij_2st_eyring(
        state_enthalpy,
        0.0,
        transition_enthalpy,
        0.0,
        kelvin - constants.zero_Celsius,
    )

    assert 0.0 < rates["kab"] < sys.float_info.min
    actual = math.log(rates["kab"]) - math.log(rates["kba"])
    # Subnormal quantization limits the log ratio here; 2e-10 is below one
    # relative spacing of the smaller returned rate.
    assert actual == pytest.approx(-10.0, rel=0.0, abs=2.0e-10)


@pytest.mark.parametrize(
    ("calculator", "transition_terms", "edges"),
    (
        (
            calculate_kij_3st_eyring_linear,
            (65_000.0, 20.0, 70_000.0, 5.0),
            ("ab", "bc"),
        ),
        (calculate_kij_3st_eyring_fork, (65_000.0, 20.0, 72_000.0, 15.0), ("ab", "ac")),
    ),
)
def test_three_state_topologies_obey_log_detailed_balance(
    calculator: RateCalculator,
    transition_terms: tuple[float, ...],
    edges: tuple[str, ...],
) -> None:
    rates = calculator(8_000.0, 10.0, 12_000.0, -5.0, *transition_terms, 25.0)

    _assert_log_detailed_balance(
        rates,
        {"a": (0.0, 0.0), "b": (8_000.0, 10.0), "c": (12_000.0, -5.0)},
        edges,
        25.0,
    )


def _four_state_parameters(transition_shift: float = 0.0) -> dict[str, float]:
    return {
        "dh_b": 8_000.0,
        "ds_b": 10.0,
        "dh_c": 12_000.0,
        "ds_c": -5.0,
        "dh_d": 15_000.0,
        "ds_d": 15.0,
        "dh_ab": 75_000.0 + transition_shift,
        "ds_ab": 50.0,
        "dh_ac": 80_000.0 + transition_shift,
        "ds_ac": 30.0,
        "dh_ad": 85_000.0 + transition_shift,
        "ds_ad": 70.0,
        "dh_bc": 70_000.0 + transition_shift,
        "ds_bc": 20.0,
        "dh_bd": 77_000.0 + transition_shift,
        "ds_bd": 40.0,
        "dh_cd": 72_000.0 + transition_shift,
        "ds_cd": 25.0,
        "temperature": 25.0,
    }


def test_four_state_rates_obey_pairwise_balance_and_cycle_closure() -> None:
    parameters = _four_state_parameters()
    rates = calculate_kij_4st_eyring(**parameters)
    states = {
        "a": (0.0, 0.0),
        "b": (parameters["dh_b"], parameters["ds_b"]),
        "c": (parameters["dh_c"], parameters["ds_c"]),
        "d": (parameters["dh_d"], parameters["ds_d"]),
    }

    _assert_log_detailed_balance(
        rates,
        states,
        ("ab", "ac", "ad", "bc", "bd", "cd"),
        parameters["temperature"],
    )
    cycle_affinity = math.fsum(
        (
            math.log(rates["kab"]) - math.log(rates["kba"]),
            math.log(rates["kbc"]) - math.log(rates["kcb"]),
            math.log(rates["kca"]) - math.log(rates["kac"]),
        )
    )
    assert cycle_affinity == pytest.approx(0.0, abs=3.0e-14)


def test_legacy_four_state_positional_wrapper_maps_all_six_barriers() -> None:
    """Protect the long retained callable surface against positional edge swaps."""
    temperature = "25.0"
    states = {
        "a": ("0", "0"),
        "b": ("6100", "4"),
        "c": ("11200", "-8"),
        "d": ("-3200", "13"),
    }
    transition_states = {
        "ab": ("-25000", "-10"),
        "ac": ("42000", "3"),
        "ad": ("48000", "-7"),
        "bc": ("55000", "19"),
        "bd": ("60000", "-13"),
        "cd": ("67000", "27"),
    }

    # Deliberately positional: this is compatibility coverage for the historical
    # 19-argument order, not a test of the cleaner keyword-facing declaration.
    rates = calculate_kij_4st_eyring(
        6100.0,
        4.0,
        11200.0,
        -8.0,
        -3200.0,
        13.0,
        -25000.0,
        -10.0,
        42000.0,
        3.0,
        48000.0,
        -7.0,
        55000.0,
        19.0,
        60000.0,
        -13.0,
        67000.0,
        27.0,
        25.0,
    )

    for edge, transition_state in transition_states.items():
        for initial, final in (edge, edge[::-1]):
            expected = _decimal_directional_rate(
                temperature,
                states[initial],
                transition_state,
            )
            actual = rates[f"k{initial}{final}"]
            # Thirty-two ulps covers the independent Decimal-to-libm evaluation
            # path while remaining far smaller than any deliberate edge swap.
            assert abs(actual - expected) <= 32 * math.ulp(expected), edge

    assert rates["kab"] > 1.0e16
    assert rates["kab"] != 1.0e16
    assert rates["kba"] == pytest.approx(
        _decimal_directional_rate(temperature, states["b"], transition_states["ab"]),
        rel=3.0e-15,
    )


def test_rate_is_not_saturated_at_former_limit() -> None:
    rates = calculate_kij_2st_eyring(0.0, 0.0, -50_000.0, 0.0, 25.0)

    assert rates["kab"] > 1.0e16
    assert rates["kba"] > 1.0e16
    assert rates["kab"] == rates["kba"]


def test_two_state_population_ratio_comes_directly_from_state_coordinates() -> None:
    populations = calculate_populations_2st_eyring(8_000.0, 10.0, 25.0)
    kelvin = 298.15
    expected_ratio = math.exp(-(8_000.0 - kelvin * 10.0) / (constants.R * kelvin))

    assert populations["pb"] / populations["pa"] == pytest.approx(
        expected_ratio, rel=2.0e-15
    )


def test_population_boundary_distinguishes_balance_from_unrepresentability() -> None:
    celsius = -273.14
    representable_enthalpy = 61.72329710480702
    rates = calculate_kij_2st_eyring(
        representable_enthalpy,
        0.0,
        representable_enthalpy,
        0.0,
        celsius,
    )
    populations = calculate_populations_2st_eyring(
        representable_enthalpy,
        0.0,
        celsius,
    )

    assert populations["pb"] > 0.0
    assert math.log(rates["kab"]) - math.log(rates["kba"]) == pytest.approx(
        math.log(populations["pb"]) - math.log(populations["pa"]),
        abs=6.0e-10,
    )

    unrepresentable_enthalpy = 61.959375430421595
    still_representable_rates = calculate_kij_2st_eyring(
        unrepresentable_enthalpy,
        0.0,
        unrepresentable_enthalpy,
        0.0,
        celsius,
    )
    assert all(rate > 0.0 for rate in still_representable_rates.values())
    with pytest.raises(ValueError, match="population.*below binary64 representability"):
        calculate_populations_2st_eyring(
            unrepresentable_enthalpy,
            0.0,
            celsius,
        )


def test_three_state_population_is_independent_of_topology() -> None:
    expected = calculate_populations_3st_eyring(
        8_000.0,
        10.0,
        12_000.0,
        -5.0,
        25.0,
    )
    linear = make_settings_3st_eyring_linear(Conditions(temperature=25.0))
    fork = make_settings_3st_eyring_fork(Conditions(temperature=25.0))

    assert expected.keys() == {"pa", "pb", "pc"}
    assert all("pop_3st_eyring" in linear[name].expr for name in expected)
    assert all(linear[name].expr == fork[name].expr for name in expected)
    assert all("pop_3st(" not in linear[name].expr for name in expected)


def test_four_state_populations_are_invariant_to_uniform_kinetic_slowdown() -> None:
    ordinary_parameters = _four_state_parameters()
    slow_parameters = _four_state_parameters(transition_shift=100_000.0)
    ordinary_rates = calculate_kij_4st_eyring(**ordinary_parameters)
    slow_rates = calculate_kij_4st_eyring(**slow_parameters)
    populations = calculate_populations_4st_eyring(
        ordinary_parameters["dh_b"],
        ordinary_parameters["ds_b"],
        ordinary_parameters["dh_c"],
        ordinary_parameters["ds_c"],
        ordinary_parameters["dh_d"],
        ordinary_parameters["ds_d"],
        ordinary_parameters["temperature"],
    )
    slow_populations = calculate_populations_4st_eyring(
        slow_parameters["dh_b"],
        slow_parameters["ds_b"],
        slow_parameters["dh_c"],
        slow_parameters["ds_c"],
        slow_parameters["dh_d"],
        slow_parameters["ds_d"],
        slow_parameters["temperature"],
    )

    assert slow_populations == populations
    rate_scales = {
        name: slow_rates[name] / ordinary_rate
        for name, ordinary_rate in ordinary_rates.items()
    }
    assert max(rate_scales.values()) < 1.0e-15
    assert max(rate_scales.values()) == pytest.approx(
        min(rate_scales.values()), rel=8.0e-14
    )
    settings = make_settings_4st_eyring(Conditions(temperature=25.0))
    assert all("pop_4st_eyring" in settings[name].expr for name in populations)
    assert all("pop_4st(" not in settings[name].expr for name in populations)


def test_settings_expose_each_directional_barrier_without_tuple_ordering() -> None:
    settings = make_settings_4st_eyring(Conditions(temperature=25.0))

    assert settings["kab"].expr == ("eyring_rate(0.0,0.0,{dh_ab},{ds_ab},25.0)['rate']")
    assert settings["kba"].expr == (
        "eyring_rate({dh_b},{ds_b},{dh_ab},{ds_ab},25.0)['rate']"
    )
    assert settings["kbc"].expr == (
        "eyring_rate({dh_b},{ds_b},{dh_bc},{ds_bc},25.0)['rate']"
    )
    assert settings["kcb"].expr == (
        "eyring_rate({dh_c},{ds_c},{dh_bc},{ds_bc},25.0)['rate']"
    )


def test_two_state_settings_use_thermodynamic_population_authority() -> None:
    settings = make_settings_2st_eyring(Conditions(temperature=25.0))

    assert settings["pa"].expr == "pop_2st_eyring({dh_b},{ds_b},25.0)['pa']"
    assert settings["pb"].expr == "pop_2st_eyring({dh_b},{ds_b},25.0)['pb']"


def test_directional_rate_registers_exact_thermodynamic_partials() -> None:
    register_2st_eyring()
    capability = next(
        item
        for item in function_linearization_registry.get("2st_eyring")
        if isinstance(item, AnalyticFunctionLinearization)
    )
    arguments = (8_000.0, 10.0, 65_000.0, 20.0, 25.0)
    rate = calculate_rate_from_coordinates(*arguments)
    kelvin = 298.15
    expected = (
        rate / (constants.R * kelvin),
        -rate / constants.R,
        -rate / (constants.R * kelvin),
        rate / constants.R,
    )

    assert (capability.function_id, capability.component) == ("eyring_rate", "rate")
    assert tuple(partial(*arguments) for partial in capability.partials[:4]) == (
        pytest.approx(expected, rel=2.0e-15)
    )


def test_eyring_parameter_scoping_preserves_current_contract() -> None:
    settings = make_settings_2st_eyring(Conditions(temperature=25.0))
    spin_g = SpinSystem.from_name("G23N-HN")
    spin_a = SpinSystem.from_name("A24N-HN")
    base = Conditions(
        temperature=0.0,
        h_larmor_frq=600.0,
        p_total=1.0e-3,
        l_total=2.0e-3,
    )
    changed_temperature_and_field = Conditions(
        temperature=25.0,
        h_larmor_frq=800.0,
        p_total=1.0e-3,
        l_total=2.0e-3,
    )
    changed_concentration = Conditions(
        temperature=0.0,
        h_larmor_frq=600.0,
        p_total=2.0e-3,
        l_total=2.0e-3,
    )

    dh_base = settings["dh_b"].name_setting.get_param_name(spin_g, base)
    dh_other = settings["dh_b"].name_setting.get_param_name(
        spin_a, changed_temperature_and_field
    )
    kab_base = settings["kab"].name_setting.get_param_name(spin_g, base)
    kab_temperature = settings["kab"].name_setting.get_param_name(
        spin_a, changed_temperature_and_field
    )
    kab_concentration = settings["kab"].name_setting.get_param_name(
        spin_g, changed_concentration
    )

    assert dh_base.id_ == dh_other.id_
    assert kab_base.id_ != kab_temperature.id_
    assert dh_base.id_ == "__DH_B__1_000E_03M_2_000E_03M"
    assert kab_base.id_ == "__KAB__0_0C_1_000E_03M_2_000E_03M"
    assert "600_0MHZ" not in kab_base.id_
    assert "G23" not in kab_base.id_
    assert kab_base.id_ != kab_concentration.id_
    assert (
        calculate_kij_2st_eyring(8_000.0, 10.0, 65_000.0, 20.0, 0.0)["kab"]
        != calculate_kij_2st_eyring(8_000.0, 10.0, 65_000.0, 20.0, 25.0)["kab"]
    )


@pytest.mark.parametrize(
    ("model_name", "legacy_function"),
    (
        ("2st_eyring", "kij_2st_eyring"),
        ("2st_eyring", "pop_2st"),
        ("3st_eyring", "kij_3st_eyring"),
        ("3st_eyring", "pop_3st"),
        ("3st_eyring_linear", "kij_3st_eyring"),
        ("3st_eyring_fork", "kij_3st_eyring_fork"),
        ("4st_eyring", "kij_4st_eyring"),
        ("4st_eyring", "pop_4st"),
    ),
)
def test_legacy_eyring_calculator_bindings_remain_available(
    model_name: str,
    legacy_function: str,
) -> None:
    session = AnalysisSession.create()
    session.set_model(model_name)

    assert legacy_function in user_function_registry.get(model_name)


@pytest.mark.parametrize(
    ("model_name", "population_names"),
    (
        ("2st_eyring", ("pa", "pb")),
        ("3st_eyring", ("pa", "pb", "pc")),
        ("3st_eyring_linear", ("pa", "pb", "pc")),
        ("3st_eyring_fork", ("pa", "pb", "pc")),
        ("4st_eyring", ("pa", "pb", "pc", "pd")),
    ),
)
def test_all_eyring_models_compile_rate_and_population_uncertainty_capabilities(
    model_name: str,
    population_names: tuple[str, ...],
) -> None:
    session, local_ids = _prepare_eyring_model(model_name)
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    output_ids = (local_ids["kab"], *(local_ids[name] for name in population_names))
    parameterization = session.compile_parameterization(Method(), set(output_ids))

    compiled = compile_model_constraint_linearization_capabilities(
        parameterization,
        output_ids,
    )

    assert any(
        isinstance(capability, FunctionAnalyticPartialCapability)
        and (capability.function_id, capability.component) == ("eyring_rate", "rate")
        for capability in compiled.capabilities
    )
    population_capabilities = tuple(
        capability
        for capability in compiled.capabilities
        if isinstance(capability, FunctionAnalyticPartialCapability)
        and capability.component in population_names
    )
    assert {capability.component for capability in population_capabilities} == set(
        population_names
    )
    assert all(
        capability.implementation_identity
        == "eyring-thermodynamic-population-partials-v1"
        for capability in population_capabilities
    )


@pytest.mark.parametrize(
    ("barrier", "message"),
    (
        (200_000.0, "below binary64 representability"),
        (-200_000.0, "exceeds maximum finite binary64"),
    ),
)
def test_unrepresentable_eyring_rate_is_a_typed_constraint_domain_failure(
    barrier: float,
    message: str,
) -> None:
    session, _ = _prepare_eyring_model(
        "2st_eyring",
        temperature=-273.149999,
        defaults={"dh_ab": barrier},
    )

    assert not session.try_build_analysis_values()
    error = session.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert message in str(error.__cause__)


def test_unrepresentable_population_is_typed_constraint_domain_failure() -> None:
    enthalpy = 61.959375430421595
    session, _ = _prepare_eyring_model(
        "2st_eyring",
        temperature=-273.14,
        defaults={"dh_b": enthalpy, "dh_ab": enthalpy},
    )

    assert not session.try_build_analysis_values()
    error = session.parameter_factory.native_construction_error
    assert isinstance(error, ConstraintDomainError)
    assert isinstance(error.__cause__, ValueError)
    assert "population" in str(error.__cause__)
    assert "below binary64 representability" in str(error.__cause__)

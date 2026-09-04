"""Unit tests for the three-state Eyring topology variants."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest
from scipy import constants

from chemex.configuration.conditions import Conditions
from chemex.models.kinetic.settings_3st_eyring import (
    MAX_RATE_CONSTANT,
    calculate_kij_3st_eyring_fork,
    calculate_kij_3st_eyring_linear,
    make_settings_3st_eyring_fork,
    make_settings_3st_eyring_linear,
)
from chemex.parameters.setting import ParamLocalSetting

RateCalculator = Callable[..., dict[str, float]]

STATE_TERMS = (8000.0, 10.0, 12000.0, -5.0)
AB_TERMS = (75000.0, 50.0)
BC_TERMS = (70000.0, 20.0)
AC_TERMS = (80000.0, 30.0)


@pytest.mark.parametrize(
    ("calculator", "arguments", "expected_rates"),
    [
        (
            calculate_kij_3st_eyring_linear,
            (*STATE_TERMS, *AB_TERMS, *BC_TERMS, 25.0),
            {"kab", "kba", "kbc", "kcb"},
        ),
        (
            calculate_kij_3st_eyring_fork,
            (*STATE_TERMS, *AB_TERMS, *AC_TERMS, 25.0),
            {"kab", "kba", "kac", "kca"},
        ),
    ],
)
def test_calculators_return_only_active_finite_rates(
    calculator: RateCalculator,
    arguments: tuple[float, ...],
    expected_rates: set[str],
) -> None:
    rates = calculator(*arguments)

    assert rates.keys() == expected_rates
    assert all(0.0 < rate <= MAX_RATE_CONSTANT for rate in rates.values())
    assert all(np.isfinite(rate) for rate in rates.values())


@pytest.mark.parametrize(
    ("calculator", "arguments", "equilibria"),
    [
        (
            calculate_kij_3st_eyring_linear,
            (*STATE_TERMS, *AB_TERMS, *BC_TERMS, 25.0),
            (("kab", "kba", 8000.0, 10.0), ("kbc", "kcb", 4000.0, -15.0)),
        ),
        (
            calculate_kij_3st_eyring_fork,
            (*STATE_TERMS, *AB_TERMS, *AC_TERMS, 25.0),
            (("kab", "kba", 8000.0, 10.0), ("kac", "kca", 12000.0, -5.0)),
        ),
    ],
)
def test_active_edges_obey_detailed_balance(
    calculator: RateCalculator,
    arguments: tuple[float, ...],
    equilibria: tuple[tuple[str, str, float, float], ...],
) -> None:
    rates = calculator(*arguments)
    kelvin = 25.0 + constants.zero_Celsius

    for forward, reverse, dh, ds in equilibria:
        expected_ratio = np.exp(-(dh - kelvin * ds) / (constants.R * kelvin))
        assert rates[forward] / rates[reverse] == pytest.approx(
            expected_ratio,
            rel=1e-12,
        )


@pytest.mark.parametrize(
    ("calculator", "transition_terms", "rate_names"),
    [
        (
            calculate_kij_3st_eyring_linear,
            (*AB_TERMS, *BC_TERMS),
            {"kab", "kba", "kbc", "kcb"},
        ),
        (
            calculate_kij_3st_eyring_fork,
            (*AB_TERMS, *AC_TERMS),
            {"kab", "kba", "kac", "kca"},
        ),
    ],
)
def test_temperature_dependence_is_preserved(
    calculator: RateCalculator,
    transition_terms: tuple[float, ...],
    rate_names: set[str],
) -> None:
    rates_25 = calculator(*STATE_TERMS, *transition_terms, 25.0)
    rates_50 = calculator(*STATE_TERMS, *transition_terms, 50.0)

    assert all(rates_50[name] > rates_25[name] for name in rate_names)


@pytest.mark.parametrize(
    "calculator",
    [calculate_kij_3st_eyring_linear, calculate_kij_3st_eyring_fork],
)
def test_rate_clipping_policy_is_preserved(calculator: RateCalculator) -> None:
    rates = calculator(0.0, 0.0, 0.0, 0.0, -2.0e5, 0.0, -2.0e5, 0.0, 25.0)

    assert set(rates.values()) == {MAX_RATE_CONSTANT}


def test_calculators_cache_repeated_inputs() -> None:
    arguments = (*STATE_TERMS, *AB_TERMS, *BC_TERMS, 25.0)

    assert calculate_kij_3st_eyring_linear(
        *arguments,
    ) is calculate_kij_3st_eyring_linear(*arguments)


@pytest.mark.parametrize(
    ("maker", "present_edges", "absent_edge"),
    [
        (make_settings_3st_eyring_linear, ("ab", "bc"), "ac"),
        (make_settings_3st_eyring_fork, ("ab", "ac"), "bc"),
    ],
)
def test_settings_encode_topology_without_absent_edge_parameters(
    maker: Callable[[Conditions], dict[str, ParamLocalSetting]],
    present_edges: tuple[str, str],
    absent_edge: str,
) -> None:
    settings = maker(Conditions(temperature=37.0))

    assert {"dh_b", "ds_b", "dh_c", "ds_c", "pa", "pb", "pc"} <= settings.keys()
    for edge in present_edges:
        assert {f"dh_{edge}", f"ds_{edge}", f"k{edge}", f"k{edge[::-1]}"} <= (
            settings.keys()
        )
    assert {
        f"dh_{absent_edge}",
        f"ds_{absent_edge}",
        f"k{absent_edge}",
        f"k{absent_edge[::-1]}",
    }.isdisjoint(
        settings,
    )
    assert "37.0" in settings["kab"].expr


@pytest.mark.parametrize(
    "maker",
    [make_settings_3st_eyring_linear, make_settings_3st_eyring_fork],
)
def test_settings_require_temperature(
    maker: Callable[[Conditions], dict[str, ParamLocalSetting]],
) -> None:
    with pytest.raises(ValueError, match="The 'temperature' is None"):
        maker(Conditions(temperature=None))

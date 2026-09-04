"""Public topology behavior for the three-state Eyring models."""

from __future__ import annotations

import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_3st
from chemex.models.factory import model_factory
from chemex.models.kinetic.settings_3st_eyring import (
    calculate_kij_3st_eyring,
    calculate_kij_3st_eyring_fork,
    calculate_kij_3st_eyring_linear,
    register,
)
from chemex.models.loader import register_kinetic_settings
from chemex.models.model import ModelSpec


def _setting_signature(name: str) -> dict[str, tuple[object, ...]]:
    settings = model_factory.create(name, Conditions(temperature=25.0))
    return {
        key: (
            setting.name_setting,
            setting.value,
            setting.min,
            setting.max,
            setting.vary,
            setting.expr,
        )
        for key, setting in settings.items()
    }


def test_registered_models_expose_only_their_active_topology() -> None:
    register()

    assert "3st_eyring_triangle" not in model_factory.set
    linear = _setting_signature("3st_eyring_linear")
    compatibility = _setting_signature("3st_eyring")
    fork = _setting_signature("3st_eyring_fork")

    assert compatibility == linear
    assert {"dh_ab", "ds_ab", "dh_bc", "ds_bc", "kab", "kba", "kbc", "kcb"} <= (
        linear.keys()
    )
    assert {"dh_ac", "ds_ac", "kac", "kca"}.isdisjoint(linear)

    assert {"dh_ab", "ds_ab", "dh_ac", "ds_ac", "kab", "kba", "kac", "kca"} <= (
        fork.keys()
    )
    assert {"dh_bc", "ds_bc", "kbc", "kcb"}.isdisjoint(fork)
    assert all(setting[1] != 1.0e10 for setting in (*linear.values(), *fork.values()))


def test_linear_rates_and_populations_match_the_scientific_oracle() -> None:
    rates = calculate_kij_3st_eyring_linear(
        8000.0,
        10.0,
        12000.0,
        -5.0,
        75000.0,
        50.0,
        70000.0,
        20.0,
        25.0,
    )

    assert calculate_kij_3st_eyring is calculate_kij_3st_eyring_linear
    assert rates == pytest.approx(
        {
            "kab": 184.29414943522417,
            "kba": 1395.4513810501355,
            "kbc": 284.22888268817394,
            "kcb": 8668.464819445371,
        },
        rel=1e-12,
    )
    populations = pop_3st(
        rates["kab"],
        rates["kba"],
        0.0,
        0.0,
        rates["kbc"],
        rates["kcb"],
    )
    assert populations == pytest.approx(
        {
            "pa": 0.8799732992396345,
            "pb": 0.11621610964835893,
            "pc": 0.003810591112006477,
        },
        rel=1e-12,
    )


def test_fork_rates_and_populations_match_the_scientific_oracle() -> None:
    rates = calculate_kij_3st_eyring_fork(
        8000.0,
        10.0,
        12000.0,
        -5.0,
        75000.0,
        50.0,
        80000.0,
        30.0,
        25.0,
    )

    assert rates == pytest.approx(
        {
            "kab": 184.29414943522417,
            "kba": 1395.4513810501355,
            "kac": 2.2124682361893853,
            "kca": 510.921512184319,
        },
        rel=1e-12,
    )
    populations = pop_3st(
        rates["kab"],
        rates["kba"],
        rates["kac"],
        rates["kca"],
        0.0,
        0.0,
    )
    assert populations == pytest.approx(
        {
            "pa": 0.8799732992396347,
            "pb": 0.11621610964835895,
            "pc": 0.0038105911120064587,
        },
        rel=1e-12,
    )


def test_linear_temperature_dependence_matches_the_scientific_oracle() -> None:
    rates = calculate_kij_3st_eyring_linear(
        8000.0,
        10.0,
        12000.0,
        -5.0,
        75000.0,
        50.0,
        70000.0,
        20.0,
        50.0,
    )

    assert rates == pytest.approx(
        {
            "kab": 2074.878616024399,
            "kba": 12239.619896043598,
            "kbc": 2132.8188755332462,
            "kcb": 57413.49646390413,
        },
        rel=1e-12,
    )


@pytest.mark.parametrize(
    ("model_name", "present", "absent"),
    [
        (
            "3st_eyring_linear.rs",
            {"dh_ab", "ds_ab", "dh_bc", "ds_bc", "kab", "kba", "kbc", "kcb"},
            {"dh_ac", "ds_ac", "kac", "kca"},
        ),
        (
            "3st_eyring_fork.rs",
            {"dh_ab", "ds_ab", "dh_ac", "ds_ac", "kab", "kba", "kac", "kca"},
            {"dh_bc", "ds_bc", "kbc", "kcb"},
        ),
    ],
)
def test_residue_specific_models_preserve_topology_and_scope(
    model_name: str,
    present: set[str],
    absent: set[str],
) -> None:
    register()
    settings = model_factory.create_for_model(
        ModelSpec.from_name(model_name),
        Conditions(temperature=25.0),
    )

    assert present <= settings.keys()
    assert absent.isdisjoint(settings)
    assert all(
        setting.name_setting.spin_system_part == "g" for setting in settings.values()
    )


@pytest.mark.parametrize(
    "model_name",
    [
        "2st_eyring",
        "3st_eyring",
        "3st_eyring_linear",
        "3st_eyring_fork",
        "4st_eyring",
    ],
)
def test_real_eyring_thermodynamic_coordinates_have_approved_bounds(
    model_name: str,
) -> None:
    register_kinetic_settings()
    settings = model_factory.create(model_name, Conditions(temperature=25.0))

    thermodynamic_settings = {
        name: setting
        for name, setting in settings.items()
        if name.startswith(("dh_", "ds_"))
    }
    assert thermodynamic_settings
    for name, setting in thermodynamic_settings.items():
        expected = (-2.0e5, 2.0e5) if name.startswith("dh_") else (-5.0e2, 5.0e2)
        assert (setting.min, setting.max) == expected
        assert setting.value is not None
        assert setting.min <= setting.value <= setting.max
        assert not setting.expr

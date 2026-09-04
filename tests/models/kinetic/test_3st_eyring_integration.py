"""Registration and expression integration for three-state Eyring models."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemex.configuration.conditions import Conditions
from chemex.configuration.parameters import DefaultSetting
from chemex.models.factory import model_factory
from chemex.models.kinetic.settings_3st_eyring import register
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.parameters.name import ParamName
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession


@pytest.mark.parametrize(
    ("model_name", "active_rates", "absent_rates"),
    [
        (
            "3st_eyring",
            {"kab", "kba", "kbc", "kcb"},
            {"kac", "kca"},
        ),
        (
            "3st_eyring_linear",
            {"kab", "kba", "kbc", "kcb"},
            {"kac", "kca"},
        ),
        (
            "3st_eyring_fork",
            {"kab", "kba", "kac", "kca"},
            {"kbc", "kcb"},
        ),
    ],
)
def test_registered_topology_builds_a_profile_parameterization(
    model_name: str,
    active_rates: set[str],
    absent_rates: set[str],
) -> None:
    register()
    session = AnalysisSession.create()
    session.set_model(model_name)
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    config = SimpleNamespace(
        conditions=Conditions(temperature=25.0),
        to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
    )

    name_map = session.parameter_factory.create_parameters(
        config,  # ty: ignore[invalid-argument-type]
        basis=basis,
        spin_system=SpinSystem.from_name("G23N-HN"),
    )

    assert active_rates <= name_map.keys()
    assert absent_rates.isdisjoint(name_map)
    assert {"pa", "pb", "pc"} <= name_map.keys()


def test_compatibility_name_and_explicit_linear_name_have_same_model_behavior() -> None:
    register()
    conditions = Conditions(temperature=25.0)
    compatibility = model_factory.create_for_model(
        ModelSpec.from_name("3st_eyring.rs"),
        conditions,
    )
    explicit = model_factory.create_for_model(
        ModelSpec.from_name("3st_eyring_linear.rs"),
        conditions,
    )

    assert compatibility.keys() == explicit.keys()
    for name in compatibility:
        left = compatibility[name]
        right = explicit[name]
        assert left.name_setting == right.name_setting
        assert left.value == right.value
        assert left.min == right.min
        assert left.max == right.max
        assert left.vary == right.vary
        assert left.expr == right.expr


@pytest.mark.parametrize(
    ("model_name", "arguments"),
    [
        ("3st_eyring_linear", "{kab},{kba},0.0,0.0,{kbc},{kcb}"),
        ("3st_eyring_fork", "{kab},{kba},{kac},{kca},0.0,0.0"),
    ],
)
def test_population_constraints_encode_the_absent_edge_as_exact_zero(
    model_name: str,
    arguments: str,
) -> None:
    register()
    settings = model_factory.create(model_name, Conditions(temperature=25.0))

    for population in ("pa", "pb", "pc"):
        expression = settings[population].expr
        assert "pop_3st" in expression
        assert arguments in expression


def test_parameter_file_bounds_override_eyring_defaults() -> None:
    register()
    session = AnalysisSession.create()
    session.set_model("3st_eyring")
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    config = SimpleNamespace(
        conditions=Conditions(temperature=25.0),
        to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
    )
    session.parameter_factory.create_parameters(
        config,  # ty: ignore[invalid-argument-type]
        basis=basis,
        spin_system=SpinSystem.from_name("G23N-HN"),
    )
    (dh_b_id,) = session.parameters.database.get_matching_ids(
        ParamName.from_section("dh_b"),
    )
    dh_b = session.parameters.get_parameters([dh_b_id])[dh_b_id]

    session.parameters.set_defaults(
        [(dh_b.param_name, DefaultSetting(1000.0, -1234.0, 5678.0))],
    )

    updated = session.parameters.get_parameters([dh_b_id])[dh_b_id]
    assert (updated.value, updated.min, updated.max) == (1000.0, -1234.0, 5678.0)

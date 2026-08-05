from __future__ import annotations

import pytest

from chemex.parameters.definitions import (
    ParameterDefinition,
    ParameterDefinitionConflictError,
    ParameterDefinitions,
    ParameterDefinitionsBuilder,
    ParameterDefinitionsCollector,
    SealedParameterModelError,
)
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting
from chemex.parameters.spin_system import SpinSystem


def _make_setting(
    name: str,
    *,
    spin_system: str = "",
    value: float = 1.0,
    min: float = 0.0,
    max: float = 100.0,
    vary: bool = False,
    expr: str = "",
) -> ParamSetting:
    param_name = ParamName(name, SpinSystem.from_name(spin_system))
    return ParamSetting(
        param_name=param_name,
        value=value,
        min=min,
        max=max,
        vary=vary,
        expr=expr,
    )


def test_builder_seals_deterministically_ordered_definitions() -> None:
    builder = ParameterDefinitionsBuilder()
    setting_b = _make_setting("kb", spin_system="G2N")
    setting_a = _make_setting("ka", spin_system="G2N")

    builder.add(setting_b.id_, setting_b, source="experiment-b")
    builder.add(setting_a.id_, setting_a, source="experiment-a")

    definitions = builder.seal()

    assert isinstance(definitions, ParameterDefinitions)
    assert list(definitions.ids) == sorted(definitions.ids)
    assert set(definitions) == {setting_a.id_, setting_b.id_}


def test_identical_contributions_do_not_conflict() -> None:
    builder = ParameterDefinitionsBuilder()
    first = _make_setting("r2", spin_system="G2N", value=10.0, expr="")
    second = _make_setting("r2", spin_system="G2N", value=10.0, expr="")

    builder.add(first.id_, first, source="profile-1")
    builder.add(second.id_, second, source="profile-2")

    definitions = builder.seal()

    assert definitions[first.id_].value == 10.0


def test_conflicting_contributions_raise() -> None:
    builder = ParameterDefinitionsBuilder()
    first = _make_setting("r2", spin_system="G2N", value=10.0)
    second = _make_setting("r2", spin_system="G2N", value=12.0)

    builder.add(first.id_, first, source="profile-1")

    with pytest.raises(ParameterDefinitionConflictError) as exc_info:
        builder.add(second.id_, second, source="profile-2")

    error = exc_info.value
    assert error.param_id == first.id_
    assert error.field_name == "value"
    assert error.values == (10.0, 12.0)
    assert error.contributors == ("profile-1", "profile-2")


def test_conflicting_expr_raises_with_correct_field() -> None:
    builder = ParameterDefinitionsBuilder()
    first = _make_setting("r2a", spin_system="G2N", expr="{r2}-{r1}")
    second = _make_setting("r2a", spin_system="G2N", expr="{r2}+{r1}")

    builder.add(first.id_, first, source="profile-1")

    with pytest.raises(ParameterDefinitionConflictError) as exc_info:
        builder.add(second.id_, second, source="profile-2")

    assert exc_info.value.field_name == "expr"


def test_sealed_definitions_cannot_be_mutated() -> None:
    builder = ParameterDefinitionsBuilder()
    setting = _make_setting("pb", spin_system="")
    builder.add(setting.id_, setting, source="profile-1")

    definitions = builder.seal()

    with pytest.raises(SealedParameterModelError):
        definitions._definitions = {}  # type: ignore[attr-defined]


def test_builder_cannot_be_used_after_sealing() -> None:
    builder = ParameterDefinitionsBuilder()
    setting = _make_setting("pb", spin_system="")
    builder.add(setting.id_, setting, source="profile-1")
    builder.seal()

    with pytest.raises(SealedParameterModelError):
        builder.add(setting.id_, setting, source="profile-2")

    with pytest.raises(SealedParameterModelError):
        builder.seal()


def test_parameter_definition_from_setting_preserves_scientific_fields() -> None:
    setting = _make_setting(
        "kex",
        spin_system="G2N",
        value=400.0,
        min=0.0,
        max=1.0e4,
        vary=True,
        expr="",
    )

    definition = ParameterDefinition.from_setting(setting.id_, setting)

    assert definition.param_id == setting.id_
    assert definition.param_name == setting.param_name
    assert definition.value == setting.value
    assert definition.min == setting.min
    assert definition.max == setting.max
    assert definition.vary == setting.vary
    assert definition.expr == setting.expr


def test_collector_seals_ordinary_and_model_free_independently() -> None:
    collector = ParameterDefinitionsCollector()
    ordinary = {"__PB": _make_setting("pb")}
    model_free = {"__TAUC_A": _make_setting("tauc_a", value=4.0)}

    collector.collect(ordinary, model_free, source="profile-1")

    ordinary_defs, model_free_defs = collector.seal()

    assert set(ordinary_defs) == {"__PB"}
    assert set(model_free_defs) == {"__TAUC_A"}

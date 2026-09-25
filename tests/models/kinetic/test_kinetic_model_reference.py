"""Keep the public kinetic reference aligned with runtime registration/settings."""

from __future__ import annotations

import re
from pathlib import Path
from types import SimpleNamespace

import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.models.loader import register_kinetic_settings
from chemex.models.model import ModelSelectionError, ModelSpec
from chemex.nmr.basis import Basis
from chemex.parameters.setting import ParamLocalSetting
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

PAGE = Path(__file__).parents[3] / "website/docs/user_guide/fitting/kinetic_models.md"
CONDITIONS = Conditions(temperature=25.0, p_total=1e-3, l_total=2e-3, d2o=0.2)
ALIASES = {
    "3st_triangle": "3st",
    "3st_eyring": "3st_eyring_linear",
}
DIMENSIONS = {"temperature": "T", "p_total": "P", "l_total": "L", "d2o": "D"}


def setup_module() -> None:
    register_kinetic_settings()


def _block(page: str, marker: str) -> str:
    start = f"<!-- {marker}:start -->"
    end = f"<!-- {marker}:end -->"
    assert page.count(start) == page.count(end) == 1
    return page.split(start, 1)[1].split(end, 1)[0].strip()


def _signature(setting: ParamLocalSetting) -> tuple[object, ...]:
    return (
        setting.name_setting,
        setting.value,
        setting.min,
        setting.max,
        setting.vary,
        setting.supports_estimation,
        setting.report_only,
        setting.expr,
    )


def _scope(setting: ParamLocalSetting) -> str:
    scope = ",".join(DIMENSIONS[part] for part in setting.name_setting.conditions_part)
    return f"{scope or 'all'}; {'group' if setting.name_setting.spin_system_part == 'g' else 'global'}"


def _input_cell(setting: ParamLocalSetting) -> str:
    role = "fit" if setting.vary else "fixed"
    if setting.supports_estimation:
        role += ", estimate"
    return (
        f"`{setting.name_setting.name.upper()}` {setting.value} "
        f"[{setting.min}, {setting.max}] ({role}; {_scope(setting)})"
    )


def _output_cell(settings: dict[str, ParamLocalSetting], *, report_only: bool) -> str:
    by_scope: dict[str, list[str]] = {}
    for setting in settings.values():
        if setting.expr and setting.report_only == report_only:
            by_scope.setdefault(_scope(setting), []).append(
                f"`{setting.name_setting.name.upper()}`"
            )
    return (
        "; ".join(f"{', '.join(names)} ({scope})" for scope, names in by_scope.items())
        or "—"
    )


def render_runtime_parameter_inventory() -> str:
    """Project parameter facts directly from registered model settings."""
    lines = [
        "| Model | Independent input: default [bounds] (role; scope) | Derived | Report-only |",
        "| --- | --- | --- | --- |",
    ]
    for name in sorted(model_factory.set):
        settings = model_factory.create(name, CONDITIONS)
        inputs = (
            "; ".join(
                _input_cell(setting)
                for setting in settings.values()
                if not setting.expr
            )
            or "—"
        )
        derived = _output_cell(settings, report_only=False)
        report_only = _output_cell(settings, report_only=True)
        lines.append(f"| `{name}` | {inputs} | {derived} | {report_only} |")
    return "\n".join(lines)


def test_public_names_and_aliases_match_runtime() -> None:
    table = _block(PAGE.read_text(encoding="utf-8"), "kinetic-model-names")
    rows = re.findall(r"^\| `([^`]+)` \| .+ \| (.+) \|$", table, re.MULTILINE)
    assert len(model_factory.set) == 34
    assert "2st_rs" not in model_factory.set
    with pytest.raises(ModelSelectionError):
        ModelSpec.from_name("2st_rs")
    assert len(rows) == len(model_factory.set)
    assert {name for name, _ in rows} == model_factory.set
    documented_aliases = {
        name: canonical.strip("`") for name, canonical in rows if canonical != "—"
    }
    assert documented_aliases == ALIASES

    for alias, canonical in ALIASES.items():
        selected = ModelSpec.from_name(alias)
        assert (
            selected.name == alias
        )  # selected spelling remains in identity/provenance
        assert selected.identity.startswith(f"{alias}|")
        actual = model_factory.create_for_model(selected, CONDITIONS)
        expected = model_factory.create_for_model(
            ModelSpec.from_name(canonical), CONDITIONS
        )
        assert actual.keys() == expected.keys()
        assert {key: _signature(value) for key, value in actual.items()} == {
            key: _signature(value) for key, value in expected.items()
        }


def test_public_parameter_inventory_matches_runtime_settings() -> None:
    actual = _block(PAGE.read_text(encoding="utf-8"), "kinetic-parameters")
    assert actual == render_runtime_parameter_inventory()


def test_modifier_order_keeps_the_same_model_settings() -> None:
    forward = ModelSpec.from_name("2st.rs.mf.tc")
    reversed_order = ModelSpec.from_name("2st.tc.rs.mf")
    assert forward == reversed_order
    assert {
        key: _signature(setting)
        for key, setting in model_factory.create_for_model(forward, CONDITIONS).items()
    } == {
        key: _signature(setting)
        for key, setting in model_factory.create_for_model(
            reversed_order, CONDITIONS
        ).items()
    }


@pytest.mark.parametrize(
    "model_name", ["2st", "2st_binding", "3st_monomer_dimer_trimer"]
)
def test_representative_inventory_facts_survive_sealed_construction(
    model_name: str,
) -> None:
    settings = model_factory.create(model_name, CONDITIONS)
    session = AnalysisSession.create()
    session.set_model(model_name)
    basis = Basis(type="iz", spin_system="nh", model=session.model.spec)
    config = SimpleNamespace(
        conditions=CONDITIONS,
        to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
    )
    spin = SpinSystem.from_name("G23N-HN")
    session.parameter_factory.create_parameters(
        config,  # ty: ignore[invalid-argument-type]
        basis=basis,
        spin_system=spin,
    )
    assert session.parameter_factory.try_seal_definitions()
    definitions = session.parameter_factory.sealed_definitions
    assert definitions is not None
    for setting in settings.values():
        param_name = setting.name_setting.get_param_name(spin, CONDITIONS)
        definition = definitions[param_name.id_]
        assert definition.name == setting.name_setting.name.upper()
        assert definition.default_value == setting.value
        assert definition.lower_bound == setting.min
        assert definition.upper_bound == setting.max
        assert definition.spin_system_name == param_name.spin_system.name


if __name__ == "__main__":
    # Maintainer aid: print the checked appendix after a kinetic settings change.
    register_kinetic_settings()
    print(render_runtime_parameter_inventory())

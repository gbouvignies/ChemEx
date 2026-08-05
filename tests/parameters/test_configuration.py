from __future__ import annotations

import pytest

from chemex.configuration.parameters import DefaultSetting
from chemex.parameters.configuration import (
    ParameterConfiguration,
    ParameterConfigurationBuilder,
    build_parameter_configuration,
)
from chemex.parameters.database import ParameterCatalog
from chemex.parameters.definitions import (
    ParameterDefinitionsBuilder,
    SealedParameterModelError,
)
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting
from chemex.parameters.spin_system import SpinSystem


def _param_name(name: str, spin_system: str = "") -> ParamName:
    return ParamName(name, SpinSystem.from_name(spin_system))


def _setting(
    name: str, spin_system: str, value: float, *, min_: float = 0.0
) -> ParamSetting:
    return ParamSetting(
        param_name=_param_name(name, spin_system),
        value=value,
        min=min_,
        max=1.0e4,
    )


def _build_definitions(settings: dict[str, ParamSetting]) -> tuple:
    builder = ParameterDefinitionsBuilder()
    builder.add_many(settings, source="test")
    return builder.seal()


def test_configuration_defaults_mirror_definitions_when_no_toml_applied() -> None:
    settings = {
        "__PB": _setting("pb", "", 0.05),
        "__KEX_AB": _setting("kex_ab", "", 200.0),
    }
    definitions = _build_definitions(settings)

    configuration = ParameterConfigurationBuilder(definitions).seal()

    assert isinstance(configuration, ParameterConfiguration)
    assert list(configuration.ids) == sorted(configuration.ids)
    for param_id, definition in definitions.items():
        entry = configuration[param_id]
        assert entry.value == definition.value
        assert entry.min == definition.min
        assert entry.max == definition.max
        assert entry.overridden is False
        assert entry.brute_step is None


def test_apply_defaults_matches_legacy_catalog_precedence() -> None:
    settings = {
        "__PB": _setting("pb", "", 0.05),
        "__CS_A_10N": _setting("cs_a", "10N", 100.0),
        "__CS_A_1N": _setting("cs_a", "1N", 110.0),
    }
    definitions = _build_definitions(settings)

    defaults = [
        (_param_name("cs_a"), DefaultSetting(value=999.0, min=-1.0, max=1.0)),
        (
            _param_name("cs_a", "10N"),
            DefaultSetting(value=123.456, min=0.0, brute_step=0.5),
        ),
    ]

    # Legacy-authoritative path: apply the same defaults to a real catalog
    # built from fresh copies of the same settings.
    legacy_catalog = ParameterCatalog()
    legacy_catalog.add_multiple(
        {
            "__PB": _setting("pb", "", 0.05),
            "__CS_A_10N": _setting("cs_a", "10N", 100.0),
            "__CS_A_1N": _setting("cs_a", "1N", 110.0),
        },
    )
    legacy_catalog.set_defaults(defaults)

    native_builder = ParameterConfigurationBuilder(definitions)
    native_builder.apply_defaults(defaults)
    configuration = native_builder.seal()

    for param_id in definitions:
        legacy_setting = legacy_catalog.get_parameters([param_id])[param_id]
        entry = configuration[param_id]
        assert entry.value == legacy_setting.value
        assert entry.min == legacy_setting.min
        assert entry.max == legacy_setting.max

    assert configuration["__CS_A_10N"].overridden is True
    assert configuration["__CS_A_10N"].brute_step == pytest.approx(0.5)
    assert configuration["__CS_A_1N"].overridden is True
    assert configuration["__PB"].overridden is False


def test_apply_defaults_reproduces_legacy_brute_step_quirk() -> None:
    """`brute_step` is only ever applied together with an explicit `min`.

    This mirrors a legacy quirk in `ParamSetting.set` (see
    `chemex.parameters.setting`): `brute_step` is gated on `default_setting.min
    is not None`, not on `brute_step` itself. Preserving this quirk keeps the
    native configuration bit-for-bit compatible with the legacy-authoritative
    path.
    """
    definitions = _build_definitions({"__PB": _setting("pb", "", 0.05)})
    defaults = [(_param_name("pb"), DefaultSetting(value=0.1, brute_step=0.01))]

    builder = ParameterConfigurationBuilder(definitions)
    builder.apply_defaults(defaults)
    configuration = builder.seal()

    assert configuration["__PB"].value == pytest.approx(0.1)
    assert configuration["__PB"].brute_step is None


def test_sealed_configuration_cannot_be_mutated() -> None:
    definitions = _build_definitions({"__PB": _setting("pb", "", 0.05)})
    configuration = ParameterConfigurationBuilder(definitions).seal()

    with pytest.raises(SealedParameterModelError):
        configuration._entries = {}  # type: ignore[attr-defined]


def test_configuration_builder_cannot_be_used_after_sealing() -> None:
    definitions = _build_definitions({"__PB": _setting("pb", "", 0.05)})
    builder = ParameterConfigurationBuilder(definitions)
    builder.seal()

    with pytest.raises(SealedParameterModelError):
        builder.apply_defaults([])

    with pytest.raises(SealedParameterModelError):
        builder.overwrite_values({})

    with pytest.raises(SealedParameterModelError):
        builder.seal()


def test_build_parameter_configuration_reads_effective_values_from_legacy_catalog() -> (
    None
):
    settings = {
        "__PB": _setting("pb", "", 0.05),
        "__KEX_AB": _setting("kex_ab", "", 200.0),
    }
    definitions = _build_definitions(settings)

    catalog = ParameterCatalog()
    catalog.add_multiple(
        {
            "__PB": _setting("pb", "", 0.05),
            "__KEX_AB": _setting("kex_ab", "", 200.0),
        },
    )
    defaults = [(_param_name("kex_ab"), DefaultSetting(value=555.0))]
    catalog.set_defaults(defaults)

    configuration = build_parameter_configuration(definitions, catalog, defaults)

    assert configuration["__KEX_AB"].value == pytest.approx(555.0)
    assert configuration["__KEX_AB"].value == catalog.get_value("__KEX_AB")
    assert configuration["__PB"].value == catalog.get_value("__PB")
    assert configuration["__KEX_AB"].overridden is True
    assert configuration["__PB"].overridden is False

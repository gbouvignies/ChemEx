from __future__ import annotations

"""Factories for creating parameter settings."""

from builtins import set as builtins_set
from collections.abc import Callable
from dataclasses import replace
from typing import TYPE_CHECKING, ClassVar

from chemex.configuration.conditions import Conditions
from chemex.parameters.setting import ParamLocalSetting

if TYPE_CHECKING:
    from chemex.models.model import ModelSpec

SettingsType = dict[str, ParamLocalSetting]
SettingMakerType = Callable[[Conditions], SettingsType]


def make_residue_specific(settings: SettingsType) -> SettingsType:
    for setting in settings.values():
        name_setting = setting.name_setting
        if not name_setting.allow_residue_specific:
            continue
        if name_setting.spin_system_part == "":
            setting.name_setting = replace(name_setting, spin_system_part="g")
    return settings


class Factory:
    """Factory for creating all parts of an experiment."""

    setting_makers_registry: ClassVar[dict[str, SettingMakerType]] = {}

    def register(self, name: str, setting_maker: SettingMakerType) -> None:
        """Register a new setting maker."""
        self.setting_makers_registry[name] = setting_maker

    def create(self, name: str, conditions: Conditions) -> SettingsType:
        try:
            return self.setting_makers_registry[name](conditions)
        except KeyError:
            msg = f"Unknown type {name!r}"
            raise ValueError(msg) from None

    def create_for_model(
        self,
        model: ModelSpec,
        conditions: Conditions,
    ) -> SettingsType:
        settings = self.create(model.name, conditions)
        if model.residue_specific:
            return make_residue_specific(settings)
        return settings

    @property
    def set(self) -> builtins_set[str]:
        return set(self.setting_makers_registry)


model_factory = Factory()

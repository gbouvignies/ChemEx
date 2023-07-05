"""Factories for creating parameter settings."""
from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

if TYPE_CHECKING:
    from collections.abc import Callable

    from chemex.configuration.conditions import Conditions
    from chemex.parameters.setting import ParamLocalSetting

    SettingsType = dict[str, ParamLocalSetting]
    SettingMakerType = Callable[[Conditions], SettingsType]


class Factory:
    """Factory for creating all parts of an experiment."""

    setting_makers_registry: ClassVar[dict[str, SettingMakerType]] = {}

    def register(self, name: str, setting_maker: SettingMakerType):
        """Register a new setting maker."""
        self.setting_makers_registry[name] = setting_maker

    def create(self, name: str, conditions: Conditions) -> SettingsType:
        try:
            return self.setting_makers_registry[name](conditions)
        except KeyError:
            msg = f"Unknown type {name!r}"
            raise ValueError(msg) from None

    @property
    def set(self) -> set[str]:
        return set(self.setting_makers_registry)


model_factory = Factory()

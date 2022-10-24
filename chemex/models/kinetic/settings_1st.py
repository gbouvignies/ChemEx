from __future__ import annotations

from typing import TYPE_CHECKING

from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting
from chemex.parameters.setting import ParamLocalSetting

if TYPE_CHECKING:
    from chemex.configuration.conditions import Conditions

NAME = "1st"


def make_settings_1st(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {"pa": ParamLocalSetting(name_setting=NameSetting("pa"), expr="1.0")}


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_1st)

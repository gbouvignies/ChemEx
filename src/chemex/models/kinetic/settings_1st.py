from __future__ import annotations

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting

NAME = "1st"


def make_settings_1st(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {"pa": ParamLocalSetting(name_setting=NameSetting("pa"), expr="1.0")}


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_1st)

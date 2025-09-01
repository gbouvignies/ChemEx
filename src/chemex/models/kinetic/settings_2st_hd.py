from __future__ import annotations

from chemex.configuration.conditions import Conditions
from chemex.models.constraints import pop_2st
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting
from chemex.parameters.userfunctions import user_function_registry

NAME = "2st_hd"


def make_settings_2st_hd(conditions: Conditions) -> dict[str, ParamLocalSetting]:
    d2o = conditions.d2o if conditions.d2o is not None else 0.1
    return {
        "d2o": ParamLocalSetting(
            name_setting=NameSetting("d2o", "", ("d2o",)),
            value=d2o,
            min=0.0,
            max=1.0,
        ),
        "kdh": ParamLocalSetting(
            name_setting=NameSetting("kdh", "g", ("temperature",)),
            value=1.0,
            min=0.0,
            vary=True,
        ),
        "phi": ParamLocalSetting(
            name_setting=NameSetting("phi", "g", ("temperature",)),
            value=1.1,
            min=0.75,
            max=1.50,
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "g", ("temperature", "d2o")),
            min=0.0,
            expr="{d2o} * {kdh} * {phi}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "g", ("temperature", "d2o")),
            expr="(1.0 - {d2o}) * {kdh}",
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "g", ("temperature", "d2o")),
            expr="pop_2st({kab}, {kba})['pa']",
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "g", ("temperature", "d2o")),
            expr="pop_2st({kab}, {kba})['pb']",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_2st_hd)
    user_function_registry.register(name=NAME, user_functions={"pop_2st": pop_2st})

from __future__ import annotations

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting

NAME = "2st"

TPL = ("temperature", "p_total", "l_total")


def make_settings_2st(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting("kex_ab", "", TPL),
            value=200.0,
            min=0.0,
            vary=True,
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "", TPL),
            value=0.05,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "", TPL),
            expr="1.0 - {pb}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "", TPL),
            expr="{kex_ab} * {pb}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "", TPL),
            expr="{kex_ab} * {pa}",
        ),
    }


def register() -> None:
    model_factory.register(name=NAME, setting_maker=make_settings_2st)

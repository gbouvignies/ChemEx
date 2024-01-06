from __future__ import annotations

from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.parameters.setting import NameSetting, ParamLocalSetting

TPL = ("temperature", "p_total", "l_total")


def make_settings_2st_rs(_conditions: Conditions) -> dict[str, ParamLocalSetting]:
    return {
        "kex_ab": ParamLocalSetting(
            name_setting=NameSetting("kex_ab", "g", TPL),
            value=200.0,
            min=0.0,
            vary=True,
        ),
        "pb": ParamLocalSetting(
            name_setting=NameSetting("pb", "g", TPL),
            value=0.05,
            min=0.0,
            max=1.0,
            vary=True,
        ),
        "pa": ParamLocalSetting(
            name_setting=NameSetting("pa", "g", TPL),
            expr="1.0 - {pb}",
        ),
        "kab": ParamLocalSetting(
            name_setting=NameSetting("kab", "g", TPL),
            expr="{kex_ab} * {pb}",
        ),
        "kba": ParamLocalSetting(
            name_setting=NameSetting("kba", "g", TPL),
            expr="{kex_ab} * {pa}",
        ),
    }


def register() -> None:
    model_factory.register(name="2st_rs", setting_maker=make_settings_2st_rs)

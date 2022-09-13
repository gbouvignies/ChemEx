from __future__ import annotations

from copy import deepcopy

from chemex.configuration.conditions import Conditions
from chemex.model import model
from chemex.nmr.constants import J_COUPLINGS
from chemex.nmr.liouvillian import Basis
from chemex.nmr.rates import get_model_free_expressions
from chemex.parameters.setting import LocalSettings
from chemex.parameters.setting import NameSetting
from chemex.parameters.setting import ParamLocalSetting


def _update_expr_for_proton_exchange(
    settings: dict[str, ParamLocalSetting], state: str, basis: Basis
) -> None:
    # As 1HN sites exchange with water, the "expr" should be modified accordingly
    # for spin system including an amide proton
    if basis.spin_system == "hn":
        expr = f"{{r2_s_{state}}} + {{r1a_is_{state}}} - {{r1_s_{state}}}"
        settings[f"r2a_s_{state}"].expr = expr
    elif basis.spin_system == "nh":
        expr = f"{{r2_i_{state}}} + {{r1a_is_{state}}} - {{r1_i_{state}}}"
        settings[f"r2a_i_{state}"].expr = expr
        settings[f"r1a_is_{state}"].expr = ""


def _add_dw_param_settings(settings: dict[str, ParamLocalSetting], state: str) -> None:
    settings[f"dw_i_a{state}"] = ParamLocalSetting(
        name_setting=NameSetting(f"dw_a{state}", "i"),
        value=0.0,
        vary=True,
    )
    settings[f"dw_s_a{state}"] = ParamLocalSetting(
        name_setting=NameSetting(f"dw_a{state}", "s"),
        value=0.0,
        vary=True,
    )
    settings[f"cs_i_{state}"].expr = f"{{cs_i_a}} + {{dw_i_a{state}}}"
    settings[f"cs_s_{state}"].expr = f"{{cs_s_a}} + {{dw_s_a{state}}}"


def _set_equal_to_a(settings: dict[str, ParamLocalSetting]) -> None:
    for name, setting in settings.items():
        if not setting.expr:
            setting.expr = f"{{{name[:-1]}a}}"


def create_base_param_settings(
    basis: Basis, state: str
) -> dict[str, ParamLocalSetting]:
    ext = basis.extension
    settings = {
        f"tauc_{state}": ParamLocalSetting(
            name_setting=NameSetting(f"tauc_{state}", "g", ("temperature",)),
            value=4.0,
            min=0.0,
        ),
        f"s2_{state}": ParamLocalSetting(
            name_setting=NameSetting(f"s2_{state}", "g", ("temperature",)),
            value=0.9,
            min=0.0,
        ),
        f"khh_{state}": ParamLocalSetting(
            name_setting=NameSetting(f"khh_{state}", "g", ("temperature",)), value=0.0
        ),
        f"r2_i_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"r2{ext}_{state}", "i", ("temperature", "h_larmor_frq")
            ),
            value=10.0,
            min=0.0,
        ),
        f"r2_s_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"r2_{state}", "s", ("temperature", "h_larmor_frq")
            ),
            value=10.0,
            min=0.0,
        ),
        f"r1_i_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"r1{ext}_{state}", "i", ("temperature", "h_larmor_frq")
            ),
            value=1.5,
            min=0.0,
        ),
        f"r1_s_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"r1_{state}", "s", ("temperature", "h_larmor_frq")
            ),
            value=1.5,
            min=0.0,
        ),
        f"r2a_i_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"r2a{ext}_{state}", "i", ("temperature", "h_larmor_frq")
            ),
            value=12.0,
            min=0.0,
            expr=f"{{r2_i_{state}}} - {{r1_s_{state}}}",
        ),
        f"r2a_s_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"r2a_{state}", "s", ("temperature", "h_larmor_frq")
            ),
            value=12.0,
            min=0.0,
            expr=f"{{r2_s_{state}}} - {{r1_i_{state}}}",
        ),
        f"r2mq_is_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"r2mq_{state}", "is", ("temperature", "h_larmor_frq")
            ),
            value=15.0,
            min=0.0,
        ),
        f"r1a_is_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"r1a{ext}_{state}", "is", ("temperature", "h_larmor_frq")
            ),
            value=4.0,
            min=0.0,
            expr=f"{{r1_i_{state}}} + {{r1_s_{state}}}",
        ),
        f"etaxy_i_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"etaxy{ext}_{state}", "i", ("temperature", "h_larmor_frq")
            ),
            value=0.0,
        ),
        f"etaxy_s_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"etaxy_{state}", "s", ("temperature", "h_larmor_frq")
            ),
            value=0.0,
        ),
        f"etaz_i_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"etaz{ext}_{state}", "i", ("temperature", "h_larmor_frq")
            ),
            value=0.0,
        ),
        f"etaz_s_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"etaz_{state}", "s", ("temperature", "h_larmor_frq")
            ),
            value=0.0,
        ),
        f"sigma_is_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"sigma_{state}", "is", ("temperature", "h_larmor_frq")
            ),
            value=0.0,
        ),
        f"mu_is_{state}": ParamLocalSetting(
            name_setting=NameSetting(
                f"mu_{state}", "is", ("temperature", "h_larmor_frq")
            ),
            value=0.0,
        ),
        f"cs_i_{state}": ParamLocalSetting(
            name_setting=NameSetting(f"cs_{state}", "i"), value=0.0
        ),
        f"cs_s_{state}": ParamLocalSetting(
            name_setting=NameSetting(f"cs_{state}", "s"), value=0.0
        ),
        f"j_is_{state}": ParamLocalSetting(
            name_setting=NameSetting(f"j_{state}", "is"),
            value=J_COUPLINGS.get(basis.spin_system, 0.0),
        ),
        f"d_{state}": ParamLocalSetting(
            name_setting=NameSetting(f"d_{state}", "", ("temperature",)), value=0.0
        ),
    }

    _update_expr_for_proton_exchange(settings, state, basis)

    if state != "a":
        _set_equal_to_a(settings)
        _add_dw_param_settings(settings, state)

    return settings


def _build_model_free_settings(
    settings: LocalSettings, basis: Basis, conditions: Conditions
) -> LocalSettings:
    model_free_expressions = get_model_free_expressions(basis, conditions)
    settings_mf = deepcopy(settings)
    for name, setting in settings_mf.items():
        if name in model_free_expressions:
            setting.expr = model_free_expressions[name]

    return settings_mf


def _select_relevant_settings(
    all_settings: LocalSettings, basis: Basis
) -> LocalSettings:

    selection = set()
    for name in set(all_settings) & set(basis.matrices):
        selection.update({name, *all_settings[name].dependencies})

    return {name: all_settings[name] for name in selection}


def build_spin_param_settings(
    basis: Basis, conditions: Conditions
) -> tuple[LocalSettings, LocalSettings]:

    all_settings: LocalSettings = {}
    for state in model.states:
        all_settings.update(create_base_param_settings(basis, state))

    all_settings_mf = _build_model_free_settings(all_settings, basis, conditions)

    settings = _select_relevant_settings(all_settings, basis)
    settings_mf = _select_relevant_settings(all_settings_mf, basis)

    return settings, settings_mf

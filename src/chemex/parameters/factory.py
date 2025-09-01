from __future__ import annotations

from functools import cache
from typing import Any

from chemex.configuration.base import ExperimentConfiguration
from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.models.model import model
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.parameters import database
from chemex.parameters.setting import LocalSettings, Parameters, ParamSetting
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spins import build_spin_param_settings

ConfigConditionsType = ExperimentConfiguration[Any, Conditions, Any]


@cache
def _build_settings(
    basis: Basis,
    conditions: Conditions,
) -> tuple[LocalSettings, LocalSettings]:
    settings_spins, settings_spins_mf = build_spin_param_settings(basis, conditions)
    settings_kinetics = model_factory.create(model.name, conditions)
    settings = settings_kinetics | settings_spins
    settings_mf = settings_kinetics | settings_spins_mf
    return settings, settings_mf


def _set_to_fit(
    parameters: Parameters,
    name_map: dict[str, str],
    fitted: list[str],
) -> None:
    pool = set(name_map)
    for fitted_name in fitted:
        selection = {name for name in pool if name.startswith(fitted_name)}
        for name in selection:
            param_id = name_map[name]
            parameters[param_id].vary = True
            parameters[param_id].expr = ""
        pool -= selection


def _build_parameters(
    settings: LocalSettings,
    spin_system: SpinSystem,
    conditions: Conditions,
) -> tuple[dict[str, str], Parameters]:
    param_names = {
        local_name: setting.name_setting.get_param_name(spin_system, conditions)
        for local_name, setting in settings.items()
    }

    name_map = {
        local_name: param_name.id_ for local_name, param_name in param_names.items()
    }

    parameters: Parameters = {}
    for local_name, setting in settings.items():
        param_id = name_map[local_name]
        param_name = param_names[local_name]
        expression = setting.expr.format_map(name_map)
        parameters[param_id] = ParamSetting(
            param_name=param_name,
            value=setting.value,
            min=setting.min,
            max=setting.max,
            vary=setting.vary,
            expr=expression,
        )

    return name_map, parameters


def create_parameters(
    config: ConfigConditionsType,
    liouvillian: LiouvillianIS,
) -> dict[str, str]:
    # A copy is done because the output of '_build_settings' is cached
    settings, settings_mf = _build_settings(liouvillian.basis, config.conditions)

    name_map, parameters = _build_parameters(
        settings,
        liouvillian.spin_system,
        config.conditions,
    )

    name_map_mf, parameters_mf = _build_parameters(
        settings_mf,
        liouvillian.spin_system,
        config.conditions,
    )

    _set_to_fit(parameters, name_map, config.to_be_fitted.rates)
    _set_to_fit(parameters_mf, name_map_mf, config.to_be_fitted.model_free)

    database.add_parameters(parameters)
    database.add_parameters_mf(parameters_mf)

    selection = set(name_map) & liouvillian.basis.required_names

    return {local_name: name_map[local_name] for local_name in selection}

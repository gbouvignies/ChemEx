from __future__ import annotations

import dataclasses as dc
import itertools as it
import re
from typing import TYPE_CHECKING

import lmfit as lf

import chemex.nmr.rates as cnr
import chemex.parameters.kinetics as cpk
import chemex.parameters.liouvillian as cpl
import chemex.parameters.name as cpn
import chemex.parameters.settings as cps

if TYPE_CHECKING:
    from chemex.containers.experiment import Experiments
    from chemex.parameters.name import ParamName


def merge(params_list):
    params_ = {}
    for params in params_list:
        for name, param in params.items():
            if name in params_ and params_[name].vary:
                continue
            params_[name] = param
    params = lf.Parameters(usersyms=cnr.rate_functions)
    params.add_many(*params_.values())
    return params


def create_params(experiments: Experiments, defaults: list[tuple[ParamName, dict]]):
    params_mf = experiments.params_mf
    params = experiments.params

    cps.set_values(params_mf, defaults)

    # Set parameters values using the model-free values
    for pname in set(params) & set(params_mf):
        params[pname].value = params_mf[pname].value

    cps.set_values(params, defaults)

    return params


def create_profile_params(config, propagator):
    basis = config["basis"]
    model = config["model"]
    conditions = config["conditions"]
    spin_system = config["spin_system"]
    observed_state = config["experiment"]["observed_state"]

    # Create settings for the parameters of the models
    settings_k = cpk.make_settings[model.name](conditions, spin_system)
    settings_l, settings_l_mf = cpl.make_settings(basis, model, conditions)

    if model.model_free:
        settings_l = settings_l_mf
        fitted = config["fit"]["model_free"]
    else:
        fitted = config["fit"]["rates"]

    settings = settings_k | settings_l
    _set_to_fit(settings, model, observed_state, fitted)
    settings_min, settings_max = _get_settings(settings, propagator)
    pnames = cpn.get_pnames(settings_min, conditions, spin_system)
    fnames = {name: pname.full for name, pname in pnames.items()}

    # Create standard parameters from settings
    params = _settings_to_params(settings_max, conditions, spin_system)

    # Create the model free parameters
    params_mf = _settings_to_params(settings_l_mf, conditions, spin_system)

    return fnames, params, params_mf


def _settings_to_params(settings, conditions, spin_system):
    pnames = cpn.get_pnames(settings, conditions, spin_system)
    fnames = {name: pname.full for name, pname in pnames.items()}
    conditions_ = dc.asdict(conditions)
    parameter_list = [
        lf.Parameter(
            name=fnames[name],
            value=setting.get("value"),
            min=setting.get("min"),
            max=setting.get("max"),
            vary=setting.get("vary"),
            expr=setting.get("expr", "").format_map(conditions_ | fnames),
            user_data={
                "pname": pnames[name],
                "print_expr": setting.get("expr", "").format_map(conditions_ | pnames),
            },
        )
        for name, setting in settings.items()
    ]
    params = lf.Parameters(usersyms=cnr.rate_functions)
    params.add_many(*parameter_list)
    return params


def _get_settings(settings_full, propagator):
    settings_min = {k: v for k, v in settings_full.items() if k in propagator.snames}
    # Find the additional parameters that are used in the profile param constraints
    remaining_names = set(settings_full) - set(settings_min)
    settings_expr = {
        name: settings_full[name]
        for setting in settings_min.values()
        for name in re.findall(r"\{(.+?)\}", setting.get("expr", ""))
        if name in remaining_names
    }
    settings_max = settings_min | settings_expr
    return settings_min, settings_max


def _set_to_fit(settings, model, observed_state, fitted):
    for sname, state in it.product(fitted, model.states):
        sname_ = sname.format(states=state, observed_state=observed_state)
        if sname_ in settings:
            settings[sname_]["vary"] = True
            settings[sname_]["expr"] = ""

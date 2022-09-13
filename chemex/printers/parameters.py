from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from re import compile
from typing import DefaultDict

from chemex.containers.experiments import Experiments
from chemex.parameters import database
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting

Parameters = dict[ParamName, ParamSetting]

RE_GROUPNAME = compile(r"^[A-Za-z0-9_-]+$")


@dataclass
class GlobalLocalParameters:
    global_: Parameters
    local: Parameters

    def __bool__(self) -> bool:
        return bool(self.global_) or bool(self.local)


@dataclass
class ClassifiedParameters:
    fitted: GlobalLocalParameters
    fixed: GlobalLocalParameters
    constrained: GlobalLocalParameters


def _format_fitted(param: ParamSetting) -> str:
    error = f"±{param.stderr:.5e}" if param.stderr else "(error not calculated)"
    return f"{param.value: .5e} # {error}"


def _format_constrained(param: ParamSetting) -> str:
    if param.value is None:
        return ""

    error = f"±{param.stderr:.5e} " if param.stderr else ""
    constraint = param.expr
    parameters = database.get_parameters(param.dependencies)
    for param_id, parameter in parameters.items():
        constraint = constraint.replace(param_id, str(parameter.param_name))
    return f"{param.value: .5e} # {error}({constraint})"


def _format_fixed(param: ParamSetting) -> str:
    return f"{param.value: .5e} # (fixed)"


_format_param = {
    "fitted": _format_fitted,
    "constrained": _format_constrained,
    "fixed": _format_fixed,
}


def _params_to_strings(
    parameters: GlobalLocalParameters, status: str
) -> dict[str, dict[str, str]]:

    result: DefaultDict[str, dict[str, str]] = defaultdict(dict)

    for pname, param in parameters.global_.items():
        result["GLOBAL"][pname.section_res] = _format_param[status](param)

    for pname, param in parameters.local.items():
        result[pname.section][str(pname.spin_system)] = _format_param[status](param)

    return result


def _quote(text: str) -> str:
    text = text.strip(" ,")
    return text if RE_GROUPNAME.match(text) else f'"{text}"'


def _format_strings(par_strings: dict[str, dict[str, str]]) -> str:
    result = []
    for section, key_values in par_strings.items():
        result.append(f"[{_quote(section)}]")
        width = len(max(key_values, key=len))
        result.extend(
            f"{_quote(key):<{width}} = {value}" for key, value in key_values.items()
        )
        result.append("")
    return "\n".join(result)


def write_file(parameters: GlobalLocalParameters, status: str, path: Path) -> None:

    if not parameters:
        return

    par_strings = _params_to_strings(parameters, status)
    formatted_strings = _format_strings(par_strings)
    filename = path / f"{status}.toml"
    filename.write_text(formatted_strings)


def classify_global(parameters: Parameters) -> GlobalLocalParameters:
    local = {}
    global_ = {}

    for pname, parameter in parameters.items():
        if parameter.param_name.spin_system:
            local[pname] = parameter
        else:
            global_[pname] = parameter

    return GlobalLocalParameters(global_, local)


def classify_parameters(experiments: Experiments) -> ClassifiedParameters:

    param_ids = experiments.param_ids
    parameters = {
        param.param_name: param for param in database.get_parameters(param_ids).values()
    }

    constrained = {pname: param for pname, param in parameters.items() if param.expr}
    fitted = {
        pname: param
        for pname, param in parameters.items()
        if param.vary and pname not in constrained
    }

    fixed_ids = set(parameters) - set(fitted) - set(constrained)
    fixed = {
        pname: parameter
        for pname, parameter in parameters.items()
        if pname in fixed_ids
    }

    return ClassifiedParameters(
        classify_global(fitted), classify_global(fixed), classify_global(constrained)
    )


def write_parameters(experiments: Experiments, path: Path):
    """Write the model parameter values and their uncertainties to a file"""
    path_par = path / "Parameters"
    path_par.mkdir(parents=True, exist_ok=True)
    classified_parameters = classify_parameters(experiments)

    write_file(classified_parameters.fitted, "fitted", path_par)
    write_file(classified_parameters.fixed, "fixed", path_par)
    write_file(classified_parameters.constrained, "constrained", path_par)

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from chemex.optimize.deterministic_uncertainty import (
    DerivationDisposition,
    DeterministicUncertainty,
)
from chemex.optimize.uncertainty import UncertaintyUnavailableKind
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ParameterRole,
    SealedParameterModel,
)
from chemex.parameters.sealed import parameter_name_from_definition

RE_GROUPNAME = re.compile(r"^[A-Za-z0-9_-]+$")


_UNAVAILABLE_REASON_TEXT = {
    UncertaintyUnavailableKind.RANK_DEFICIENT: "rank deficient",
    UncertaintyUnavailableKind.INSUFFICIENT_INFORMATION: "insufficient information",
    UncertaintyUnavailableKind.BOUNDARY_LIMITED: (
        "boundary limited (active relaxation-PSD boundary)"
    ),
    UncertaintyUnavailableKind.NORMALIZATION_INVALID: "normalization invalid",
    UncertaintyUnavailableKind.JACOBIAN_UNAVAILABLE: "Jacobian unavailable",
    UncertaintyUnavailableKind.COVARIANCE_NUMERICAL_FAILURE: (
        "covariance numerical failure"
    ),
    UncertaintyUnavailableKind.UNSUPPORTED_CONSTRAINED_DERIVATIVE: (
        "unsupported constrained derivative"
    ),
    UncertaintyUnavailableKind.DERIVATION_STOPPED: ("derivation interrupted/cancelled"),
    UncertaintyUnavailableKind.CONSTRAINED_PROPAGATION_UNAVAILABLE: (
        "constrained propagation unavailable"
    ),
}

BOUNDARY_WARNING_TEXT = "boundary may make uncertainty asymmetric"
GRID_UNCERTAINTY_WITHHELD_TEXT = (
    "withheld: GRID coordinates were selected on a discrete profiled surface, "
    "not by a full continuous Direct TRF"
)

if set(_UNAVAILABLE_REASON_TEXT) != set(UncertaintyUnavailableKind):
    raise RuntimeError("Uncertainty unavailability vocabulary is not exhaustive")


def uncertainty_unavailable_reason(kind: UncertaintyUnavailableKind) -> str:
    """Return the one human-readable phrase for a typed claim classification."""
    return _UNAVAILABLE_REASON_TEXT[kind]


def uncertainty_progress_status(uncertainty: DeterministicUncertainty) -> str:
    """Render the established terminal progress text from scientific facts."""
    if uncertainty.boundary_warning:
        return "covariance available with boundary warnings"
    if uncertainty.root_covariance_available:
        return "covariance available"
    if uncertainty.block_covariance_available:
        return "covariance partially available"
    controlled_ids = uncertainty.accepted_anchor.controlled_ids
    reason = next(
        (
            uncertainty_unavailable_reason(item.unavailable_kind)
            for param_id in controlled_ids
            if (item := uncertainty.parameter(param_id)) is not None
            and item.unavailable_kind is not None
        ),
        "insufficient information",
    )
    return f"uncertainty unavailable: {reason}"


@dataclass(frozen=True, slots=True)
class ReportParameter:
    """One disposable presentation row joined from sealed metadata and values."""

    param_id: str
    param_name: ParamName
    value: float


type Parameters = dict[ParamName, ReportParameter]


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


def _parameter_uncertainty_rendering(
    param_id: str,
    uncertainty: DeterministicUncertainty | None,
) -> tuple[float | None, str | None, str | None]:
    """Translate one predetermined conclusion into presentation values."""
    conclusion = None if uncertainty is None else uncertainty.parameter(param_id)
    standard_error = None if conclusion is None else conclusion.standard_error
    warning = (
        BOUNDARY_WARNING_TEXT
        if conclusion is not None and conclusion.boundary_warning
        else None
    )
    reason = (
        GRID_UNCERTAINTY_WITHHELD_TEXT
        if uncertainty is not None
        and uncertainty.disposition is DerivationDisposition.WITHHELD
        and conclusion is not None
        else uncertainty_unavailable_reason(conclusion.unavailable_kind)
        if conclusion is not None and conclusion.unavailable_kind is not None
        else None
    )
    return standard_error, warning, reason


def _format_fitted(
    param: ReportParameter,
    uncertainty: DeterministicUncertainty | None,
) -> str:
    standard_error, warning, reason = _parameter_uncertainty_rendering(
        param.param_id,
        uncertainty,
    )
    error = (
        f"±{standard_error:.5e}" + (f" ({warning})" if warning is not None else "")
        if standard_error is not None
        else f"(error unavailable: {reason})"
        if reason is not None
        else "(error not calculated)"
    )
    return f"{param.value: .5e} # {error}"


def _format_constrained(
    param: ReportParameter,
    uncertainty: DeterministicUncertainty | None,
    constraint_expression: str,
) -> str:
    standard_error, warning, reason = _parameter_uncertainty_rendering(
        param.param_id,
        uncertainty,
    )
    error = (
        f"±{standard_error:.5e}" + (f" ({warning}); " if warning is not None else " ")
        if standard_error is not None
        else f"error unavailable: {reason}; "
        if reason is not None
        else ""
    )
    return f"{param.value: .5e} # {error}({constraint_expression})"


def _format_fixed(param: ReportParameter) -> str:
    return f"{param.value: .5e} # (fixed)"


def _params_to_strings(
    parameters: GlobalLocalParameters,
    status: str,
    uncertainty: DeterministicUncertainty | None,
    constraint_expressions: Mapping[str, str],
) -> dict[str, dict[str, str]]:
    result: defaultdict[str, dict[str, str]] = defaultdict(dict)

    for pname, param in parameters.global_.items():
        result["GLOBAL"][pname.section_res] = _format_parameter(
            param,
            status,
            uncertainty,
            constraint_expressions.get(param.param_id, ""),
        )

    for pname, param in parameters.local.items():
        result[pname.section][str(pname.spin_system)] = _format_parameter(
            param,
            status,
            uncertainty,
            constraint_expressions.get(param.param_id, ""),
        )

    return result


def _format_parameter(
    param: ReportParameter,
    status: str,
    uncertainty: DeterministicUncertainty | None,
    constraint_expression: str = "",
) -> str:
    if status == "fitted":
        return _format_fitted(param, uncertainty)
    if status == "fixed":
        return _format_fixed(param)
    if status == "constrained":
        return _format_constrained(
            param,
            uncertainty,
            constraint_expression,
        )
    msg = (
        f"Unknown parameter status: {status!r}. Expected 'fitted', 'fixed', "
        "or 'constrained'."
    )
    raise ValueError(msg)


def _quote(text: str) -> str:
    text = text.strip(" ,")
    return text if RE_GROUPNAME.match(text) else f'"{text}"'


def _format_strings(par_strings: dict[str, dict[str, str]]) -> str:
    result: list[str] = []
    for section, key_values in par_strings.items():
        result.append(f"[{_quote(section)}]")
        quoted_keys = [_quote(key) for key in key_values]
        width = len(max(quoted_keys, key=len))
        result.extend(
            f"{_quote(key):<{width}} = {value}" for key, value in key_values.items()
        )
        result.append("")
    return "\n".join(result)


def write_file(
    parameters: GlobalLocalParameters,
    status: str,
    path: Path,
    uncertainty: DeterministicUncertainty | None = None,
    constraint_expressions: Mapping[str, str] | None = None,
) -> None:
    if not parameters:
        return

    par_strings = _params_to_strings(
        parameters,
        status,
        uncertainty,
        {} if constraint_expressions is None else constraint_expressions,
    )
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


def classify_parameters(
    parameter_model: SealedParameterModel,
    parameter_values: Mapping[str, float],
    parameterization: ActiveParameterization,
    fitted_ids: tuple[str, ...],
) -> ClassifiedParameters:
    parameters = {
        parameter_name_from_definition(definition): ReportParameter(
            definition.param_id,
            parameter_name_from_definition(definition),
            parameter_values[definition.param_id],
        )
        for definition in parameter_model.definitions
        if definition.param_id in parameterization.scope_ids
    }
    fitted_id_set = set(fitted_ids)
    constrained = {
        pname: param
        for pname, param in parameters.items()
        if parameterization.role(param.param_id) is ParameterRole.DERIVED
    }
    fitted = {
        pname: param
        for pname, param in parameters.items()
        if param.param_id in fitted_id_set
    }

    fixed_ids = set(parameters) - set(fitted) - set(constrained)
    fixed = {
        pname: parameter
        for pname, parameter in parameters.items()
        if pname in fixed_ids
    }

    return ClassifiedParameters(
        classify_global(fitted),
        classify_global(fixed),
        classify_global(constrained),
    )


def write_parameters(
    path: Path,
    *,
    parameter_model: SealedParameterModel,
    parameter_values: Mapping[str, float],
    parameterization: ActiveParameterization,
    fitted_ids: tuple[str, ...] = (),
    deterministic_uncertainty: DeterministicUncertainty | None = None,
) -> None:
    """Write the model parameter values and their uncertainties to a file."""
    path_par = path / "Parameters"
    path_par.mkdir(parents=True, exist_ok=True)
    classified_parameters = classify_parameters(
        parameter_model,
        parameter_values,
        parameterization,
        fitted_ids,
    )
    names = {
        definition.param_id: parameter_name_from_definition(definition)
        for definition in parameter_model.definitions
    }
    constraint_expressions = {
        item.target_id: _replace_parameter_ids(item.expression_text, names)
        for item in parameterization.program.constraints
    }

    write_file(
        classified_parameters.fitted,
        "fitted",
        path_par,
        deterministic_uncertainty,
    )
    write_file(
        classified_parameters.fixed,
        "fixed",
        path_par,
    )
    write_file(
        classified_parameters.constrained,
        "constrained",
        path_par,
        deterministic_uncertainty,
        constraint_expressions,
    )


def _replace_parameter_ids(
    expression: str,
    parameter_names: Mapping[str, ParamName],
) -> str:
    """Render stable IDs in one constraint using sealed presentation names."""
    if not parameter_names:
        return expression
    pattern = re.compile(
        "|".join(
            re.escape(param_id)
            for param_id in sorted(parameter_names, key=len, reverse=True)
        )
    )
    return pattern.sub(
        lambda match: str(parameter_names[match.group(0)]),
        expression,
    )

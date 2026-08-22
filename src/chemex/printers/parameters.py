from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from chemex.containers.experiments import Experiments
from chemex.optimize.uncertainty import (
    UncertaintyEvidence,
    UncertaintyUnavailableKind,
)
from chemex.parameters.database import ParameterStore
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting

Parameters = dict[ParamName, ParamSetting]

RE_GROUPNAME = re.compile(r"^[A-Za-z0-9_-]+$")


@dataclass(frozen=True, slots=True)
class ParameterUncertaintyView:
    """Immutable report-only covariance-derived uncertainty selection."""

    standard_errors: tuple[tuple[str, float], ...] = ()
    warnings: tuple[tuple[str, str], ...] = ()
    unavailable_reasons: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        error_ids = {param_id for param_id, _value in self.standard_errors}
        warning_ids = {param_id for param_id, _warning in self.warnings}
        unavailable_ids = {param_id for param_id, _reason in self.unavailable_reasons}
        if warning_ids - error_ids:
            raise ValueError("Uncertainty warnings require an available standard error")
        if error_ids & unavailable_ids:
            raise ValueError(
                "A parameter uncertainty cannot be available and unavailable"
            )

    def standard_error(self, param_id: str) -> float | None:
        return dict(self.standard_errors).get(param_id)

    def warning(self, param_id: str) -> str | None:
        return dict(self.warnings).get(param_id)

    def unavailable_reason(self, param_id: str) -> str | None:
        return dict(self.unavailable_reasons).get(param_id)


_UNAVAILABLE_REASON_TEXT = {
    UncertaintyUnavailableKind.RANK_DEFICIENT: "rank deficient",
    UncertaintyUnavailableKind.INSUFFICIENT_INFORMATION: "insufficient information",
    UncertaintyUnavailableKind.BOUNDARY_LIMITED: "boundary limited",
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

BOUNDARY_WARNING_TEXT = "boundary may make local error asymmetric"

if set(_UNAVAILABLE_REASON_TEXT) != set(UncertaintyUnavailableKind):
    raise RuntimeError("Uncertainty unavailability vocabulary is not exhaustive")


def uncertainty_unavailable_reason(kind: UncertaintyUnavailableKind) -> str:
    """Return the one human-readable phrase for a typed claim classification."""
    return _UNAVAILABLE_REASON_TEXT[kind]


def _controlled_unavailability_kind(
    evidence: UncertaintyEvidence,
) -> UncertaintyUnavailableKind:
    categories = {failure.category for failure in evidence.failures}
    if categories & {"cancelled", "interrupted"}:
        return UncertaintyUnavailableKind.DERIVATION_STOPPED
    if "rank_deficient" in categories:
        return UncertaintyUnavailableKind.RANK_DEFICIENT
    if categories & {
        "insufficient_effective_observations",
        "non_positive_nominal_residual_degrees_of_freedom",
    }:
        return UncertaintyUnavailableKind.INSUFFICIENT_INFORMATION
    covariance = evidence.covariance
    if covariance is not None:
        claims = {item.name: item.state.value for item in covariance.claims}
        if claims.get("PROFILED_NORMALIZATION_REGULARITY") == "violated":
            return UncertaintyUnavailableKind.NORMALIZATION_INVALID
    if any(failure.stage == "residual_linearization" for failure in evidence.failures):
        return UncertaintyUnavailableKind.JACOBIAN_UNAVAILABLE
    return UncertaintyUnavailableKind.COVARIANCE_NUMERICAL_FAILURE


def parameter_uncertainty_view(
    evidence: UncertaintyEvidence,
    unsupported_constrained_ids: tuple[str, ...] = (),
) -> ParameterUncertaintyView:
    """Select only scientifically reportable covariance-derived errors."""
    standard_errors: list[tuple[str, float]] = []
    warnings: list[tuple[str, str]] = []
    unavailable_reasons: list[tuple[str, str]] = []
    marginal = evidence.marginal_errors
    if marginal is not None and marginal.scope_reportable:
        standard_errors.extend(
            (entry.param_id, entry.value)
            for entry in marginal.entries
            if entry.value is not None
        )
        if evidence.covariance is not None and evidence.covariance.boundary_warning:
            warnings.extend(
                (entry.param_id, BOUNDARY_WARNING_TEXT)
                for entry in marginal.entries
                if entry.value is not None
            )
    else:
        reason = uncertainty_unavailable_reason(
            _controlled_unavailability_kind(evidence)
        )
        unavailable_reasons.extend(
            (param_id, reason) for param_id in evidence.accepted_anchor.controlled_ids
        )
    constrained = evidence.constrained_marginal_errors
    if constrained is not None and constrained.scope_reportable:
        standard_errors.extend(
            (entry.param_id, entry.value)
            for entry in constrained.entries
            if entry.value is not None
        )
        if evidence.covariance is not None and evidence.covariance.boundary_warning:
            warnings.extend(
                (entry.param_id, BOUNDARY_WARNING_TEXT)
                for entry in constrained.entries
                if entry.value is not None
            )
    elif evidence.requested_output_scope:
        unavailable_reasons.extend(
            (
                param_id,
                uncertainty_unavailable_reason(
                    UncertaintyUnavailableKind.CONSTRAINED_PROPAGATION_UNAVAILABLE
                ),
            )
            for param_id in evidence.requested_output_scope
        )
    unavailable_reasons.extend(
        (
            param_id,
            uncertainty_unavailable_reason(
                UncertaintyUnavailableKind.UNSUPPORTED_CONSTRAINED_DERIVATIVE
            ),
        )
        for param_id in unsupported_constrained_ids
    )
    return ParameterUncertaintyView(
        standard_errors=tuple(standard_errors),
        warnings=tuple(warnings),
        unavailable_reasons=tuple(unavailable_reasons),
    )


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


def _format_fitted(
    param_id: str,
    param: ParamSetting,
    uncertainty: ParameterUncertaintyView | None,
) -> str:
    standard_error = (
        None if uncertainty is None else uncertainty.standard_error(param_id)
    )
    warning = None if uncertainty is None else uncertainty.warning(param_id)
    reason = None if uncertainty is None else uncertainty.unavailable_reason(param_id)
    error = (
        f"±{standard_error:.5e}" + (f" ({warning})" if warning is not None else "")
        if standard_error is not None
        else f"(error unavailable: {reason})"
        if reason is not None
        else "(error not calculated)"
    )
    return f"{param.value: .5e} # {error}"


def _format_constrained(
    param_id: str,
    param: ParamSetting,
    parameter_store: ParameterStore,
    uncertainty: ParameterUncertaintyView | None,
) -> str:
    if param.value is None:
        return ""

    standard_error = (
        None if uncertainty is None else uncertainty.standard_error(param_id)
    )
    warning = None if uncertainty is None else uncertainty.warning(param_id)
    reason = None if uncertainty is None else uncertainty.unavailable_reason(param_id)
    error = (
        f"±{standard_error:.5e}" + (f" ({warning}); " if warning is not None else " ")
        if standard_error is not None
        else f"error unavailable: {reason}; "
        if reason is not None
        else ""
    )
    constraint = param.expr
    parameters = parameter_store.get_parameters(param.dependencies)
    for param_id, parameter in parameters.items():
        constraint = constraint.replace(param_id, str(parameter.param_name))
    return f"{param.value: .5e} # {error}({constraint})"


def _format_fixed(param: ParamSetting) -> str:
    return f"{param.value: .5e} # (fixed)"


def _params_to_strings(
    parameters: GlobalLocalParameters,
    status: str,
    parameter_store: ParameterStore,
    parameter_ids: dict[ParamName, str],
    uncertainty: ParameterUncertaintyView | None,
) -> dict[str, dict[str, str]]:
    result: defaultdict[str, dict[str, str]] = defaultdict(dict)

    for pname, param in parameters.global_.items():
        result["GLOBAL"][pname.section_res] = _format_parameter(
            parameter_ids[pname],
            param,
            status,
            parameter_store,
            uncertainty,
        )

    for pname, param in parameters.local.items():
        result[pname.section][str(pname.spin_system)] = _format_parameter(
            parameter_ids[pname],
            param,
            status,
            parameter_store,
            uncertainty,
        )

    return result


def _format_parameter(
    param_id: str,
    param: ParamSetting,
    status: str,
    parameter_store: ParameterStore,
    uncertainty: ParameterUncertaintyView | None,
) -> str:
    if status == "fitted":
        return _format_fitted(param_id, param, uncertainty)
    if status == "fixed":
        return _format_fixed(param)
    if status == "constrained":
        return _format_constrained(param_id, param, parameter_store, uncertainty)
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
    parameter_store: ParameterStore,
    parameter_ids: dict[ParamName, str] | None = None,
    uncertainty: ParameterUncertaintyView | None = None,
) -> None:
    if not parameters:
        return

    ids = (
        {pname: parameter.id_ for pname, parameter in parameters.global_.items()}
        | {pname: parameter.id_ for pname, parameter in parameters.local.items()}
        if parameter_ids is None
        else parameter_ids
    )
    par_strings = _params_to_strings(
        parameters,
        status,
        parameter_store,
        ids,
        uncertainty,
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
    experiments: Experiments,
    parameter_store: ParameterStore,
) -> ClassifiedParameters:
    param_ids = experiments.param_ids
    parameters = {
        param.param_name: param
        for param in parameter_store.get_parameters(param_ids).values()
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
        classify_global(fitted),
        classify_global(fixed),
        classify_global(constrained),
    )


def write_parameters(
    experiments: Experiments,
    path: Path,
    parameter_store: ParameterStore,
    uncertainty: ParameterUncertaintyView | None = None,
) -> None:
    """Write the model parameter values and their uncertainties to a file."""
    path_par = path / "Parameters"
    path_par.mkdir(parents=True, exist_ok=True)
    classified_parameters = classify_parameters(experiments, parameter_store)
    parameter_ids = {
        parameter.param_name: param_id
        for param_id, parameter in parameter_store.get_parameters(
            experiments.param_ids
        ).items()
    }

    write_file(
        classified_parameters.fitted,
        "fitted",
        path_par,
        parameter_store,
        parameter_ids,
        uncertainty,
    )
    write_file(
        classified_parameters.fixed,
        "fixed",
        path_par,
        parameter_store,
        parameter_ids,
    )
    write_file(
        classified_parameters.constrained,
        "constrained",
        path_par,
        parameter_store,
        parameter_ids,
        uncertainty,
    )

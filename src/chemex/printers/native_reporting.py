"""Typed artifact rendering for native step-root publications."""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Protocol, cast

from scipy import stats

from chemex.evaluation.native import EvaluationPlan, EvaluationResult
from chemex.native_provenance import serialize_independent_parameters
from chemex.optimize.direct_trf import canonical_chi_square
from chemex.optimize.native_mcmc import (
    McmcEvidence,
    PosteriorSampleEvidence,
    PosteriorScalarEstimate,
    PosteriorSummary,
)
from chemex.optimize.native_resampling import (
    ResamplingEvidence,
    ResamplingSummaryOutcome,
    SummaryFailure,
)
from chemex.optimize.uncertainty import (
    RootAnchoredBlockCovarianceEvidence,
    UncertaintyEvidence,
)
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ParameterRole,
    SealedParameterModel,
)
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.printers.parameters import (
    parameter_uncertainty_view,
    uncertainty_unavailable_reason,
)


class ComponentDiagnosticRecord(Protocol):
    """Fields needed to render a non-authoritative component diagnostic."""

    @property
    def identity(self) -> str: ...

    @property
    def disposition(self) -> str: ...

    @property
    def controlled_ids(self) -> tuple[str, ...]: ...

    @property
    def local_chi_square(self) -> float | None: ...


def _toml_float(value: float) -> str:
    scalar = float(value)
    if not math.isfinite(scalar):
        raise ValueError("Published statistics must be finite")
    return repr(0.0 if scalar == 0.0 else scalar)


def _toml_key(value: str) -> str:
    return json.dumps(value, ensure_ascii=True)


def write_json(path: Path, record: object) -> None:
    """Write one deterministic, strict JSON artifact."""
    path.write_text(
        json.dumps(
            record,
            allow_nan=False,
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def write_parameter_reports(
    path: Path,
    parameterization: ActiveParameterization,
    result: EvaluationResult,
    uncertainty: UncertaintyEvidence | None = None,
) -> None:
    """Separate central parameter values by their authoritative role."""
    path.mkdir()
    constraints = {
        item.target_id: item.expression_text
        for item in parameterization.program.constraints
    }
    uncertainty_view = (
        None if uncertainty is None else parameter_uncertainty_view(uncertainty)
    )
    by_role: dict[ParameterRole, list[tuple[str, float]]] = {
        role: [] for role in ParameterRole
    }
    for param_id in parameterization.scope_ids:
        by_role[parameterization.role(param_id)].append(
            (param_id, result.resolved_values[param_id])
        )
    filenames = {
        ParameterRole.FIT: "fitted.toml",
        ParameterRole.FIX: "fixed.toml",
        ParameterRole.DERIVED: "constrained.toml",
    }
    for role, items in by_role.items():
        if not items:
            continue
        lines = ["[parameters]"]
        for param_id, value in items:
            comment = ""
            if role is ParameterRole.FIX:
                comment = " # fixed"
            elif role is ParameterRole.DERIVED:
                expression = constraints.get(param_id, "model-derived")
                comment = f" # constrained: {expression}"
            if role in {ParameterRole.FIT, ParameterRole.DERIVED} and (
                uncertainty_view is not None
            ):
                standard_error = uncertainty_view.standard_error(param_id)
                reason = uncertainty_view.unavailable_reason(param_id)
                if standard_error is not None:
                    prefix = f" # ±{_toml_float(standard_error)}"
                    comment = (
                        prefix
                        if role is ParameterRole.FIT
                        else f"{prefix};{comment.removeprefix(' #')}"
                    )
                elif reason is not None:
                    prefix = f" # error unavailable: {reason}"
                    comment = (
                        prefix
                        if role is ParameterRole.FIT
                        else f"{prefix};{comment.removeprefix(' #')}"
                    )
            lines.append(f"{_toml_key(param_id)} = {_toml_float(value)}{comment}")
        (path / filenames[role]).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_evaluated_parameters(
    path: Path,
    parameterization: ActiveParameterization,
    result: EvaluationResult,
) -> None:
    """Write evaluate-only values without claiming fitted or restart authority."""
    path.mkdir()
    lines = ["[parameters]"]
    for param_id in parameterization.scope_ids:
        role = parameterization.role(param_id).value
        value = result.resolved_values[param_id]
        lines.append(f"{_toml_key(param_id)} = {_toml_float(value)} # {role}")
    (path / "evaluated.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_restart_parameters(
    path: Path,
    parameter_model: SealedParameterModel,
    parameterization: ActiveParameterization,
    snapshot: AnalysisValuesSnapshot,
) -> None:
    """Write complete committed independent state in parameter-input form."""
    path.write_text(
        serialize_independent_parameters(
            parameter_model,
            parameterization.independent_ids,
            snapshot,
            state_kind="committed",
        ),
        encoding="utf-8",
    )


def write_data(
    path: Path,
    plan: EvaluationPlan,
    result: EvaluationResult,
) -> None:
    """Write aggregate calculated data and retained residual membership."""
    path.mkdir()
    for ordinal, (descriptor, profile) in enumerate(
        zip(plan.profiles, result.profiles, strict=True)
    ):
        start = descriptor.observation_offset
        retained = set(profile.retained_observation_indices)
        residual_by_index = dict(
            zip(
                profile.retained_observation_indices,
                result.residuals[
                    profile.residual_offset : profile.residual_offset
                    + profile.residual_count
                ],
                strict=True,
            )
        )
        lines = [
            (
                "index\texperimental\tuncertainty\tretained\t"
                "calculated_unscaled\tcalculated\tresidual\tnormalization"
            )
        ]
        for index in range(descriptor.observation_count):
            residual = (
                _toml_float(float(residual_by_index[index]))
                if index in retained
                else ""
            )
            lines.append(
                "\t".join(
                    (
                        str(index),
                        _toml_float(descriptor.experimental_values[index]),
                        _toml_float(descriptor.uncertainties[index]),
                        "true" if index in retained else "false",
                        _toml_float(float(result.unscaled_calculations[start + index])),
                        _toml_float(
                            float(result.normalized_calculations[start + index])
                        ),
                        residual,
                        _toml_float(profile.normalization_factor),
                    )
                )
            )
        (path / f"profile_{ordinal:04d}.tsv").write_text(
            "\n".join(lines) + "\n", encoding="utf-8"
        )


def write_statistics(
    path: Path,
    plan: EvaluationPlan,
    result: EvaluationResult,
    controlled_count: int,
) -> None:
    """Write aggregate statistics from retained residuals and approved dof."""
    residual_count = len(result.residuals)
    normalization_count = sum(profile.is_scaled for profile in plan.profiles)
    parameter_count = controlled_count + normalization_count
    degrees_of_freedom = residual_count - parameter_count
    chi_square = canonical_chi_square(result.residuals)
    reduced_chi_square: float | str = (
        chi_square / degrees_of_freedom if degrees_of_freedom > 0 else "unavailable"
    )
    chi_square_p_value: float | str = (
        float(stats.chi2.sf(chi_square, degrees_of_freedom))
        if degrees_of_freedom > 0
        else "unavailable"
    )
    ks_p_value = float(stats.kstest(result.residuals, "norm").pvalue)
    aic = chi_square + 2.0 * parameter_count
    bic = chi_square + math.log(residual_count) * parameter_count

    def value(item: int | float | str) -> str:
        if isinstance(item, str):
            return _toml_key(item)
        if isinstance(item, int):
            return str(item)
        return _toml_float(item)

    items: tuple[tuple[str, int | float | str], ...] = (
        ("number of data points", residual_count),
        ("number of variables", parameter_count),
        ("retained residual count", residual_count),
        ("controlled parameter count", controlled_count),
        ("profiled normalization count", normalization_count),
        ("nominal residual degrees of freedom", degrees_of_freedom),
        ("chi-square", chi_square),
        ("reduced-chi-square", reduced_chi_square),
        ("chi-squared test", chi_square_p_value),
        ("Kolmogorov-Smirnov test", ks_p_value),
        ("Akaike Information Criterion (AIC)", aic),
        ("Bayesian Information Criterion (BIC)", bic),
    )
    lines = [f"{_toml_key(key)} = {value(item)}" for key, item in items]
    (path / "statistics.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def component_diagnostic_records(
    components: Sequence[ComponentDiagnosticRecord],
) -> list[dict[str, object]]:
    """Return the shared non-authoritative component record representation."""
    return [
        {
            "ordinal": ordinal,
            "identity": item.identity,
            "disposition": item.disposition,
            "controlled_ids": list(item.controlled_ids),
            "local_chi_square": item.local_chi_square,
        }
        for ordinal, item in enumerate(components)
    ]


def write_components(
    path: Path,
    components: Sequence[ComponentDiagnosticRecord],
) -> None:
    """Write optional component execution diagnostics without authority copies."""
    if not components:
        return
    path.mkdir()
    write_json(
        path / "index.json",
        {
            "schema_version": 1,
            "authority": "diagnostic_only",
            "components": component_diagnostic_records(components),
        },
    )


def write_uncertainty(path: Path, uncertainty: UncertaintyEvidence | None) -> None:
    """Write covariance and constrained evidence in distinct typed families."""
    if uncertainty is None:
        return
    record = uncertainty.to_record()
    payload = record.get("payload")
    if not isinstance(payload, dict):
        raise TypeError("Uncertainty evidence has no typed payload")
    common = {
        "schema_version": record["schema_version"],
        "accepted_result_identity": uncertainty.accepted_result_identity,
        "accepted_occurrence_identity": uncertainty.accepted_occurrence_identity,
        "request_identity": uncertainty.request_identity,
        "policy_identity": uncertainty.policy_identity,
        "bundle_identity": uncertainty.identity,
    }
    covariance_path = path / "Covariance"
    covariance_path.mkdir()
    write_json(
        covariance_path / "evidence.json",
        common
        | {
            "artifact_type": "native_covariance_evidence",
            "residual_jacobian": payload.get("residual_jacobian"),
            "rank_diagnostic": payload.get("rank_diagnostic"),
            "covariance": payload.get("covariance"),
            "marginal_standard_errors": payload.get("marginal_errors"),
            "correlations": payload.get("correlations"),
            "operations": payload.get("operations"),
            "failures": payload.get("failures"),
        },
    )
    if uncertainty.requested_output_scope:
        constrained_path = path / "Constrained"
        constrained_path.mkdir()
        write_json(
            constrained_path / "evidence.json",
            common
            | {
                "artifact_type": "native_constrained_evidence",
                "requested_output_scope": list(uncertainty.requested_output_scope),
                "constraint_jacobian": payload.get("constraint_jacobian"),
                "constrained_propagation": payload.get("constrained_propagation"),
                "marginal_standard_errors": payload.get("constrained_marginal_errors"),
                "correlations": payload.get("constrained_correlations"),
                "operations": payload.get("operations"),
                "failures": payload.get("failures"),
            },
        )


def write_block_uncertainty(
    path: Path,
    evidence: RootAnchoredBlockCovarianceEvidence | None,
) -> None:
    """Write root-authoritative failure-isolated covariance blocks."""
    if evidence is None:
        return
    covariance_path = path / "Covariance"
    covariance_path.mkdir(exist_ok=True)
    record = evidence.to_record()
    blocks = record["blocks"]
    if not isinstance(blocks, list):
        raise TypeError("Block covariance record has an invalid payload")
    for block, block_record in zip(evidence.blocks, blocks, strict=True):
        if not isinstance(block_record, dict):
            raise TypeError("Block covariance record has an invalid block")
        typed_block_record = cast("dict[str, object]", block_record)
        typed_block_record["unavailable_reason"] = (
            ""
            if block.unavailable_kind is None
            else uncertainty_unavailable_reason(block.unavailable_kind)
        )
    write_json(covariance_path / "blocks.json", record)


def resampling_evidence_record(evidence: ResamplingEvidence) -> dict[str, object]:
    """Serialize one primary resampling execution artifact."""
    return {
        "artifact_type": "native_resampling_evidence",
        "schema_version": 1,
        "identity": evidence.identity,
        "plan_identity": evidence.plan.identity,
        "accepted_result_identity": evidence.plan.accepted_result_identity,
        "accepted_occurrence_identity": evidence.plan.accepted_occurrence_identity,
        "scheme": evidence.plan.scheme.value,
        "lifecycle": evidence.lifecycle.value,
        "intended_count": evidence.plan.replicate_count,
        "completed_count": evidence.completed_count,
        "successful_count": evidence.successful_count,
        "failed_count": evidence.failed_count,
        "coverage_satisfied": evidence.coverage_satisfied,
        "claims": [
            {
                "name": claim.name,
                "state": claim.state.value,
                "policy_identity": claim.policy_identity,
                "details": list(claim.details),
            }
            for claim in evidence.claims
        ],
        "outcomes": [
            {
                "ordinal": outcome.ordinal,
                "seed": outcome.seed,
                "stage": outcome.stage.value,
                "disposition": outcome.disposition.value,
                "identity": outcome.identity,
                "draw_identity": outcome.draw_identity,
                "success": (
                    None if outcome.success is None else outcome.success.to_record()
                ),
                "failure": (
                    None
                    if outcome.failure is None
                    else {
                        "category": outcome.failure.category,
                        "message": outcome.failure.message,
                        "identity": outcome.failure.identity,
                    }
                ),
            }
            for outcome in evidence.outcomes
        ],
    }


def write_resampling(
    path: Path,
    publications: Sequence[tuple[ResamplingEvidence, ResamplingSummaryOutcome]],
) -> None:
    """Write resampling evidence and summaries under scheme-specific families."""
    if not publications:
        return
    root = path / "Resampling"
    root.mkdir()
    for evidence, outcome in publications:
        family = root / evidence.plan.scheme.value.upper()
        family.mkdir()
        write_json(family / "evidence.json", resampling_evidence_record(evidence))
        if outcome.summary is not None:
            record = outcome.summary.to_record()
            write_json(
                family / "summary.json",
                {
                    "artifact_type": "native_resampling_summary",
                    "schema_version": record["schema_version"],
                    "scheme": evidence.plan.scheme.value,
                    **record,
                },
            )
        elif outcome.failure is not None:
            write_json(
                family / "summary.json",
                {
                    "artifact_type": "native_resampling_summary_unavailable",
                    "schema_version": 1,
                    "scheme": evidence.plan.scheme.value,
                    "source_evidence_identity": outcome.failure.source_evidence_identity,
                    "category": outcome.failure.category,
                    "message": outcome.failure.message,
                    "identity": outcome.failure.identity,
                },
            )


def _posterior_estimate_record(
    estimate: PosteriorScalarEstimate,
) -> dict[str, str | None]:
    return {
        "status": estimate.status.value,
        "reason": None if estimate.reason is None else estimate.reason.value,
        "value": None if estimate.value is None else float(estimate.value).hex(),
    }


def _posterior_sample_record(samples: PosteriorSampleEvidence) -> dict[str, object]:
    return {
        "artifact_type": "native_posterior_sample_evidence",
        "schema_version": 1,
        "identity": samples.identity,
        "selection_identity": samples.selection_identity,
        "output_scope": list(samples.output_scope),
        "successful_labels": [list(label) for label in samples.successful_labels],
        "outcomes": [
            {
                "state_ordinal": item.state_ordinal,
                "walker_ordinal": item.walker_ordinal,
                "disposition": item.disposition.value,
                "independent_items": [
                    [param_id, float(value).hex()]
                    for param_id, value in item.independent_items
                ],
                "resolved_items": (
                    None
                    if item.resolved_items is None
                    else [
                        [param_id, float(value).hex()]
                        for param_id, value in item.resolved_items
                    ]
                ),
                "failure": (
                    None
                    if item.failure is None
                    else {
                        "category": item.failure.category,
                        "message": item.failure.message,
                    }
                ),
                "identity": item.identity,
            }
            for item in samples.outcomes
        ],
    }


def _posterior_summary_record(summary: PosteriorSummary) -> dict[str, object]:
    return {
        "artifact_type": "native_posterior_summary",
        "schema_version": 1,
        "identity": summary.identity,
        "source_identity": summary.source_identity,
        "policy_identity": summary.policy_identity,
        "included_labels": [list(label) for label in summary.included_labels],
        "excluded_labels": [list(label) for label in summary.excluded_labels],
        "parameter_summaries": [
            {
                "parameter_id": item.parameter_id,
                "unit": item.unit.value,
                "quantiles": [
                    [float(probability).hex(), float(value).hex()]
                    for probability, value in item.quantiles
                ],
                "credible_interval_name": item.credible_interval_name,
                "credible_interval": [
                    float(value).hex() for value in item.credible_interval
                ],
                "posterior_standard_deviation": float(
                    item.posterior_standard_deviation
                ).hex(),
                "autocorrelation_time": _posterior_estimate_record(
                    item.autocorrelation_time
                ),
                "effective_sample_size": _posterior_estimate_record(
                    item.effective_sample_size
                ),
                "monte_carlo_standard_error": _posterior_estimate_record(
                    item.monte_carlo_standard_error
                ),
            }
            for item in summary.parameter_summaries
        ],
        "covariance": [
            {
                "parameter_a": item.parameter_a,
                "parameter_b": item.parameter_b,
                "value": float(item.value).hex(),
            }
            for item in summary.covariance
        ],
        "correlations": [
            {
                "parameter_a": item.parameter_a,
                "parameter_b": item.parameter_b,
                "estimate": _posterior_estimate_record(item.estimate),
            }
            for item in summary.correlations
        ],
        "acceptance_diagnostics": summary.acceptance.to_record(),
    }


def write_mcmc(
    path: Path,
    evidence: McmcEvidence | None,
    posterior_samples: PosteriorSampleEvidence | None,
    summary: PosteriorSummary | None,
) -> None:
    """Write primary MCMC and explicitly derived posterior artifacts."""
    if evidence is None:
        return
    if posterior_samples is None or summary is None:
        raise ValueError("Completed MCMC publication requires derived evidence")
    root = path / "MCMC"
    root.mkdir()
    write_json(
        root / "evidence.json",
        {"artifact_type": "native_mcmc_evidence", **evidence.to_record()},
    )
    write_json(
        root / "posterior-samples.json", _posterior_sample_record(posterior_samples)
    )
    write_json(root / "posterior-summary.json", _posterior_summary_record(summary))


def write_suppressed_outcome(
    path: Path,
    *,
    lifecycle: str,
    operation_record: dict[str, object],
    accepted_result_identity: str | None,
    accepted_occurrence_identity: str | None,
    components: Sequence[ComponentDiagnosticRecord],
) -> None:
    """Write minimal classified provenance for a suppressed workflow."""
    write_json(
        path,
        {
            "artifact_type": "native_suppressed_publication",
            "schema_version": 1,
            "lifecycle": lifecycle,
            "operation": operation_record,
            "accepted_result_identity": accepted_result_identity,
            "accepted_occurrence_identity": accepted_occurrence_identity,
            "components": component_diagnostic_records(components),
        },
    )


def write_partial_resampling(path: Path, evidence: ResamplingEvidence) -> None:
    """Write genuine partial primary resampling evidence."""
    write_json(path, resampling_evidence_record(evidence))


def write_resampling_summary_failure(path: Path, failure: SummaryFailure) -> None:
    """Write the typed reason completed resampling lacked a usable summary."""
    write_json(
        path,
        {
            "artifact_type": "native_resampling_summary_failure",
            "schema_version": 1,
            "identity": failure.identity,
            "source_evidence_identity": failure.source_evidence_identity,
            "category": failure.category,
            "message": failure.message,
        },
    )


def write_partial_mcmc(path: Path, evidence: McmcEvidence) -> None:
    """Write genuine partial primary MCMC evidence."""
    write_json(
        path,
        {"artifact_type": "native_partial_mcmc_evidence", **evidence.to_record()},
    )

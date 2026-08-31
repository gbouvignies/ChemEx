"""Authoritative interpretation of uncertainty for one deterministic fit.

Low-level numerical derivation and Evidence construction remain in
``chemex.optimize.uncertainty``.  This module owns the production policy,
scope, Evidence selection, recovery precedence, and final typed scientific
conclusions associated with one exact accepted deterministic fit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import StrEnum

from chemex.evaluation.native import EvaluationEngine
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    OptimizationProblem,
    accepted_occurrence_is_authoritative,
)
from chemex.optimize.grouped_direct_trf import FitPartitionProof
from chemex.optimize.uncertainty import (
    CompiledConstraintLinearizationCapabilities,
    MissingFunctionLinearizationCapability,
    OperationTerminal,
    ParameterUnit,
    ResidualVarianceScaling,
    RootAnchoredBlockCovarianceEvidence,
    UncertaintyConstructionError,
    UncertaintyEvidence,
    UncertaintyPolicy,
    UncertaintyUnavailableKind,
    compile_constraint_linearization_capabilities,
    derive_root_anchored_block_covariance,
    derive_uncertainty_evidence,
)
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ParameterRole,
)


class InterpretationCompleteness(StrEnum):
    """Whether the complete production interpretation finished."""

    COMPLETE = "complete"
    INCOMPLETE = "incomplete"


class DerivationDisposition(StrEnum):
    """Whether uncertainty derivation was evaluated or withheld by policy."""

    EVALUATED = "evaluated"
    WITHHELD = "withheld"


@dataclass(frozen=True, slots=True)
class ContinuousTrfBasis:
    """Existing partition proof for continuous Direct-TRF uncertainty."""

    partition_proof: FitPartitionProof = field(repr=False, compare=False)


@dataclass(frozen=True, slots=True)
class ProfiledGridBasis:
    """Marker for the existing profiled-GRID withholding policy."""


type DeterministicUncertaintyBasis = ContinuousTrfBasis | ProfiledGridBasis


@dataclass(frozen=True, slots=True)
class AcceptedDeterministicFitFacts:
    """Narrow authoritative inputs required for uncertainty interpretation."""

    accepted: AcceptedFitResult = field(repr=False, compare=False)
    problem: OptimizationProblem = field(repr=False, compare=False)
    parameterization: ActiveParameterization = field(repr=False, compare=False)
    engine: EvaluationEngine = field(repr=False, compare=False)
    basis: DeterministicUncertaintyBasis
    resolved_environment_identity: str


@dataclass(frozen=True, slots=True)
class ParameterUncertaintyConclusion:
    """Presentation-neutral reportability facts for one parameter."""

    param_id: str
    standard_error: float | None = None
    unavailable_kind: UncertaintyUnavailableKind | None = None
    boundary_warning: bool = False

    def __post_init__(self) -> None:
        if not self.param_id:
            raise UncertaintyConstructionError(
                "Parameter uncertainty requires a parameter ID"
            )
        if self.standard_error is not None and (
            not math.isfinite(self.standard_error) or self.standard_error < 0.0
        ):
            raise UncertaintyConstructionError(
                "Reportable parameter uncertainty must be finite and non-negative"
            )
        if self.standard_error is not None and self.unavailable_kind is not None:
            raise UncertaintyConstructionError(
                "Parameter uncertainty cannot be reportable and unavailable"
            )
        if self.boundary_warning and self.standard_error is None:
            raise UncertaintyConstructionError(
                "A boundary warning requires a reportable standard error"
            )

    @property
    def reportable(self) -> bool:
        """Whether this parameter has a reportable standard error."""
        return self.standard_error is not None


@dataclass(frozen=True, slots=True)
class DeterministicUncertainty:
    """Complete interpreted uncertainty for one accepted deterministic fit."""

    accepted_anchor: AcceptedFitResult = field(repr=False, compare=False)
    policy: UncertaintyPolicy = field(repr=False, compare=False)
    completeness: InterpretationCompleteness
    disposition: DerivationDisposition
    parameters: tuple[ParameterUncertaintyConclusion, ...]
    root_evidence: UncertaintyEvidence | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    block_evidence: RootAnchoredBlockCovarianceEvidence | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    incomplete_terminal: OperationTerminal | None = None

    def __post_init__(self) -> None:
        if not accepted_occurrence_is_authoritative(self.accepted_anchor):
            raise UncertaintyConstructionError(
                "Deterministic uncertainty requires an authoritative accepted result"
            )
        param_ids = tuple(item.param_id for item in self.parameters)
        if len(param_ids) != len(set(param_ids)):
            raise UncertaintyConstructionError(
                "Deterministic uncertainty parameter IDs must be unique"
            )
        self._validate_evidence_lineage()
        self._validate_interpretation_state()

    def _validate_evidence_lineage(self) -> None:
        if self.block_evidence is not None and self.root_evidence is None:
            raise UncertaintyConstructionError(
                "Block uncertainty Evidence requires root Evidence"
            )
        if self.root_evidence is not None:
            root = self.root_evidence
            if (
                root.accepted_anchor is not self.accepted_anchor
                or root.source_policy is not self.policy
            ):
                raise UncertaintyConstructionError(
                    "Root Evidence does not belong to this deterministic uncertainty"
                )
        if self.block_evidence is not None:
            block = self.block_evidence
            if block.source_bundle is not self.root_evidence:
                raise UncertaintyConstructionError(
                    "Block Evidence does not belong to this deterministic uncertainty"
                )

    def _validate_interpretation_state(self) -> None:
        if self.disposition is DerivationDisposition.WITHHELD:
            if (
                self.completeness is not InterpretationCompleteness.COMPLETE
                or self.root_evidence is not None
                or self.block_evidence is not None
                or self.incomplete_terminal is not None
                or any(
                    item.reportable
                    or item.unavailable_kind is not None
                    or item.boundary_warning
                    for item in self.parameters
                )
            ):
                raise UncertaintyConstructionError(
                    "Withheld uncertainty must be complete and carry no Evidence claims"
                )
            return
        if self.completeness is InterpretationCompleteness.INCOMPLETE:
            if self.incomplete_terminal is not OperationTerminal.INTERRUPTED:
                raise UncertaintyConstructionError(
                    "Incomplete evaluated uncertainty requires interruption"
                )
        elif self.incomplete_terminal is not None:
            raise UncertaintyConstructionError(
                "Complete uncertainty cannot carry an incomplete terminal"
            )
        if (
            self.completeness is InterpretationCompleteness.COMPLETE
            and self.root_evidence is None
        ):
            raise UncertaintyConstructionError(
                "Complete evaluated uncertainty requires root Evidence"
            )
        if any(
            not item.reportable and item.unavailable_kind is None
            for item in self.parameters
        ):
            raise UncertaintyConstructionError(
                "Evaluated unavailable uncertainty requires a typed classification"
            )

    def parameter(self, param_id: str) -> ParameterUncertaintyConclusion | None:
        """Return the predetermined conclusion for one parameter, if relevant."""
        return next(
            (item for item in self.parameters if item.param_id == param_id), None
        )

    @property
    def boundary_warning(self) -> bool:
        """Whether interpreted root or recovered-block Evidence has a caveat."""
        root_warning = (
            self.root_evidence is not None
            and self.root_evidence.covariance is not None
            and self.root_evidence.covariance.boundary_warning
        )
        block_warning = self.block_evidence is not None and any(
            block.boundary_warning for block in self.block_evidence.blocks
        )
        return root_warning or block_warning

    @property
    def root_covariance_available(self) -> bool:
        """Whether the joint root covariance is available for all controls."""
        return (
            self.root_evidence is not None
            and self.root_evidence.covariance is not None
            and all(
                (conclusion := self.parameter(param_id)) is not None
                and conclusion.reportable
                for param_id in self.accepted_anchor.controlled_ids
            )
        )

    @property
    def block_covariance_available(self) -> bool:
        """Whether any root-anchored recovery block has covariance."""
        return self.block_evidence is not None and any(
            block.covariance is not None for block in self.block_evidence.blocks
        )


@dataclass(frozen=True, slots=True)
class _ResolvedUncertaintyInputs:
    """Private resolved production policy and constrained-output scope."""

    policy: UncertaintyPolicy
    constrained_scope: tuple[str, ...]
    compiled_capabilities: CompiledConstraintLinearizationCapabilities
    unsupported_constrained_ids: tuple[str, ...]

    @property
    def parameter_ids(self) -> tuple[str, ...]:
        return self.constrained_scope + self.unsupported_constrained_ids


def _validate_facts(facts: AcceptedDeterministicFitFacts) -> None:
    accepted = facts.accepted
    problem = facts.problem
    parameterization = facts.parameterization
    engine = facts.engine
    if not facts.resolved_environment_identity:
        raise UncertaintyConstructionError(
            "Deterministic uncertainty requires a resolved environment identity"
        )
    if not accepted_occurrence_is_authoritative(accepted):
        raise UncertaintyConstructionError(
            "Deterministic uncertainty requires an exact authoritative occurrence"
        )
    if (
        accepted.problem_identity != problem.identity
        or accepted.parameterization_identity != problem.parameterization_identity
        or accepted.evaluator_parameterization_identity
        != problem.evaluator_parameterization_identity
        or accepted.controlled_ids != problem.controlled_ids
        or accepted.source_occurrence_identity
        != problem.source_snapshot.occurrence_identity
        or accepted.source_revision != problem.source_snapshot.revision
        or accepted.evaluation_result.plan_identity != problem.evaluation_plan_identity
    ):
        raise UncertaintyConstructionError(
            "Accepted result and uncertainty problem lineage differ"
        )
    try:
        problem.validate_parameterization(parameterization)
    except ValueError as error:
        raise UncertaintyConstructionError(
            "Uncertainty parameterization lineage differs from the problem"
        ) from error
    if (
        engine.plan.identity != problem.evaluation_plan_identity
        or engine.plan.parameterization_identity
        != problem.evaluator_parameterization_identity
        or engine.plan.constraint_program_identity
        != problem.constraint_program_identity
    ):
        raise UncertaintyConstructionError(
            "Uncertainty evaluation engine lineage differs from the problem"
        )
    if isinstance(facts.basis, ContinuousTrfBasis):
        proof = facts.basis.partition_proof
        if (
            proof.root_plan_identity != engine.plan.identity
            or proof.constraint_program_identity != parameterization.program.fingerprint
            or proof.controlled_ids != problem.controlled_ids
        ):
            raise UncertaintyConstructionError(
                "Continuous uncertainty requires the exact fit partition proof"
            )


def _resolve_inputs(
    facts: AcceptedDeterministicFitFacts,
) -> _ResolvedUncertaintyInputs:
    problem = facts.problem
    parameterization = facts.parameterization
    policy = UncertaintyPolicy(
        calibration_identity="native-product-local-covariance-numerics-v2",
        numerical_compatibility_requirement=(
            "binary64-retained-scipy-2point-gesdd-column-equilibrated-v1"
        ),
        coordinate_scales=tuple((param_id, 1.0) for param_id in problem.controlled_ids),
        coordinate_units=tuple(
            (param_id, ParameterUnit.UNSPECIFIED) for param_id in problem.controlled_ids
        ),
        residual_variance_scaling=(
            ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES
        ),
        relative_step_tolerance=1.0e-4,
        roundoff_multiplier=64.0,
        smaller_step_extent=8,
        larger_step_extent=8,
        svd_driver="gesdd",
        rank_absolute_tolerance=0.0,
        rank_relative_tolerance=0.0,
        weak_relative_tolerance=1.0e-6,
        singular_value_cluster_relative_tolerance=1.0e-10,
        conditioning_limit=1.0e12,
        correlation_roundoff_multiplier=64.0,
        affine_feasibility_policy="canonical-root-affine-halfspace-zeta-gt-3-v1",
        residual_jacobian_strategy="retained-backend-or-accepted-2-point",
    )
    constraints = {
        item.target_id: item for item in parameterization.program.constraints
    }
    controlled = frozenset(problem.controlled_ids)

    def depends_on_controlled(
        param_id: str,
        visiting: frozenset[str] = frozenset(),
    ) -> bool:
        if param_id in controlled:
            return True
        if param_id in visiting or param_id not in constraints:
            return False
        constraint = constraints[param_id]
        return any(
            depends_on_controlled(dependency, visiting | {param_id})
            for dependency in constraint.dependencies
        )

    propagation_candidates = tuple(
        param_id
        for param_id in parameterization.scope_ids
        if parameterization.role(param_id) is ParameterRole.DERIVED
        and depends_on_controlled(param_id)
    )
    supported: list[str] = []
    unsupported: list[str] = []
    for param_id in propagation_candidates:
        try:
            compile_constraint_linearization_capabilities(
                parameterization,
                (param_id,),
                (),
            )
        except MissingFunctionLinearizationCapability:
            unsupported.append(param_id)
        else:
            supported.append(param_id)
    constrained_scope = tuple(supported)
    compiled_capabilities = compile_constraint_linearization_capabilities(
        parameterization,
        constrained_scope,
        (),
    )
    return _ResolvedUncertaintyInputs(
        policy,
        constrained_scope,
        compiled_capabilities,
        tuple(unsupported),
    )


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
    if "active_relaxation_feasibility_boundary" in categories:
        return UncertaintyUnavailableKind.BOUNDARY_LIMITED
    covariance = evidence.covariance
    if covariance is not None:
        claims = {item.name: item.state.value for item in covariance.claims}
        if claims.get("PROFILED_NORMALIZATION_REGULARITY") == "violated":
            return UncertaintyUnavailableKind.NORMALIZATION_INVALID
    if any(failure.stage == "residual_linearization" for failure in evidence.failures):
        return UncertaintyUnavailableKind.JACOBIAN_UNAVAILABLE
    return UncertaintyUnavailableKind.COVARIANCE_NUMERICAL_FAILURE


def _constrained_boundary_warning_ids(
    evidence: UncertaintyEvidence,
    controlled_warning_ids: frozenset[str],
) -> frozenset[str]:
    constraint_jacobian = evidence.constraint_jacobian
    if constraint_jacobian is None or not controlled_warning_ids:
        return frozenset()
    return frozenset(
        output_id
        for output_id, dependencies in zip(
            constraint_jacobian.output_ids,
            constraint_jacobian.structural_dependencies,
            strict=True,
        )
        if controlled_warning_ids.intersection(dependencies)
    )


def _interpret_evidence(
    evidence: UncertaintyEvidence,
    block_evidence: RootAnchoredBlockCovarianceEvidence | None,
    inputs: _ResolvedUncertaintyInputs,
) -> tuple[ParameterUncertaintyConclusion, ...]:
    standard_errors: dict[str, float] = {}
    unavailable: dict[str, UncertaintyUnavailableKind] = {}
    warnings: set[str] = set()
    controlled_warning_ids = (
        frozenset()
        if evidence.covariance is None
        else evidence.covariance.simple_bound_warning_ids
    )
    marginal = evidence.marginal_errors
    if marginal is not None and marginal.scope_reportable:
        standard_errors.update(
            (entry.param_id, entry.value)
            for entry in marginal.entries
            if entry.value is not None
        )
        warnings.update(
            entry.param_id
            for entry in marginal.entries
            if entry.value is not None and entry.param_id in controlled_warning_ids
        )
    else:
        kind = _controlled_unavailability_kind(evidence)
        unavailable.update(
            (param_id, kind) for param_id in evidence.accepted_anchor.controlled_ids
        )
    constrained = evidence.constrained_marginal_errors
    if constrained is not None and constrained.scope_reportable:
        standard_errors.update(
            (entry.param_id, entry.value)
            for entry in constrained.entries
            if entry.value is not None
        )
        constrained_warning_ids = _constrained_boundary_warning_ids(
            evidence,
            controlled_warning_ids,
        )
        warnings.update(
            entry.param_id
            for entry in constrained.entries
            if entry.value is not None and entry.param_id in constrained_warning_ids
        )
    elif evidence.requested_output_scope:
        unavailable.update(
            (
                param_id,
                UncertaintyUnavailableKind.CONSTRAINED_PROPAGATION_UNAVAILABLE,
            )
            for param_id in evidence.requested_output_scope
        )
    unavailable.update(
        (
            param_id,
            UncertaintyUnavailableKind.UNSUPPORTED_CONSTRAINED_DERIVATIVE,
        )
        for param_id in inputs.unsupported_constrained_ids
    )
    if block_evidence is not None:
        block_controlled_warning_ids = block_evidence.simple_bound_warning_ids
        block_constrained_warning_ids = _constrained_boundary_warning_ids(
            block_evidence.source_bundle,
            block_controlled_warning_ids,
        )
        for block in block_evidence.blocks:
            recovered_ids: set[str] = set()
            for entry in (*block.marginal_errors, *block.constrained_marginal_errors):
                if entry.value is not None:
                    standard_errors[entry.param_id] = entry.value
                    recovered_ids.add(entry.param_id)
            for param_id in recovered_ids:
                unavailable.pop(param_id, None)
            if block.unavailable_kind is not None:
                unavailable.update(
                    (param_id, block.unavailable_kind)
                    for param_id in block.controlled_ids
                )
            warnings.update(
                entry.param_id
                for entry in block.marginal_errors
                if entry.value is not None
                and entry.param_id in block_controlled_warning_ids
            )
            warnings.update(
                entry.param_id
                for entry in block.constrained_marginal_errors
                if entry.value is not None
                and entry.param_id in block_constrained_warning_ids
            )
    parameter_ids = evidence.accepted_anchor.controlled_ids + inputs.parameter_ids
    return tuple(
        ParameterUncertaintyConclusion(
            param_id,
            standard_errors.get(param_id),
            unavailable.get(param_id),
            param_id in warnings,
        )
        for param_id in parameter_ids
    )


def _interrupted_conclusions(
    facts: AcceptedDeterministicFitFacts,
    inputs: _ResolvedUncertaintyInputs,
) -> tuple[ParameterUncertaintyConclusion, ...]:
    return tuple(
        ParameterUncertaintyConclusion(
            param_id,
            unavailable_kind=(
                UncertaintyUnavailableKind.UNSUPPORTED_CONSTRAINED_DERIVATIVE
                if param_id in inputs.unsupported_constrained_ids
                else UncertaintyUnavailableKind.DERIVATION_STOPPED
            ),
        )
        for param_id in facts.problem.controlled_ids + inputs.parameter_ids
    )


def derive_deterministic_uncertainty(
    facts: AcceptedDeterministicFitFacts,
) -> DeterministicUncertainty:
    """Interpret uncertainty for one exact accepted deterministic fit."""
    _validate_facts(facts)
    inputs = _resolve_inputs(facts)
    parameter_ids = facts.problem.controlled_ids + inputs.parameter_ids
    if isinstance(facts.basis, ProfiledGridBasis):
        return DeterministicUncertainty(
            facts.accepted,
            inputs.policy,
            InterpretationCompleteness.COMPLETE,
            DerivationDisposition.WITHHELD,
            tuple(
                ParameterUncertaintyConclusion(param_id) for param_id in parameter_ids
            ),
        )
    try:
        root_evidence = derive_uncertainty_evidence(
            facts.accepted,
            problem=facts.problem,
            parameterization=facts.parameterization,
            engine=facts.engine,
            policy=inputs.policy,
            constrained_scope=inputs.constrained_scope,
            constrained_units=tuple(
                (param_id, ParameterUnit.UNSPECIFIED)
                for param_id in inputs.constrained_scope
            ),
            constrained_scales=tuple(
                (param_id, 1.0) for param_id in inputs.constrained_scope
            ),
            compiled_constraint_linearization=inputs.compiled_capabilities,
            resolved_environment_identity=facts.resolved_environment_identity,
        )
    except KeyboardInterrupt:
        return DeterministicUncertainty(
            facts.accepted,
            inputs.policy,
            InterpretationCompleteness.INCOMPLETE,
            DerivationDisposition.EVALUATED,
            _interrupted_conclusions(facts, inputs),
            incomplete_terminal=OperationTerminal.INTERRUPTED,
        )
    try:
        block_evidence = derive_root_anchored_block_covariance(
            root_evidence,
            facts.basis.partition_proof,
        )
    except KeyboardInterrupt:
        return DeterministicUncertainty(
            facts.accepted,
            inputs.policy,
            InterpretationCompleteness.INCOMPLETE,
            DerivationDisposition.EVALUATED,
            _interpret_evidence(root_evidence, None, inputs),
            root_evidence,
            incomplete_terminal=OperationTerminal.INTERRUPTED,
        )
    return DeterministicUncertainty(
        facts.accepted,
        inputs.policy,
        InterpretationCompleteness.COMPLETE,
        DerivationDisposition.EVALUATED,
        _interpret_evidence(root_evidence, block_evidence, inputs),
        root_evidence,
        block_evidence,
    )

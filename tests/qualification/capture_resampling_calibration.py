"""Frozen prospective acquisition specification for issue #603."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
from pydantic import BaseModel

from chemex.containers.data import Data
from chemex.containers.profile import Profile, PulseSequence
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.migration_core import migration_core_authority_selection
from chemex.nmr.spectrometer import Spectrometer
from chemex.numerical_lanes import (
    LaneAttestation,
    NumericalLane,
    RuntimeEnvironment,
    canonical_lanes,
)
from chemex.optimize.direct_trf import AcceptedFitResult, OptimizationProblem
from chemex.optimize.native_resampling import (
    OperationTerminal,
    OptimizationStrategy,
    ReplicateDisposition,
    ResamplingDatasetManifest,
    ResamplingLifecycle,
    ResamplingOperation,
    ResamplingPlan,
    ResamplingScheme,
    ResamplingSummaryOutcome,
    ResamplingSummaryPolicy,
    SummaryTerminal,
    execute_resampling_evidence,
    summarize_resampling_evidence,
)
from chemex.parameters.parameterization import (
    ActiveParameterization,
    ConstraintProgram,
    ParameterRole,
    ScientificFunctionBinder,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.printers.data import Printer
from chemex.runtime import ExecutionSettings

# Frozen qualification specifications are intentionally kept compact and inspectable.
# fmt: off

SPECIFICATION_ID = "chemex-resampling-calibration-v2"
CANDIDATE_COUNTS = (64, 128, 256)
OUTPUT_SCOPE = ("A", "B")
THETA = (1.0, 2.0)
TRUTH_FIELDS = ("mean", "bias", "spread", "q2.5", "median", "q97.5", "failure_prevalence")
X = (-2, -1, 0, 1, 2, 3)
CALCULATED = (-3, -1, 1, 3, 5, 7)
RESIDUAL = (1, -2, 1, -1, 2, -1)
OBSERVED = (-2, -3, 2, 2, 7, 6)
ERRORS = (1, 1, 1, 1, 1, 1)
REFERENCES = (True, False, False, True, False, False)
NUCLEUS_GROUPS = ("N1", "N2", "N3", "N1", "N2", "N3")
PROFILE_INDICES = ((0, 1, 2), (3, 4, 5))
SUMMARY_POLICY = ResamplingSummaryPolicy()
SUMMARY_POLICY_IDENTITY = "83b0a16dbfd6df67337bd7d17ed9fc7cd8e82283688fa0946fce7d88edf7ae44"
THRESHOLDS = {
    "bias": 0.50,
    "spread": 0.35,
    "q2.5": 0.90,
    "median": 0.50,
    "q97.5": 0.90,
    "failure_prevalence": 0.0,
}
FAMILY_SCHEMES = {
    "mc": ResamplingScheme.MONTE_CARLO,
    "bs": ResamplingScheme.BOOTSTRAP,
    "bsn": ResamplingScheme.NUCLEUS_BOOTSTRAP,
}


def _roots(literals: str) -> tuple[int, ...]:
    return tuple(int(value) for value in literals.split())


ROOTS = {
    "mc": {
        "calibration": _roots("3352135948555585357 2472045721331367299 18179771078494673657 7172209355565492995 4428942634759812093 16846731500777651301 1276328454428800799 1007434282609839600 7623382162802148034 17314847138863603676 17123956007879464799 9852929745021454887 18118882682188821729 1998968888814372653 11174296571166983796 9342205682851433915"),
        "holdout": _roots("1683597629452963050 16208179182483689898 8138106942819703130 15515473093383081872 2976836775184601562 14083146987549980469 5591399098166460148 13265399871182262342"),
    },
    "bs": {
        "calibration": _roots("17119401643894300767 16464330990378071471 3287997670289142773 13004064220656459317 9175610863457614675 13493705660742239536 7566518655171357174 9129588113907940364 15581839316159386495 7621495334857864309 9939444681066633173 15450513739920280325 13405509690366816759 8661374386535711639 7796930234651015581 15997482312548809078"),
        "holdout": _roots("8841647693099355132 14673849892314593488 62288681685728296 9553479069377181116 1616470759781790097 3438610658747297316 796585017052296001 10339959776139241745"),
    },
    "bsn": {
        "calibration": _roots("7639148882328525895 12072455617509212606 8535324653541211058 9476519509736332819 5478984577077683145 2525069006632874260 1105367432024468405 3898413442510559254 3384662453840932065 14315773297000690354 16176071534482763393 15507941143551554954 9258971447573535519 5287528659740383302 1416814305191906398 14303118813658220336"),
        "holdout": _roots("10916365485830788432 16652394355300541660 2840841312473151217 11985423680230601622 4282993254360459241 16406857506768168855 12066811090466048655 677902911216112496"),
    },
}
V1_HOLDOUT_ROOTS = {
    "mc": _roots("11227737360748962646 17470780368220950770 882017822349353614 6006908154484889383 682688700279365841 8190374274261472580 18196953354430856719 15474585266757114297"),
    "bs": _roots("10049477911408417482 12926280384734704413 3894811774098425923 14835158499454950328 14745417758411211654 16213361473081760634 9514719903799318634 2156741863032581007"),
    "bsn": _roots("12671600559231941323 10264161786616261730 9189865121253494806 7938868261418550067 6646078350059218959 8516678745947789619 18038196089002735289 9413147880661201059"),
}
REFERENCE_TRUTH = {
    "mc": {
        "A": (1.0, 0.0, 0.4253849796976627, 0.16626076022827885, 1.0, 1.833739239771721, 0.0),
        "B": (2.0, 0.0, 0.23904572186687872, 1.5314789944825389, 2.0, 2.468521005517461, 0.0),
    },
    "bs": {
        "A": (1.0318488108646653, 0.0318488108646653, 0.4937978568774975, 0.04255319148936169, 1.0, 1.9130434782608696, 0.0),
        "B": (2.0273454890306235, 0.0273454890306235, 0.28735337973018993, 1.5221238938053097, 2.0, 2.6067415730337076, 0.0),
    },
    "bsn": {
        "A": (0.9856816147956649, -0.0143183852043351, 0.3244618767900137, 0.3360174781523096, 1.0, 1.989970501474927, 0.0),
        "B": (2.0286367704086707, 0.0286367704086707, 0.48010532709585624, 1.3333333333333333, 2.0, 3.3212234706616734, 0.0),
    },
}


def derive_root(family: str, phase: str, index: int) -> int:
    version = {"calibration": "v1", "holdout": "v2"}[phase]
    label = f"chemex-issue-603-resampling-calibration-{version}|{family}|{phase}|{index:02d}"
    return int.from_bytes(hashlib.sha256(label.encode("ascii")).digest()[:8], "big")


def type8_quantile(values: tuple[Fraction, ...], probability: Fraction) -> Fraction:
    ordered = tuple(sorted(values))
    h = (Fraction(len(ordered)) + Fraction(1, 3)) * probability + Fraction(1, 3)
    if h <= 1:
        return ordered[0]
    if h >= len(ordered):
        return ordered[-1]
    lower = h.numerator // h.denominator
    return ordered[lower - 1] + (h - lower) * (ordered[lower] - ordered[lower - 1])


def _ols(indices: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    xs = tuple(Fraction(X[index]) for index in indices)
    ys = tuple(Fraction(OBSERVED[index]) for index in indices)
    size, sum_x, sum_y = Fraction(len(indices)), sum(xs), sum(ys)
    sum_xx = sum(value * value for value in xs)
    sum_xy = sum(x * y for x, y in zip(xs, ys, strict=True))
    determinant = size * sum_xx - sum_x * sum_x
    return (
        (sum_y * sum_xx - sum_x * sum_xy) / determinant,
        (size * sum_xy - sum_x * sum_y) / determinant,
    )


def _finite_truth(
    samples: tuple[tuple[Fraction, Fraction], ...],
) -> dict[str, tuple[float, ...]]:
    result = {}
    for parameter_index, (parameter_id, theta) in enumerate(zip(OUTPUT_SCOPE, THETA, strict=True)):
        values = tuple(sample[parameter_index] for sample in samples)
        mean = sum(values, Fraction()) / len(values)
        variance = sum((value - mean) ** 2 for value in values) / len(values)
        quantiles = tuple(float(type8_quantile(values, probability)) for probability in (Fraction(1, 40), Fraction(1, 2), Fraction(39, 40)))
        result[parameter_id] = (float(mean), float(mean - Fraction(theta)), math.sqrt(float(variance)), *quantiles, 0.0)
    return result


def truth_estimands() -> dict[str, dict[str, tuple[float, ...]]]:
    inverse_diagonal = (Fraction(19, 105), Fraction(2, 35))
    z975 = 1.959963984540054
    mc = {}
    for parameter_id, theta, variance in zip(OUTPUT_SCOPE, THETA, inverse_diagonal, strict=True):
        spread = math.sqrt(float(variance))
        mc[parameter_id] = (theta, 0.0, spread, theta - z975 * spread, theta, theta + z975 * spread, 0.0)
    bs = tuple(_ols((0, *left, 3, *right)) for left in product((1, 2), repeat=2) for right in product((4, 5), repeat=2))
    groups = ((0, 3), (1, 4), (2, 5))
    bsn = tuple(_ols(tuple(index for group in draw for index in groups[group])) for draw in product(range(3), repeat=3))
    return {"mc": mc, "bs": _finite_truth(bs), "bsn": _finite_truth(bsn)}


def select_candidate(passing_counts: set[int]) -> tuple[str, int | None]:
    table = {
        frozenset(): ("NO_ADEQUATE_CANDIDATE", None),
        frozenset({256}): ("GRID_SATURATED", None),
        frozenset({128, 256}): ("SELECTED", 128),
        frozenset({64, 128, 256}): ("SELECTED", 64),
    }
    if not passing_counts <= set(CANDIDATE_COUNTS):
        return "NON_MONOTONE_ADEQUACY", None
    return table.get(frozenset(passing_counts), ("NON_MONOTONE_ADEQUACY", None))


def holdout_decision(selected: int, passed: tuple[bool, ...]) -> tuple[str, int | None]:
    if len(passed) == 8 and all(passed):
        return "SELECTED", selected
    return "HOLDOUT_FAILED_UNSUPPORTED", None


def replay_signature(operation: ResamplingOperation, outcome: ResamplingSummaryOutcome) -> dict[str, object]:
    evidence = operation.evidence
    summary = outcome.summary
    if evidence is None or summary is None:
        return {"unsupported": True}
    return {
        "operation_terminal": operation.terminal.value,
        "unstarted_ordinals": operation.unstarted_ordinals,
        "plan_identity": evidence.plan.identity,
        "evidence_identity": evidence.identity,
        "population_identity": evidence.population_identity,
        "lifecycle": evidence.lifecycle.value,
        "counts": (evidence.completed_count, evidence.successful_count, evidence.failed_count),
        "outcomes": tuple((item.ordinal, item.request_identity, item.seed, item.draw_identity, item.stage.value, item.disposition.value, item.success.identity if item.success else item.failure.identity, item.identity) for item in evidence.outcomes),
        "summary_terminal": outcome.terminal.value,
        "summary_record": summary.to_record(),
    }


def replay_equal(left: tuple[ResamplingOperation, ResamplingSummaryOutcome], right: tuple[ResamplingOperation, ResamplingSummaryOutcome]) -> bool:
    left_signature, right_signature = replay_signature(*left), replay_signature(*right)
    return "unsupported" not in left_signature and left_signature == right_signature


class _KernelSettings(BaseModel):
    kind: str = "issue-603-linear-kernel"


class _LinearSpectrometer:
    def __init__(self, name: str) -> None:
        self.spin_system = SpinSystem.from_name(name)
        self.values = {"a": 0.0, "b": 0.0}

    def update(self, values: dict[str, float]) -> None:
        self.values = dict(values)

    def new_native_workspace(self) -> _LinearSpectrometer:
        return _LinearSpectrometer(str(self.spin_system))

    def native_kernel_descriptor(self) -> dict[str, object]:
        return {"kind": "issue-603-linear-spectrometer"}


class _LinearPulseSequence:
    settings = _KernelSettings()

    def calculate(self, spectrometer: _LinearSpectrometer, data: Data) -> np.ndarray:
        return spectrometer.values["a"] + spectrometer.values["b"] * np.asarray(data.metadata, dtype=np.float64)

    def is_reference(self, metadata: np.ndarray) -> np.ndarray:
        return metadata == metadata[0]


def _profile(name: str, indices: tuple[int, ...]) -> Profile:
    data = Data(exp=np.asarray([OBSERVED[index] for index in indices], dtype=np.float64), err=np.ones(len(indices)), metadata=np.asarray([X[index] for index in indices], dtype=np.float64))
    return Profile(data, cast("Spectrometer", _LinearSpectrometer(name)), cast("PulseSequence", _LinearPulseSequence()), {"a": "A", "b": "B"}, cast("Printer", None), is_scaled=False)


@dataclass(frozen=True, slots=True)
class NativeFixture:
    accepted: AcceptedFitResult
    dataset: ResamplingDatasetManifest
    problem: OptimizationProblem
    parameterization: ActiveParameterization
    engine: EvaluationEngine


def native_fixture() -> NativeFixture:
    binder = ScientificFunctionBinder("issue-603", {})
    program = ConstraintProgram("parameter-model", "model", "definitions", "configuration", binder.identity, OUTPUT_SCOPE, OUTPUT_SCOPE, (), (), ())
    parameterization = ActiveParameterization(program, binder, "source-occurrence", 603, (("A", ParameterRole.FIT), ("B", ParameterRole.FIT)))
    snapshot = AnalysisValuesSnapshot("source-occurrence", "model", "definitions", "configuration", 603, tuple(zip(OUTPUT_SCOPE, THETA, strict=True)))
    profiles = tuple(_profile(f"P{ordinal}", indices) for ordinal, indices in enumerate(PROFILE_INDICES, 1))
    engine = EvaluationEngine.from_experiments(cast("Any", (SimpleNamespace(profiles=profiles),)), parameterization)
    evaluated = engine.new_evaluator().evaluate(EvaluationFrame.from_lifecycle_frame(parameterization, parameterization.frame_from_snapshot(snapshot)))
    assert isinstance(evaluated, EvaluationResult)
    dataset = ResamplingDatasetManifest(engine.plan, tuple(float(value) for value in evaluated.normalized_calculations), REFERENCES, NUCLEUS_GROUPS, tuple(f"index={index}" for index in range(len(X))))
    problem = OptimizationProblem(engine.plan.identity, parameterization.identity, parameterization.evaluator_identity, program.fingerprint, "configuration", snapshot, tuple(zip(OUTPUT_SCOPE, THETA, strict=True)), OUTPUT_SCOPE, (), THETA, (-100.0, -100.0), (100.0, 100.0), OUTPUT_SCOPE)
    accepted = AcceptedFitResult.for_qualification(
        occurrence_identity="issue-603-accepted",
        problem_identity=problem.identity,
        invocation_identity="issue-603-invocation",
        execution_identity="issue-603-execution",
        materialization_identity="issue-603-materialization",
        parameterization_identity=parameterization.identity,
        evaluator_parameterization_identity=parameterization.evaluator_identity,
        source_occurrence_identity=snapshot.occurrence_identity,
        source_revision=snapshot.revision,
        controlled_ids=problem.controlled_ids,
        vector=problem.start,
        chi_square=float(np.dot(evaluated.residuals, evaluated.residuals)),
        evaluation_result=evaluated,
        commit_scope=problem.commit_scope,
        commit_items=problem.independent_items,
        origin_context_identity="issue-603-origin",
    )
    assert dataset.calculated == tuple(float(value) for value in CALCULATED)
    return NativeFixture(accepted, dataset, problem, parameterization, engine)


def make_plan(fixture: NativeFixture, family: str, count: int, root: int) -> ResamplingPlan:
    return ResamplingPlan.for_accepted(fixture.accepted, dataset=fixture.dataset, source_problem=fixture.problem, parameterization=fixture.parameterization, source_engine=fixture.engine, scheme=FAMILY_SCHEMES[family], replicate_count=count, replicate_structural_identities=tuple(f"replicate-{index:03d}" for index in range(count)), replicate_component_identities=tuple((f"component-{index:03d}",) for index in range(count)), root_seed=root, output_scope=OUTPUT_SCOPE, output_units=("arbitrary", "arbitrary"), minimum_successful_count=count, strategy=OptimizationStrategy.DIRECT_TRF)


def assess_run(family: str, count: int, root: int, operation: ResamplingOperation, outcome: ResamplingSummaryOutcome) -> dict[str, object]:
    evidence, summary = operation.evidence, outcome.summary
    reasons = []
    if evidence is None or summary is None:
        reasons.append("missing_evidence_or_summary")
    else:
        plan = evidence.plan
        expected_ordinals = tuple(range(count))
        if operation.terminal is not OperationTerminal.COMPLETED or evidence.lifecycle is not ResamplingLifecycle.COMPLETED or evidence.completed_count != count or evidence.successful_count != count or evidence.failed_count != 0:
            reasons.append("incomplete_or_unsuccessful_population")
        if plan.scheme is not FAMILY_SCHEMES[family] or plan.replicate_count != count or plan.root_seed != root or plan.output_scope != OUTPUT_SCOPE or plan.minimum_successful_count != count or plan.strategy is not OptimizationStrategy.DIRECT_TRF or plan.strategy_settings:
            reasons.append("plan_mismatch")
        if outcome.terminal is not SummaryTerminal.COMPLETED or summary.policy_identity != SUMMARY_POLICY_IDENTITY or summary.output_scope != OUTPUT_SCOPE or summary.included_ordinals != expected_ordinals or summary.unstarted_ordinals or summary.exclusions or tuple(item.ordinal for item in summary.samples) != expected_ordinals or any(item.disposition is not ReplicateDisposition.SUCCEEDED for item in evidence.outcomes):
            reasons.append("summary_or_policy_mismatch")
        try:
            evidence.validate_integrity()
            summary.validate_integrity()
        except Exception as error:  # noqa: BLE001 - acquisition records typed disqualifiers
            reasons.append(f"integrity:{type(error).__name__}")
    if reasons:
        return {"status": "DISQUALIFIED", "reasons": reasons}
    truth = truth_estimands()[family]
    distributions = {item.parameter_id: item for item in summary.distributions}
    failure = abs((count - evidence.successful_count) / count - truth["A"][6])
    errors_by_parameter = {}
    for parameter_id in OUTPUT_SCOPE:
        exact, distribution = truth[parameter_id], distributions[parameter_id]
        scale = exact[2]
        errors_by_parameter[parameter_id] = {"bias": abs(distribution.mean - exact[0]) / scale, "spread": abs(distribution.standard_deviation - scale) / scale, "q2.5": abs(distribution.percentile_95_lower - exact[3]) / scale, "median": abs(distribution.median - exact[4]) / scale, "q97.5": abs(distribution.percentile_95_upper - exact[5]) / scale, "failure_prevalence": failure}
    maxima = {name: max(values[name] for values in errors_by_parameter.values()) for name in THRESHOLDS}
    if not all(math.isfinite(value) for value in maxima.values()):
        return {"status": "DISQUALIFIED", "reasons": ["non_finite_metric"], "errors": errors_by_parameter, "maxima": maxima}
    violations = tuple(name for name, value in maxima.items() if value > THRESHOLDS[name])
    return {"status": "PASS" if not violations else "DISQUALIFIED", "reasons": [f"threshold:{name}" for name in violations], "errors": errors_by_parameter, "maxima": maxima}


def candidate_maxima(runs: list[dict[str, object]]) -> dict[str, float]:
    return {name: max(cast("dict[str, float]", run["maxima"])[name] for run in runs) for name in THRESHOLDS}


def _run(fixture: NativeFixture, family: str, count: int, root: int, workers: int) -> tuple[dict[str, object], ResamplingOperation, ResamplingSummaryOutcome]:
    plan = make_plan(fixture, family, count, root)
    operation = execute_resampling_evidence(fixture.accepted, plan, execution=ExecutionSettings(workers=workers, native_threads=1))
    assert operation.evidence is not None
    outcome = summarize_resampling_evidence(operation.evidence, SUMMARY_POLICY)
    record = assess_run(family, count, root, operation, outcome)
    record.update({"count": count, "root": root, "workers": workers, "evidence_identity": operation.evidence.identity, "summary_identity": None if outcome.summary is None else outcome.summary.identity})
    return record, operation, outcome


def validate_canonical_lane_records(records: Mapping[str, object]) -> None:
    names = {"numerical_lane", "lane_attestation", "runtime_environment"}
    if set(records) != names or any(not isinstance(records[name], Mapping) for name in names):
        raise RuntimeError("canonical acquisition requires complete typed lane records")
    lane = NumericalLane.from_record(cast("Mapping[str, object]", records["numerical_lane"]))
    attestation = LaneAttestation.from_record(cast("Mapping[str, object]", records["lane_attestation"]))
    environment = RuntimeEnvironment.from_record(cast("Mapping[str, object]", records["runtime_environment"]))
    expected_lane, selection = canonical_lanes()[0], migration_core_authority_selection()
    if (
        lane != expected_lane
        or lane.identity != selection.lane_identity
        or lane.role != selection.lane_role
        or lane.semantics.image_digest != selection.image_digest
        or attestation.identity != selection.attestation_identity
        or attestation.lane_identity != lane.identity
        or attestation.environment_identity != selection.environment_identity
        or attestation.workers != selection.workers
        or attestation.native_threads != selection.native_threads
        or environment.identity != selection.environment_identity
        or environment.semantics != lane.semantics
    ):
        raise RuntimeError("canonical acquisition lane records do not match #588")


def attest_canonical_lane(image_digest: str) -> dict[str, object]:
    lane = canonical_lanes()[0]
    authority = lane.attest_current_process(image_digest)
    records = {
        "numerical_lane": lane.to_record(),
        "lane_attestation": authority.to_record(),
        "runtime_environment": RuntimeEnvironment(lane.semantics).to_record(),
    }
    validate_canonical_lane_records(records)
    return records


def acquire(image_digest: str) -> dict[str, object]:
    lane_records = attest_canonical_lane(image_digest)
    if SUMMARY_POLICY.identity != SUMMARY_POLICY_IDENTITY:
        raise RuntimeError("frozen summary policy identity changed")
    fixture, families = native_fixture(), {}
    for family in FAMILY_SCHEMES:
        calibration, serial = {}, {}
        passing = set()
        for count in CANDIDATE_COUNTS:
            runs = []
            for root in ROOTS[family]["calibration"]:
                record, operation, outcome = _run(fixture, family, count, root, 1)
                runs.append(record)
                if root == ROOTS[family]["calibration"][0]:
                    serial[count] = (operation, outcome)
            aggregate = None if any("maxima" not in run for run in runs) else candidate_maxima(runs)
            calibration[str(count)] = {"runs": runs, "maxima": aggregate}
            if all(run["status"] == "PASS" for run in runs) and aggregate is not None and all(value <= THRESHOLDS[name] for name, value in aggregate.items()):
                passing.add(count)
        status, selected = select_candidate(passing)
        family_record: dict[str, object] = {"calibration": calibration, "selection_status": status, "selected_count": selected}
        if selected is not None:
            replay_record, replay_operation, replay_outcome = _run(fixture, family, selected, ROOTS[family]["calibration"][0], 2)
            replay_ok = replay_equal(serial[selected], (replay_operation, replay_outcome))
            family_record["replay"] = {"passed": replay_ok, "run": replay_record}
            if not replay_ok:
                family_record.update(selection_status="REPLAY_MISMATCH_UNSUPPORTED", selected_count=None)
            else:
                holdouts = [_run(fixture, family, selected, root, 1)[0] for root in ROOTS[family]["holdout"]]
                holdout_status, final = holdout_decision(selected, tuple(run["status"] == "PASS" for run in holdouts))
                family_record.update(holdout=holdouts, selection_status=holdout_status, selected_count=final)
        families[family] = family_record
    return {"specification_id": SPECIFICATION_ID, "policy_identity": SUMMARY_POLICY_IDENTITY, "truth_fields": TRUTH_FIELDS, "truth": truth_estimands(), "families": families, "canonical_lane": lane_records}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--image-digest", required=True)
    arguments = parser.parse_args()
    arguments.output.write_text(json.dumps(acquire(arguments.image_digest), allow_nan=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

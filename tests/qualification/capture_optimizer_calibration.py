"""Frozen one-round prospective optimizer calibration for issue #602.

This is deliberately a direct acquisition script, not a reusable framework.  The
constants below are the reviewed prospective specification.  Authoritative work
may start only from the exact commit that first contains this file and its tests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, cast
from unittest.mock import patch

from chemex.numerical_lanes import RuntimeEnvironment, canonical_lanes

ROOT = Path(__file__).parents[2]
SPECIFICATION_ID = "chemex-optimizer-calibration-v1"

# fmt: off
SPECIFICATION: dict[str, object] = {
    "truth_cases": {
        "numerical_probes": (
            "trf-routine-quadratic-v1",
            "trf-difficult-rosenbrock-v1",
            "grid-27-seed-coverage-v1",
            "grid-candidate-ordering-v1",
            "de-bounded-search-v1",
        ),
        "routine_trf": "RELAXATION_HZNZ:G2N-HN",
        "difficult_trf": "DCEST_15N_3States:K19N-HN:3st_fork",
        "coupled_grid": "RELAXATION_HZNZ:G2N-HN",
        "decomposed_grid": "RELAXATION_HZNZ:all-five-profiles",
        "de_trf": "DCEST_15N_3States:K19N-HN:3st_fork",
    },
    "difficult_starts": {
        "dce_st": {
            "dw_ab": 2.51,
            "dw_ac": -1.97,
            "kex_ab": 10.0,
            "kex_ac": 10.0,
            "pb": 0.03,
            "pc": 0.03,
        }
    },
    "candidates": {
        "routine_trf_budgets": (8, 16, 32),
        "difficult_trf": (
            ("uniform", 2048, (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
            ("physical", 256, (5.0, 5.0, 100.0, 100.0, 0.1, 0.1)),
            ("physical", 512, (5.0, 5.0, 100.0, 100.0, 0.1, 0.1)),
        ),
        "coupled_grid_budgets_per_seed": (16, 32, 64, 128),
        "decomposed_grid_budgets_per_component": (16, 32, 64, 128),
        "de": (
            ("compact-30x16", 5, 16, 510, 512),
            ("compact-48x32", 8, 32, 1584, 512),
        ),
    },
    "grid_coordinates": {
        "coupled_pb": (0.0, 0.05, 0.1),
        "decomposed_first_r1_factor": (0.75, 1.0, 1.25),
        "immutable_grid_seed_count": 27,
    },
    "de_coordinates": (
        ("dw_ab", -5.0, 5.0, "linear"),
        ("dw_ac", -5.0, 5.0, "linear"),
        ("kex_ab", 1.0, 1000.0, "log"),
        ("kex_ac", 1.0, 1000.0, "log"),
        ("pb", 0.001, 0.15, "linear"),
        ("pc", 0.001, 0.15, "linear"),
    ),
    "de_roots": {
        "calibration_count": 2,
        "replay_calibration_index": 0,
        "holdout_count": 1,
        "derivation": "u64be(sha256(ascii(chemex-602-v1|de|{phase}|{index:02d}))[0:8])",
    },
    "disqualifiers": (
        "typed_terminal_not_accepted",
        "non_finite_objective_or_vector",
        "objective_request_budget_exceeded",
        "counter_ordering_violation",
        "solver_backend_counter_missing",
        "grid_seed_or_component_incomplete",
        "candidate_ordering_mismatch",
        "transition_count_not_exactly_one",
        "fresh_root_materialization_count_not_exactly_one",
        "de_objective_not_below_90_percent_of_direct_reference",
        "replay_identity_mismatch",
        "holdout_failure",
        "resource_ceiling_exceeded",
    ),
    "selection": {
        "budget": "smallest-passing-with-all-larger-candidates-passing",
        "difficult_trf": "objective-requests,budget,policy-ordinal",
        "de": "maximum-objective-ratio,total-objective-requests,topology-ordinal",
        "tie_break": "frozen-candidate-declaration-order",
        "edge_only": "grid_saturated_unsupported",
        "runner_up_after_holdout": "forbidden",
    },
    "accounting": {
        "authoritative_unit": "chemex-objective-request",
        "nonfungible": (
            "solver_requests_received",
            "objective_requests_accepted",
            "objective_evaluations_completed",
            "backend_nfev",
            "backend_njev",
            "grid_seed_allocations",
            "grouped_component_allocations",
            "de_population_size",
            "de_generation_limit",
            "de_stage_allocation",
            "trf_polish_allocation",
            "search_to_polish_transfers",
            "root_materializations",
            "cache_hits",
            "cache_misses",
            "elapsed_seconds_diagnostic_only",
        ),
        "retries": 0,
        "borrowing": False,
        "redistribution": False,
    },
    "resource_limits": {
        "routine_trf_objective_requests": 88,
        "difficult_trf_objective_requests": 3328,
        "coupled_grid_objective_requests": 1104,
        "decomposed_grid_objective_requests": 5520,
        "de_trf_objective_requests": 10428,
        "wall_clock_scientific_budget": None,
        "workers": 1,
        "native_threads": 1,
    },
    "holdouts": {
        "routine_trf": "RELAXATION_HZNZ:Y3N-HN",
        "difficult_trf": "DCEST_15N_3States:D20N-HN:3st_fork",
        "coupled_grid": "RELAXATION_HZNZ:Y3N-HN",
        "decomposed_grid": "same-case-selected-policy-exact-replay",
        "de_trf": "one-unopened-derived-root-on-selected-K19N-HN-policy",
    },
    "unsupported": {
        "de_beyond_selected_six_coordinates": True,
        "adaptive_grid_refinement": True,
        "adaptive_de_topology_or_window_extension": True,
        "wall_clock_budgets": True,
        "hidden_retries": True,
    },
}
# fmt: on

RELAXATION_EXPERIMENT = (
    ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
)
RELAXATION_PARAMETERS = (
    ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
)
DCEST_DIRECTORY = ROOT / "examples/Experiments/DCEST_15N_3States"
DCEST_EXPERIMENTS = tuple(sorted((DCEST_DIRECTORY / "Experiments").glob("*.toml")))
DCEST_PARAMETERS = DCEST_DIRECTORY / "Parameters/parameters.toml"


def _canonical_json(value: object) -> bytes:
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


def _identity(kind: str, value: object) -> str:
    return hashlib.sha256(
        _canonical_json({"kind": kind, "schema_version": 1, "record": value})
    ).hexdigest()


def specification_identity() -> str:
    return _identity(SPECIFICATION_ID, SPECIFICATION)


def derive_de_root(phase: str, index: int) -> int:
    if phase not in {"calibration", "holdout"} or index < 0:
        raise ValueError("unknown DE root phase or ordinal")
    label = f"chemex-602-v1|de|{phase}|{index:02d}"
    return int.from_bytes(hashlib.sha256(label.encode("ascii")).digest()[:8], "big")


def select_monotone_budget(
    candidates: Sequence[int], passing: set[int]
) -> tuple[str, int | None]:
    ordered = tuple(candidates)
    if (
        not ordered
        or tuple(sorted(set(ordered))) != ordered
        or not passing <= set(ordered)
    ):
        return "NON_MONOTONE_ADEQUACY", None
    if not passing:
        return "NO_ADEQUATE_CANDIDATE", None
    selected_index = min(ordered.index(value) for value in passing)
    expected = set(ordered[selected_index:])
    if passing != expected:
        return "NON_MONOTONE_ADEQUACY", None
    if selected_index == len(ordered) - 1:
        return "GRID_SATURATED", None
    return "SELECTED", ordered[selected_index]


def candidate_ordering_key(
    record: Mapping[str, object],
) -> tuple[float, tuple[float, ...], int]:
    return (
        cast("float", record["objective"]),
        tuple(cast("Sequence[float]", record["vector"])),
        int(cast("int", record["ordinal"])),
    )


def holdout_decision(selected: str, passed: Sequence[bool]) -> tuple[str, str | None]:
    if passed and all(passed):
        return "SELECTED", selected
    return "HOLDOUT_FAILED_UNSUPPORTED", None


def validate_specification_commit(expected: str) -> str:
    head = (ROOT / ".git/HEAD").read_text(encoding="ascii").strip()
    if (
        len(expected) != 40
        or any(character not in "0123456789abcdef" for character in expected)
        or expected != head
    ):
        raise RuntimeError(
            "authoritative acquisition requires the frozen specification commit "
            "as a detached HEAD"
        )
    return head


def _file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _attest(image_digest: str) -> tuple[object, dict[str, object]]:
    lane = canonical_lanes()[0]
    authority = lane.attest_current_process(image_digest)
    records = {
        "numerical_lane": lane.to_record(),
        "lane_attestation": authority.to_record(),
        "runtime_environment": RuntimeEnvironment(lane.semantics).to_record(),
    }
    return authority, records


def _counter_record(counters: Any) -> dict[str, int]:
    return {
        "solver_requests_received": int(counters.solver_requests_received),
        "objective_requests_accepted": int(counters.objective_requests_accepted),
        "objective_evaluations_completed": int(
            counters.objective_evaluations_completed
        ),
    }


def _counter_reasons(counters: Mapping[str, int], budget: int) -> list[str]:
    received = counters["solver_requests_received"]
    accepted = counters["objective_requests_accepted"]
    completed = counters["objective_evaluations_completed"]
    reasons = []
    if not completed <= accepted <= received:
        reasons.append("counter_ordering_violation")
    if accepted > budget:
        reasons.append("objective_request_budget_exceeded")
    return reasons


def _backend_record(backend: Any | None) -> dict[str, object] | None:
    if backend is None:
        return None
    return {"nfev": int(backend.nfev), "njev": backend.njev}


def _direct_record(outcome: Any, budget: int, ordinal: int) -> dict[str, object]:
    execution = outcome.execution
    accepted = outcome.accepted_result
    counters = _counter_record(execution.counters)
    backend = _backend_record(execution.backend)
    reasons = _counter_reasons(counters, budget)
    if execution.terminal.value != "converged" or accepted is None:
        reasons.append("typed_terminal_not_accepted")
    vector = () if accepted is None else tuple(accepted.vector)
    objective = math.inf if accepted is None else float(accepted.chi_square)
    if not math.isfinite(objective) or any(
        not math.isfinite(value) for value in vector
    ):
        reasons.append("non_finite_objective_or_vector")
    if backend is None:
        reasons.append("solver_backend_counter_missing")
    return {
        "ordinal": ordinal,
        "budget": budget,
        "status": "PASS" if not reasons else "DISQUALIFIED",
        "reasons": tuple(dict.fromkeys(reasons)),
        "terminal": execution.terminal.value,
        "objective": objective,
        "vector": vector,
        "counters": counters,
        "backend": backend,
        "execution_identity": execution.identity,
        "accepted_identity": None if accepted is None else accepted.identity,
    }


def _build_relaxation(
    spin_system: str | None, *, grouped: bool = False
) -> tuple[Any, Any, Any, Any, Any]:
    from chemex.configuration.methods import Method, Selection, read_methods
    from chemex.configuration.parameters import read_defaults
    from chemex.evaluation.native import EvaluationEngine
    from chemex.experiments.builder import build_experiments
    from chemex.optimize.direct_trf import OptimizationProblem
    from chemex.parameters.spin_system import SpinSystem
    from chemex.runtime import AnalysisSession

    session = AnalysisSession.create()
    session.set_model("2st")
    include = None if spin_system is None else [SpinSystem.from_name(spin_system)]
    experiments = build_experiments(
        [RELAXATION_EXPERIMENT],
        Selection(include=include, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([RELAXATION_PARAMETERS]))
    if not session.try_build_analysis_values():
        raise RuntimeError("RELAXATION_HZNZ native parameter construction failed")
    method = (
        read_methods(
            [ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"]
        )["DEFAULT"]
        if grouped
        else Method(fit=["PB"], fix=["KEX_AB"])
    )
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    if configuration is None:
        raise RuntimeError("RELAXATION_HZNZ configuration did not seal")
    problem = OptimizationProblem.from_native(
        engine.plan, parameterization, configuration, session.analysis_values.snapshot()
    )
    return session, experiments, parameterization, engine, problem


def _build_dcest(spin_system: str) -> tuple[Any, Any, Any, Any, Any]:
    from chemex.configuration.methods import Method, Selection
    from chemex.configuration.parameters import read_defaults
    from chemex.evaluation.native import EvaluationEngine
    from chemex.experiments.builder import build_experiments
    from chemex.optimize.direct_trf import OptimizationProblem
    from chemex.parameters.spin_system import SpinSystem
    from chemex.runtime import AnalysisSession

    session = AnalysisSession.create()
    session.set_model("3st_fork")
    experiments = build_experiments(
        list(DCEST_EXPERIMENTS),
        Selection(include=[SpinSystem.from_name(spin_system)], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([DCEST_PARAMETERS]))
    if not session.try_build_analysis_values():
        raise RuntimeError("DCEST native parameter construction failed")
    method = Method(fix=["R1_A", "R2_A", "R2_B", "R2_C"])
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    if configuration is None:
        raise RuntimeError("DCEST configuration did not seal")
    problem = OptimizationProblem.from_native(
        engine.plan, parameterization, configuration, session.analysis_values.snapshot()
    )
    expected = ("DW_AB", "DW_AC", "KEX_AB", "KEX_AC", "PB", "PC")
    observed = tuple(
        _semantic_parameter_name(param_id) for param_id in problem.controlled_ids
    )
    if observed != expected:
        raise RuntimeError(f"frozen DCEST controlled coordinates changed: {observed!r}")
    return session, experiments, parameterization, engine, problem


def _semantic_parameter_name(param_id: str) -> str:
    value = param_id.removeprefix("__")
    for semantic in ("DW_AB", "DW_AC", "KEX_AB", "KEX_AC", "PB", "PC"):
        if value == semantic or value.startswith(f"{semantic}_"):
            return semantic
    return value


def _run_direct(
    problem: Any,
    parameterization: Any,
    engine: Any,
    budget: int,
    scale: Sequence[float],
) -> object:
    from chemex.optimize.direct_trf import DirectTrfInvocation, execute_direct_trf

    invocation = DirectTrfInvocation.for_problem(
        problem, objective_request_budget=budget, x_scale=tuple(scale)
    )
    return execute_direct_trf(problem, invocation, parameterization, engine)


def _direct_policy_signature(record: Mapping[str, object]) -> dict[str, object]:
    return {
        "terminal": record["terminal"],
        "objective": record["objective"],
        "vector": record["vector"],
        "counters": record["counters"],
        "backend": record["backend"],
        "execution_identity": record["execution_identity"],
        "accepted_identity": record["accepted_identity"],
    }


def _calibrate_routine_direct() -> dict[str, object]:
    _session, _experiments, parameterization, engine, problem = _build_relaxation(
        "G2N-HN"
    )
    budgets = cast(
        "tuple[int, ...]",
        cast("Mapping[str, object]", SPECIFICATION["candidates"])[
            "routine_trf_budgets"
        ],
    )
    records = []
    passing: set[int] = set()
    for ordinal, budget in enumerate(budgets):
        record = _direct_record(
            _run_direct(
                problem,
                parameterization,
                engine,
                budget,
                (1.0,) * len(problem.controlled_ids),
            ),
            budget,
            ordinal,
        )
        records.append(record)
        if record["status"] == "PASS":
            passing.add(budget)
    status, selected = select_monotone_budget(budgets, passing)
    result: dict[str, object] = {
        "candidates": records,
        "selection": {"status": status, "budget": selected},
    }
    if selected is None:
        result.update(replay=None, holdout=None)
        return result
    selected_record = next(item for item in records if item["budget"] == selected)
    replay = _direct_record(
        _run_direct(
            problem,
            parameterization,
            engine,
            selected,
            (1.0,) * len(problem.controlled_ids),
        ),
        selected,
        -1,
    )
    replay_passed = _direct_policy_signature(replay) == _direct_policy_signature(
        selected_record
    )
    result["replay"] = {
        "status": "PASS" if replay_passed else "DISQUALIFIED",
        "record": replay,
    }
    if not replay_passed:
        result["selection"] = {"status": "REPLAY_MISMATCH_UNSUPPORTED", "budget": None}
        result["holdout"] = None
        return result
    _hs, _he, hp, heng, hproblem = _build_relaxation("Y3N-HN")
    holdout = _direct_record(
        _run_direct(
            hproblem, hp, heng, selected, (1.0,) * len(hproblem.controlled_ids)
        ),
        selected,
        0,
    )
    holdout_status, final = holdout_decision(
        str(selected), (holdout["status"] == "PASS",)
    )
    result["holdout"] = holdout
    result["selection"] = {
        "status": holdout_status,
        "budget": None if final is None else int(final),
    }
    return result


def _calibrate_difficult_direct() -> tuple[dict[str, object], float | None]:
    _session, _experiments, parameterization, engine, problem = _build_dcest("K19N-HN")
    candidates = cast(
        "tuple[tuple[str, int, tuple[float, ...]], ...]",
        cast("Mapping[str, object]", SPECIFICATION["candidates"])["difficult_trf"],
    )
    records = []
    for ordinal, (name, budget, scale) in enumerate(candidates):
        record = _direct_record(
            _run_direct(problem, parameterization, engine, budget, scale),
            budget,
            ordinal,
        )
        record["policy"] = name
        records.append(record)
    passing = [item for item in records if item["status"] == "PASS"]
    result: dict[str, object] = {"candidates": records}
    if not passing:
        result.update(
            selection={
                "status": "NO_ADEQUATE_CANDIDATE",
                "policy": None,
                "budget": None,
            },
            replay=None,
            holdout=None,
        )
        return result, None
    selected_record = min(
        passing,
        key=lambda item: (
            cast("Mapping[str, int]", item["counters"])["objective_requests_accepted"],
            int(cast("int", item["budget"])),
            int(cast("int", item["ordinal"])),
        ),
    )
    name, budget = (
        str(selected_record["policy"]),
        int(cast("int", selected_record["budget"])),
    )
    scale = candidates[int(cast("int", selected_record["ordinal"]))][2]
    replay = _direct_record(
        _run_direct(problem, parameterization, engine, budget, scale), budget, -1
    )
    replay_passed = _direct_policy_signature(replay) == _direct_policy_signature(
        selected_record
    )
    result["replay"] = {
        "status": "PASS" if replay_passed else "DISQUALIFIED",
        "record": replay,
    }
    if not replay_passed:
        result.update(
            selection={
                "status": "REPLAY_MISMATCH_UNSUPPORTED",
                "policy": None,
                "budget": None,
            },
            holdout=None,
        )
        return result, None
    _hs, _he, hp, heng, hproblem = _build_dcest("D20N-HN")
    holdout = _direct_record(_run_direct(hproblem, hp, heng, budget, scale), budget, 0)
    holdout_status, final = holdout_decision(
        f"{name}-{budget}", (holdout["status"] == "PASS",)
    )
    result["holdout"] = holdout
    result["selection"] = {
        "status": holdout_status,
        "policy": name if final else None,
        "budget": budget if final else None,
        "x_scale": scale if final else None,
    }
    return result, cast("float", selected_record["objective"]) if final else None


def _candidate_execution_record(outcome: Any, budget: int) -> dict[str, object]:
    execution = outcome.execution
    counters = _counter_record(execution.counters)
    backend = _backend_record(execution.backend)
    reasons = _counter_reasons(counters, budget)
    if backend is None:
        reasons.append("solver_backend_counter_missing")
    return {
        "terminal": outcome.terminal.value,
        "execution_terminal": execution.terminal.value,
        "execution_identity": execution.identity,
        "counters": counters,
        "backend": backend,
        "status": "PASS" if not reasons else "DISQUALIFIED",
        "reasons": tuple(reasons),
    }


def _grid_signature(record: Mapping[str, object]) -> dict[str, object]:
    return {
        "terminal": record["terminal"],
        "attempts": record["attempts"],
        "selection_identity": record["selection_identity"],
        "accepted_identity": record["accepted_identity"],
        "root_materialization": record["root_materialization"],
        "local_executions": record["local_executions"],
    }


def _run_coupled_grid(
    problem: Any,
    parameterization: Any,
    engine: Any,
    budget: int,
) -> dict[str, object]:
    import chemex.optimize.grid_direct_trf as grid_module
    from chemex.optimize.grid_direct_trf import (
        GridDirectTrfInvocation,
        execute_grid_direct_trf,
    )

    pb_id = next(
        param_id
        for param_id in problem.controlled_ids
        if _semantic_parameter_name(param_id) == "PB"
    )
    axes = cast(
        "tuple[float, ...]",
        cast("Mapping[str, object]", SPECIFICATION["grid_coordinates"])["coupled_pb"],
    )
    invocation = GridDirectTrfInvocation.for_problem(
        problem, axes=((pb_id, axes),), objective_request_budget=budget
    )
    local_outcomes: list[object] = []
    original = grid_module.execute_direct_trf_candidate

    def tracked(*args: Any, **kwargs: Any) -> object:
        outcome = original(*args, **kwargs)
        local_outcomes.append(outcome)
        return outcome

    started = time.perf_counter()
    with patch.object(grid_module, "execute_direct_trf_candidate", tracked):
        outcome = execute_grid_direct_trf(problem, invocation, parameterization, engine)
    elapsed = time.perf_counter() - started
    selected_ordinal = (
        None if outcome.selection is None else outcome.selection.selected_seed_ordinal
    )
    attempts = tuple(
        {
            "seed_identity": item.seed_identity,
            "seed_ordinal": item.seed_ordinal,
            "coordinates": tuple(item.axis_items),
            "disposition": item.disposition.value,
            "objective": item.objective,
            "candidate_identity": (
                None if item.candidate is None else item.candidate.identity
            ),
            "selected": item.seed_ordinal == selected_ordinal,
        }
        for item in outcome.attempts
    )
    locals_record = tuple(
        _candidate_execution_record(item, budget) for item in local_outcomes
    )
    reasons = []
    if outcome.terminal.value != "accepted":
        reasons.append("typed_terminal_not_accepted")
    if len(attempts) != len(invocation.seeds) or any(
        item["disposition"] != "eligible" for item in attempts
    ):
        reasons.append("grid_seed_or_component_incomplete")
    if len(locals_record) != len(invocation.seeds) or any(
        item["status"] != "PASS" for item in locals_record
    ):
        reasons.append("grid_seed_or_component_incomplete")
    selected = outcome.selection
    candidates = tuple(
        {
            "objective": float(item.objective),
            "vector": tuple(item.candidate.candidate.vector),
            "ordinal": item.seed_ordinal,
            "identity": item.candidate.identity,
        }
        for item in outcome.attempts
        if item.candidate is not None and item.objective is not None
    )
    expected_order = tuple(
        item["identity"] for item in sorted(candidates, key=candidate_ordering_key)
    )
    if selected is None or selected.eligible_candidate_identities != expected_order:
        reasons.append("candidate_ordering_mismatch")
    materialization = outcome.materialization
    root_materialization = (
        None
        if materialization is None
        else {
            "identity": materialization.identity,
            "evaluation_count": materialization.evaluation_count,
            "cache_hits": materialization.cache_hits,
            "cache_misses": materialization.cache_misses,
        }
    )
    if materialization is None or materialization.evaluation_count != 1:
        reasons.append("fresh_root_materialization_count_not_exactly_one")
    return {
        "budget_per_seed": budget,
        "seed_count": len(invocation.seeds),
        "allocation": budget * len(invocation.seeds),
        "terminal": outcome.terminal.value,
        "status": "PASS" if not reasons else "DISQUALIFIED",
        "reasons": tuple(dict.fromkeys(reasons)),
        "attempts": attempts,
        "selection_identity": None if selected is None else selected.identity,
        "accepted_identity": None
        if outcome.accepted_result is None
        else outcome.accepted_result.identity,
        "root_materialization": root_materialization,
        "local_executions": locals_record,
        "elapsed_seconds_diagnostic_only": elapsed,
    }


def _calibrate_coupled_grid() -> dict[str, object]:
    _s, _e, parameterization, engine, problem = _build_relaxation("G2N-HN")
    budgets = cast(
        "tuple[int, ...]",
        cast("Mapping[str, object]", SPECIFICATION["candidates"])[
            "coupled_grid_budgets_per_seed"
        ],
    )
    records = [
        _run_coupled_grid(problem, parameterization, engine, budget)
        for budget in budgets
    ]
    passing = {
        cast("int", item["budget_per_seed"])
        for item in records
        if item["status"] == "PASS"
    }
    status, selected = select_monotone_budget(budgets, passing)
    result: dict[str, object] = {
        "candidates": records,
        "selection": {"status": status, "budget_per_seed": selected},
    }
    if selected is None:
        result.update(replay=None, holdout=None)
        return result
    chosen = next(item for item in records if item["budget_per_seed"] == selected)
    replay = _run_coupled_grid(problem, parameterization, engine, selected)
    replay_passed = _grid_signature(replay) == _grid_signature(chosen)
    result["replay"] = {
        "status": "PASS" if replay_passed else "DISQUALIFIED",
        "record": replay,
    }
    if not replay_passed:
        result["selection"] = {
            "status": "REPLAY_MISMATCH_UNSUPPORTED",
            "budget_per_seed": None,
        }
        result["holdout"] = None
        return result
    _hs, _he, hp, heng, hproblem = _build_relaxation("Y3N-HN")
    holdout = _run_coupled_grid(hproblem, hp, heng, selected)
    holdout_status, final = holdout_decision(
        str(selected), (holdout["status"] == "PASS",)
    )
    result["holdout"] = holdout
    result["selection"] = {
        "status": holdout_status,
        "budget_per_seed": None if final is None else int(final),
    }
    return result


def _run_grouped_grid(
    problem: Any,
    parameterization: Any,
    engine: Any,
    decomposition: Any,
    budget: int,
) -> dict[str, object]:
    import chemex.optimize.grouped_direct_trf as grouped_module
    from chemex.optimize.grid_direct_trf import GridDirectTrfInvocation
    from chemex.optimize.grouped_grid_direct_trf import execute_grouped_grid_direct_trf

    axis_id = problem.controlled_ids[0]
    start = problem.start[0]
    factors = cast(
        "tuple[float, ...]",
        cast("Mapping[str, object]", SPECIFICATION["grid_coordinates"])[
            "decomposed_first_r1_factor"
        ],
    )
    invocation = GridDirectTrfInvocation.for_problem(
        problem,
        axes=((axis_id, tuple(start * factor for factor in factors)),),
        objective_request_budget=budget,
    )
    local_outcomes: list[object] = []
    original = grouped_module.execute_direct_trf_candidate

    def tracked(*args: Any, **kwargs: Any) -> object:
        outcome = original(*args, **kwargs)
        local_outcomes.append(outcome)
        return outcome

    started = time.perf_counter()
    with patch.object(grouped_module, "execute_direct_trf_candidate", tracked):
        outcome = execute_grouped_grid_direct_trf(
            problem, decomposition, invocation, parameterization, engine
        )
    elapsed = time.perf_counter() - started
    attempts = tuple(
        {
            "seed_identity": item.seed_identity,
            "seed_ordinal": item.seed_ordinal,
            "disposition": item.disposition.value,
            "component_dispositions": tuple(
                component.disposition.value for component in item.components
            ),
            "component_identities": tuple(
                component.identity for component in item.components
            ),
            "aggregate_identity": None
            if item.candidate is None
            else item.candidate.identity,
            "objective": item.objective,
        }
        for item in outcome.attempts
    )
    locals_record = tuple(
        _candidate_execution_record(item, budget) for item in local_outcomes
    )
    expected_local_count = len(invocation.seeds) * len(decomposition.components)
    reasons = []
    if outcome.terminal.value != "accepted":
        reasons.append("typed_terminal_not_accepted")
    if len(attempts) != len(invocation.seeds) or any(
        item["disposition"] != "eligible"
        or any(
            value != "succeeded"
            for value in cast("Sequence[str]", item["component_dispositions"])
        )
        for item in attempts
    ):
        reasons.append("grid_seed_or_component_incomplete")
    if len(locals_record) != expected_local_count or any(
        item["status"] != "PASS" for item in locals_record
    ):
        reasons.append("grid_seed_or_component_incomplete")
    aggregates = tuple(
        {
            "objective": float(item.objective),
            "vector": tuple(item.candidate.vector),
            "ordinal": item.seed_ordinal,
            "identity": item.candidate.identity,
        }
        for item in outcome.attempts
        if item.candidate is not None and item.objective is not None
    )
    selection = outcome.selection
    expected_order = tuple(
        item["identity"] for item in sorted(aggregates, key=candidate_ordering_key)
    )
    if selection is None or selection.eligible_candidate_identities != expected_order:
        reasons.append("candidate_ordering_mismatch")
    accepted = outcome.accepted_result
    materialization = (
        None if accepted is None else accepted.fresh_candidate.materialization
    )
    root_materialization = (
        None
        if materialization is None
        else {
            "identity": materialization.identity,
            "evaluation_count": materialization.evaluation_count,
            "cache_hits": materialization.cache_hits,
            "cache_misses": materialization.cache_misses,
        }
    )
    if materialization is None or materialization.evaluation_count != 1:
        reasons.append("fresh_root_materialization_count_not_exactly_one")
    return {
        "budget_per_component": budget,
        "seed_count": len(invocation.seeds),
        "component_count": len(decomposition.components),
        "allocation": budget * expected_local_count,
        "terminal": outcome.terminal.value,
        "status": "PASS" if not reasons else "DISQUALIFIED",
        "reasons": tuple(dict.fromkeys(reasons)),
        "decomposition_identity": decomposition.identity,
        "attempts": attempts,
        "selection_identity": None if selection is None else selection.identity,
        "accepted_identity": None if accepted is None else accepted.identity,
        "root_materialization": root_materialization,
        "local_executions": locals_record,
        "elapsed_seconds_diagnostic_only": elapsed,
    }


def _calibrate_grouped_grid() -> dict[str, object]:
    from chemex.optimize.grouped_direct_trf import FitDecomposition

    _s, _e, parameterization, engine, problem = _build_relaxation(None, grouped=True)
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    budgets = cast(
        "tuple[int, ...]",
        cast("Mapping[str, object]", SPECIFICATION["candidates"])[
            "decomposed_grid_budgets_per_component"
        ],
    )
    records = [
        _run_grouped_grid(problem, parameterization, engine, decomposition, budget)
        for budget in budgets
    ]
    passing = {
        cast("int", item["budget_per_component"])
        for item in records
        if item["status"] == "PASS"
    }
    status, selected = select_monotone_budget(budgets, passing)
    result: dict[str, object] = {
        "candidates": records,
        "selection": {"status": status, "budget_per_component": selected},
        "holdout": {"kind": "selected-policy-exact-replay"},
    }
    if selected is None:
        result.update(replay=None)
        result["holdout"] = None
        return result
    chosen = next(item for item in records if item["budget_per_component"] == selected)
    replay = _run_grouped_grid(
        problem, parameterization, engine, decomposition, selected
    )
    replay_passed = _grid_signature(replay) == _grid_signature(chosen)
    result["replay"] = {
        "status": "PASS" if replay_passed else "DISQUALIFIED",
        "record": replay,
    }
    if not replay_passed:
        result["selection"] = {
            "status": "REPLAY_MISMATCH_UNSUPPORTED",
            "budget_per_component": None,
        }
        result["holdout"] = {"status": "HOLDOUT_FAILED_UNSUPPORTED"}
    else:
        result["holdout"] = {"status": "PASS", "exact_replay": True}
    return result


def _de_signature(record: Mapping[str, object]) -> dict[str, object]:
    return {
        "terminal": record["terminal"],
        "search_identity": record["search_identity"],
        "accounting_identity": record["accounting_identity"],
        "accepted_identity": record["accepted_identity"],
        "objective": record["objective"],
        "vector": record["vector"],
        "de_counters": record["de_counters"],
        "polish_counters": record["polish_counters"],
        "backend": record["backend"],
    }


def _run_de(
    problem: Any,
    parameterization: Any,
    engine: Any,
    topology: tuple[str, int, int, int, int],
    root_seed: int,
    direct_reference: float,
    ordinal: int,
) -> dict[str, object]:
    from chemex.optimize.de_direct_trf import (
        DeDirectTrfInvocation,
        execute_de_direct_trf,
    )

    label, multiplier, generations, de_budget, polish_budget = topology
    ids = {
        _semantic_parameter_name(param_id): param_id
        for param_id in problem.controlled_ids
    }
    coordinate_specs = cast(
        "tuple[tuple[str, float, float, str], ...]", SPECIFICATION["de_coordinates"]
    )
    search_coordinates = tuple(
        (ids[name.upper()], lower, upper, semantics)
        for name, lower, upper, semantics in coordinate_specs
    )
    scales = tuple(
        5.0
        if _semantic_parameter_name(param_id).startswith("DW_")
        else 100.0
        if _semantic_parameter_name(param_id).startswith("KEX_")
        else 0.1
        for param_id in problem.controlled_ids
    )
    invocation = DeDirectTrfInvocation.for_problem(
        problem,
        search_coordinates=search_coordinates,
        root_seed=root_seed,
        de_objective_request_budget=de_budget,
        polish_objective_request_budget=polish_budget,
        population_multiplier=multiplier,
        maximum_generations=generations,
        tol=0.0,
        atol=0.0,
        polish_x_scale=scales,
    )
    started = time.perf_counter()
    outcome = execute_de_direct_trf(problem, invocation, parameterization, engine)
    elapsed = time.perf_counter() - started
    accounting = outcome.accounting
    de_counters = _counter_record(accounting.de_counters)
    polish_counters = (
        None
        if accounting.polish_counters is None
        else _counter_record(accounting.polish_counters)
    )
    accepted = outcome.accepted_result
    objective = math.inf if accepted is None else float(accepted.chi_square)
    vector = () if accepted is None else tuple(accepted.vector)
    reasons = _counter_reasons(de_counters, de_budget)
    if polish_counters is None:
        reasons.append("typed_terminal_not_accepted")
    else:
        reasons.extend(_counter_reasons(polish_counters, polish_budget))
    if outcome.terminal.value != "accepted" or accepted is None:
        reasons.append("typed_terminal_not_accepted")
    if not math.isfinite(objective) or any(
        not math.isfinite(value) for value in vector
    ):
        reasons.append("non_finite_objective_or_vector")
    if accounting.search_to_polish_transfers != 1:
        reasons.append("transition_count_not_exactly_one")
    if accounting.root_materializations != 1:
        reasons.append("fresh_root_materialization_count_not_exactly_one")
    if not objective <= 0.9 * direct_reference:
        reasons.append("de_objective_not_below_90_percent_of_direct_reference")
    backend = (
        None
        if outcome.search.backend is None
        else {
            "generation_count": outcome.search.backend.generation_count,
            "backend_evaluation_count": outcome.search.backend.backend_evaluation_count,
            "population_size": outcome.search.backend.population_size,
        }
    )
    if backend is None:
        reasons.append("solver_backend_counter_missing")
    materialization = outcome.root_materialization
    return {
        "ordinal": ordinal,
        "topology": label,
        "root_seed": root_seed,
        "population": {
            "multiplier": multiplier,
            "size": invocation.population.size,
            "maximum_generations": generations,
            "identity": invocation.population.identity,
        },
        "budgets": {"de": de_budget, "polish": polish_budget},
        "allocation": de_budget + polish_budget,
        "terminal": outcome.terminal.value,
        "search_terminal": outcome.search.terminal.value,
        "status": "PASS" if not reasons else "DISQUALIFIED",
        "reasons": tuple(dict.fromkeys(reasons)),
        "objective": objective,
        "objective_ratio": objective / direct_reference,
        "vector": vector,
        "de_counters": de_counters,
        "polish_counters": polish_counters,
        "backend": backend,
        "candidate_ordering_policy": outcome.search.candidate_ordering_policy,
        "search_identity": outcome.search.identity,
        "accounting_identity": accounting.identity,
        "accepted_identity": None if accepted is None else accepted.identity,
        "materialization_identity": None
        if materialization is None
        else materialization.identity,
        "cache": None
        if materialization is None
        else {
            "hits": materialization.cache_hits,
            "misses": materialization.cache_misses,
        },
        "elapsed_seconds_diagnostic_only": elapsed,
    }


def _calibrate_de(direct_reference: float | None) -> dict[str, object]:
    if direct_reference is None:
        return {
            "candidates": (),
            "selection": {
                "status": "DIRECT_REFERENCE_UNQUALIFIED",
                "topology": None,
            },
            "replay": None,
            "holdout": None,
        }
    _s, _e, parameterization, engine, problem = _build_dcest("K19N-HN")
    topologies = cast(
        "tuple[tuple[str, int, int, int, int], ...]",
        cast("Mapping[str, object]", SPECIFICATION["candidates"])["de"],
    )
    roots = tuple(derive_de_root("calibration", index) for index in range(2))
    candidates: list[dict[str, Any]] = []
    adequate_ordinals: set[int] = set()
    for ordinal, topology in enumerate(topologies):
        runs = tuple(
            _run_de(
                problem,
                parameterization,
                engine,
                topology,
                root_seed,
                direct_reference,
                ordinal,
            )
            for root_seed in roots
        )
        candidates.append({"ordinal": ordinal, "topology": topology[0], "runs": runs})
        if all(run["status"] == "PASS" for run in runs):
            adequate_ordinals.add(ordinal)
    if not adequate_ordinals:
        status, selected_ordinal = "NO_ADEQUATE_CANDIDATE", None
    elif adequate_ordinals == {len(topologies) - 1}:
        status, selected_ordinal = "GRID_SATURATED", None
    elif adequate_ordinals != set(range(min(adequate_ordinals), len(topologies))):
        status, selected_ordinal = "NON_MONOTONE_ADEQUACY", None
    else:
        selected_ordinal = min(
            adequate_ordinals,
            key=lambda value: (
                max(
                    cast("float", run["objective_ratio"])
                    for run in cast(
                        "Sequence[Mapping[str, object]]",
                        candidates[value]["runs"],
                    )
                ),
                sum(
                    cast("Mapping[str, int]", run["de_counters"])[
                        "objective_requests_accepted"
                    ]
                    + cast("Mapping[str, int]", run["polish_counters"])[
                        "objective_requests_accepted"
                    ]
                    for run in cast(
                        "Sequence[Mapping[str, object]]",
                        candidates[value]["runs"],
                    )
                ),
                value,
            ),
        )
        status = "SELECTED"
    result: dict[str, object] = {
        "direct_reference_objective": direct_reference,
        "candidates": tuple(candidates),
        "selection": {
            "status": status,
            "topology": None
            if selected_ordinal is None
            else topologies[selected_ordinal][0],
            "ordinal": selected_ordinal,
        },
    }
    if selected_ordinal is None:
        result.update(replay=None, holdout=None)
        return result
    topology = topologies[selected_ordinal]
    original = cast(
        "Sequence[Mapping[str, object]]", candidates[selected_ordinal]["runs"]
    )[0]
    replay = _run_de(
        problem,
        parameterization,
        engine,
        topology,
        roots[0],
        direct_reference,
        selected_ordinal,
    )
    replay_passed = _de_signature(replay) == _de_signature(original)
    result["replay"] = {
        "status": "PASS" if replay_passed else "DISQUALIFIED",
        "record": replay,
    }
    if not replay_passed:
        result["selection"] = {
            "status": "REPLAY_MISMATCH_UNSUPPORTED",
            "topology": None,
            "ordinal": None,
        }
        result["holdout"] = None
        return result
    holdout = _run_de(
        problem,
        parameterization,
        engine,
        topology,
        derive_de_root("holdout", 0),
        direct_reference,
        selected_ordinal,
    )
    holdout_status, final = holdout_decision(
        topology[0], (holdout["status"] == "PASS",)
    )
    result["holdout"] = holdout
    result["selection"] = {
        "status": holdout_status,
        "topology": final,
        "ordinal": selected_ordinal if final else None,
        "population_multiplier": topology[1] if final else None,
        "maximum_generations": topology[2] if final else None,
        "de_budget": topology[3] if final else None,
        "polish_budget": topology[4] if final else None,
    }
    return result


def _qualification_probes(authority: object) -> dict[str, object]:
    from chemex.optimize.numerical_probes import (
        CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
        run_numerical_probe_baseline,
    )

    baseline = run_numerical_probe_baseline(
        cast("Any", authority),
        expected_manifest_identity=CANONICAL_NUMERICAL_PROBE_MANIFEST_IDENTITY,
    )
    required = cast(
        "tuple[str, ...]",
        cast("Mapping[str, object]", SPECIFICATION["truth_cases"])["numerical_probes"],
    )
    selected = tuple(
        artifact.to_record()
        for artifact in baseline.artifacts
        if artifact.probe_id in required
    )
    if tuple(item["probe_id"] for item in selected) != required:
        raise RuntimeError("frozen numerical probe set changed")
    return {
        "baseline_identity": baseline.identity,
        "manifest_identity": baseline.manifest_identity,
        "qualification": baseline.historical_qualification,
        "artifacts": selected,
    }


def _guard(name: str, function: Any) -> dict[str, object]:
    try:
        value = function()
    except Exception as error:  # noqa: BLE001 - one round preserves typed failure
        return {
            "status": "ARCHITECTURAL_FAILURE",
            "stage": name,
            "error_type": type(error).__name__,
            "message": str(error),
        }
    if not isinstance(value, dict):
        raise TypeError(f"{name} did not return a record")
    return value


def acquire(image_digest: str, specification_commit: str) -> dict[str, object]:
    validate_specification_commit(specification_commit)
    authority, lane_records = _attest(image_digest)
    started = time.perf_counter()
    probes = _guard("truth_probes", lambda: _qualification_probes(authority))
    routine = _guard("routine_direct_trf", _calibrate_routine_direct)
    try:
        difficult, direct_reference = _calibrate_difficult_direct()
    except Exception as error:  # noqa: BLE001 - preserve the sole round
        difficult = {
            "status": "ARCHITECTURAL_FAILURE",
            "stage": "difficult_direct_trf",
            "error_type": type(error).__name__,
            "message": str(error),
        }
        direct_reference = None
    coupled = _guard("coupled_grid", _calibrate_coupled_grid)
    grouped = _guard("decomposed_grid", _calibrate_grouped_grid)
    de = _guard("de_trf", lambda: _calibrate_de(direct_reference))
    record: dict[str, object] = {
        "schema_version": 1,
        "record_version": "canonical-optimizer-calibration-v1",
        "specification": {
            "id": SPECIFICATION_ID,
            "identity": specification_identity(),
            "record": SPECIFICATION,
        },
        "source": {
            "specification_commit": specification_commit,
            "qualification_script_sha256": _file_hash(Path(__file__)),
            "qualification_test_sha256": _file_hash(
                ROOT / "tests/test_optimizer_calibration.py"
            ),
            "dependency_lock_sha256": _file_hash(ROOT / "uv.lock"),
            "input_sha256": {
                str(path.relative_to(ROOT)): _file_hash(path)
                for path in (
                    RELAXATION_EXPERIMENT,
                    RELAXATION_PARAMETERS,
                    *DCEST_EXPERIMENTS,
                    DCEST_PARAMETERS,
                )
            },
        },
        "canonical_lane": lane_records,
        "authoritative_acquisition": {
            "round": 1,
            "retry_count": 0,
            "retuned": False,
            "adaptive_extension": False,
        },
        "truth_probes": probes,
        "strata": {
            "routine_direct_trf": routine,
            "difficult_direct_trf": difficult,
            "coupled_grid_trf": coupled,
            "decomposed_grouped_grid_trf": grouped,
            "de_trf": de,
        },
        "operational": {
            "elapsed_seconds_diagnostic_only": time.perf_counter() - started,
            "workers": 1,
            "native_threads": 1,
        },
    }
    record["identity"] = _identity("canonical-optimizer-calibration", record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--image-digest", required=True)
    parser.add_argument("--specification-commit", required=True)
    arguments = parser.parse_args()
    record = acquire(arguments.image_digest, arguments.specification_commit)
    arguments.output.write_text(
        json.dumps(record, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

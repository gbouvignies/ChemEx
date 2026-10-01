"""Issue #712 experiment; no live commits or public search policy.

Run with the thread environment in benchmarks/README.md. JSON contains every
outer point and converged endpoint. Singleton GRID supplies all profiling,
factor proof, reconstruction and validation; the temporary patch adds only
factor-local alternate starts. GRID's internal acceptance is never committed
or used as a final fit. Only complete grouped TRF endpoints are ranked.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import platform
import sys
from contextlib import redirect_stdout
from dataclasses import dataclass, replace
from itertools import product
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

import numpy as np
import scipy
from scipy.stats import qmc

import chemex.optimize.grouped_direct_trf as grouped
import chemex.optimize.profiled_grid as grid
from chemex import __file__ as chemex_source
from chemex.configuration.method_plan import ProfileSelection
from chemex.configuration.methods import Method, Selection, read_method_plan
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import (
    BoundEvaluator,
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
)
from chemex.experiments.builder import build_experiments
from chemex.optimize.de_direct_trf import DeSearchInvocation, execute_de_search
from chemex.optimize.direct_trf import OptimizationProblem, canonical_chi_square
from chemex.optimize.method_compiler import (
    DeSearchInstruction,
    FitStep,
    GridSearchInstruction,
    compile_method_plan,
)
from chemex.optimize.native_deterministic import _build_invocation
from chemex.parameters.parameterization import ActiveParameterization, ParameterRole
from chemex.parameters.sealed import SealedConfiguration
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples/Experiments"


@dataclass
class Context:
    problem: OptimizationProblem
    parameterization: ActiveParameterization
    engine: EvaluationEngine
    configuration: SealedConfiguration
    coordinates: tuple[tuple[str, float, float, str], ...]

    def at(self, vector: tuple[float, ...]) -> OptimizationProblem:
        """Fresh research snapshot, without changing the session's live values.

        GRID derives nuisance starts from independent_items, so restart_from
        alone (which preserves those items) is insufficient for this experiment.
        Build through the native constructor rather than modifying child records.
        """
        updates = dict(zip(self.problem.controlled_ids, vector, strict=True))
        snapshot = replace(
            self.problem.source_snapshot,
            _items=tuple(
                (key, updates.get(key, value))
                for key, value in self.problem.source_snapshot.items()
            ),
        )
        return OptimizationProblem.from_native(
            self.engine.plan, self.parameterization, self.configuration, snapshot
        )


def build_context(case: str) -> Context:
    session = AnalysisSession.create()
    session.set_model("3st" if case == "dcest" else "2st")
    example = (
        EXAMPLES
        / {
            "cpmg": "CPMG_15N_IP",
            "cest": "CEST_13C_LABEL_CN",
            "dcest": "DCEST_15N_3States",
        }[case]
    )
    with redirect_stdout(io.StringIO()):
        experiments = build_experiments(
            sorted((example / "Experiments").glob("*.toml")),
            Selection(include=None, exclude=None),
            session=session,
        )
    session.parameters.set_defaults(
        read_defaults([example / "Parameters/parameters.toml"])
    )
    if not session.try_build_analysis_values():
        raise RuntimeError(session.parameter_factory.native_construction_error)
    configuration = session.parameter_factory.sealed_configuration
    model = session.parameter_factory.sealed_parameter_model
    if configuration is None or model is None:
        raise RuntimeError("Native configuration was not sealed")
    if case == "cest":
        # Small real-data scope, retaining the historical ambiguous residue.
        method = Method(
            include=["L18CD1", "L18CB", "K25CA"],
            constraints=["[R2_B] = [R2_A]"],
            fit=["PB", "KEX_AB", "CS_A"],
        )
        experiments.select_profiles(
            ProfileSelection(include=("L18CD1", "L18CB", "K25CA"))
        )
        parameterization = session.compile_parameterization(
            method, experiments.param_ids
        )
        ranges = (("KEX_AB", 100.0, 600.0), ("PB", 0.03, 0.15))
        # Same coupling search box as shipped CPMG GRID; explicitly an
        # exploratory CEST range around its supplied PB=.0583, KEX=116.
        coordinates = tuple(
            (key, lo, hi, "log")
            for name, lo, hi in ranges
            for key in parameterization.independent_ids
            if parameterization.role(key) is ParameterRole.FIT
            and model.definitions[key].name == name
        )
    else:
        plan = read_method_plan(
            [
                example
                / "Methods"
                / ("method_de.toml" if case == "dcest" else "method_grid.toml")
            ]
        )
        step = compile_method_plan(plan, model, experiments).steps[0]
        if not isinstance(step, FitStep):
            raise RuntimeError("Benchmark step has no profiles")
        for binding in step.bindings:
            binding.activate()
        parameterization = step.parameterization.bind(
            session.analysis_values.snapshot()
        )
        if case == "dcest":
            if not isinstance(step.search, DeSearchInstruction):
                raise RuntimeError("Expected shipped DE search")
            coordinates = tuple(
                (key, lo, hi, str(scale))
                for key, lo, hi, scale in step.search.coordinates
            )
        else:
            if not isinstance(step.search, GridSearchInstruction):
                raise RuntimeError("Expected shipped GRID search")
            coordinates = tuple(
                (axis.param_id, axis.values[0], axis.values[-1], "log")
                for axis in step.search.axes
                if model.definitions[axis.param_id].name in {"PB", "KEX_AB"}
            )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    problem = OptimizationProblem.from_native(
        engine.plan, parameterization, configuration, session.analysis_values.snapshot()
    )
    return Context(problem, parameterization, engine, configuration, coordinates)


def outer_points(coordinates, start, exponent: int, seed: int):
    """Keep the supplied start separately from the balanced Sobol design."""
    points = qmc.Sobol(len(coordinates), rng=np.random.default_rng(seed)).random_base2(
        exponent
    )
    bounds = np.asarray(
        [
            np.log((lo, hi)) if scale == "log" else (lo, hi)
            for _, lo, hi, scale in coordinates
        ]
    )
    mapped = bounds[:, 0] + points * (bounds[:, 1] - bounds[:, 0])
    for index, (_, _, _, scale) in enumerate(coordinates):
        if scale == "log":
            mapped[:, index] = np.exp(mapped[:, index])
    return [tuple(start[key] for key, *_ in coordinates), *map(tuple, mapped)]


def branch_vectors(problem, factor, normal):
    """Enumerate DW signs only inside this proven independent factor.

    Never invent a magnitude: use the supplied nonzero start, or a successful
    normal endpoint for a zero start. If both are zero, no branch is available.
    Multi-DW combinations are intentionally bounded by the factor's DW count.
    """
    fitted = dict(normal.nuisance_items)
    starts = dict(zip(problem.controlled_ids, problem.start, strict=True))
    magnitudes = {
        key: abs(starts[key] or fitted.get(key, 0.0))
        for key in factor.nuisance_ids
        if key.startswith("__DW_")
    }
    active = tuple(key for key, magnitude in magnitudes.items() if magnitude != 0.0)
    for signs in product((-1.0, 1.0), repeat=len(active)):
        updates = dict(
            zip(
                active,
                (
                    sign * magnitudes[key]
                    for key, sign in zip(active, signs, strict=True)
                ),
                strict=True,
            )
        )
        vector = tuple(
            updates.get(key, value)
            for key, value in zip(problem.controlled_ids, problem.start, strict=True)
        )
        if vector != problem.start:
            yield vector


def measured(
    context, operation, *, branches: bool = False, endpoint_mirror: bool = False
):
    calls = []
    model_calculations = 0
    original_candidate = grid.execute_direct_trf_candidate
    original_factor = grid._fit_factor_point
    original_calculate = BoundEvaluator._calculate_unscaled

    def count_calculations(*args, **kwargs):
        nonlocal model_calculations
        model_calculations += 1
        return original_calculate(*args, **kwargs)

    def capture(problem, invocation, parameterization, engine, **kwargs):
        outcome = original_candidate(
            problem, invocation, parameterization, engine, **kwargs
        )
        counters = outcome.execution.counters
        calls.append(
            {
                "requests": counters.objective_requests_accepted,
                "evaluations": counters.objective_evaluations_completed,
                # Work proxy, not cache-adjusted model calculations.
                "profile_evaluations_upper": counters.objective_evaluations_completed
                * len(engine.plan.profiles),
                "terminal": outcome.terminal.value,
                "ids": problem.controlled_ids,
                "start": problem.start,
                "chi2": None
                if outcome.candidate is None
                else outcome.candidate.chi_square,
            }
        )
        return outcome

    def branch_factor(
        problem, factor, engine, parameterization, ordinal, axes, **kwargs
    ):
        normal = original_factor(
            problem, factor, engine, parameterization, ordinal, axes, **kwargs
        )
        attempted = [normal]
        branch_problem = problem
        if endpoint_mirror and normal.status is grid.ProfiledGridPointStatus.SUCCESS:
            values = dict(zip(problem.controlled_ids, problem.start, strict=True))
            values.update(normal.nuisance_items)
            branch_problem = context.at(
                tuple(values[key] for key in problem.controlled_ids)
            )
        for vector in branch_vectors(branch_problem, factor, normal):
            alternate = context.at(vector)
            attempted.append(
                original_factor(
                    alternate, factor, engine, parameterization, ordinal, axes, **kwargs
                )
            )
        successful = [
            point
            for point in attempted
            if point.status is grid.ProfiledGridPointStatus.SUCCESS
        ]
        return (
            min(successful, key=lambda point: point.chi_square)
            if successful
            else normal
        )

    start = perf_counter()
    with (
        patch.object(grid, "execute_direct_trf_candidate", capture),
        patch.object(grouped, "execute_direct_trf_candidate", capture),
        patch.object(BoundEvaluator, "_calculate_unscaled", count_calculations),
        patch.object(
            grid, "_fit_factor_point", branch_factor if branches else original_factor
        ),
    ):
        outcome = operation()
    return outcome, {
        "seconds": perf_counter() - start,
        "model_calculations": model_calculations,
        "requests": sum(call["requests"] for call in calls),
        "evaluations": sum(call["evaluations"] for call in calls),
        "profile_evaluations_upper": sum(
            call["profile_evaluations_upper"] for call in calls
        ),
        "local_refinements": len(calls),
        "calls": calls,
    }


def refine(context, vector):
    problem = context.at(vector)
    decomposition = grouped.FitDecomposition.from_root(
        problem, context.parameterization, context.engine
    )
    outcome, trace = measured(
        context,
        lambda: grouped.execute_grouped_direct_trf(
            problem,
            decomposition,
            _build_invocation(problem, decomposition),
            context.parameterization,
            context.engine,
        ),
    )
    trace.update(
        {
            "terminal": outcome.terminal.value,
            "chi2": None
            if outcome.accepted_result is None
            else outcome.accepted_result.chi_square,
            "vector": None
            if outcome.accepted_result is None
            else outcome.accepted_result.vector,
        }
    )
    return trace


def selected_indices(records, points, coordinates, count: int, distance: float):
    ordered = sorted(
        (index for index, record in enumerate(records) if record["vector"] is not None),
        key=lambda index: records[index]["chi2"],
    )
    normalized = np.asarray(points, dtype=float).copy()
    for axis, (_, lo, hi, scale) in enumerate(coordinates):
        if scale == "log":
            normalized[:, axis] = np.log(normalized[:, axis])
            lo, hi = np.log((lo, hi))
        normalized[:, axis] = (normalized[:, axis] - lo) / (hi - lo)
    selected = []
    for index in ordered:
        if all(
            np.linalg.norm(normalized[index] - normalized[previous]) >= distance
            for previous in selected
        ):
            selected.append(index)
        if len(selected) == count:
            break
    return selected


def run_case(case: str, exponent: int, seed: int, nuisance_budget: int):
    context = build_context(case)
    root = context.problem
    # DCEST's shipped DE includes both local DW coordinates. Profile only the
    # four coupling coordinates; compare DE exactly as currently declared.
    outer = tuple(
        coordinate
        for coordinate in context.coordinates
        if not coordinate[0].startswith("__DW_")
    )
    points = outer_points(outer, dict(root.independent_items), exponent, seed)
    result = {
        "case": case,
        "nuisance_budget": nuisance_budget,
        "seed": seed,
        "ids": root.controlled_ids,
        "start": root.start,
        "coordinates": context.coordinates,
        "outer_coordinates": outer,
        "points": points,
        "profiles": len(context.engine.plan.profiles),
        "residuals": context.engine.plan.retained_observation_count,
        "normal": refine(context, root.start),
    }
    print(case, "normal", result["normal"]["chi2"], flush=True)
    result["multistart"] = []
    for point in points:
        updates = dict(zip((key for key, *_ in outer), point, strict=True))
        vector = tuple(
            updates.get(key, value)
            for key, value in zip(root.controlled_ids, root.start, strict=True)
        )
        result["multistart"].append(refine(context, vector))
        print(
            case,
            "multistart",
            len(result["multistart"]),
            result["multistart"][-1]["chi2"],
            flush=True,
        )
    polished = {}
    for branch_aware in (False, True):
        mode = "profiled_branches" if branch_aware else "profiled"
        records = []
        for point in points:
            axes = {
                coordinate[0]: (value,)
                for coordinate, value in zip(outer, point, strict=True)
            }
            outcome, trace = measured(
                context,
                lambda axes=axes: grid.execute_profiled_grid(
                    root,
                    axes,
                    context.parameterization,
                    context.engine,
                    objective_request_budget=nuisance_budget,
                ),
                branches=branch_aware,
            )
            trace.update(
                {
                    "terminal": outcome.terminal.value,
                    "chi2": None
                    if outcome.accepted_result is None
                    else outcome.accepted_result.chi_square,
                    "vector": None
                    if outcome.accepted_result is None
                    else outcome.accepted_result.vector,
                    "factor_sizes": [
                        len(item.factor.nuisance_ids) for item in outcome.factors
                    ],
                }
            )
            records.append(trace)
            print(case, mode, len(records), trace["chi2"], flush=True)
        selections = {}
        # Prefixes answer sample-count questions without repeating the search.
        for size in (3, 5, len(points)):
            for count, distance in ((1, 0.0), (2, 0.0), (4, 0.0), (4, 0.25)):
                indices = selected_indices(
                    records[:size], points[:size], outer, count, distance
                )
                key = f"n{size}_top{count}_distance{distance}"
                selections[key] = indices
                for index in indices:
                    polish_key = f"{mode}:{index}"
                    if polish_key not in polished:
                        polished[polish_key] = refine(context, records[index]["vector"])
        result[mode] = {"records": records, "selections": selections}
    result["polished"] = polished
    invocation = DeSearchInvocation.for_product_problem(
        root, search_coordinates=context.coordinates, root_seed=seed
    )
    de, de_trace = measured(
        context,
        lambda: execute_de_search(
            root, invocation, context.parameterization, context.engine
        ),
    )
    de_refinement = (
        refine(context, de.best_candidate.full_vector)
        if de.restart_eligible and de.best_candidate is not None
        else None
    )
    result["de"] = {
        "seconds": de_trace["seconds"],
        "model_calculations": de_trace["model_calculations"],
        "terminal": de.terminal.value,
        "requests": de.counters.objective_requests_accepted,
        "evaluations": de.counters.objective_evaluations_completed,
        "profile_evaluations_upper": de.counters.objective_evaluations_completed
        * len(context.engine.plan.profiles),
        "search_chi2": None
        if de.best_candidate is None
        else de.best_candidate.chi_square,
        "refinement": de_refinement,
    }
    print(
        case, "DE", None if de_refinement is None else de_refinement["chi2"], flush=True
    )
    return result


def run_diagnostics(case: str, seed: int):
    """Targeted alternative-start checks, separate from the outer-search runs."""
    context = build_context(case)
    root = context.problem
    outer = tuple(
        item for item in context.coordinates if not item[0].startswith("__DW_")
    )
    axes = {key: (dict(root.independent_items)[key],) for key, *_ in outer}

    def profile(local_context, *, endpoint_mirror):
        outcome, trace = measured(
            local_context,
            lambda: grid.execute_profiled_grid(
                local_context.problem,
                axes,
                local_context.parameterization,
                local_context.engine,
                objective_request_budget=4000,
            ),
            branches=True,
            endpoint_mirror=endpoint_mirror,
        )
        trace["terminal"] = outcome.terminal.value
        trace["chi2"] = (
            None
            if outcome.accepted_result is None
            else outcome.accepted_result.chi_square
        )
        if outcome.accepted_result is not None:
            trace["polish"] = refine(local_context, outcome.accepted_result.vector)
        return trace

    result = {"endpoint_mirror": profile(context, endpoint_mirror=True)}
    if case != "cpmg":
        return result
    dws = tuple(key for key in root.controlled_ids if key.startswith("__DW_"))
    starts = dict(zip(root.controlled_ids, root.start, strict=True))
    result["start_sensitivity"] = {}
    for name, sign in (("mirrored", -1), ("zero", 0)):
        vector = tuple(
            value * sign if key in dws else value for key, value in starts.items()
        )
        local = replace(context, problem=context.at(vector))
        result["start_sensitivity"][name] = {
            "full": refine(local, vector),
            "profile": profile(local, endpoint_mirror=False),
        }
    points = outer_points(outer, dict(root.independent_items), 3, seed)
    signs = np.where(
        qmc.Sobol(len(dws), rng=np.random.default_rng(seed)).random_base2(3) < 0.5,
        -1,
        1,
    )
    patterns = [(-1,) * len(dws), *map(tuple, signs)]
    result["signed_multistart"] = []
    for point, pattern in zip(points, patterns, strict=True):
        updates = dict(zip((key for key, *_ in outer), point, strict=True))
        updates.update(
            (key, float(sign) * abs(starts[key]))
            for key, sign in zip(dws, pattern, strict=True)
        )
        result["signed_multistart"].append(
            refine(
                context, tuple(updates.get(key, value) for key, value in starts.items())
            )
        )
    normal = refine(context, root.start)["vector"]
    result["same_vector_signs"] = []
    for key in dws:
        mirror = tuple(
            -value if param_id == key else value
            for param_id, value in zip(root.controlled_ids, normal, strict=True)
        )
        chi_squares = []
        for vector in (normal, mirror):
            problem = context.at(vector)
            frame = EvaluationFrame.from_lifecycle_frame(
                context.parameterization,
                problem.lifecycle_frame(vector, context.parameterization),
            )
            residuals = context.engine.new_evaluator().evaluate_residuals(frame)
            if isinstance(residuals, EvaluationFailure):
                raise RuntimeError(  # noqa: TRY004 - scientific failure, not invalid type
                    f"Scientific evaluation failed: {residuals}"
                )
            chi_squares.append(canonical_chi_square(residuals))
        result["same_vector_signs"].append({"id": key, "chi2": chi_squares})
    return result


def main():
    if Path(chemex_source).resolve() != ROOT / "src/chemex/__init__.py":
        raise RuntimeError("Benchmark must import this checkout's source")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("cpmg", "cest", "dcest"))
    parser.add_argument("--exponent", type=int, default=3)
    parser.add_argument("--seed", type=int, default=597)
    parser.add_argument("--nuisance-budget", type=int, default=2000)
    parser.add_argument("--diagnostics", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = (
        run_diagnostics(args.case, args.seed)
        if args.diagnostics
        else run_case(args.case, args.exponent, args.seed, args.nuisance_budget)
    )
    result["environment"] = {
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "logical_cpus": os.cpu_count(),
        "threads": {
            name: os.environ.get(name)
            for name in (
                "OPENBLAS_NUM_THREADS",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()

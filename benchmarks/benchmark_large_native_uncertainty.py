#!/usr/bin/env python3
"""Phase benchmark for the shipped 96-profile native covariance workload."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import statistics
import sys
import tempfile
from collections import Counter, defaultdict
from collections.abc import Callable
from dataclasses import asdict, dataclass
from functools import wraps
from pathlib import Path
from time import perf_counter
from typing import Any, cast
from unittest.mock import patch

import numpy as np
import scipy

import chemex.optimize.grouped_direct_trf as grouped_owner
import chemex.optimize.method_step as method_step_owner
import chemex.optimize.native_deterministic as product_owner
from chemex.chemex import run
from chemex.cli import build_parser
from chemex.containers.data import Data
from chemex.experiments.catalog.cpmg_15n_ip_0013 import Cpmg15N0013IpSequence
from chemex.nmr.spectrometer import Spectrometer
from chemex.optimize.direct_trf import DirectTrfCandidateOutcome
from chemex.optimize.grouped_direct_trf import (
    ComponentDisposition,
    GroupedDirectTrfOutcome,
)
from chemex.optimize.method_step import MethodStepOutcome
from chemex.optimize.uncertainty import UncertaintyEvidence
from chemex.runtime import AnalysisSession
from chemex.typing import Array

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/Experiments/CPMG_15N_IP_0013"
EXPECTED_VECTOR_SHA256 = (
    "df17c7858383672b876bd588f0207db45946e9e0f3db82476e0470f556e8f26d"
)


@dataclass(frozen=True, slots=True)
class LargeFitMeasurement:
    """One exact no-plot workload measurement and deterministic acceptance record."""

    component_seconds: float
    root_materialization_composition_seconds: float
    covariance_seconds: float
    output_seconds: float
    total_seconds: float
    profiles: int
    residuals: int
    controls: int
    constrained_outputs: int
    components: int
    successful_components: int
    objective_requests: int
    fit_acceptance_kernels: int
    covariance_residual_evaluations: int
    covariance_kernels: int
    retained_matrix_conversions: int
    chi_square: float
    rank: int
    vector_sha256: str
    covariance_sha256: str


def _timed[**P, R](
    phase: str,
    function: Callable[P, R],
    timings: dict[str, list[float]],
    active_phase: list[str],
) -> Callable[P, R]:
    @wraps(function)
    def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
        previous = active_phase[0]
        active_phase[0] = phase
        started = perf_counter()
        try:
            return function(*args, **kwargs)
        finally:
            timings[phase].append(perf_counter() - started)
            active_phase[0] = previous

    return wrapped


def _matrix_digest(matrix: tuple[tuple[float, ...], ...]) -> str:
    return hashlib.sha256(np.asarray(matrix, dtype=np.float64).tobytes()).hexdigest()


def _vector_digest(vector: tuple[float, ...]) -> str:
    return hashlib.sha256(np.asarray(vector, dtype=np.float64).tobytes()).hexdigest()


def _environment_context() -> dict[str, object]:
    """Record enough runtime context to interpret exploratory phase timings."""
    numpy_config = cast(dict[str, Any], np.show_config(mode="dicts"))
    dependencies = cast(
        dict[str, object],
        numpy_config.get("Build Dependencies", {}),
    )
    thread_variables = {
        name: os.environ.get(name)
        for name in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
        )
    }
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "blas": dependencies.get("blas"),
        "lapack": dependencies.get("lapack"),
        "thread_environment": thread_variables,
    }


def _benchmark_arguments(output: Path) -> argparse.Namespace:
    return build_parser().parse_args(
        [
            "fit",
            "-e",
            *(str(path) for path in sorted((EXAMPLE / "Experiments").glob("*.toml"))),
            "-p",
            str(EXAMPLE / "Parameters/parameters.toml"),
            "-o",
            str(output),
            "--plot",
            "nothing",
            "--workers",
            "1",
        ]
    )


def measure_large_native_uncertainty() -> LargeFitMeasurement:  # noqa: C901
    """Run and validate the exact retained-Jacobian large-fit acceptance case."""
    timings: dict[str, list[float]] = defaultdict(list)
    active_phase = ["setup"]
    kernel_calls: Counter[str] = Counter()
    retained_matrices: list[tuple[tuple[float, ...], ...]] = []
    retained_matrix_conversions = 0
    candidate_outcomes: list[DirectTrfCandidateOutcome] = []
    method_outcomes: list[MethodStepOutcome] = []
    original_numpy = grouped_owner.np
    original_calculate = Cpmg15N0013IpSequence.calculate

    class NumpyConversionProbe:
        def __getattr__(self, name: str) -> object:
            return getattr(original_numpy, name)

        def asarray(
            self,
            value: object,
            dtype: type[np.float64] | np.dtype[np.float64] | None = None,
        ) -> Array:
            nonlocal retained_matrix_conversions
            if any(value is matrix for matrix in retained_matrices):
                retained_matrix_conversions += 1
            return original_numpy.asarray(value, dtype=dtype)

    def count_kernel(
        self: Cpmg15N0013IpSequence,
        spectrometer: Spectrometer,
        data: Data,
    ) -> Array:
        kernel_calls[active_phase[0]] += 1
        return original_calculate(self, spectrometer, data)

    original_candidate = grouped_owner.execute_direct_trf_candidate

    @wraps(original_candidate)
    def capture_candidate(*args: Any, **kwargs: Any) -> DirectTrfCandidateOutcome:
        outcome = original_candidate(*args, **kwargs)
        candidate_outcomes.append(outcome)
        return outcome

    original_components = grouped_owner.execute_direct_trf_components

    @wraps(original_components)
    def capture_components(*args: Any, **kwargs: Any):
        outcomes = original_components(*args, **kwargs)
        retained_matrices.extend(
            item.final_residual_jacobian.matrix
            for item in outcomes
            if item.final_residual_jacobian is not None
        )
        return outcomes

    original_method_step = product_owner.execute_method_step

    @wraps(original_method_step)
    def capture_method_step(*args: Any, **kwargs: Any) -> MethodStepOutcome:
        outcome = original_method_step(*args, **kwargs)
        method_outcomes.append(outcome)
        return outcome

    started = perf_counter()
    with tempfile.TemporaryDirectory() as temporary:
        output = Path(temporary) / "Output"
        with (
            patch.object(Cpmg15N0013IpSequence, "calculate", count_kernel),
            patch.object(grouped_owner, "np", NumpyConversionProbe()),
            patch.object(
                grouped_owner, "execute_direct_trf_candidate", capture_candidate
            ),
            patch.object(
                grouped_owner,
                "execute_direct_trf_components",
                _timed("component", capture_components, timings, active_phase),
            ),
            patch.object(
                grouped_owner,
                "materialize_grouped_direct_trf",
                _timed(
                    "root",
                    grouped_owner.materialize_grouped_direct_trf,
                    timings,
                    active_phase,
                ),
            ),
            patch.object(
                method_step_owner,
                "_execute_uncertainty",
                _timed(
                    "covariance",
                    method_step_owner._execute_uncertainty,
                    timings,
                    active_phase,
                ),
            ),
            patch.object(product_owner, "execute_method_step", capture_method_step),
            patch.object(
                product_owner,
                "execute_post_fit",
                _timed(
                    "output",
                    product_owner.execute_post_fit,
                    timings,
                    active_phase,
                ),
            ),
        ):
            run(_benchmark_arguments(output), session=AnalysisSession.create())
    total_seconds = perf_counter() - started

    if len(method_outcomes) != 1 or len(candidate_outcomes) != 1:
        raise AssertionError(
            "Large-fit benchmark expected one method and one component"
        )
    method_outcome = method_outcomes[0]
    primary = cast(GroupedDirectTrfOutcome, method_outcome.primary_execution)
    accepted = method_outcome.accepted_result
    candidate = candidate_outcomes[0]
    if accepted is None or candidate.execution.backend is None:
        raise AssertionError(
            "Large-fit benchmark did not retain its accepted TRF result"
        )
    uncertainty = next(
        (
            artifact
            for derivation in method_outcome.derivations
            for artifact in derivation.artifacts
            if isinstance(artifact, UncertaintyEvidence)
        ),
        None,
    )
    if (
        uncertainty is None
        or uncertainty.residual_jacobian is None
        or uncertainty.rank_diagnostic is None
        or uncertainty.covariance is None
        or uncertainty.constraint_jacobian is None
    ):
        raise AssertionError("Large-fit benchmark covariance evidence is incomplete")

    successful_components = sum(
        item.disposition is ComponentDisposition.SUCCEEDED
        for item in primary.components
    )
    objective_requests = sum(
        item.execution.counters.objective_requests_accepted
        for item in candidate_outcomes
    )
    fit_acceptance_kernels = kernel_calls["component"] + kernel_calls["root"]
    covariance_kernels = kernel_calls["covariance"]
    measurement = LargeFitMeasurement(
        component_seconds=sum(timings["component"]),
        root_materialization_composition_seconds=sum(timings["root"]),
        covariance_seconds=sum(timings["covariance"]),
        output_seconds=sum(timings["output"]),
        total_seconds=total_seconds,
        profiles=len(accepted.evaluation_result.profiles),
        residuals=accepted.evaluation_result.residuals.size,
        controls=len(accepted.vector),
        constrained_outputs=len(uncertainty.constraint_jacobian.output_ids),
        components=len(primary.components),
        successful_components=successful_components,
        objective_requests=objective_requests,
        fit_acceptance_kernels=fit_acceptance_kernels,
        covariance_residual_evaluations=uncertainty.residual_jacobian.evaluation_count,
        covariance_kernels=covariance_kernels,
        retained_matrix_conversions=retained_matrix_conversions,
        chi_square=accepted.chi_square,
        rank=uncertainty.rank_diagnostic.rank,
        vector_sha256=_vector_digest(accepted.vector),
        covariance_sha256=_matrix_digest(uncertainty.covariance.covariance),
    )
    if (
        (measurement.profiles, measurement.residuals, measurement.controls)
        != (96, 2784, 146)
        or measurement.constrained_outputs != 147
        or (measurement.components, measurement.successful_components) != (1, 1)
        or measurement.objective_requests != 2058
        or measurement.fit_acceptance_kernels != 6912
        or measurement.covariance_residual_evaluations != 0
        or measurement.covariance_kernels != 0
        or measurement.retained_matrix_conversions > successful_components
        or not math.isclose(
            measurement.chi_square,
            12965.326784811383,
            rel_tol=1.0e-13,
            abs_tol=1.0e-9,
        )
        or measurement.rank != 146
        or measurement.vector_sha256 != EXPECTED_VECTOR_SHA256
    ):
        raise AssertionError(f"Large-fit acceptance contract failed: {measurement!r}")
    return measurement


def run_benchmark(runs: int = 3) -> None:
    """Run the median phase gate in the existing benchmark lane."""
    if runs < 1:
        raise ValueError("runs must be positive")
    measurements = tuple(measure_large_native_uncertainty() for _ in range(runs))
    component_median = statistics.median(
        item.component_seconds for item in measurements
    )
    post_fit_median = statistics.median(
        item.root_materialization_composition_seconds + item.covariance_seconds
        for item in measurements
    )
    ratio = post_fit_median / component_median
    if ratio > 0.5:
        raise AssertionError(
            "Root composition plus covariance exceeded the 0.5x component-time gate: "
            f"{ratio:.3f}"
        )
    print(
        json.dumps(
            {
                "workload": "CPMG_15N_IP_0013",
                "plot": "nothing",
                "environment": _environment_context(),
                "runs": [asdict(item) for item in measurements],
                "median_component_seconds": component_median,
                "median_root_plus_covariance_seconds": post_fit_median,
                "root_plus_covariance_to_component_ratio": ratio,
            },
            indent=2,
            sort_keys=True,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=int, default=3)
    arguments = parser.parse_args()
    if arguments.runs < 1:
        parser.error("--runs must be positive")
    run_benchmark(arguments.runs)


if __name__ == "__main__":
    main()

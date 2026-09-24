#!/usr/bin/env python3
"""Measure a few representative ChemEx workloads from shipped inputs."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import tomllib
from contextlib import redirect_stdout
from pathlib import Path
from statistics import median
from time import perf_counter

import numpy as np
import scipy

from chemex import __file__ as chemex_source
from chemex import __version__
from chemex.configuration.methods import Method, Selection
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import (
    EvaluationEngine,
    EvaluationFailure,
    EvaluationFrame,
)
from chemex.experiments.builder import build_experiments
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples" / "Experiments"
CPMG = EXAMPLES / "CPMG_15N_IP_0013"
DCEST = EXAMPLES / "DCEST_15N_HD_EXCH"
THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "MPLCONFIGDIR",
    "XDG_CACHE_HOME",
)


def _revision() -> str | None:
    git = shutil.which("git")
    if git is None:
        return None
    try:
        result = subprocess.run(  # noqa: S603 - fixed git command
            [git, "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def _source_tree_sha256() -> str:
    """Identify the Python source contents, including uncommitted edits."""
    digest = hashlib.sha256()
    for path in sorted((ROOT / "src" / "chemex").rglob("*.py")):
        digest.update(path.relative_to(ROOT).as_posix().encode())
        digest.update(b"\0")
        digest.update(hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest()


def _assert_source_checkout(*, child: bool) -> None:
    expected = (ROOT / "src" / "chemex" / "__init__.py").resolve()
    if Path(chemex_source).resolve() != expected:
        raise RuntimeError(
            f"Imported ChemEx source is not this checkout: {chemex_source}"
        )
    if child:
        completed = subprocess.run(
            [sys.executable, "-c", "import chemex; print(chemex.__file__)"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
        if (
            completed.returncode != 0
            or Path(completed.stdout.strip()).resolve() != expected
        ):
            raise RuntimeError(
                "Fit subprocess would import a different ChemEx source: "
                f"{completed.stdout.strip()}\n{completed.stderr.strip()}"
            )


def _environment(source_tree_sha256: str, runner_sha256: str) -> dict[str, object]:
    source_version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"][
        "version"
    ]
    return {
        "revision": _revision(),
        "python": sys.version.split()[0],
        "source_version": source_version,
        "installed_distribution_version": __version__,
        "chemex_source": str(chemex_source),
        "source_tree_sha256": source_tree_sha256,
        "runner_sha256": runner_sha256,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "logical_cpus": os.cpu_count(),
        "thread_environment": {name: os.environ.get(name) for name in THREAD_VARIABLES},
    }


def _residual_workload(name: str, runs: int) -> dict[str, object]:
    if name == "cpmg-residual":
        example = CPMG
        experiment = example / "Experiments" / "600mhz.toml"
        model = "2st"
        spin_name = "L3N"
    else:
        example = DCEST
        experiment = example / "Experiments" / "3hz.toml"
        model = "2st_hd"
        spin_name = "1N"

    parameters = example / "Parameters" / "parameters.toml"
    session = AnalysisSession.create()
    session.set_model(model)
    with redirect_stdout(io.StringIO()):
        experiments = build_experiments(
            [experiment],
            Selection(include=[SpinSystem.from_name(spin_name)], exclude=None),
            session=session,
        )
    session.parameters.set_defaults(read_defaults([parameters]))
    if not session.try_build_analysis_values():
        raise RuntimeError(
            "Cannot initialize benchmark parameters: "
            f"{session.parameter_factory.native_construction_error!r}"
        )
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        parameterization.frame_from_snapshot(session.analysis_values.snapshot()),
    )

    cold_times: list[float] = []
    warm_times: list[float] = []
    reference: np.ndarray | None = None
    for _ in range(runs):
        evaluator = engine.new_evaluator()
        start = perf_counter()
        cold = evaluator.evaluate_residuals(frame)
        cold_times.append(perf_counter() - start)
        if isinstance(cold, EvaluationFailure):
            raise RuntimeError(  # noqa: TRY004 - evaluation failure, not a type error
                f"Native residual evaluation failed: {cold}"
            )
        start = perf_counter()
        warm = evaluator.evaluate_residuals(frame)
        warm_times.append(perf_counter() - start)
        if isinstance(warm, EvaluationFailure) or not np.array_equal(cold, warm):
            raise RuntimeError("Cold and cached residual evaluations disagree")
        if reference is not None and not np.array_equal(reference, cold):
            raise RuntimeError("Repeated residual evaluations disagree")
        reference = cold

    if reference is None or not np.isfinite(reference).all():
        raise RuntimeError("Benchmark residuals are missing or non-finite")
    return {
        "workload": name,
        "kind": "native residual evaluation; construction excluded",
        "experiments": [str(experiment.relative_to(ROOT))],
        "parameters": str(parameters.relative_to(ROOT)),
        "model": model,
        "spin_system": spin_name,
        "profiles": len(engine.plan.profiles),
        "residuals": int(reference.size),
        "runs": runs,
        "cold_median_seconds": median(cold_times),
        "cached_median_seconds": median(warm_times),
        "residual_sha256": hashlib.sha256(reference.tobytes()).hexdigest(),
    }


def _fit_workload(runs: int, timeout: int) -> dict[str, object]:
    experiment = CPMG / "Experiments" / "600mhz.toml"
    parameters = CPMG / "Parameters" / "parameters.toml"
    times: list[float] = []
    for _ in range(runs):
        with tempfile.TemporaryDirectory(prefix="chemex-benchmark-") as temporary:
            output = Path(temporary) / "Output"
            command = [
                sys.executable,
                "-m",
                "chemex",
                "fit",
                "-e",
                str(experiment),
                "-p",
                str(parameters),
                "-d",
                "2st",
                "--include",
                "L3N",
                "--plot",
                "nothing",
                "--workers",
                "1",
                "--native-threads",
                "1",
                "-o",
                str(output),
            ]
            start = perf_counter()
            completed = subprocess.run(  # noqa: S603 - fixed ChemEx command
                command,
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
                timeout=timeout,
            )
            times.append(perf_counter() - start)
            if completed.returncode != 0:
                raise RuntimeError(
                    f"ChemEx fit exited {completed.returncode}:\n"
                    f"{completed.stdout}\n{completed.stderr}"
                )
            outcome = tomllib.loads(
                (output / "run_info" / "outcome.toml").read_text(encoding="utf-8")
            )
            if outcome.get("status") != "complete":
                raise RuntimeError(f"ChemEx fit did not complete: {outcome}")
    return {
        "workload": "cpmg-fit",
        "kind": "end-to-end child-process fit, including startup and output",
        "experiments": [str(experiment.relative_to(ROOT))],
        "parameters": str(parameters.relative_to(ROOT)),
        "model": "2st",
        "spin_system": "L3N",
        "plot": "nothing",
        "workers": 1,
        "native_threads": 1,
        "runs": runs,
        "median_seconds": median(times),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "workload",
        choices=("cpmg-residual", "dcest-residual", "cpmg-fit"),
    )
    parser.add_argument("--runs", type=int, help="Repeat count (residual: 3, fit: 1)")
    parser.add_argument(
        "--timeout", type=int, default=180, help="Fit timeout in seconds"
    )
    args = parser.parse_args()
    runs = (
        args.runs
        if args.runs is not None
        else (1 if args.workload == "cpmg-fit" else 3)
    )
    if runs < 1 or args.timeout < 1:
        parser.error("--runs and --timeout must be positive")
    _assert_source_checkout(child=args.workload == "cpmg-fit")
    source_tree_sha256 = _source_tree_sha256()
    runner_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    measurement = (
        _fit_workload(runs, args.timeout)
        if args.workload == "cpmg-fit"
        else _residual_workload(args.workload, runs)
    )
    if (
        _source_tree_sha256() != source_tree_sha256
        or hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != runner_sha256
    ):
        raise RuntimeError(
            "ChemEx source or benchmark runner changed during measurement"
        )
    print(
        json.dumps(
            {
                "environment": _environment(source_tree_sha256, runner_sha256),
                "measurement": measurement,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

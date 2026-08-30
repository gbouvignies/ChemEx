"""Canonical execution semantics for an authoritative Method Plan."""

from __future__ import annotations

import secrets
from pathlib import Path

from chemex.configuration.method_plan import (
    MethodPlan,
    ResamplingKind,
    ResamplingRequest,
    StatisticsPlan,
)
from chemex.containers.experiments import Experiments
from chemex.messages import (
    print_fitmethod,
    print_mcmc_no_vary_warning,
    print_no_data,
    print_running_statistics,
    print_step_name,
)
from chemex.optimize.mcmc import run_native_mcmc
from chemex.optimize.native_deterministic import (
    NativeDeterministicFit,
    run_native_deterministic,
)
from chemex.optimize.resampling import run_native_resampling
from chemex.run_info import RunInfo, mark_failure_stage
from chemex.runtime import AnalysisSession


def _run_resampling(
    experiments: Experiments,
    path: Path,
    kind: ResamplingKind,
    request: ResamplingRequest,
    fit: NativeDeterministicFit,
    *,
    session: AnalysisSession,
    run_info: RunInfo | None,
    step_name: str,
) -> None:
    root_seed = secrets.randbits(64) if request.seed is None else request.seed
    if run_info is not None:
        run_info.record_stochastic_operation(step_name, kind, root_seed)
    run_native_resampling(
        experiments,
        path,
        kind,
        request,
        fit,
        execution=session.execution,
        root_seed=root_seed,
    )


def _run_native_statistics(
    experiments: Experiments,
    path: Path,
    statistics: StatisticsPlan,
    fit: NativeDeterministicFit,
    *,
    session: AnalysisSession,
    run_info: RunInfo | None = None,
    step_name: str = "DEFAULT",
) -> None:
    committed = session.analysis_values.snapshot()
    try:
        for kind, request in (
            ("mc", statistics.mc),
            ("bs", statistics.bs),
            ("bsn", statistics.bsn),
        ):
            if request is not None:
                _run_resampling(
                    experiments,
                    path,
                    kind,
                    request,
                    fit,
                    session=session,
                    run_info=run_info,
                    step_name=step_name,
                )
        if statistics.mcmc is not None:
            print_running_statistics("MCMC")
            run_native_mcmc(
                experiments,
                fit,
                statistics.mcmc,
                path,
                execution=session.execution,
                seed_recorder=(
                    None
                    if run_info is None
                    else lambda seed: run_info.record_stochastic_operation(
                        step_name,
                        "mcmc",
                        seed,
                    )
                ),
            )
    finally:
        if session.analysis_values.snapshot() != committed:
            raise RuntimeError("Native statistics mutated the committed central fit")


def _requests_only_mcmc(statistics: StatisticsPlan) -> bool:
    return (
        statistics.mcmc is not None
        and statistics.mc is None
        and statistics.bs is None
        and statistics.bsn is None
    )


def _run_requested_native_statistics(
    experiments: Experiments,
    path: Path,
    statistics: StatisticsPlan | None,
    fit: NativeDeterministicFit | None,
    *,
    session: AnalysisSession,
    run_info: RunInfo | None = None,
    step_name: str = "DEFAULT",
) -> None:
    if statistics is None:
        return
    if fit is None:
        if _requests_only_mcmc(statistics):
            print_mcmc_no_vary_warning()
            return
        raise RuntimeError("Native resampling requires a committed deterministic fit")
    try:
        _run_native_statistics(
            experiments,
            path,
            statistics,
            fit,
            session=session,
            run_info=run_info,
            step_name=step_name,
        )
    except (Exception, KeyboardInterrupt) as error:
        mark_failure_stage(error, "statistics")
        raise


def execute_method_plan(
    experiments: Experiments,
    plan: MethodPlan,
    path: Path,
    plot_level: str,
    *,
    session: AnalysisSession,
    run_info: RunInfo | None = None,
) -> None:
    """Execute a validated Method Plan afresh against committed analysis state."""
    effective_actions = plan.effective_role_actions()
    for index, step in enumerate(plan.steps, start=1):
        section = step.name
        if section:
            print_step_name(section, index, len(plan.steps))

        experiments.select_profiles(step.selection)
        if not experiments:
            print_no_data()
            continue

        print_fitmethod("trf")
        parameterization = session.compile_parameterization_from_actions(
            effective_actions[section],
            experiments.param_ids,
        )
        step_path = path / section if len(plan.steps) > 1 else path
        step_name = section or "DEFAULT"
        fit = run_native_deterministic(
            experiments,
            step_path,
            plot_level,
            session=session,
            parameterization=parameterization,
            search=step.search,
            run_info=run_info,
            step_name=step_name,
        )
        _run_requested_native_statistics(
            experiments,
            step_path,
            step.statistics,
            fit,
            session=session,
            run_info=run_info,
            step_name=step_name,
        )

"""Behavioral tests for native minimization progress reporting."""

from __future__ import annotations

from io import StringIO
from unittest.mock import patch

import pytest
from rich.console import Console

from chemex.messages import MinimizationProgressReporter, UncertaintyProgressReporter
from chemex.optimize.progress import (
    FitProgressContext,
    ProgressEvent,
    ProgressPhase,
    ProgressRateLimiter,
)


def _event(
    phase: ProgressPhase,
    *,
    elapsed: float,
    completed: int = 0,
    current: float | None = None,
    best: float | None = None,
    terminal: str | None = None,
    observations: int = 8,
    parameters: int = 2,
) -> ProgressEvent:
    return ProgressEvent(
        phase=phase,
        solver_requests_received=completed,
        objective_requests_accepted=completed,
        objective_evaluations_completed=completed,
        current_vector=(1.0, 2.0),
        current_chi_square=current,
        best_chi_square=best,
        retained_observation_count=observations,
        controlled_parameter_count=parameters,
        objective_request_budget=100,
        elapsed_seconds=elapsed,
        terminal_status=terminal,
    )


def test_progress_policy_reports_start_improvement_heartbeat_and_terminal() -> None:
    limiter = ProgressRateLimiter(interactive=True)

    started = limiter.observe(_event(ProgressPhase.STARTED, elapsed=0.0))
    first = limiter.observe(
        _event(
            ProgressPhase.EVALUATED,
            elapsed=0.1,
            completed=1,
            current=100.0,
            best=100.0,
        )
    )
    below_threshold = limiter.observe(
        _event(
            ProgressPhase.EVALUATED,
            elapsed=1.0,
            completed=2,
            current=99.95,
            best=99.95,
        )
    )
    improved = limiter.observe(
        _event(
            ProgressPhase.EVALUATED,
            elapsed=2.0,
            completed=3,
            current=99.8,
            best=99.8,
        )
    )
    heartbeat = limiter.observe(
        _event(
            ProgressPhase.EVALUATED,
            elapsed=7.0,
            completed=4,
            current=100.5,
            best=99.8,
        )
    )
    terminal = limiter.observe(
        _event(
            ProgressPhase.TERMINATED,
            elapsed=7.1,
            completed=4,
            current=99.8,
            best=99.8,
            terminal="converged",
        )
    )

    assert started is not None
    assert first is not None
    assert first.relative_best_change is None
    assert below_threshold is None
    assert improved is not None
    assert improved.relative_best_change == pytest.approx(-0.002)
    assert heartbeat is not None
    assert heartbeat.relative_best_change == pytest.approx(0.0)
    assert terminal is not None
    assert terminal.event.terminal_status == "converged"


def test_noninteractive_progress_uses_ten_second_heartbeat() -> None:
    limiter = ProgressRateLimiter(interactive=False)
    limiter.observe(_event(ProgressPhase.STARTED, elapsed=0.0))
    early_improvement = limiter.observe(
        _event(
            ProgressPhase.EVALUATED,
            elapsed=0.1,
            completed=1,
            current=100.0,
            best=100.0,
        )
    )

    before_heartbeat = limiter.observe(
        _event(
            ProgressPhase.EVALUATED,
            elapsed=9.9,
            completed=2,
            current=90.0,
            best=90.0,
        )
    )
    heartbeat = limiter.observe(
        _event(
            ProgressPhase.EVALUATED,
            elapsed=10.0,
            completed=3,
            current=89.0,
            best=89.0,
        )
    )

    assert early_improvement is None
    assert before_heartbeat is None
    assert heartbeat is not None


@pytest.mark.parametrize(
    ("observations", "parameters", "expected"),
    [(8, 2, 2.0), (2, 2, 12.0), (1, 2, 12.0)],
)
def test_progress_reduced_chi_square_uses_existing_nonpositive_dof_convention(
    observations: int,
    parameters: int,
    expected: float,
) -> None:
    event = _event(
        ProgressPhase.EVALUATED,
        elapsed=1.0,
        current=12.0,
        best=12.0,
        observations=observations,
        parameters=parameters,
    )

    assert event.reduced_chi_square == expected


def _capturing_console(*, interactive: bool = False) -> tuple[Console, StringIO]:
    stream = StringIO()
    return (
        Console(
            file=stream,
            force_terminal=interactive,
            color_system=None,
            width=180,
        ),
        stream,
    )


def test_noninteractive_reporter_renders_truthful_context_and_final_result() -> None:
    output_console, stream = _capturing_console()
    context = FitProgressContext(2, 3, ("__R1A_A_L18CD1",))
    times = iter((0.0, 0.0, 12.4, 12.5, 12.6))
    reporter = MinimizationProgressReporter(
        output_console,
        interactive=False,
        retained_observation_count=8,
        controlled_parameter_count=2,
        step_name="STEP2",
        group_labels={frozenset(context.controlled_ids): "L18CD1"},
        clock=lambda: next(times),
    )

    with reporter:
        reporter.observe(context, _event(ProgressPhase.STARTED, elapsed=0.0))
        reporter.observe(
            context,
            _event(
                ProgressPhase.EVALUATED,
                elapsed=0.1,
                completed=1,
                current=12.0,
                best=12.0,
            ),
        )
        reporter.observe(
            context,
            _event(
                ProgressPhase.TERMINATED,
                elapsed=0.2,
                completed=1,
                current=12.0,
                best=12.0,
                terminal="success",
            ),
        )
        reporter.finish(final_chi_square=12.0, terminal_status="committed")

    rendered = stream.getvalue()
    assert "STEP2" in rendered
    assert "group 2/3" in rendered
    assert "L18CD1" in rendered
    assert "eval 1/100" in rendered
    assert "best χ² 12" in rendered
    assert "red. χ² 2" in rendered
    assert "eval 1 total" in rendered
    assert "final χ² 12" in rendered
    assert "committed" in rendered
    assert "Iteration" not in rendered


def test_uncertainty_reporter_separates_post_fit_elapsed_time_and_status() -> None:
    output_console, stream = _capturing_console()
    times = iter((20.0, 22.25))
    reporter = UncertaintyProgressReporter(
        output_console,
        clock=lambda: next(times),
    )

    reporter.start()
    reporter.finish("covariance available")

    assert stream.getvalue().splitlines() == [
        "  • Estimating parameter uncertainties...",
        "    covariance available · 2.2 s",
    ]


def test_noninteractive_reporter_rate_limits_group_transitions_fit_wide() -> None:
    output_console, stream = _capturing_console()
    times = iter((0.0, 0.0, 0.1, 0.2, 0.3, 10.1, 10.2, 10.3))
    reporter = MinimizationProgressReporter(
        output_console,
        interactive=False,
        retained_observation_count=8,
        controlled_parameter_count=2,
        clock=lambda: next(times),
    )

    with reporter:
        for ordinal in range(1, 4):
            context = FitProgressContext(ordinal, 3, (f"__P{ordinal}",))
            reporter.observe(context, _event(ProgressPhase.STARTED, elapsed=0.0))
            reporter.observe(
                context,
                _event(
                    ProgressPhase.TERMINATED,
                    elapsed=0.1,
                    completed=1,
                    current=12.0 - ordinal,
                    best=12.0 - ordinal,
                    terminal="success",
                ),
            )
        reporter.finish(final_chi_square=9.0, terminal_status="committed")

    rendered = stream.getvalue()
    assert rendered.count("group ") == 2
    assert "group 1/3" in rendered
    assert "group 3/3" in rendered
    assert "eval 3 total" in rendered


def test_disabled_reporter_emits_no_output() -> None:
    output_console, stream = _capturing_console()
    reporter = MinimizationProgressReporter(
        output_console,
        interactive=False,
        retained_observation_count=8,
        controlled_parameter_count=2,
        enabled=False,
    )

    with reporter:
        reporter.observe(
            FitProgressContext(1, 1, ("__PB",)),
            _event(ProgressPhase.STARTED, elapsed=0.0),
        )
        reporter.finish(final_chi_square=12.0, terminal_status="committed")

    assert stream.getvalue() == ""


def test_grid_reporter_aggregates_local_attempts_instead_of_streaming_each_seed() -> (
    None
):
    output_console, stream = _capturing_console()
    reporter = MinimizationProgressReporter(
        output_console,
        interactive=False,
        retained_observation_count=8,
        controlled_parameter_count=2,
        grid=True,
        clock=lambda: 0.0,
    )

    with reporter:
        for seed in range(1, 4):
            context = FitProgressContext(
                1,
                1,
                ("__PB",),
                grid_seed_ordinal=seed,
                grid_seed_total=3,
            )
            reporter.observe(
                context,
                _event(ProgressPhase.STARTED, elapsed=0.0),
            )
            reporter.observe(
                context,
                _event(
                    ProgressPhase.EVALUATED,
                    elapsed=0.1,
                    completed=1,
                    current=12.0,
                    best=12.0,
                ),
            )
            reporter.observe(
                context,
                _event(
                    ProgressPhase.TERMINATED,
                    elapsed=0.2,
                    completed=1,
                    current=12.0,
                    best=12.0,
                    terminal="success",
                ),
            )
        reporter.finish(final_chi_square=10.0, terminal_status="committed")

    rendered = stream.getvalue()
    assert rendered.count("GRID seed") == 1
    assert "eval 3 total" in rendered


def test_keyboard_interrupt_stops_live_reporting_and_prints_a_terminal_line() -> None:
    output_console, stream = _capturing_console(interactive=True)
    times = iter((0.0, 0.0, 2.0))
    reporter = MinimizationProgressReporter(
        output_console,
        interactive=True,
        retained_observation_count=8,
        controlled_parameter_count=2,
        clock=lambda: next(times),
    )

    with pytest.raises(KeyboardInterrupt), reporter:
        reporter.observe(
            FitProgressContext(1, 1, ("__PB",)),
            _event(ProgressPhase.STARTED, elapsed=0.0),
        )
        raise KeyboardInterrupt

    assert not reporter.is_active
    assert "interrupted" in stream.getvalue()


def test_exception_fallback_does_not_label_a_local_best_as_final() -> None:
    output_console, stream = _capturing_console()
    times = iter((0.0, 0.0, 1.0))
    reporter = MinimizationProgressReporter(
        output_console,
        interactive=False,
        retained_observation_count=8,
        controlled_parameter_count=2,
        clock=lambda: next(times),
    )

    with pytest.raises(RuntimeError, match="scientific failure"), reporter:
        reporter.observe(
            FitProgressContext(1, 2, ("__P1",)),
            _event(
                ProgressPhase.EVALUATED,
                elapsed=0.1,
                completed=1,
                current=12.0,
                best=12.0,
            ),
        )
        raise RuntimeError("scientific failure")

    rendered = stream.getvalue()
    assert "eval 1 total" in rendered
    assert "failed" in rendered
    assert "final χ²" not in rendered


@pytest.mark.parametrize("grid", [False, True])
def test_live_elapsed_time_is_fit_wide_across_local_attempts(grid: bool) -> None:
    output_console, stream = _capturing_console(interactive=True)
    times = iter((0.0, 0.0, 1.0, 6.0, 6.1))
    reporter = MinimizationProgressReporter(
        output_console,
        interactive=True,
        retained_observation_count=8,
        controlled_parameter_count=2,
        grid=grid,
        clock=lambda: next(times),
    )
    first = FitProgressContext(
        1,
        2,
        ("__P1",),
        grid_seed_ordinal=1 if grid else None,
        grid_seed_total=2 if grid else None,
    )
    second = FitProgressContext(
        2,
        2,
        ("__P2",),
        grid_seed_ordinal=2 if grid else None,
        grid_seed_total=2 if grid else None,
    )
    rendered_updates: list[str] = []

    with (
        patch.object(reporter, "_render", side_effect=rendered_updates.append),
        reporter,
    ):
        reporter.observe(first, _event(ProgressPhase.STARTED, elapsed=0.0))
        reporter.observe(second, _event(ProgressPhase.STARTED, elapsed=0.0))
        reporter.observe(
            second,
            _event(
                ProgressPhase.EVALUATED,
                elapsed=0.1,
                completed=1,
                current=12.0,
                best=12.0,
            ),
        )
        reporter.finish(final_chi_square=12.0, terminal_status="committed")

    assert any("6.0 s" in line for line in rendered_updates)

"""Semantic reporting tests for native MCMC sampling progress."""

from __future__ import annotations

from io import StringIO

import pytest
from rich.console import Console

from chemex.messages import McmcProgressReporter
from chemex.optimize.progress import McmcProgressEvent, McmcProgressPhase


def test_mcmc_progress_tracks_complete_ensemble_transitions() -> None:
    output = StringIO()
    events: list[McmcProgressEvent] = []
    reporter = McmcProgressReporter(
        Console(file=output, force_terminal=False),
        requested_steps=3,
        interactive=False,
        observer=events.append,
    )

    reporter.start()
    reporter.observe(1)
    reporter.observe(2)
    reporter.observe(3)
    reporter.finish("completed")

    assert [event.phase for event in events] == [
        McmcProgressPhase.STARTED,
        McmcProgressPhase.ADVANCED,
        McmcProgressPhase.ADVANCED,
        McmcProgressPhase.ADVANCED,
        McmcProgressPhase.TERMINATED,
    ]
    assert [event.completed_steps for event in events] == [0, 1, 2, 3, 3]
    assert all(event.requested_steps == 3 for event in events)
    assert events[-1].terminal_status == "completed"
    rendered = output.getvalue()
    assert "MCMC sampling 0/3" in rendered
    assert "MCMC sampling 3/3" in rendered
    assert "complete" in rendered


@pytest.mark.parametrize("terminal", ["interrupted", "failed"])
def test_mcmc_progress_closes_partial_operations(terminal: str) -> None:
    output = StringIO()
    events: list[McmcProgressEvent] = []
    reporter = McmcProgressReporter(
        Console(file=output, force_terminal=False),
        requested_steps=4,
        interactive=False,
        observer=events.append,
    )

    reporter.start()
    reporter.observe(1)
    reporter.finish(terminal)
    reporter.observe(2)

    assert events[-1].phase is McmcProgressPhase.TERMINATED
    assert events[-1].completed_steps == 1
    assert events[-1].terminal_status == terminal
    assert f"1/4 -> {terminal}" in output.getvalue()


def test_mcmc_progress_ignores_duplicate_or_out_of_range_steps() -> None:
    events: list[McmcProgressEvent] = []
    reporter = McmcProgressReporter(
        Console(file=StringIO(), force_terminal=False),
        requested_steps=2,
        interactive=False,
        observer=events.append,
    )

    reporter.start()
    reporter.observe(1)
    reporter.observe(1)
    reporter.observe(3)
    reporter.finish("interrupted")

    assert [event.completed_steps for event in events] == [0, 1, 1]


def test_mcmc_progress_observer_failure_cannot_escape_into_sampling() -> None:
    def fail_reporting(_event: McmcProgressEvent) -> None:
        raise RuntimeError("terminal reporter failed")

    output = StringIO()
    reporter = McmcProgressReporter(
        Console(file=output, force_terminal=False),
        requested_steps=2,
        interactive=False,
        observer=fail_reporting,
    )

    reporter.start()
    reporter.observe(1)
    reporter.finish("completed")

    assert "MCMC sampling 0/2" in output.getvalue()
    assert "1/2" not in output.getvalue()

"""Portable observation data and rate limiting for native minimization."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum

_RELATIVE_IMPROVEMENT_THRESHOLD = 1.0e-3
_INTERACTIVE_HEARTBEAT_SECONDS = 5.0
_NONINTERACTIVE_HEARTBEAT_SECONDS = 10.0


class ProgressPhase(StrEnum):
    """Lifecycle phase for one native Direct-TRF attempt."""

    STARTED = "started"
    EVALUATED = "evaluated"
    TERMINATED = "terminated"


@dataclass(frozen=True, slots=True)
class ProgressEvent:
    """Already-computed scalar state at one objective-observation boundary."""

    phase: ProgressPhase
    solver_requests_received: int
    objective_requests_accepted: int
    objective_evaluations_completed: int
    current_vector: tuple[float, ...]
    current_chi_square: float | None
    best_chi_square: float | None
    retained_observation_count: int
    controlled_parameter_count: int
    objective_request_budget: int
    elapsed_seconds: float
    terminal_status: str | None = None

    @property
    def reduced_chi_square(self) -> float | None:
        """Return reduced chi-square using ChemEx's established dof convention."""
        if self.best_chi_square is None:
            return None
        degrees_of_freedom = max(
            1,
            self.retained_observation_count - self.controlled_parameter_count,
        )
        return self.best_chi_square / degrees_of_freedom


type ProgressObserver = Callable[[ProgressEvent], None]


@dataclass(frozen=True, slots=True)
class FitProgressContext:
    """Orchestration context for one local Direct-TRF attempt."""

    group_ordinal: int
    group_total: int
    controlled_ids: tuple[str, ...]
    grid_seed_ordinal: int | None = None
    grid_seed_total: int | None = None


type ContextualProgressObserver = Callable[[FitProgressContext, ProgressEvent], None]


@dataclass(frozen=True, slots=True)
class ProgressUpdate:
    """One progress event selected for user-visible reporting."""

    event: ProgressEvent
    relative_best_change: float | None


class ProgressRateLimiter:
    """Select useful progress updates without changing objective execution."""

    def __init__(self, *, interactive: bool) -> None:
        self._interactive = interactive
        self._heartbeat_seconds = (
            _INTERACTIVE_HEARTBEAT_SECONDS
            if interactive
            else _NONINTERACTIVE_HEARTBEAT_SECONDS
        )
        self._last_reported_elapsed: float | None = None
        self._last_reported_best: float | None = None

    def observe(self, event: ProgressEvent) -> ProgressUpdate | None:
        """Return a visible update when the event satisfies the UX policy."""
        relative_change = self._relative_change(event.best_chi_square)
        report = event.phase in {ProgressPhase.STARTED, ProgressPhase.TERMINATED}

        if event.phase is ProgressPhase.EVALUATED:
            if self._interactive:
                report = (
                    self._last_reported_best is None
                    and event.best_chi_square is not None
                )
                if relative_change is not None:
                    report = (
                        report or relative_change <= -_RELATIVE_IMPROVEMENT_THRESHOLD
                    )
            if self._last_reported_elapsed is not None:
                report = report or (
                    event.elapsed_seconds - self._last_reported_elapsed
                    >= self._heartbeat_seconds
                )

        if not report:
            return None

        update = ProgressUpdate(event, relative_change)
        self._last_reported_elapsed = event.elapsed_seconds
        if event.best_chi_square is not None:
            self._last_reported_best = event.best_chi_square
        return update

    def _relative_change(self, best_chi_square: float | None) -> float | None:
        previous = self._last_reported_best
        if best_chi_square is None or previous is None:
            return None
        if previous == 0.0:
            return 0.0
        return (best_chi_square - previous) / abs(previous)


__all__ = [
    "ContextualProgressObserver",
    "FitProgressContext",
    "ProgressEvent",
    "ProgressPhase",
    "ProgressObserver",
    "ProgressRateLimiter",
    "ProgressUpdate",
]

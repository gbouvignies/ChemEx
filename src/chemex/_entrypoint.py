"""Lightweight terminal boundary shared by ChemEx command-line entry paths."""

from __future__ import annotations

import sys
from collections.abc import Callable
from contextlib import suppress
from importlib import import_module
from typing import cast

type DiagnosticRenderer = Callable[[str, tuple[str, ...], bool], None]

_BASE_EXCEPTION_DICT = BaseException.__dict__["__dict__"]


def _plain_diagnostic(
    message: str,
    notes: tuple[str, ...],
    interrupted: bool,
) -> None:
    prefix = "" if interrupted else "ERROR: "
    sys.stderr.write(f"{prefix}{message}\n")
    for note in notes:
        sys.stderr.write(f"NOTE: {note}\n")


def _instance_data(error: BaseException) -> dict[str, object]:
    """Read built-in exception storage without invoking subclass descriptors."""
    try:
        data = _BASE_EXCEPTION_DICT.__get__(error, type(error))
    except BaseException:  # noqa: BLE001 - diagnostics must be fail-safe
        return {}
    return data if type(data) is dict else {}


def _audited_notes(error: BaseException) -> tuple[str, ...]:
    notes = _instance_data(error).get("__notes__")
    if type(notes) not in {list, tuple}:
        return ()
    trusted_notes = cast("list[object] | tuple[object, ...]", notes)
    return tuple(
        note
        for note in trusted_notes
        if type(note) is str and str.startswith(note, "ChemEx ")
    )


def _is_interrupted(error: Exception) -> bool:
    terminal = _instance_data(error).get("terminal")
    return isinstance(terminal, str) and str.__eq__(terminal, "interrupted") is True


def _error_message(error: Exception) -> str:
    try:
        message = str(error)
    except BaseException:  # noqa: BLE001 - diagnostics must be fail-safe
        message = ""
    if message:
        return message
    try:
        name = type.__getattribute__(type(error), "__name__")
    except BaseException:  # noqa: BLE001 - diagnostics must be fail-safe
        name = "Exception"
    if type(name) is not str or not name:
        name = "Exception"
    return f"{name} occurred without an error message."


def _load_renderer() -> DiagnosticRenderer:
    """Load optional application rendering without making it a failure source."""
    try:
        messages = import_module("chemex.messages")
        renderer = messages.print_cli_diagnostic
    except BaseException:  # noqa: BLE001 - renderer acquisition is optional
        return _plain_diagnostic
    else:
        return renderer


def _render(
    renderer: DiagnosticRenderer,
    message: str,
    notes: tuple[str, ...],
    *,
    interrupted: bool,
) -> None:
    try:
        renderer(message, notes, interrupted)
    except BaseException:  # noqa: BLE001 - diagnostics cannot replace exit semantics
        with suppress(BaseException):
            _plain_diagnostic(message, notes, interrupted)


def main() -> None:
    """Run ChemEx and translate terminal failures without exposing tracebacks."""
    renderer = _plain_diagnostic
    try:
        renderer = _load_renderer()
        application = import_module("chemex.chemex")
        application.main()
    except KeyboardInterrupt as error:
        _render(
            renderer,
            "ChemEx interrupted.",
            _audited_notes(error),
            interrupted=True,
        )
        raise SystemExit(130) from None
    except Exception as error:  # noqa: BLE001 - terminal CLI translation seam
        interrupted = _is_interrupted(error)
        message = "ChemEx interrupted." if interrupted else _error_message(error)
        _render(
            renderer,
            message,
            _audited_notes(error),
            interrupted=interrupted,
        )
        raise SystemExit(130 if interrupted else 1) from None

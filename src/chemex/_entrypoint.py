"""Lightweight terminal boundary shared by ChemEx command-line entry paths."""

from __future__ import annotations

import sys
from collections.abc import Callable
from contextlib import suppress
from importlib import import_module

type DiagnosticRenderer = Callable[[str], None]

_BASE_EXCEPTION_DICT = BaseException.__dict__["__dict__"]


def _plain_diagnostic(message: str) -> None:
    sys.stderr.write(f"{message}\n")


def _instance_data(error: BaseException) -> dict[str, object]:
    """Read built-in exception storage without invoking subclass descriptors."""
    try:
        data = _BASE_EXCEPTION_DICT.__get__(error, type(error))
    except BaseException:  # noqa: BLE001 - diagnostics must be fail-safe
        return {}
    return data if type(data) is dict else {}


def _is_interrupted(error: Exception) -> bool:
    terminal = _instance_data(error).get("terminal")
    return isinstance(terminal, str) and str.__eq__(terminal, "interrupted") is True


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
) -> None:
    try:
        renderer(message)
    except BaseException:  # noqa: BLE001 - diagnostics cannot replace exit semantics
        with suppress(BaseException):
            _plain_diagnostic(message)


def main() -> None:
    """Run ChemEx and translate terminal failures without exposing tracebacks."""
    renderer = _plain_diagnostic
    try:
        renderer = _load_renderer()
        application = import_module("chemex.chemex")
        application.main()
    except KeyboardInterrupt:
        _render(renderer, "ChemEx interrupted.")
        raise SystemExit(130) from None
    except Exception as error:  # noqa: BLE001 - terminal CLI translation seam
        interrupted = _is_interrupted(error)
        _render(
            renderer,
            (
                "ChemEx interrupted."
                if interrupted
                else "ChemEx encountered an unexpected error."
            ),
        )
        raise SystemExit(130 if interrupted else 1) from None

"""Lightweight terminal boundary shared by ChemEx command-line entry paths."""

from __future__ import annotations

import sys
from collections.abc import Callable
from contextlib import suppress
from importlib import import_module
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from chemex.exceptions import ChemExError

type DiagnosticRenderer = Callable[[str], None]

_UNEXPECTED_MESSAGE = (
    "ChemEx encountered an unexpected internal error.\n"
    "Rerun with --debug to show the traceback and chained causes."
)

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


def _debug_requested(arguments: list[str]) -> bool:
    """Recognize explicit debug opt-in without scanning past end-of-options."""
    try:
        end_of_options = arguments.index("--")
    except ValueError:
        return "--debug" in arguments
    return "--debug" in arguments[:end_of_options]


def _known_failure_message(error: Exception) -> str:
    try:
        terminal = import_module("chemex.terminal")
        message = terminal.format_known_failure(cast("ChemExError", error))
    except BaseException:  # noqa: BLE001 - diagnostics must remain fail-safe
        return "ChemEx could not complete the requested operation."
    return message or "ChemEx could not complete the requested operation."


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
    """Run ChemEx and apply the three terminal failure policies."""
    renderer = _plain_diagnostic
    chemex_error_type: type[Exception] | None = None
    debug = _debug_requested(sys.argv[1:])
    try:
        renderer = _load_renderer()
        exceptions = import_module("chemex.exceptions")
        chemex_error_type = exceptions.ChemExError
        application = import_module("chemex.chemex")
        application.main()
    except KeyboardInterrupt as error:
        try:
            terminal = import_module("chemex.terminal")
            message = terminal.format_interruption(error)
        except BaseException:  # noqa: BLE001 - diagnostics must remain fail-safe
            message = "ChemEx interrupted."
        _render(renderer, message)
        raise SystemExit(130) from None
    except Exception as error:
        known = chemex_error_type is not None and isinstance(error, chemex_error_type)
        if known and _is_interrupted(error):
            try:
                terminal = import_module("chemex.terminal")
                message = terminal.format_interruption(error)
            except BaseException:  # noqa: BLE001 - diagnostics must remain fail-safe
                message = "ChemEx interrupted."
            _render(renderer, message)
            raise SystemExit(130) from None
        if known:
            _render(renderer, _known_failure_message(error))
            raise SystemExit(1) from None
        if debug:
            raise
        _render(renderer, _UNEXPECTED_MESSAGE)
        raise SystemExit(1) from None

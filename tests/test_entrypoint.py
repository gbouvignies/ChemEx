from __future__ import annotations

import os
import subprocess
import sys
import textwrap
import tomllib
from collections.abc import Sequence
from configparser import Error as ConfigParserError
from pathlib import Path

import pytest

from chemex import chemex as application_module

TRACEBACK_HEADER = "Traceback (most recent call last)"
SECRET_CANARY = "SECRET_CANARY"  # noqa: S105 - security regression canary
UNEXPECTED_MESSAGE = "ChemEx encountered an unexpected error."
INTERRUPTED_MESSAGE = "ChemEx interrupted."
ROOT = Path(__file__).parent.parent
EXAMPLE = ROOT / "examples/Experiments/RELAXATION_HZNZ"
EXPERIMENT = EXAMPLE / "Experiments/800mhz.toml"
PARAMETERS = EXAMPLE / "Parameters/parameters.toml"


def _entry_commands() -> tuple[tuple[str, ...], ...]:
    return (
        (str(Path(sys.executable).with_name("chemex")),),
        (sys.executable, "-m", "chemex"),
    )


def _run(
    command: Sequence[str],
    *arguments: str,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(  # noqa: S603 - fixed local test entry paths
        (*command, *arguments),
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )


def _run_bootstrap(
    failure_source: str,
    *,
    setup_source: str = "",
    messages_source: str = "",
) -> subprocess.CompletedProcess[str]:
    failure = textwrap.indent(textwrap.dedent(failure_source).strip(), "    ")
    source = f"""
import sys
import types
{textwrap.dedent(setup_source)}
{textwrap.dedent(messages_source)}
application = types.ModuleType("chemex.chemex")

def fail():
{failure}

application.main = fail
sys.modules["chemex.chemex"] = application

from chemex._entrypoint import main
main()
"""
    return _run((sys.executable, "-c"), source)


def _assert_bootstrap_result(
    completed: subprocess.CompletedProcess[str],
    status: int,
    message: str,
    *,
    absent: Sequence[str] = (),
) -> None:
    combined = completed.stdout + completed.stderr
    assert completed.returncode == status
    assert completed.stdout == ""
    assert completed.stderr == f"{message}\n"
    assert all(item not in combined for item in absent)
    assert TRACEBACK_HEADER not in combined


def _application_import_failure_env(
    tmp_path: Path,
    failure_source: str,
) -> dict[str, str]:
    (tmp_path / "sitecustomize.py").write_text(
        f"""
import importlib.abc
import sys


class BlockApplicationImport(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "chemex.chemex":
            raise {failure_source}
        return None


sys.meta_path.insert(0, BlockApplicationImport())
""".lstrip(),
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, (str(tmp_path), env.get("PYTHONPATH", "")))
    )
    return env


@pytest.mark.parametrize("command", _entry_commands())
def test_cli_entrypoints_translate_application_import_failure(
    command: tuple[str, ...],
    tmp_path: Path,
) -> None:
    env = _application_import_failure_env(
        tmp_path,
        f"ImportError('{SECRET_CANARY}')",
    )

    completed = _run(command, env=env)

    _assert_bootstrap_result(
        completed,
        1,
        UNEXPECTED_MESSAGE,
        absent=(SECRET_CANARY,),
    )


@pytest.mark.parametrize("command", _entry_commands())
@pytest.mark.parametrize(
    ("argument", "status", "stream", "message"),
    (
        ("--help", 0, "stdout", "usage:"),
        ("--version", 0, "stdout", "chemex"),
        (
            "--not-a-real-option",
            2,
            "stderr",
            "unrecognized arguments: --not-a-real-option",
        ),
    ),
)
def test_cli_entrypoints_preserve_argparse_behavior(
    command: tuple[str, ...],
    argument: str,
    status: int,
    stream: str,
    message: str,
) -> None:
    completed = _run(command, argument)

    combined = completed.stdout + completed.stderr
    assert completed.returncode == status
    assert message in getattr(completed, stream).lower()
    assert TRACEBACK_HEADER not in combined


@pytest.mark.parametrize("command", _entry_commands())
@pytest.mark.parametrize("model", ("missing-model", "2st.unsupported"))
def test_cli_entrypoints_report_invalid_model_as_unsuccessful_error(
    command: tuple[str, ...],
    model: str,
) -> None:
    completed = _run(
        command,
        "simulate",
        "-e",
        str(EXPERIMENT),
        "-p",
        str(PARAMETERS),
        "-d",
        model,
        "--plot",
        "nothing",
    )

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 1
    assert "ERROR:" in completed.stderr
    assert model in completed.stderr
    assert TRACEBACK_HEADER not in combined


@pytest.mark.parametrize("command", _entry_commands())
def test_cli_entrypoints_preserve_no_selected_profiles_success(
    command: tuple[str, ...],
) -> None:
    completed = _run(
        command,
        "simulate",
        "-e",
        str(EXPERIMENT),
        "-p",
        str(PARAMETERS),
        "--include",
        "DOES-NOT-EXIST",
        "--plot",
        "nothing",
    )

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 0
    assert "No data to fit" in completed.stdout
    assert TRACEBACK_HEADER not in combined


@pytest.mark.parametrize("command", _entry_commands())
def test_cli_entrypoints_translate_uncaught_runtime_failure(
    command: tuple[str, ...],
) -> None:
    completed = _run(
        command,
        "plot_param",
        "-p",
        str(ROOT / "pyproject.toml"),
        "-n",
        "BAD",
    )

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 1
    assert completed.stderr == f"{UNEXPECTED_MESSAGE}\n"
    assert "parsing errors" not in combined
    assert TRACEBACK_HEADER not in combined


def test_programmatic_main_retains_ordinary_exception_semantics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        ["chemex", "plot_param", "-p", "pyproject.toml", "-n", "BAD"],
    )

    with pytest.raises(ConfigParserError):
        application_module.main()


@pytest.mark.parametrize(
    ("failure_source", "expected_status", "expected_message"),
    [
        ("raise KeyboardInterrupt", 130, INTERRUPTED_MESSAGE),
        (
            "; ".join(
                (
                    f"error = RuntimeError('{SECRET_CANARY}')",
                    f"error.add_note('ChemEx {SECRET_CANARY}')",
                    "error.terminal = 'interrupted'",
                    "raise error",
                )
            ),
            130,
            INTERRUPTED_MESSAGE,
        ),
        (f"raise RuntimeError('{SECRET_CANARY}')", 1, UNEXPECTED_MESSAGE),
        ("raise RuntimeError('[bold]literal[/]')", 1, UNEXPECTED_MESSAGE),
        ("raise RuntimeError", 1, UNEXPECTED_MESSAGE),
    ],
)
def test_cli_bootstrap_translates_terminal_failures(
    failure_source: str,
    expected_status: int,
    expected_message: str,
) -> None:
    completed = _run_bootstrap(failure_source)

    _assert_bootstrap_result(
        completed,
        expected_status,
        expected_message,
        absent=(SECRET_CANARY,),
    )


def test_cli_bootstrap_ignores_all_exception_notes() -> None:
    completed = _run_bootstrap(
        f"""
        error = RuntimeError("{SECRET_CANARY}")
        error.add_note("ChemEx {SECRET_CANARY}")
        error.add_note("third-party {SECRET_CANARY}")
        raise error
        """
    )

    _assert_bootstrap_result(
        completed,
        1,
        UNEXPECTED_MESSAGE,
        absent=(SECRET_CANARY, "NOTE:"),
    )


@pytest.mark.parametrize(
    ("failure_source", "setup_source", "expected_status", "message", "absent"),
    (
        pytest.param(
            "error = RuntimeError('application failed'); "
            "error.__notes__ = None; raise error",
            "",
            1,
            UNEXPECTED_MESSAGE,
            (),
            id="notes-none",
        ),
        pytest.param(
            "error = RuntimeError('application failed'); "
            "error.__notes__ = object(); raise error",
            "",
            1,
            UNEXPECTED_MESSAGE,
            (),
            id="notes-non-sequence",
        ),
        pytest.param(
            """
            HostileNote = type('HostileNote', (), {
                '__str__': lambda self: (_ for _ in ()).throw(
                    RuntimeError('note stringification leaked')
                )
            })
            error = RuntimeError('SECRET_CANARY')
            error.__notes__ = [HostileNote(), 'ChemEx SECRET_CANARY']
            raise error
            """,
            "",
            1,
            UNEXPECTED_MESSAGE,
            ("note stringification leaked", SECRET_CANARY),
            id="hostile-note",
        ),
        pytest.param(
            """
            class HostileError(RuntimeError):
                def __str__(self):
                    raise RuntimeError("stringification leaked")
            raise HostileError("SECRET_CANARY")
            """,
            "",
            1,
            UNEXPECTED_MESSAGE,
            ("stringification leaked", SECRET_CANARY),
            id="hostile-message",
        ),
        pytest.param(
            """
            class HostileTerminalError(RuntimeError):
                @property
                def terminal(self):
                    raise RuntimeError("terminal property leaked")
            raise HostileTerminalError("SECRET_CANARY")
            """,
            "",
            1,
            UNEXPECTED_MESSAGE,
            ("terminal property leaked", SECRET_CANARY),
            id="hostile-terminal-property",
        ),
        pytest.param(
            """
            error = RuntimeError("SECRET_CANARY")
            error.add_note("ChemEx SECRET_CANARY")
            error.terminal = Terminal.INTERRUPTED
            raise error
            """,
            """
            from enum import StrEnum
            class Terminal(StrEnum):
                INTERRUPTED = "interrupted"
            """,
            130,
            INTERRUPTED_MESSAGE,
            (SECRET_CANARY,),
            id="strenum-terminal",
        ),
    ),
)
def test_cli_bootstrap_handles_hostile_exception_metadata(
    failure_source: str,
    setup_source: str,
    expected_status: int,
    message: str,
    absent: tuple[str, ...],
) -> None:
    completed = _run_bootstrap(failure_source, setup_source=setup_source)

    _assert_bootstrap_result(
        completed,
        expected_status,
        message,
        absent=absent,
    )


@pytest.mark.parametrize(
    ("phase", "renderer_failure"),
    (
        *(
            ("acquisition", failure)
            for failure in (
                "ImportError('renderer unavailable')",
                "SystemExit(71)",
                "KeyboardInterrupt()",
            )
        ),
        *(
            ("execution", failure)
            for failure in (
                "RuntimeError('renderer failed')",
                "SystemExit(77)",
                "KeyboardInterrupt()",
            )
        ),
    ),
)
@pytest.mark.parametrize(
    ("failure_source", "expected_status", "expected_message"),
    (
        (f"raise RuntimeError('{SECRET_CANARY}')", 1, UNEXPECTED_MESSAGE),
        ("raise KeyboardInterrupt", 130, INTERRUPTED_MESSAGE),
    ),
)
def test_cli_bootstrap_renderer_failures_use_plain_stderr(
    phase: str,
    renderer_failure: str,
    failure_source: str,
    expected_status: int,
    expected_message: str,
) -> None:
    messages_source = (
        f"""
        class BrokenMessages(types.ModuleType):
            def __getattribute__(self, name):
                if name == "print_cli_diagnostic":
                    raise {renderer_failure}
                return super().__getattribute__(name)
        sys.modules["chemex.messages"] = BrokenMessages("chemex.messages")
        """
        if phase == "acquisition"
        else f"""
        messages = types.ModuleType("chemex.messages")
        def broken_renderer(*args, **kwargs):
            raise {renderer_failure}
        messages.print_cli_diagnostic = broken_renderer
        sys.modules["chemex.messages"] = messages
        """
    )
    completed = _run_bootstrap(
        failure_source,
        messages_source=messages_source,
    )

    _assert_bootstrap_result(
        completed,
        expected_status,
        expected_message,
        absent=("renderer unavailable", "renderer failed", SECRET_CANARY),
    )


@pytest.mark.parametrize(
    ("import_failure", "expected_status"),
    ((f"SyntaxError('{SECRET_CANARY}')", 1), ("SystemExit(23)", 23)),
)
def test_cli_bootstrap_handles_protected_application_import_baseexceptions(
    import_failure: str,
    expected_status: int,
    tmp_path: Path,
) -> None:
    env = _application_import_failure_env(tmp_path, import_failure)

    completed = _run((sys.executable, "-m", "chemex"), env=env)

    if expected_status == 1:
        _assert_bootstrap_result(
            completed,
            1,
            UNEXPECTED_MESSAGE,
            absent=(SECRET_CANARY,),
        )
    else:
        assert completed.returncode == expected_status
        assert completed.stdout == completed.stderr == ""


def test_real_fit_interruption_flows_through_run_info_and_shared_bootstrap(
    tmp_path: Path,
) -> None:
    (tmp_path / "sitecustomize.py").write_text(
        """
import chemex.optimize.direct_trf as direct_trf

original = direct_trf.least_squares
calls = 0

def interrupt_third_fit(*args, **kwargs):
    global calls
    calls += 1
    if calls == 3:
        raise KeyboardInterrupt
    return original(*args, **kwargs)

direct_trf.least_squares = interrupt_third_fit
""".lstrip(),
        encoding="utf-8",
    )
    method = tmp_path / "method.toml"
    method.write_text(
        '[STEP1]\nFIX = ["PB", "KEX_AB"]\n\n[STEP2]\nFIX = ["PB", "KEX_AB"]\n',
        encoding="utf-8",
    )
    output = tmp_path / "Output"
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, (str(tmp_path), env.get("PYTHONPATH", "")))
    )

    completed = _run(
        (sys.executable, "-m", "chemex"),
        "fit",
        "-e",
        str(EXPERIMENT),
        "-p",
        str(PARAMETERS),
        "-m",
        str(method),
        "-o",
        str(output),
        "-d",
        "2st",
        "--include",
        "G2N-HN",
        "H3N-HN",
        "--plot",
        "nothing",
        "--workers",
        "1",
        env=env,
    )

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 130
    assert completed.stderr.count("ChemEx interrupted.") == 1
    assert TRACEBACK_HEADER not in combined
    assert (output / "STEP1" / "Parameters" / "fitted.toml").is_file()
    assert not (output / "STEP2" / "Parameters" / "fitted.toml").exists()
    outcome = tomllib.loads(
        (output / "run_info" / "outcome.toml").read_text(encoding="utf-8")
    )
    assert outcome["status"] == "incomplete"
    assert outcome["terminal"] == "interrupted"
    assert outcome["latest_committed_revision"] == 1

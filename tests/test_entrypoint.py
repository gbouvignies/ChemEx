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
from chemex.configuration.parameters import ParameterConfigurationError
from chemex.exceptions import InputFileReadError

TRACEBACK_HEADER = "Traceback (most recent call last)"
SECRET_CANARY = "SECRET_CANARY"  # noqa: S105 - security regression canary
UNEXPECTED_MESSAGE = (
    "ChemEx encountered an unexpected internal error.\n"
    "Rerun with --debug to show the traceback and chained causes."
)
INTERRUPTED_MESSAGE = "Analysis interrupted by user."
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
    arguments: Sequence[str] = (),
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
    return _run((sys.executable, "-c"), source, *arguments)


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
@pytest.mark.parametrize("subcommand", ("fit", "simulate", "pick_cest", "plot_param"))
@pytest.mark.parametrize("debug_first", (False, True))
def test_global_debug_option_is_discoverable_for_every_command(
    command: tuple[str, ...],
    subcommand: str,
    debug_first: bool,
) -> None:
    arguments = (
        ("--debug", subcommand, "--help")
        if debug_first
        else (subcommand, "--debug", "--help")
    )
    completed = _run(command, *arguments)

    assert completed.returncode == 0
    assert "--debug" in completed.stdout
    assert TRACEBACK_HEADER not in completed.stderr


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
    assert "Exchange model selection is invalid" in completed.stderr
    assert model in completed.stderr
    assert TRACEBACK_HEADER not in combined


@pytest.mark.parametrize("command", _entry_commands())
@pytest.mark.parametrize("input_kind", ("experiment", "parameters"))
def test_cli_entrypoints_report_malformed_input_with_source(
    command: tuple[str, ...],
    input_kind: str,
    tmp_path: Path,
) -> None:
    malformed = tmp_path / f"malformed-{input_kind}.toml"
    malformed.write_text("[broken\n", encoding="utf-8")
    experiment = malformed if input_kind == "experiment" else EXPERIMENT
    parameters = malformed if input_kind == "parameters" else PARAMETERS

    completed = _run(
        command,
        "simulate",
        "-e",
        str(experiment),
        "-p",
        str(parameters),
        "--plot",
        "nothing",
    )

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 1
    assert str(malformed) in completed.stderr
    assert "TOML" in completed.stderr
    assert "unexpected internal error" not in completed.stderr
    assert TRACEBACK_HEADER not in combined


@pytest.mark.parametrize("command", _entry_commands())
def test_cli_entrypoints_report_invalid_method_plan_with_source(
    command: tuple[str, ...],
    tmp_path: Path,
) -> None:
    method = tmp_path / "invalid-method.toml"
    method.write_text("[DEFAULT]\nFIT = 42\n", encoding="utf-8")

    completed = _run(
        command,
        "fit",
        "-e",
        str(EXPERIMENT),
        "-p",
        str(PARAMETERS),
        "-m",
        str(method),
        "--plot",
        "nothing",
    )

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 1
    assert "Method Plan is invalid" in completed.stderr
    assert str(method) in completed.stderr
    assert "unexpected internal error" not in completed.stderr
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
def test_cli_entrypoints_report_invalid_plot_parameter_input(
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
    assert "Parameter configuration is invalid" in completed.stderr
    assert str(ROOT / "pyproject.toml") in completed.stderr
    assert "parsing errors" not in combined
    assert "unexpected internal error" not in combined
    assert TRACEBACK_HEADER not in combined


@pytest.mark.parametrize("command", _entry_commands())
def test_plot_param_multiple_inputs_is_a_known_cli_failure(
    command: tuple[str, ...],
) -> None:
    completed = _run(
        command,
        "plot_param",
        "-p",
        str(PARAMETERS),
        str(PARAMETERS),
        "-n",
        "PB",
    )

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 1
    assert completed.stderr.count("Multiple parameter files were given") == 1
    assert "unexpected internal error" not in combined
    assert TRACEBACK_HEADER not in combined


@pytest.mark.parametrize("command", _entry_commands())
def test_plot_param_missing_input_is_a_known_cli_failure(
    command: tuple[str, ...],
    tmp_path: Path,
) -> None:
    missing = tmp_path / "missing-fitted.toml"
    completed = _run(
        command,
        "plot_param",
        "-p",
        str(missing),
        "-n",
        "PB",
    )

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 1
    assert "Input file could not be read" in completed.stderr
    assert str(missing) in completed.stderr
    assert "unexpected internal error" not in combined
    assert TRACEBACK_HEADER not in combined


def test_plot_param_missing_input_retains_original_oserror(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    missing = tmp_path / "missing-fitted.toml"
    monkeypatch.setattr(
        sys,
        "argv",
        ["chemex", "plot_param", "-p", str(missing), "-n", "PB"],
    )

    with pytest.raises(InputFileReadError) as error_info:
        application_module.main()

    assert error_info.value.path == missing
    assert isinstance(error_info.value.__cause__, FileNotFoundError)
    assert error_info.value.error is error_info.value.__cause__


def test_programmatic_main_retains_ordinary_exception_semantics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        ["chemex", "plot_param", "-p", "pyproject.toml", "-n", "BAD"],
    )

    with pytest.raises(ParameterConfigurationError) as error_info:
        application_module.main()

    assert isinstance(error_info.value.__cause__, ConfigParserError)


@pytest.mark.parametrize(
    ("failure_source", "expected_status", "expected_message"),
    [
        ("raise KeyboardInterrupt", 130, INTERRUPTED_MESSAGE),
        (
            "; ".join(
                (
                    "from chemex.exceptions import ChemExError",
                    f"error = ChemExError('{SECRET_CANARY}')",
                    f"error.add_note('ChemEx {SECRET_CANARY}')",
                    "error.terminal = 'interrupted'",
                    "raise error",
                )
            ),
            130,
            INTERRUPTED_MESSAGE,
        ),
        (f"raise RuntimeError('{SECRET_CANARY}')", 1, UNEXPECTED_MESSAGE),
        (f"raise AssertionError('{SECRET_CANARY}')", 1, UNEXPECTED_MESSAGE),
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


def test_cli_bootstrap_renders_known_chemex_failure_once() -> None:
    completed = _run_bootstrap(
        "from chemex.exceptions import ChemExError; "
        "raise ChemExError('Known ChemEx failure.')"
    )

    _assert_bootstrap_result(completed, 1, "Known ChemEx failure.")


def test_cli_bootstrap_keeps_known_failure_concise_in_debug_mode() -> None:
    completed = _run_bootstrap(
        "from chemex.exceptions import ChemExError; "
        "raise ChemExError('Known ChemEx failure.')",
        arguments=("--debug",),
    )

    _assert_bootstrap_result(completed, 1, "Known ChemEx failure.")


def test_cli_bootstrap_does_not_trust_interrupted_attribute_on_unknown_error() -> None:
    completed = _run_bootstrap(
        f"error = RuntimeError('{SECRET_CANARY}'); "
        "error.terminal = 'interrupted'; raise error"
    )

    _assert_bootstrap_result(
        completed,
        1,
        UNEXPECTED_MESSAGE,
        absent=(SECRET_CANARY,),
    )


def test_cli_bootstrap_debug_reraises_unexpected_failure_with_cause() -> None:
    completed = _run_bootstrap(
        f"""
        try:
            raise ValueError("ROOT_CAUSE")
        except ValueError as cause:
            error = RuntimeError("{SECRET_CANARY}")
            error.add_note("DEBUG_NOTE")
            raise error from cause
        """,
        arguments=("--debug",),
    )

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 1
    assert TRACEBACK_HEADER in combined
    assert "ValueError: ROOT_CAUSE" in combined
    assert f"RuntimeError: {SECRET_CANARY}" in combined
    assert "DEBUG_NOTE" in combined
    assert "The above exception was the direct cause" in combined


def test_debug_token_after_end_of_options_does_not_expose_internal_failure() -> None:
    completed = _run_bootstrap(
        f"raise RuntimeError('{SECRET_CANARY}')",
        arguments=("--", "--debug"),
    )

    _assert_bootstrap_result(
        completed,
        1,
        UNEXPECTED_MESSAGE,
        absent=(SECRET_CANARY,),
    )


def test_cli_bootstrap_reports_incomplete_resampling_artifacts_once() -> None:
    completed = _run_bootstrap(
        """
        from pathlib import Path
        from chemex.optimize.resampling import NativeResamplingIncompleteError
        error = NativeResamplingIncompleteError("completed 1 of 2 requested replicates")
        error.diagnostics_path = Path("Output/Statistics/MonteCarlo/diagnostics.toml")
        error.samples_path = Path("Output/Statistics/MonteCarlo/samples.tsv")
        error.failures_path = Path("Output/Statistics/MonteCarlo/failures.tsv")
        raise error
        """
    )

    rendered = " ".join(completed.stderr.split())
    assert completed.returncode == 1
    assert rendered.count("Resampling analysis is incomplete") == 1
    assert "completed 1 of 2 requested replicates" in rendered
    assert "diagnostics.toml" in rendered
    assert "samples.tsv" in rendered
    assert "failures.tsv" in rendered
    assert "Inspect the resampling diagnostics and failed replicates" in rendered
    assert "unexpected internal error" not in rendered


def test_incomplete_statistics_without_artifacts_does_not_invent_advice() -> None:
    completed = _run_bootstrap(
        """
        from chemex.optimize.mcmc import NativeMcmcIncompleteError
        raise NativeMcmcIncompleteError("scientific interpretation was withheld")
        """
    )

    rendered = " ".join(completed.stderr.split())
    assert completed.returncode == 1
    assert "MCMC analysis is incomplete" in rendered
    assert "Next:" not in rendered
    assert "Inspect the MCMC diagnostics" not in rendered


def test_cli_bootstrap_reports_recognized_deterministic_nonconvergence() -> None:
    completed = _run_bootstrap(
        """
        from types import SimpleNamespace
        from chemex.optimize.direct_trf import TerminalFailure
        from chemex.optimize.native_deterministic import NativeDeterministicAnalysisError
        failure = TerminalFailure("non_converged", "maximum evaluations reached")
        raise NativeDeterministicAnalysisError(
            "Native deterministic fitting did not produce an accepted fit.",
            reason="The optimizer did not converge.",
            outcome=SimpleNamespace(),
            failures=(failure,),
        )
        """
    )

    rendered = " ".join(completed.stderr.split())
    assert completed.returncode == 1
    assert rendered.count("Deterministic analysis did not produce") == 1
    assert "optimizer did not converge" in rendered
    assert "maximum evaluations reached" in rendered
    assert "unexpected internal error" not in rendered
    assert TRACEBACK_HEADER not in rendered


def test_cli_bootstrap_reports_contextualized_publication_failure() -> None:
    completed = _run_bootstrap(
        """
        from pathlib import Path
        from chemex.exceptions import ArtifactPublicationError
        cause = OSError(28, "No space left on device")
        error = ArtifactPublicationError(
            "publish the restart state",
            Path("Output/run_info/restart.toml"),
            cause,
        )
        error.outcome_path = Path("Output/run_info/outcome.toml")
        raise error from cause
        """
    )

    rendered = " ".join(completed.stderr.split())
    assert completed.returncode == 1
    assert "could not publish an output artifact" in rendered
    assert "publish the restart state" in rendered
    assert "Output/run_info/restart.toml" in rendered
    assert "No space left on device" in rendered
    assert "Output/run_info/outcome.toml" in rendered
    assert "unexpected internal error" not in rendered
    assert TRACEBACK_HEADER not in rendered


def test_interrupted_publication_failure_keeps_exit_130_and_integrity_details() -> None:
    completed = _run_bootstrap(
        """
        from pathlib import Path
        from chemex.exceptions import ArtifactPublicationError
        cause = OSError(28, "No space left on device")
        error = ArtifactPublicationError(
            "publish the interrupted run outcome",
            Path("Output/run_info/outcome.toml"),
            cause,
        )
        error.terminal = "interrupted"
        error.restart_path = Path("Output/run_info/restart.toml")
        raise error from cause
        """
    )

    rendered = " ".join(completed.stderr.split())
    assert completed.returncode == 130
    assert "Analysis interrupted by user" in rendered
    assert "Interruption finalization failed" in rendered
    assert "publish the interrupted run outcome" in rendered
    assert "Output/run_info/outcome.toml" in rendered
    assert "No space left on device" in rendered
    assert "Output/run_info/restart.toml" in rendered
    assert TRACEBACK_HEADER not in rendered


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
            from chemex.exceptions import ChemExError
            error = ChemExError("SECRET_CANARY")
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
        (
            (
                "from chemex.exceptions import ChemExError; "
                "raise ChemExError('Known ChemEx failure.')"
            ),
            1,
            "Known ChemEx failure.",
        ),
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


def test_debug_reraises_application_import_failure(tmp_path: Path) -> None:
    env = _application_import_failure_env(
        tmp_path,
        f"RuntimeError('{SECRET_CANARY}')",
    )

    completed = _run((sys.executable, "-m", "chemex"), "--debug", env=env)

    combined = completed.stdout + completed.stderr
    assert completed.returncode == 1
    assert TRACEBACK_HEADER in combined
    assert f"RuntimeError: {SECRET_CANARY}" in combined
    assert UNEXPECTED_MESSAGE not in combined


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
    assert completed.stderr.count("Analysis interrupted by user.") == 1
    assert "During: deterministic fit (Method Step STEP2)" in completed.stderr
    assert str(output / "run_info" / "outcome.toml") in completed.stderr
    assert str(output / "run_info" / "restart.toml") in completed.stderr
    assert TRACEBACK_HEADER not in combined
    assert (output / "STEP1" / "Parameters" / "fitted.toml").is_file()
    assert not (output / "STEP2" / "Parameters" / "fitted.toml").exists()
    outcome = tomllib.loads(
        (output / "run_info" / "outcome.toml").read_text(encoding="utf-8")
    )
    assert outcome["status"] == "incomplete"
    assert outcome["terminal"] == "interrupted"
    assert outcome["latest_committed_revision"] == 1

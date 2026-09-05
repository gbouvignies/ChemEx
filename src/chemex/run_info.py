from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import uuid
from argparse import Namespace
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from hashlib import blake2b
from pathlib import Path
from typing import Literal

from chemex.atomic import write_text_atomic
from chemex.exceptions import ArtifactPublicationError, InputFileReadError
from chemex.parameters.parameterization import (
    ParameterRole,
    SealedParameterModel,
    baseline_parameter_role,
)
from chemex.parameters.sealed import parameter_name_from_definition
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.runtime.environment import RuntimeEnvironment
from chemex.runtime.execution import ExecutionSettings

SCHEMA_VERSION = 2
_INPUT_CATEGORIES = ("experiments", "parameters", "methods")
RunOutcomeStatus = Literal["running", "complete", "incomplete"]
FailureStage = Literal[
    "deterministic_fit",
    "restart_publication",
    "simulation",
    "uncertainty",
    "statistics",
    "output",
    "run_info",
    "unknown",
]
StochasticKind = Literal["de", "mc", "bs", "bsn", "mcmc"]


def mark_failure_stage(error: BaseException, stage: FailureStage) -> BaseException:
    """Assign a product failure stage without changing the exception type."""
    if getattr(error, "failure_stage", None) is None:
        error.failure_stage = stage  # ty: ignore[unresolved-attribute]
    return error


def mark_method_step(error: BaseException, step: str) -> BaseException:
    """Assign a Method Step without changing the propagated failure type."""
    if getattr(error, "method_step", None) is None:
        error.method_step = step  # ty: ignore[unresolved-attribute]
    return error


@dataclass(frozen=True, slots=True)
class InputFile:
    category: str
    provided_path: Path
    resolved_path: Path
    content: bytes


@dataclass(frozen=True, slots=True)
class CopiedInputFile:
    provided_path: Path
    resolved_path: Path
    copied_path: Path


@dataclass(slots=True)
class RunInfo:
    """Small invocation-owned publisher for restart, seeds, and lifecycle."""

    path: Path
    parameter_model: SealedParameterModel
    starting_snapshot: AnalysisValuesSnapshot
    restart_revision: int = 0
    stochastic_operations: tuple[tuple[str, StochasticKind, int], ...] = ()

    def publish_restart(self, snapshot: AnalysisValuesSnapshot) -> None:
        """Atomically publish the latest committed continuation state."""
        try:
            if (
                snapshot.occurrence_identity
                != self.starting_snapshot.occurrence_identity
                or snapshot.model_identity != self.starting_snapshot.model_identity
                or snapshot.definitions_identity
                != self.starting_snapshot.definitions_identity
                or snapshot.configuration_identity
                != self.starting_snapshot.configuration_identity
                or snapshot.revision <= self.restart_revision
            ):
                raise ValueError("Restart snapshot does not advance this invocation")
            text = serialize_parameter_file(
                self.parameter_model,
                snapshot,
                state_kind="restart",
            )
            write_text_atomic(self.path / "restart.toml", text)
            self.restart_revision = snapshot.revision
        except OSError as error:
            failure = ArtifactPublicationError(
                "publish the restart state",
                self.path / "restart.toml",
                error,
            )
            mark_failure_stage(failure, "restart_publication")
            raise failure from error
        except (Exception, KeyboardInterrupt) as error:
            mark_failure_stage(error, "restart_publication")
            raise

    def record_stochastic_operation(
        self,
        step: str,
        kind: StochasticKind,
        seed: int,
    ) -> None:
        """Persist one resolved root seed before its operation consumes RNG state."""
        try:
            if not step or kind not in {"de", "mc", "bs", "bsn", "mcmc"}:
                raise ValueError("Invalid stochastic operation record")
            if (
                isinstance(seed, bool)
                or not isinstance(seed, int)
                or not 0 <= seed < 2**64
            ):
                raise ValueError(
                    "Stochastic root seed must be an unsigned 64-bit integer"
                )
            if any(
                recorded_step == step and recorded_kind == kind
                for recorded_step, recorded_kind, _recorded_seed in self.stochastic_operations
            ):
                raise ValueError("Stochastic operation was already recorded")
            record = (step, kind, seed)
            lines = [
                "[[stochastic_operations]]",
                f"step = {_quote_string(step)}",
                f"kind = {_quote_string(kind)}",
                f"seed = {_quote_string(str(seed))}",
                "",
            ]
            run_path = self.path / "run.toml"
            text = run_path.read_text(encoding="utf-8") + "\n".join(lines)
            write_text_atomic(run_path, text)
            self.stochastic_operations += (record,)
        except OSError as error:
            failure = ArtifactPublicationError(
                "record the stochastic operation",
                self.path / "run.toml",
                error,
            )
            mark_failure_stage(failure, "run_info")
            raise failure from error
        except (Exception, KeyboardInterrupt) as error:
            mark_failure_stage(error, "run_info")
            raise

    def write_outcome(
        self,
        status: RunOutcomeStatus,
        snapshot: AnalysisValuesSnapshot,
        *,
        failure: BaseException | None = None,
        failure_stage: FailureStage | None = None,
    ) -> None:
        """Atomically publish the invocation lifecycle and truthful revisions."""
        stage = (
            getattr(failure, "failure_stage", None) if failure is not None else None
        ) or failure_stage
        try:
            write_run_outcome(
                self.path.parent,
                status,
                failure=failure,
                latest_committed_revision=snapshot.revision,
                restart_revision=self.restart_revision,
                failure_stage=stage,
            )
        except OSError as error:
            failure = ArtifactPublicationError(
                f"publish the {status} run outcome",
                self.path / "outcome.toml",
                error,
            )
            mark_failure_stage(failure, "run_info")
            raise failure from error
        except (Exception, KeyboardInterrupt) as error:
            mark_failure_stage(error, "run_info")
            raise


def _quote_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def _format_string_list(values: Sequence[str]) -> str:
    return "[" + ", ".join(_quote_string(value) for value in values) + "]"


def _parameter_sections(
    parameter_model: SealedParameterModel,
) -> Mapping[str, Mapping[str, str]]:
    sections: defaultdict[str, dict[str, str]] = defaultdict(dict)
    for definition in parameter_model.definitions:
        declaration = parameter_model.declarations[definition.param_id]
        if baseline_parameter_role(declaration) is ParameterRole.DERIVED:
            continue
        param_name = parameter_name_from_definition(definition)
        if param_name.spin_system:
            section = param_name.section
            key = str(param_name.spin_system)
        else:
            section = "GLOBAL"
            key = param_name.section_res
        sections[section][key] = definition.param_id
    return sections


def _parameter_values(
    param_id: str,
    parameter_model: SealedParameterModel,
    snapshot: AnalysisValuesSnapshot,
) -> str:
    configuration = parameter_model.configuration[param_id]
    values = (
        snapshot[param_id],
        configuration.lower_bound,
        configuration.upper_bound,
    )
    return "[" + ", ".join(repr(float(value)) for value in values) + "]"


def serialize_parameter_file(
    parameter_model: SealedParameterModel,
    snapshot: AnalysisValuesSnapshot,
    *,
    state_kind: Literal["original", "restart"],
) -> str:
    """Serialize an ordinary parameter file from authoritative central values."""
    lines = (
        [
            "# Original invocation start used by ChemEx.",
            "# This file reproduces the original start; it is not the latest fitted state.",
            "# Each array is [value, minimum, maximum].",
            "",
        ]
        if state_kind == "original"
        else [
            "# Latest committed continuation state reached by this invocation.",
            "# Use this ordinary parameter file with: chemex fit ... -p restart.toml ...",
            "# Each array is [value, minimum, maximum].",
            "",
        ]
    )
    for section, parameters in _parameter_sections(parameter_model).items():
        lines.append(f"[{_quote_string(section)}]")
        lines.extend(
            f"{_quote_string(key)} = {_parameter_values(param_id, parameter_model, snapshot)}"
            for key, param_id in parameters.items()
        )
        lines.append("")
    return "\n".join(lines)


def _resolve_path(path: Path, working_directory: Path) -> Path:
    if path.is_absolute():
        return path.resolve()
    return (working_directory / path).resolve()


def capture_input_files(
    args: Namespace,
    working_directory: Path | None = None,
) -> tuple[InputFile, ...]:
    """Capture exact CLI TOML bytes before any consumer parses them."""
    cwd = (
        Path.cwd().resolve()
        if working_directory is None
        else working_directory.resolve()
    )
    files: list[InputFile] = []
    for category in _INPUT_CATEGORIES:
        paths = getattr(args, category if category != "methods" else "method", None)
        if paths is None:
            continue
        for path in paths:
            resolved_path = _resolve_path(path, cwd)
            try:
                content = resolved_path.read_bytes()
            except OSError as error:
                raise InputFileReadError(resolved_path, error) from error
            files.append(InputFile(category, path, resolved_path, content))
    return tuple(files)


def _copy_name(path: Path, *, collides: bool) -> str:
    if not collides:
        return path.name
    digest = blake2b(str(path).encode(), digest_size=4).hexdigest()
    return f"{path.stem}__{digest}{path.suffix}"


def _copy_inputs(
    input_files: Sequence[InputFile],
    run_info_path: Path,
) -> dict[str, list[CopiedInputFile]]:
    copied_inputs: defaultdict[str, list[CopiedInputFile]] = defaultdict(list)
    files_by_category: defaultdict[str, list[InputFile]] = defaultdict(list)
    for category in _INPUT_CATEGORIES:
        (run_info_path / "inputs" / category).mkdir(parents=True, exist_ok=True)
    for input_file in input_files:
        if input_file.resolved_path not in {
            item.resolved_path for item in files_by_category[input_file.category]
        }:
            files_by_category[input_file.category].append(input_file)
    for category, files in files_by_category.items():
        name_counts = Counter(item.resolved_path.name for item in files)
        destination_directory = run_info_path / "inputs" / category
        used_names: set[str] = set()
        for input_file in files:
            copy_name = _copy_name(
                input_file.resolved_path,
                collides=name_counts[input_file.resolved_path.name] > 1,
            )
            candidate = Path(copy_name)
            index = 2
            while copy_name in used_names:
                copy_name = f"{candidate.stem}__{index}{candidate.suffix}"
                index += 1
            used_names.add(copy_name)
            destination = destination_directory / copy_name
            destination.write_bytes(input_file.content)
            copied_inputs[category].append(
                CopiedInputFile(
                    input_file.provided_path,
                    input_file.resolved_path,
                    destination.relative_to(run_info_path),
                )
            )
    return dict(copied_inputs)


def _find_git_root(start: Path) -> Path | None:
    for directory in (start, *start.parents):
        if (directory / ".git").exists():
            return directory
    return None


def _run_git(repository: Path, *arguments: str) -> str:
    executable = shutil.which("git")
    if executable is None:
        raise FileNotFoundError
    result = subprocess.run(  # noqa: S603
        [executable, "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        text=True,
        timeout=2,
    )
    return result.stdout.strip()


def _git_metadata() -> dict[str, str | bool] | None:
    source_file = Path(__file__).resolve()
    repository = _find_git_root(source_file.parent)
    if repository is None:
        return None
    try:
        source_relative = source_file.relative_to(repository)
        _run_git(
            repository,
            "ls-files",
            "--error-unmatch",
            "--",
            source_relative.as_posix(),
        )
        commit = _run_git(repository, "rev-parse", "HEAD")
        branch = _run_git(repository, "branch", "--show-current")
        status = _run_git(repository, "status", "--porcelain")
    except (OSError, ValueError, subprocess.SubprocessError):
        return None
    metadata: dict[str, str | bool] = {"commit": commit}
    if branch:
        metadata["branch"] = branch
    metadata["working_tree_dirty"] = bool(status)
    return metadata


def _replace_run_info(staging_path: Path, run_info_path: Path) -> None:
    backup_path = run_info_path.with_name(f".run_info-backup-{uuid.uuid4().hex}")
    had_existing = run_info_path.exists() or run_info_path.is_symlink()
    if had_existing:
        run_info_path.replace(backup_path)
    try:
        staging_path.replace(run_info_path)
    except OSError:
        if had_existing:
            backup_path.replace(run_info_path)
        raise
    if not had_existing:
        return
    if backup_path.is_symlink() or backup_path.is_file():
        backup_path.unlink()
    else:
        shutil.rmtree(backup_path)


def _selection_values(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, (list, tuple)):
        return tuple(str(item) for item in value)
    return (str(value),)


def _serialize_run(
    *,
    timestamp: datetime,
    working_directory: Path,
    requested_output_path: Path,
    output_directory: Path,
    model: str,
    include: Sequence[str],
    exclude: Sequence[str],
    execution: ExecutionSettings,
    environment: RuntimeEnvironment,
    argv: Sequence[str],
    copied_inputs: Mapping[str, Sequence[CopiedInputFile]],
    git_metadata: Mapping[str, str | bool] | None,
) -> str:
    lines = [
        "# ChemEx run information.",
        "# Archived inputs are verbatim copies of the invocation TOML files.",
        "",
        f"schema_version = {SCHEMA_VERSION}",
        f"created_at_utc = {_quote_string(timestamp.astimezone(UTC).isoformat())}",
        "",
        "[run]",
        'kind = "fit"',
        f"working_directory = {_quote_string(str(working_directory))}",
        f"requested_output_path = {_quote_string(str(requested_output_path))}",
        f"resolved_output_path = {_quote_string(str(output_directory))}",
        f"model = {_quote_string(model)}",
        "",
        "[command]",
        f"arguments = {_format_string_list(argv)}",
        "",
        "[selection]",
        f"include = {_format_string_list(include)}",
        f"exclude = {_format_string_list(exclude)}",
        "",
        "[execution]",
        f"workers = {execution.workers}",
        "native_threads = "
        + (
            '"automatic"'
            if execution.native_threads is None
            else str(execution.native_threads)
        ),
        "",
        "[software]",
        f"chemex = {_quote_string(environment.chemex_version)}",
        f"python = {_quote_string(environment.python_version)}",
        f"numpy = {_quote_string(environment.numpy_version)}",
        f"scipy = {_quote_string(environment.scipy_version)}",
    ]
    if environment.emcee_version != "unavailable":
        lines.append(f"emcee = {_quote_string(environment.emcee_version)}")
    lines.append("")
    for category in _INPUT_CATEGORIES:
        for copied_input in copied_inputs.get(category, ()):
            lines.extend(
                (
                    f"[[inputs.{category}]]",
                    f"provided_path = {_quote_string(str(copied_input.provided_path))}",
                    f"resolved_path = {_quote_string(str(copied_input.resolved_path))}",
                    f"copied_path = {_quote_string(str(copied_input.copied_path))}",
                    "",
                )
            )
    if git_metadata is not None:
        lines.append("[git]")
        for key, value in git_metadata.items():
            formatted = (
                str(value).lower() if isinstance(value, bool) else _quote_string(value)
            )
            lines.append(f"{key} = {formatted}")
        lines.append("")
    return "\n".join(lines)


def _serialize_outcome(
    status: RunOutcomeStatus,
    failure: BaseException | None = None,
    *,
    latest_committed_revision: int = 0,
    restart_revision: int = 0,
    failure_stage: FailureStage | None = None,
) -> str:
    lines = [
        f"schema_version = {SCHEMA_VERSION}",
        f"status = {_quote_string(status)}",
        f"latest_committed_revision = {latest_committed_revision}",
        f"restart_revision = {restart_revision}",
    ]
    if failure is not None:
        reported_terminal = getattr(failure, "terminal", None)
        terminal = (
            str(reported_terminal)
            if reported_terminal in {"failed", "interrupted"}
            else "interrupted"
            if isinstance(failure, KeyboardInterrupt)
            else "failed"
        )
        message = str(failure).replace("\n", " ")
        lines.extend(
            (
                f"terminal = {_quote_string(terminal)}",
                f"failure_stage = {_quote_string(failure_stage or 'unknown')}",
                f"failure_type = {_quote_string(type(failure).__name__)}",
                f"failure_message = {_quote_string(message)}",
            )
        )
    return "\n".join(lines) + "\n"


def write_run_outcome(
    output_directory: Path,
    status: RunOutcomeStatus,
    *,
    failure: BaseException | None = None,
    latest_committed_revision: int = 0,
    restart_revision: int = 0,
    failure_stage: FailureStage | None = None,
) -> None:
    """Atomically replace the lifecycle marker for the current fit invocation."""
    outcome_path = Path(output_directory) / "run_info" / "outcome.toml"
    write_text_atomic(
        outcome_path,
        _serialize_outcome(
            status,
            failure,
            latest_committed_revision=latest_committed_revision,
            restart_revision=restart_revision,
            failure_stage=failure_stage,
        ),
    )


def write_run_info(
    args: Namespace,
    *,
    parameter_model: SealedParameterModel,
    starting_values: AnalysisValuesSnapshot,
    execution: ExecutionSettings | None = None,
    input_files: Sequence[InputFile] | None = None,
    argv: Sequence[str] | None = None,
    working_directory: Path | None = None,
    timestamp: datetime | None = None,
) -> RunInfo:
    """Publish the invocation record, contextualizing only concrete I/O failures."""
    output_directory = _resolve_path(
        args.output,
        Path.cwd().resolve()
        if working_directory is None
        else working_directory.resolve(),
    )
    try:
        return _write_run_info(
            args,
            parameter_model=parameter_model,
            starting_values=starting_values,
            execution=execution,
            input_files=input_files,
            argv=argv,
            working_directory=working_directory,
            timestamp=timestamp,
        )
    except InputFileReadError:
        raise
    except OSError as error:
        failure = ArtifactPublicationError(
            "publish run information",
            output_directory / "run_info",
            error,
        )
        mark_failure_stage(failure, "run_info")
        raise failure from error


def _write_run_info(
    args: Namespace,
    *,
    parameter_model: SealedParameterModel,
    starting_values: AnalysisValuesSnapshot,
    execution: ExecutionSettings | None = None,
    input_files: Sequence[InputFile] | None = None,
    argv: Sequence[str] | None = None,
    working_directory: Path | None = None,
    timestamp: datetime | None = None,
) -> RunInfo:
    """Stage and atomically publish the lightweight invocation record."""
    cwd = (
        Path.cwd().resolve()
        if working_directory is None
        else working_directory.resolve()
    )
    output_directory = _resolve_path(args.output, cwd)
    run_info_path = output_directory / "run_info"
    captured_inputs = (
        capture_input_files(args, cwd) if input_files is None else input_files
    )
    git_metadata = _git_metadata()
    resolved_execution = (
        ExecutionSettings.from_counts(
            workers=getattr(args, "workers", 1),
            native_threads=getattr(args, "native_threads", "auto"),
        )
        if execution is None
        else execution
    )
    environment = RuntimeEnvironment.from_current_process()
    output_directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=output_directory,
        prefix=".run_info-",
    ) as staging:
        staging_path = Path(staging)
        copied_inputs = _copy_inputs(captured_inputs, staging_path)
        (staging_path / "parameters_used.toml").write_text(
            serialize_parameter_file(
                parameter_model,
                starting_values,
                state_kind="original",
            ),
            encoding="utf-8",
        )
        (staging_path / "run.toml").write_text(
            _serialize_run(
                timestamp=datetime.now(UTC) if timestamp is None else timestamp,
                working_directory=cwd,
                requested_output_path=args.output,
                output_directory=output_directory,
                model=str(getattr(args, "model", parameter_model.model_name)),
                include=_selection_values(getattr(args, "include", None)),
                exclude=_selection_values(getattr(args, "exclude", None)),
                execution=resolved_execution,
                environment=environment,
                argv=tuple(sys.argv if argv is None else argv),
                copied_inputs=copied_inputs,
                git_metadata=git_metadata,
            ),
            encoding="utf-8",
        )
        (staging_path / "outcome.toml").write_text(
            _serialize_outcome("running"),
            encoding="utf-8",
        )
        _replace_run_info(staging_path, run_info_path)
    return RunInfo(run_info_path, parameter_model, starting_values)

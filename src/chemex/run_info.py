from __future__ import annotations

import json
import platform
import shutil
import subprocess
import sys
import tempfile
import tomllib
import uuid
from argparse import Namespace
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from enum import StrEnum
from hashlib import blake2b, sha256
from pathlib import Path
from typing import Literal, cast

from chemex import __version__
from chemex.atomic import publish_directory_noreplace, write_text_atomic
from chemex.containers.experiments import Experiments
from chemex.native_provenance import (
    ArtifactReference,
    ArtifactRole,
    CommittedRestartRecord,
    NativeProvenanceError,
    PublishedStepReference,
    WorkflowProvenance,
    published_step_reference_identity,
    serialize_independent_parameters,
    validate_native_step_manifest_bytes,
)
from chemex.parameters.database import ParameterStore
from chemex.parameters.parameterization import SealedParameterModel
from chemex.parameters.setting import ParamSetting
from chemex.parameters.values import AnalysisValuesSnapshot

SCHEMA_VERSION = 1
_INPUT_CATEGORIES = ("experiments", "parameters", "methods")
RunOutcomeStatus = Literal["running", "complete", "incomplete"]


@dataclass(frozen=True, slots=True)
class InputFile:
    category: str
    provided_path: Path
    resolved_path: Path


@dataclass(frozen=True, slots=True)
class CopiedInputFile:
    provided_path: Path
    resolved_path: Path
    copied_path: Path
    sha256: str
    identity: str | None = None


@dataclass(frozen=True, slots=True)
class CapturedNativeInput:
    """Immutable occurrence-start bytes plus non-authoritative source metadata."""

    category: str
    provided_path: Path
    resolved_path: Path
    archive_name: str
    content: bytes = field(repr=False)
    sha256: str = field(init=False)
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        digest = sha256(self.content).hexdigest()
        if (
            self.category not in _INPUT_CATEGORIES
            or not self.archive_name
            or Path(self.archive_name).name != self.archive_name
            or self.archive_name in {".", ".."}
        ):
            raise ValueError("Captured native input has an unsafe archive identity")
        object.__setattr__(self, "sha256", digest)
        object.__setattr__(
            self,
            "identity",
            sha256(
                json.dumps(
                    (
                        "native-input-v2",
                        self.category,
                        str(self.provided_path),
                        str(self.resolved_path),
                        self.archive_name,
                        digest,
                        len(self.content),
                    ),
                    ensure_ascii=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        )


@dataclass(frozen=True, slots=True)
class CapturedNativeInputs:
    """Canonical immutable input bundle captured before native execution begins."""

    members: tuple[CapturedNativeInput, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        members = tuple(self.members)
        keys = tuple((item.category, item.archive_name) for item in members)
        if not members or len(set(keys)) != len(keys):
            raise ValueError("Captured native inputs must have unique archive members")
        object.__setattr__(self, "members", members)
        object.__setattr__(
            self,
            "identity",
            sha256(
                json.dumps(
                    (
                        "native-input-bundle-v2",
                        tuple(item.identity for item in members),
                    ),
                    ensure_ascii=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        )


class RunInformationKind(StrEnum):
    """Supported meanings of historical and native run-information records."""

    HISTORICAL_V1 = "historical_schema_v1"
    NATIVE_V2 = "native_schema_v2"


@dataclass(frozen=True, slots=True)
class NativeRunInformation:
    """Complete invocation state and links to atomic native step outcomes."""

    invocation_identity: str
    inputs: CapturedNativeInputs
    parameter_model: SealedParameterModel
    starting_snapshot: AnalysisValuesSnapshot
    steps: tuple[PublishedStepReference, ...]
    execution_occurrence_identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not self.invocation_identity.strip():
            raise ValueError("Native invocation identity cannot be empty")
        if not self.steps:
            raise ValueError("Native run information requires at least one step")
        snapshot = self.starting_snapshot
        if (
            snapshot.revision != 0
            or snapshot.model_identity != self.parameter_model.model_identity
            or snapshot.definitions_identity
            != self.parameter_model.definitions.identity
            or snapshot.configuration_identity
            != self.parameter_model.configuration.identity
        ):
            raise ValueError(
                "Native run starting state must be revision zero of its parameter model"
            )
        environments = {step.provenance.environment.identity for step in self.steps}
        if len(environments) != 1:
            raise ValueError("Native steps must share one resolved environment")
        object.__setattr__(
            self,
            "execution_occurrence_identity",
            _execution_occurrence_identity(
                self.inputs.identity,
                snapshot.occurrence_identity,
                snapshot.revision,
                tuple(step.identity for step in self.steps),
            ),
        )


@dataclass(frozen=True, slots=True)
class HistoricalNativeRunInformation:
    """Recursively validated durable schema-v2 provenance without live authority."""

    invocation_identity: str
    input_bundle_identity: str
    record_identity: str
    published_step_identities: tuple[str, ...]
    execution_occurrence_identity: str


def classify_run_information(path: Path) -> RunInformationKind:
    """Classify schema-v1 as historical and schema-v2 as native product data."""
    record = tomllib.loads(Path(path).read_text(encoding="utf-8"))
    schema_version = record.get("schema_version")
    if schema_version == 1:
        return RunInformationKind.HISTORICAL_V1
    if (
        schema_version == 2
        and record.get("schema_kind") == "native_product_run_information"
    ):
        return RunInformationKind.NATIVE_V2
    raise ValueError("Unsupported run-information schema")


def _native_run_record_identity(record: Mapping[str, object]) -> str:
    payload = dict(record)
    payload.pop("record_identity", None)
    return sha256(
        json.dumps(
            {
                "kind": "native-run-information-v2",
                "schema_version": 1,
                "record": payload,
            },
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode()
    ).hexdigest()


def _execution_occurrence_identity(
    input_bundle_identity: str,
    starting_occurrence_identity: str,
    starting_revision: int,
    published_step_identities: tuple[str, ...],
) -> str:
    return sha256(
        json.dumps(
            (
                "native-execution-occurrence-v2",
                input_bundle_identity,
                starting_occurrence_identity,
                starting_revision,
                published_step_identities,
            ),
            ensure_ascii=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


def _quote_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def _quote_key(value: str) -> str:
    return _quote_string(value)


def _format_string_list(values: Sequence[str]) -> str:
    return "[" + ", ".join(_quote_string(value) for value in values) + "]"


def _format_float(value: float | None) -> str:
    return "nan" if value is None else repr(float(value))


def _parameter_values(parameter: ParamSetting) -> str:
    values = [
        _format_float(parameter.value),
        _format_float(parameter.min),
        _format_float(parameter.max),
    ]
    if parameter.brute_step is not None:
        values.append(_format_float(parameter.brute_step))
    return "[" + ", ".join(values) + "]"


def _parameter_sections(
    experiments: Experiments,
    parameter_store: ParameterStore,
) -> Mapping[str, Mapping[str, ParamSetting]]:
    sections: defaultdict[str, dict[str, ParamSetting]] = defaultdict(dict)
    parameters = parameter_store.get_parameters(experiments.param_ids)

    for parameter in sorted(parameters.values(), key=lambda item: item.param_name):
        if parameter.expr:
            continue
        param_name = parameter.param_name
        if param_name.spin_system:
            section = param_name.section
            key = str(param_name.spin_system)
        else:
            section = "GLOBAL"
            key = param_name.section_res
        sections[section][key] = parameter

    return sections


def _serialize_parameters(
    experiments: Experiments,
    parameter_store: ParameterStore,
) -> str:
    lines = [
        "# Starting independent parameters used by ChemEx for this fit.",
        "# Parameters defined by expressions are omitted here because they are",
        "# reconstructed from the independent parameters and model definition.",
        "# This file is intended to be usable as a starting parameter file for",
        "# re-running the fit.",
        "# Each array is [value, minimum, maximum, optional brute-force step].",
        "",
    ]

    for section, parameters in _parameter_sections(
        experiments,
        parameter_store,
    ).items():
        lines.append(f"[{_quote_key(section)}]")
        lines.extend(
            f"{_quote_key(key)} = {_parameter_values(parameter)}"
            for key, parameter in parameters.items()
        )
        lines.append("")

    return "\n".join(lines)


def _resolve_path(path: Path, working_directory: Path) -> Path:
    if path.is_absolute():
        return path.resolve()
    return (working_directory / path).resolve()


def _collect_input_files(args: Namespace, working_directory: Path) -> list[InputFile]:
    files: list[InputFile] = []
    for category in _INPUT_CATEGORIES:
        paths = getattr(args, category if category != "methods" else "method", None)
        if paths is None:
            continue
        files.extend(
            InputFile(
                category,
                path,
                _resolve_path(path, working_directory),
            )
            for path in paths
        )
    return files


def _copy_name(path: Path, *, collides: bool) -> str:
    if not collides:
        return path.name
    digest = blake2b(str(path).encode(), digest_size=4).hexdigest()
    return f"{path.stem}__{digest}{path.suffix}"


def capture_native_inputs(
    args: Namespace,
    *,
    working_directory: Path | None = None,
) -> CapturedNativeInputs:
    """Capture native occurrence input bytes exactly once before execution."""
    cwd = (
        Path.cwd().resolve()
        if working_directory is None
        else working_directory.resolve()
    )
    input_files = _collect_input_files(args, cwd)
    unique: list[tuple[InputFile, bytes]] = []
    seen: set[tuple[str, Path]] = set()
    for item in input_files:
        key = (item.category, item.resolved_path)
        if key in seen:
            continue
        seen.add(key)
        unique.append((item, item.resolved_path.read_bytes()))
    counts = Counter((item.category, item.resolved_path.name) for item, _ in unique)
    members: list[CapturedNativeInput] = []
    for item, content in unique:
        basename = item.resolved_path.name
        if not basename or basename in {".", ".."}:
            basename = "input"
        collides = counts[(item.category, item.resolved_path.name)] > 1
        if collides:
            digest = blake2b(
                (str(item.resolved_path) + sha256(content).hexdigest()).encode(),
                digest_size=4,
            ).hexdigest()
            source = Path(basename)
            basename = f"{source.stem}__{digest}{source.suffix}"
        members.append(
            CapturedNativeInput(
                item.category,
                item.provided_path,
                item.resolved_path,
                basename,
                bytes(content),
            )
        )
    return CapturedNativeInputs(tuple(members))


def _write_captured_native_inputs(
    inputs: CapturedNativeInputs,
    run_info_path: Path,
) -> dict[str, list[CopiedInputFile]]:
    copied: defaultdict[str, list[CopiedInputFile]] = defaultdict(list)
    for member in inputs.members:
        destination = run_info_path / "inputs" / member.category / member.archive_name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(member.content)
        copied[member.category].append(
            CopiedInputFile(
                member.provided_path,
                member.resolved_path,
                destination.relative_to(run_info_path),
                member.sha256,
                member.identity,
            )
        )
    return dict(copied)


def _copy_inputs(
    input_files: Sequence[InputFile],
    run_info_path: Path,
) -> dict[str, list[CopiedInputFile]]:
    copied_inputs: defaultdict[str, list[CopiedInputFile]] = defaultdict(list)
    files_by_category: defaultdict[str, list[InputFile]] = defaultdict(list)

    for input_file in input_files:
        resolved_paths = {
            item.resolved_path for item in files_by_category[input_file.category]
        }
        if input_file.resolved_path not in resolved_paths:
            files_by_category[input_file.category].append(input_file)

    for category, files in files_by_category.items():
        name_counts = Counter(item.resolved_path.name for item in files)
        destination_directory = run_info_path / "inputs" / category
        destination_directory.mkdir(parents=True, exist_ok=True)
        used_names: set[str] = set()

        for input_file in files:
            resolved_path = input_file.resolved_path
            copy_name = _copy_name(
                resolved_path,
                collides=name_counts[resolved_path.name] > 1,
            )
            candidate = Path(copy_name)
            index = 2
            while copy_name in used_names:
                copy_name = f"{candidate.stem}__{index}{candidate.suffix}"
                index += 1
            used_names.add(copy_name)
            destination = destination_directory / copy_name
            shutil.copy2(resolved_path, destination)
            copied_inputs[category].append(
                CopiedInputFile(
                    input_file.provided_path,
                    resolved_path,
                    destination.relative_to(run_info_path),
                    sha256(destination.read_bytes()).hexdigest(),
                ),
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


def _serialize_run(
    *,
    timestamp: datetime,
    working_directory: Path,
    output_directory: Path,
    argv: Sequence[str],
    copied_inputs: Mapping[str, Sequence[CopiedInputFile]],
    git_metadata: Mapping[str, str | bool] | None,
) -> str:
    lines = [
        "# ChemEx run information.",
        "# Input records distinguish the path provided by the user, its resolved",
        "# absolute path, and the archived copy relative to this run_info directory.",
        "",
        f"schema_version = {SCHEMA_VERSION}",
        f"created_at_utc = {_quote_string(timestamp.astimezone(UTC).isoformat())}",
        "",
        "[run]",
        'kind = "fit"',
        f"working_directory = {_quote_string(str(working_directory))}",
        f"output_directory = {_quote_string(str(output_directory))}",
        "",
        "[chemex]",
        f"version = {_quote_string(__version__)}",
        "",
        "[python]",
        f"version = {_quote_string(platform.python_version())}",
        f"platform = {_quote_string(platform.platform())}",
        "",
        "[command]",
        "# Exact process arguments, including the executable.",
        f"arguments = {_format_string_list(argv)}",
        "",
    ]

    for category, files in copied_inputs.items():
        for copied_input in files:
            lines.extend(
                (
                    f"[[inputs.{category}]]",
                    f"provided_path = {_quote_string(str(copied_input.provided_path))}",
                    f"resolved_path = {_quote_string(str(copied_input.resolved_path))}",
                    f"copied_path = {_quote_string(str(copied_input.copied_path))}",
                    "",
                ),
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
) -> str:
    lines = [
        f"schema_version = {SCHEMA_VERSION}",
        f"status = {_quote_string(status)}",
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
) -> None:
    """Atomically replace the lifecycle marker for the current fit invocation."""
    outcome_path = Path(output_directory) / "run_info" / "outcome.toml"
    write_text_atomic(outcome_path, _serialize_outcome(status, failure))


def write_run_info(
    args: Namespace,
    experiments: Experiments,
    *,
    argv: Sequence[str] | None = None,
    working_directory: Path | None = None,
    timestamp: datetime | None = None,
) -> None:
    """Write a lightweight description of a fit to its output directory."""
    cwd = (
        Path.cwd().resolve()
        if working_directory is None
        else working_directory.resolve()
    )
    output_directory = _resolve_path(args.output, cwd)
    run_info_path = output_directory / "run_info"
    input_files = _collect_input_files(args, cwd)
    git_metadata = _git_metadata()

    output_directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=output_directory,
        prefix=".run_info-",
    ) as staging_directory:
        staging_path = Path(staging_directory)
        copied_inputs = _copy_inputs(input_files, staging_path)
        parameters_text = _serialize_parameters(
            experiments,
            experiments.parameter_store,
        )
        (staging_path / "parameters_used.toml").write_text(
            parameters_text,
            encoding="utf-8",
        )

        run_text = _serialize_run(
            timestamp=datetime.now(UTC) if timestamp is None else timestamp,
            working_directory=cwd,
            output_directory=output_directory,
            argv=tuple(sys.argv if argv is None else argv),
            copied_inputs=copied_inputs,
            git_metadata=git_metadata,
        )
        (staging_path / "run.toml").write_text(run_text, encoding="utf-8")
        (staging_path / "outcome.toml").write_text(
            _serialize_outcome("running"),
            encoding="utf-8",
        )

        _replace_run_info(staging_path, run_info_path)


def _native_independent_ids(run: NativeRunInformation) -> tuple[str, ...]:
    requested = {param_id for step in run.steps for param_id in step.independent_ids}
    if not requested:
        raise ValueError("Native run information has no independent parameter scope")
    ordered = tuple(
        definition.param_id
        for definition in run.parameter_model.definitions
        if definition.param_id in requested
    )
    if set(ordered) != requested:
        raise ValueError("Native step independent scope is not in the parameter model")
    return ordered


def _relative_step_path(step: PublishedStepReference, output: Path) -> Path:
    try:
        return step.path.relative_to(output)
    except ValueError as error:
        raise ValueError("Native step publication is outside the run output") from error


def _validated_workflow_records(  # noqa: C901 - closed live/durable boundary
    run: NativeRunInformation,
    output: Path,
) -> dict[str, object]:
    workflows: list[dict[str, object]] = []
    current_occurrence = run.starting_snapshot.occurrence_identity
    current_revision = run.starting_snapshot.revision
    for step in run.steps:
        step.require_exact_live_publication()
        relative_step = _relative_step_path(step, output)
        manifest_path = step.path / "fit-manifest.toml"
        manifest_bytes = manifest_path.read_bytes()
        if sha256(manifest_bytes).hexdigest() != step.manifest_sha256:
            raise ValueError("Native step manifest changed before run publication")
        manifest = validate_native_step_manifest_bytes(manifest_bytes)
        manifest_workflow = manifest.get("workflow")
        manifest_method = manifest.get("method")
        if not isinstance(manifest_workflow, dict) or not isinstance(
            manifest_method, dict
        ):
            raise TypeError("Native step manifest provenance is malformed")
        if manifest.get("manifest_identity") != step.manifest_identity:
            raise ValueError(
                "Native step manifest identity changed before run publication"
            )
        if (
            manifest.get("publication_occurrence_identity")
            != step.publication_occurrence_identity
            or manifest.get("lifecycle") != step.lifecycle
            or manifest.get("authority") != step.authority
            or manifest_workflow.get("identity") != step.provenance.workflow_identity
            or manifest_workflow.get("provenance_identity") != step.provenance.identity
            or manifest_method.get("identity") != step.provenance.method_identity
            or manifest_method.get("parameterization_identity")
            != step.provenance.parameterization_identity
            or manifest_method.get("evaluation_plan_identity")
            != step.provenance.evaluation_plan_identity
        ):
            raise ValueError("Native step reference contradicts its manifest")
        manifest_artifacts = manifest.get("artifacts")
        expected_artifacts = {
            artifact.path: {
                "role": artifact.role,
                "sha256": artifact.sha256,
            }
            for artifact in step.artifacts
        }
        if manifest_artifacts != expected_artifacts:
            raise ValueError("Native step artifact references contradict its manifest")
        actual_members = {
            str(item.relative_to(step.path))
            for item in step.path.rglob("*")
            if item.is_file()
        }
        if actual_members != {"fit-manifest.toml", *expected_artifacts}:
            raise ValueError("Native step artifact catalogue is incomplete")
        artifact_records: list[dict[str, str]] = []
        for artifact in step.artifacts:
            artifact_path = step.path / artifact.path
            if sha256(artifact_path.read_bytes()).hexdigest() != artifact.sha256:
                raise ValueError("Native step artifact changed before run publication")
            artifact_records.append(
                {
                    "path": str(relative_step / artifact.path),
                    "role": artifact.role,
                    "sha256": artifact.sha256,
                }
            )
        if step.lifecycle == "committed":
            restart_record = CommittedRestartRecord.from_record(
                json.loads(
                    (step.path / "Parameters" / "restart-provenance.json").read_text(
                        encoding="utf-8"
                    )
                )
            )
            if (
                manifest.get("starting_occurrence_identity") != current_occurrence
                or manifest.get("starting_revision") != current_revision
                or restart_record.starting_occurrence_identity != current_occurrence
                or restart_record.starting_revision != current_revision
                or restart_record.parameter_model_identity
                != run.parameter_model.identity
                or restart_record.configuration_identity
                != run.parameter_model.configuration.identity
                or restart_record.workflow_identity != step.provenance.workflow_identity
                or tuple(item[0] for item in restart_record.independent_items)
                != step.independent_ids
            ):
                raise NativeProvenanceError(
                    "Committed restart belongs to another run occurrence"
                )
            current_occurrence = restart_record.committed_occurrence_identity
            current_revision = restart_record.committed_revision
        record = step.provenance.to_record()
        record["published_step_identity"] = step.identity
        record["publication_occurrence_identity"] = step.publication_occurrence_identity
        record["independent_ids"] = list(step.independent_ids)
        record["outcome"] = {
            "lifecycle": step.lifecycle,
            "authority": step.authority,
        }
        record["manifest"] = {
            "path": str(relative_step / "fit-manifest.toml"),
            "identity": step.manifest_identity,
            "sha256": step.manifest_sha256,
        }
        record["artifacts"] = artifact_records
        workflows.append(record)
    return {
        "schema_version": 1,
        "invocation_identity": run.invocation_identity,
        "execution_occurrence_identity": run.execution_occurrence_identity,
        "workflows": workflows,
    }


def _serialize_native_run(
    *,
    run: NativeRunInformation,
    timestamp: datetime,
    working_directory: Path,
    output_directory: Path,
    argv: Sequence[str],
    copied_inputs: Mapping[str, Sequence[CopiedInputFile]],
    git_metadata: Mapping[str, str | bool] | None,
    starting_sha256: str,
    normalized_methods: Sequence[tuple[str, str, str, str]],
    workflow_records_sha256: str,
    record_identity: str,
) -> str:
    environment = run.steps[0].provenance.environment
    lines = [
        "# Native ChemEx product-facing run information.",
        "# schema_version = 1 records remain historical and are never rewritten.",
        "",
        "schema_version = 2",
        'schema_kind = "native_product_run_information"',
        f"invocation_identity = {_quote_string(run.invocation_identity)}",
        f"input_bundle_identity = {_quote_string(run.inputs.identity)}",
        f"execution_occurrence_identity = {_quote_string(run.execution_occurrence_identity)}",
        f"record_identity = {_quote_string(record_identity)}",
        f"created_at_utc = {_quote_string(timestamp.astimezone(UTC).isoformat())}",
        "",
        "[run]",
        'kind = "fit"',
        f"working_directory = {_quote_string(str(working_directory))}",
        f"output_directory = {_quote_string(str(output_directory))}",
        "",
        "[command]",
        f"arguments = {_format_string_list(argv)}",
        "",
        "[starting_state]",
        'kind = "starting_independent_state"',
        'path = "parameters_used.toml"',
        f"sha256 = {_quote_string(starting_sha256)}",
        f"occurrence_identity = {_quote_string(run.starting_snapshot.occurrence_identity)}",
        f"revision = {run.starting_snapshot.revision}",
        "",
        "[workflow_records]",
        'path = "workflows.json"',
        f"sha256 = {_quote_string(workflow_records_sha256)}",
        "",
        "[environment]",
        f"identity = {_quote_string(environment.identity)}",
        f"chemex_version = {_quote_string(environment.chemex_version)}",
        f"python_version = {_quote_string(environment.python_version)}",
        f"python_implementation = {_quote_string(environment.python_implementation)}",
        f"platform = {_quote_string(environment.platform)}",
        f"numpy_version = {_quote_string(environment.numpy_version)}",
        f"scipy_version = {_quote_string(environment.scipy_version)}",
        f"emcee_version = {_quote_string(environment.emcee_version)}",
        "",
    ]
    for workflow_identity, method_identity, path, digest in normalized_methods:
        lines.extend(
            (
                "[[normalized_methods]]",
                f"workflow_identity = {_quote_string(workflow_identity)}",
                f"method_identity = {_quote_string(method_identity)}",
                f"path = {_quote_string(path)}",
                f"sha256 = {_quote_string(digest)}",
                "",
            )
        )
    for category, files in copied_inputs.items():
        for copied_input in files:
            lines.extend(
                (
                    f"[[inputs.{category}]]",
                    f"provided_path = {_quote_string(str(copied_input.provided_path))}",
                    f"resolved_path = {_quote_string(str(copied_input.resolved_path))}",
                    f"copied_path = {_quote_string(str(copied_input.copied_path))}",
                    f"sha256 = {_quote_string(copied_input.sha256)}",
                    f"identity = {_quote_string(copied_input.identity or '')}",
                    "",
                )
            )
    for step in run.steps:
        relative_step = _relative_step_path(step, output_directory)
        lines.extend(
            (
                "[[steps]]",
                f"path = {_quote_string(str(relative_step))}",
                f"workflow_identity = {_quote_string(step.provenance.workflow_identity)}",
                f"execution_identity = {_quote_string(step.provenance.execution_identity)}",
                f"method_identity = {_quote_string(step.provenance.method_identity)}",
                f"provenance_identity = {_quote_string(step.provenance.identity)}",
                f"published_step_identity = {_quote_string(step.identity)}",
                f"publication_occurrence_identity = {_quote_string(step.publication_occurrence_identity)}",
                f"lifecycle = {_quote_string(step.lifecycle)}",
                f"authority = {_quote_string(step.authority)}",
                f"manifest_path = {_quote_string(str(relative_step / 'fit-manifest.toml'))}",
                f"manifest_identity = {_quote_string(step.manifest_identity)}",
                f"manifest_sha256 = {_quote_string(step.manifest_sha256)}",
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


def _safe_link(root: Path, value: object) -> Path:
    text = str(value)
    relative = Path(text)
    if (
        relative.is_absolute()
        or ".." in relative.parts
        or not relative.parts
        or str(relative) != text
    ):
        raise NativeProvenanceError("Native run information contains an unsafe link")
    return root / relative


def read_native_run_information(  # noqa: C901 - closed durable schema
    path: Path,
) -> HistoricalNativeRunInformation:
    """Recursively validate durable schema-v2 provenance without live authority."""
    run_path = Path(path)
    if run_path.is_dir():
        run_path = run_path / "run.toml"
    root = run_path.parent
    try:
        record = tomllib.loads(run_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise NativeProvenanceError("Native run information is unreadable") from error
    required = {
        "schema_version",
        "schema_kind",
        "invocation_identity",
        "input_bundle_identity",
        "execution_occurrence_identity",
        "record_identity",
        "created_at_utc",
        "run",
        "command",
        "starting_state",
        "normalized_methods",
        "workflow_records",
        "environment",
        "inputs",
        "steps",
    }
    if (
        not required <= set(record)
        or set(record) - required - {"git"}
        or record.get("schema_version") != 2
        or record.get("schema_kind") != "native_product_run_information"
        or record.get("record_identity") != _native_run_record_identity(record)
    ):
        raise NativeProvenanceError("Native run record identity or fields are invalid")
    run_record = record.get("run")
    command = record.get("command")
    run_environment = record.get("environment")
    if (
        not isinstance(record.get("created_at_utc"), str)
        or not isinstance(record.get("invocation_identity"), str)
        or not str(record["invocation_identity"]).strip()
        or not isinstance(run_record, dict)
        or set(run_record) != {"kind", "working_directory", "output_directory"}
        or run_record.get("kind") != "fit"
        or not isinstance(command, dict)
        or set(command) != {"arguments"}
        or not isinstance(command.get("arguments"), list)
        or any(not isinstance(item, str) for item in command["arguments"])
        or not isinstance(run_environment, dict)
        or set(run_environment)
        != {
            "identity",
            "chemex_version",
            "python_version",
            "python_implementation",
            "platform",
            "numpy_version",
            "scipy_version",
            "emcee_version",
        }
    ):
        raise NativeProvenanceError("Native run metadata fields are malformed")
    git = record.get("git")
    if git is not None and (
        not isinstance(git, dict)
        or not {"commit", "working_tree_dirty"} <= set(git)
        or set(git) - {"commit", "branch", "working_tree_dirty"}
        or not isinstance(git.get("commit"), str)
        or not isinstance(git.get("working_tree_dirty"), bool)
        or ("branch" in git and not isinstance(git["branch"], str))
    ):
        raise NativeProvenanceError("Native run git metadata is malformed")
    starting = record.get("starting_state")
    normalized = record.get("normalized_methods")
    workflow_link = record.get("workflow_records")
    if (
        not isinstance(starting, dict)
        or set(starting)
        != {"kind", "path", "sha256", "occurrence_identity", "revision"}
        or starting.get("kind") != "starting_independent_state"
        or not isinstance(normalized, list)
        or not normalized
        or any(
            not isinstance(item, dict)
            or set(item) != {"workflow_identity", "method_identity", "path", "sha256"}
            for item in normalized
        )
        or not isinstance(workflow_link, dict)
        or set(workflow_link) != {"path", "sha256"}
    ):
        raise NativeProvenanceError("Native run child links are malformed")
    normalized_links = cast(list[dict[str, object]], normalized)
    for link in (starting, workflow_link, *normalized_links):
        member = _safe_link(root, link["path"])
        if sha256(member.read_bytes()).hexdigest() != link["sha256"]:
            raise NativeProvenanceError("Native run child hash does not match")
    if len({str(link["path"]) for link in normalized_links}) != len(normalized_links):
        raise NativeProvenanceError("Normalized method links are duplicated")

    inputs = record.get("inputs")
    if not isinstance(inputs, dict) or set(inputs) - set(_INPUT_CATEGORIES):
        raise NativeProvenanceError("Native input archive is malformed")
    captured: list[CapturedNativeInput] = []
    expected_root_members = {
        "run.toml",
        str(starting["path"]),
        str(workflow_link["path"]),
        *(str(link["path"]) for link in normalized_links),
    }
    for category in _INPUT_CATEGORIES:
        entries = inputs.get(category, [])
        if not isinstance(entries, list):
            raise NativeProvenanceError("Native input archive is malformed")
        for item in entries:
            if not isinstance(item, dict) or set(item) != {
                "provided_path",
                "resolved_path",
                "copied_path",
                "sha256",
                "identity",
            }:
                raise NativeProvenanceError("Native input member is malformed")
            member_path = _safe_link(root, item["copied_path"])
            content = member_path.read_bytes()
            archived = CapturedNativeInput(
                category,
                Path(str(item["provided_path"])),
                Path(str(item["resolved_path"])),
                member_path.name,
                content,
            )
            if (
                item["copied_path"] != f"inputs/{category}/{archived.archive_name}"
                or item["sha256"] != archived.sha256
                or item["identity"] != archived.identity
            ):
                raise NativeProvenanceError("Native input member identity is invalid")
            captured.append(archived)
            expected_root_members.add(str(item["copied_path"]))
    bundle = CapturedNativeInputs(tuple(captured))
    if record.get("input_bundle_identity") != bundle.identity:
        raise NativeProvenanceError("Native input bundle identity is invalid")

    try:
        workflows = json.loads(_safe_link(root, workflow_link["path"]).read_text())
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise NativeProvenanceError("Native workflow records are malformed") from error
    if (
        not isinstance(workflows, dict)
        or set(workflows)
        != {
            "schema_version",
            "invocation_identity",
            "execution_occurrence_identity",
            "workflows",
        }
        or workflows.get("schema_version") != 1
        or workflows.get("invocation_identity") != record.get("invocation_identity")
        or not isinstance(workflows.get("workflows"), list)
    ):
        raise NativeProvenanceError("Native workflow record envelope is invalid")
    steps = record.get("steps")
    workflow_items = workflows["workflows"]
    if (
        not isinstance(steps, list)
        or len(steps) != len(workflow_items)
        or len(steps) != len(normalized_links)
        or not steps
    ):
        raise NativeProvenanceError("Native run step links are incomplete")
    published_ids: list[str] = []
    output = root.parent
    current_occurrence = str(starting["occurrence_identity"])
    current_revision = int(starting["revision"])
    for step, workflow, normalized_link in zip(
        steps, workflow_items, normalized_links, strict=True
    ):
        if (
            not isinstance(step, dict)
            or set(step)
            != {
                "path",
                "workflow_identity",
                "execution_identity",
                "method_identity",
                "provenance_identity",
                "published_step_identity",
                "publication_occurrence_identity",
                "lifecycle",
                "authority",
                "manifest_path",
                "manifest_identity",
                "manifest_sha256",
            }
            or not isinstance(workflow, dict)
        ):
            raise NativeProvenanceError("Native run step record is malformed")
        step = cast(dict[str, object], step)
        workflow = cast(dict[str, object], workflow)
        step_path = _safe_link(output, step["path"])
        if step["manifest_path"] != f"{step['path']}/fit-manifest.toml":
            raise NativeProvenanceError("Native step manifest link is not canonical")
        manifest_path = _safe_link(output, step["manifest_path"])
        manifest_bytes = manifest_path.read_bytes()
        manifest = validate_native_step_manifest_bytes(manifest_bytes)
        provenance = WorkflowProvenance.from_manifest_record(manifest)
        artifacts_record = workflow.get("artifacts")
        independent_raw = workflow.get("independent_ids")
        workflow_keys = {
            "identity",
            "execution_identity",
            "workflow_type",
            "grouping_topology",
            "provenance_identity",
            "method",
            "selection",
            "policies",
            "budgets",
            "seeds",
            "execution",
            "environment",
            "baseline_references",
            "published_step_identity",
            "publication_occurrence_identity",
            "independent_ids",
            "outcome",
            "manifest",
            "artifacts",
        }
        if (
            set(workflow) != workflow_keys
            or not isinstance(artifacts_record, list)
            or not isinstance(independent_raw, list)
            or any(
                not isinstance(item, dict)
                or set(item) != {"path", "role", "sha256"}
                or not str(cast(dict[str, object], item).get("path")).startswith(
                    f"{step['path']}/"
                )
                for item in artifacts_record
            )
        ):
            raise NativeProvenanceError("Native workflow children are malformed")
        artifact_records = cast(list[dict[str, object]], artifacts_record)
        artifacts = tuple(
            ArtifactReference(
                str(item["path"]).removeprefix(f"{step['path']}/"),
                cast(ArtifactRole, str(item["role"])),
                str(item["sha256"]),
            )
            for item in artifact_records
        )
        independent_ids = tuple(str(item) for item in independent_raw)
        rebuilt_step_identity = published_step_reference_identity(
            publication_occurrence_identity=str(
                manifest["publication_occurrence_identity"]
            ),
            lifecycle=str(manifest["lifecycle"]),
            authority=str(manifest["authority"]),
            manifest_identity=str(manifest["manifest_identity"]),
            manifest_sha256=sha256(manifest_bytes).hexdigest(),
            provenance=provenance,
            independent_ids=independent_ids,
            artifacts=artifacts,
        )
        canonical_workflow = provenance.to_record()
        if (
            step["manifest_sha256"] != sha256(manifest_bytes).hexdigest()
            or step["manifest_identity"] != manifest["manifest_identity"]
            or step["publication_occurrence_identity"]
            != manifest["publication_occurrence_identity"]
            or step["lifecycle"] != manifest["lifecycle"]
            or step["authority"] != manifest["authority"]
            or step["workflow_identity"] != provenance.workflow_identity
            or step["execution_identity"] != provenance.execution_identity
            or step["method_identity"] != provenance.method_identity
            or step["provenance_identity"] != provenance.identity
            or step["published_step_identity"] != rebuilt_step_identity
            or any(
                workflow.get(key) != value for key, value in canonical_workflow.items()
            )
            or workflow.get("published_step_identity") != rebuilt_step_identity
            or workflow.get("publication_occurrence_identity")
            != manifest["publication_occurrence_identity"]
            or workflow.get("independent_ids") != list(independent_ids)
            or workflow.get("outcome")
            != {
                "lifecycle": manifest["lifecycle"],
                "authority": manifest["authority"],
            }
            or workflow.get("manifest")
            != {
                "path": step["manifest_path"],
                "identity": manifest["manifest_identity"],
                "sha256": step["manifest_sha256"],
            }
            or normalized_link["workflow_identity"] != provenance.workflow_identity
            or normalized_link["method_identity"] != provenance.method_identity
            or _safe_link(root, normalized_link["path"]).read_text(encoding="utf-8")
            != provenance.normalized_method_text
            or provenance.environment.to_record()
            != {
                **run_environment,
                "numerical_libraries": provenance.environment.to_record()[
                    "numerical_libraries"
                ],
            }
        ):
            raise NativeProvenanceError("Native run step lineage is contradictory")
        manifest_artifacts = manifest["artifacts"]
        if {
            item.path: {"role": item.role, "sha256": item.sha256} for item in artifacts
        } != manifest_artifacts:
            raise NativeProvenanceError(
                "Native run artifact catalogue is contradictory"
            )
        actual_step_members = {
            str(item.relative_to(step_path))
            for item in step_path.rglob("*")
            if item.is_file()
        }
        if actual_step_members != {
            "fit-manifest.toml",
            *(item.path for item in artifacts),
        }:
            raise NativeProvenanceError(
                "Native run step artifact catalogue is incomplete"
            )
        for artifact in artifacts:
            artifact_path = _safe_link(step_path, artifact.path)
            if sha256(artifact_path.read_bytes()).hexdigest() != artifact.sha256:
                raise NativeProvenanceError("Native run artifact hash does not match")
        if manifest["lifecycle"] == "committed":
            restart_record = CommittedRestartRecord.from_record(
                json.loads(
                    (step_path / "Parameters" / "restart-provenance.json").read_text()
                )
            )
            if (
                restart_record.identity != manifest["restart_record_identity"]
                or restart_record.starting_occurrence_identity != current_occurrence
                or restart_record.starting_revision != current_revision
                or restart_record.workflow_identity != provenance.workflow_identity
                or restart_record.accepted_result_identity
                != manifest["accepted_result_identity"]
                or restart_record.accepted_occurrence_identity
                != manifest["accepted_occurrence_identity"]
                or restart_record.commit_operation_identity
                != manifest["commit_operation_identity"]
                or restart_record.commit_occurrence_identity
                != manifest["commit_occurrence_identity"]
                or restart_record.committed_occurrence_identity
                != manifest["committed_occurrence_identity"]
                or restart_record.committed_revision != manifest["committed_revision"]
                or restart_record.problem_identity != manifest["problem_identity"]
                or restart_record.parameterization_identity
                != manifest["parameterization_identity"]
                or tuple(item[0] for item in restart_record.independent_items)
                != independent_ids
            ):
                raise NativeProvenanceError(
                    "Native committed restart lineage is invalid"
                )
            current_occurrence = restart_record.committed_occurrence_identity
            current_revision = restart_record.committed_revision
        published_ids.append(rebuilt_step_identity)

    execution_occurrence_identity = _execution_occurrence_identity(
        bundle.identity,
        str(starting["occurrence_identity"]),
        int(starting["revision"]),
        tuple(published_ids),
    )
    if (
        record.get("execution_occurrence_identity") != execution_occurrence_identity
        or workflows.get("execution_occurrence_identity")
        != execution_occurrence_identity
    ):
        raise NativeProvenanceError("Native execution occurrence identity is invalid")

    actual_root_members = {
        str(item.relative_to(root)) for item in root.rglob("*") if item.is_file()
    }
    if actual_root_members != expected_root_members:
        raise NativeProvenanceError("Native run information member catalogue changed")
    return HistoricalNativeRunInformation(
        str(record["invocation_identity"]),
        bundle.identity,
        str(record["record_identity"]),
        tuple(published_ids),
        execution_occurrence_identity,
    )


def write_native_run_info(
    args: Namespace,
    run: NativeRunInformation,
    *,
    argv: Sequence[str] | None = None,
    working_directory: Path | None = None,
    timestamp: datetime | None = None,
) -> None:
    """Atomically write schema-v2 native invocation provenance and archives."""
    cwd = (
        Path.cwd().resolve()
        if working_directory is None
        else working_directory.resolve()
    )
    output_directory = _resolve_path(args.output, cwd)
    run_info_path = output_directory / "run_info"
    git_metadata = _git_metadata()
    independent_ids = _native_independent_ids(run)
    workflow_records = _validated_workflow_records(run, output_directory)

    output_directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=output_directory,
        prefix=".run_info-",
    ) as staging_directory:
        staging_path = Path(staging_directory)
        copied_inputs = _write_captured_native_inputs(run.inputs, staging_path)
        starting_text = serialize_independent_parameters(
            run.parameter_model,
            independent_ids,
            run.starting_snapshot,
            state_kind="starting",
        )
        starting_path = staging_path / "parameters_used.toml"
        starting_path.write_text(starting_text, encoding="utf-8")
        normalized_directory = staging_path / "normalized"
        normalized_directory.mkdir()
        normalized_methods: list[tuple[str, str, str, str]] = []
        for index, step in enumerate(run.steps, start=1):
            relative_path = f"normalized/method-{index:04d}.toml"
            normalized_path = staging_path / relative_path
            normalized_path.write_text(
                step.provenance.normalized_method_text,
                encoding="utf-8",
            )
            normalized_methods.append(
                (
                    step.provenance.workflow_identity,
                    step.provenance.method_identity,
                    relative_path,
                    sha256(normalized_path.read_bytes()).hexdigest(),
                )
            )
        workflow_path = staging_path / "workflows.json"
        workflow_path.write_text(
            json.dumps(
                workflow_records,
                allow_nan=False,
                ensure_ascii=True,
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        created_at = datetime.now(UTC) if timestamp is None else timestamp
        starting_sha256 = sha256(starting_path.read_bytes()).hexdigest()
        workflow_records_sha256 = sha256(workflow_path.read_bytes()).hexdigest()
        provisional = _serialize_native_run(
            run=run,
            timestamp=created_at,
            working_directory=cwd,
            output_directory=output_directory,
            argv=tuple(sys.argv if argv is None else argv),
            copied_inputs=copied_inputs,
            git_metadata=git_metadata,
            starting_sha256=starting_sha256,
            normalized_methods=normalized_methods,
            workflow_records_sha256=workflow_records_sha256,
            record_identity="0" * 64,
        )
        provisional_record = tomllib.loads(provisional)
        run_text = _serialize_native_run(
            run=run,
            timestamp=created_at,
            working_directory=cwd,
            output_directory=output_directory,
            argv=tuple(sys.argv if argv is None else argv),
            copied_inputs=copied_inputs,
            git_metadata=git_metadata,
            starting_sha256=starting_sha256,
            normalized_methods=normalized_methods,
            workflow_records_sha256=workflow_records_sha256,
            record_identity=_native_run_record_identity(provisional_record),
        )
        (staging_path / "run.toml").write_text(run_text, encoding="utf-8")
        historical = read_native_run_information(staging_path)
        if (
            historical.invocation_identity != run.invocation_identity
            or historical.input_bundle_identity != run.inputs.identity
            or historical.published_step_identities
            != tuple(step.identity for step in run.steps)
            or historical.execution_occurrence_identity
            != run.execution_occurrence_identity
        ):
            raise NativeProvenanceError(
                "Staged native run information differs from the live occurrence"
            )
        publish_directory_noreplace(staging_path, run_info_path)

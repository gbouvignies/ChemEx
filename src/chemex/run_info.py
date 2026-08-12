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
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import StrEnum
from hashlib import blake2b, sha256
from pathlib import Path

from chemex import __version__
from chemex.atomic import publish_directory_noreplace
from chemex.containers.experiments import Experiments
from chemex.native_provenance import (
    PublishedStepReference,
    serialize_independent_parameters,
)
from chemex.parameters.database import ParameterStore
from chemex.parameters.parameterization import SealedParameterModel
from chemex.parameters.setting import ParamSetting
from chemex.parameters.values import AnalysisValuesSnapshot

SCHEMA_VERSION = 1
_INPUT_CATEGORIES = ("experiments", "parameters", "methods")


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


class RunInformationKind(StrEnum):
    """Supported meanings of historical and native run-information records."""

    HISTORICAL_V1 = "historical_schema_v1"
    NATIVE_V2 = "native_schema_v2"


@dataclass(frozen=True, slots=True)
class NativeRunInformation:
    """Complete invocation state and links to atomic native step outcomes."""

    invocation_identity: str
    parameter_model: SealedParameterModel
    starting_snapshot: AnalysisValuesSnapshot
    steps: tuple[PublishedStepReference, ...]

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
        normalized_methods = {
            step.provenance.normalized_method_text for step in self.steps
        }
        environments = {step.provenance.environment.identity for step in self.steps}
        workflow_ids = {step.provenance.workflow_identity for step in self.steps}
        if len(normalized_methods) != 1:
            raise ValueError(
                "Native steps must reference one normalized method archive"
            )
        if len(environments) != 1:
            raise ValueError("Native steps must share one resolved environment")
        if len(workflow_ids) != len(self.steps):
            raise ValueError("Native workflow identities must be unique")


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


def _validated_workflow_records(
    run: NativeRunInformation,
    output: Path,
) -> dict[str, object]:
    workflows: list[dict[str, object]] = []
    for step in run.steps:
        relative_step = _relative_step_path(step, output)
        manifest_path = step.path / "fit-manifest.toml"
        manifest_bytes = manifest_path.read_bytes()
        if sha256(manifest_bytes).hexdigest() != step.manifest_sha256:
            raise ValueError("Native step manifest changed before run publication")
        manifest = tomllib.loads(manifest_bytes.decode("utf-8"))
        if manifest.get("manifest_identity") != step.manifest_identity:
            raise ValueError(
                "Native step manifest identity changed before run publication"
            )
        if (
            manifest.get("lifecycle") != step.lifecycle
            or manifest.get("authority") != step.authority
            or manifest.get("workflow", {}).get("identity")
            != step.provenance.workflow_identity
            or manifest.get("workflow", {}).get("provenance_identity")
            != step.provenance.identity
            or manifest.get("method", {}).get("identity")
            != step.provenance.method_identity
            or manifest.get("method", {}).get("parameterization_identity")
            != step.provenance.parameterization_identity
            or manifest.get("method", {}).get("evaluation_plan_identity")
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
        record = step.provenance.to_record()
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
    normalized_method_sha256: str,
    workflow_records_sha256: str,
) -> str:
    environment = run.steps[0].provenance.environment
    lines = [
        "# Native ChemEx product-facing run information.",
        "# schema_version = 1 records remain historical and are never rewritten.",
        "",
        "schema_version = 2",
        'schema_kind = "native_product_run_information"',
        f"invocation_identity = {_quote_string(run.invocation_identity)}",
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
        "[normalized_method]",
        'path = "normalized/method.toml"',
        f"sha256 = {_quote_string(normalized_method_sha256)}",
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
    for category, files in copied_inputs.items():
        for copied_input in files:
            lines.extend(
                (
                    f"[[inputs.{category}]]",
                    f"provided_path = {_quote_string(str(copied_input.provided_path))}",
                    f"resolved_path = {_quote_string(str(copied_input.resolved_path))}",
                    f"copied_path = {_quote_string(str(copied_input.copied_path))}",
                    f"sha256 = {_quote_string(copied_input.sha256)}",
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
                f"method_identity = {_quote_string(step.provenance.method_identity)}",
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
    input_files = _collect_input_files(args, cwd)
    git_metadata = _git_metadata()
    independent_ids = _native_independent_ids(run)
    normalized_method = run.steps[0].provenance.normalized_method_text
    workflow_records = _validated_workflow_records(run, output_directory)

    output_directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=output_directory,
        prefix=".run_info-",
    ) as staging_directory:
        staging_path = Path(staging_directory)
        copied_inputs = _copy_inputs(input_files, staging_path)
        starting_text = serialize_independent_parameters(
            run.parameter_model,
            independent_ids,
            run.starting_snapshot,
            state_kind="starting",
        )
        starting_path = staging_path / "parameters_used.toml"
        starting_path.write_text(starting_text, encoding="utf-8")
        normalized_path = staging_path / "normalized" / "method.toml"
        normalized_path.parent.mkdir()
        normalized_path.write_text(normalized_method, encoding="utf-8")
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
        run_text = _serialize_native_run(
            run=run,
            timestamp=datetime.now(UTC) if timestamp is None else timestamp,
            working_directory=cwd,
            output_directory=output_directory,
            argv=tuple(sys.argv if argv is None else argv),
            copied_inputs=copied_inputs,
            git_metadata=git_metadata,
            starting_sha256=sha256(starting_path.read_bytes()).hexdigest(),
            normalized_method_sha256=sha256(normalized_path.read_bytes()).hexdigest(),
            workflow_records_sha256=sha256(workflow_path.read_bytes()).hexdigest(),
        )
        (staging_path / "run.toml").write_text(run_text, encoding="utf-8")
        publish_directory_noreplace(staging_path, run_info_path)

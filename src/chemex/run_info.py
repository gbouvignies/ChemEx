from __future__ import annotations

import json
import platform
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

from chemex import __version__
from chemex.containers.experiments import Experiments
from chemex.parameters.database import ParameterStore
from chemex.parameters.setting import ParamSetting

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

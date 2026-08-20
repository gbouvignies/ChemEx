"""Compare the four #657 native product cases with the frozen legacy outputs.

This is intentionally a narrow reader for the immutable v1 anchor archive, not
a revival of migration-core authority or gating.  Tolerances are declared in
``CASES`` before any ChemEx process is started.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import sys
import tarfile
import tomllib
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).parents[2]
ARCHIVE = ROOT / "tests/fixtures/migration_core_canonical_anchor_release_v1.tar.xz"
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
UNCERTAINTY = re.compile(rf"±\s*(?P<value>{FLOAT})")


@dataclass(frozen=True, slots=True)
class FallbackTolerance:
    population_abs: float
    shift_abs: float
    relaxation_abs: float
    rate_relative: float
    other_relative: float
    other_abs: float


@dataclass(frozen=True, slots=True)
class Case:
    slug: str
    attempt: str
    example: Path
    model: str
    experiments: tuple[str, ...]
    parameters: tuple[str, ...]
    methods: tuple[str, ...]
    fallback: FallbackTolerance


# Frozen uncertainties remain the first choice.  These fallbacks cover values
# for which the legacy report has no finite positive uncertainty: populations
# are fractions, shifts are ppm, relaxation values are s^-1, and kinetic rates
# span orders of magnitude.  CPMG/CEST use 1% rate/other relative tolerance;
# binding and the more weakly conditioned three-state D-CEST case use 2%.
CASES = (
    Case(
        "cpmg-15n-ip",
        "96fdb5d9f0425b92a892e7ea784d95481ebeebbfa62655785b79fc1ef3a3ba4b",
        ROOT / "examples/Experiments/CPMG_15N_IP",
        "2st",
        ("Experiments/*.toml",),
        ("Parameters/parameters.toml",),
        ("Methods/method.toml",),
        FallbackTolerance(0.002, 0.01, 0.02, 0.01, 0.01, 1.0e-5),
    ),
    Case(
        "cest-13c-label-cn",
        "b2aa92f23ed505402bd9419b717c8a245f06ebd525d8ff8be11fce3465f2f614",
        ROOT / "examples/Experiments/CEST_13C_LABEL_CN",
        "2st",
        ("Experiments/*.toml",),
        ("Parameters/parameters.toml",),
        ("Methods/method.toml",),
        FallbackTolerance(0.002, 0.01, 0.02, 0.01, 0.01, 1.0e-5),
    ),
    Case(
        "2st-binding",
        "2c87cd4769195e59860646f380adb0d72fb05b6844b41c51e40969f922664b1a",
        ROOT / "examples/Combinations/2stBinding",
        "2st_binding",
        ("Experiments/*.toml",),
        ("Parameters/*.toml",),
        ("Methods/*.toml",),
        FallbackTolerance(0.003, 0.015, 0.03, 0.02, 0.02, 2.0e-5),
    ),
    Case(
        "dcest-fifu-drd",
        "714f523c0fedb5a4bc422ddfa51e0a1fbf45bcb2ecc81444e1ce0a70ba87fbb5",
        ROOT / "examples/Experiments/DCEST_15N_3States",
        "3st_fork",
        ("Experiments/*.toml", "Experiments/DRD/*.toml"),
        ("Parameters/parameters.toml",),
        ("Methods/method.toml",),
        FallbackTolerance(0.003, 0.015, 0.03, 0.02, 0.02, 2.0e-5),
    ),
)


@dataclass(frozen=True, slots=True)
class ParameterValue:
    identifier: tuple[str, str]
    value: float
    uncertainty: float | None


@dataclass(frozen=True, slots=True)
class DataRow:
    profile: str
    masked: bool
    values: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class FitRow:
    profile: str
    coordinates: tuple[float, ...]


def _members(archive: tarfile.TarFile, case: Case) -> dict[str, str]:
    root = f"publisher/attempts/{case.attempt}/terminal/success"
    manifest_file = archive.extractfile(f"{root}/manifest.json")
    if manifest_file is None:
        raise AssertionError(f"missing frozen manifest for {case.slug}")
    manifest = json.load(manifest_file)
    return {
        member["role"].removeprefix("legacy-output:"): (
            f"{root}/members/{member['content_hash']}"
        )
        for member in manifest["bundle"]["members"]
        if member["role"].startswith("legacy-output:")
    }


def _read_member(archive: tarfile.TarFile, name: str) -> str:
    file = archive.extractfile(name)
    if file is None:
        raise AssertionError(f"missing frozen member {name}")
    return file.read().decode()


def _comparable_paths(paths: set[str]) -> set[str]:
    # The frozen reference intentionally omitted regenerated PDF and .exp plot
    # companions.  Every retained frozen product path (including .fit curves
    # and statistics) must otherwise match exactly; run_info is a newer extra.
    return {
        path
        for path in paths
        if "run_info" not in Path(path).parts
        and Path(path).suffix not in {".exp", ".pdf"}
    }


def _data_rows(text: str) -> tuple[tuple[str, ...], tuple[DataRow, ...]]:
    profiles: list[str] = []
    rows: list[DataRow] = []
    profile = ""
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            profile = line
            profiles.append(line)
            continue
        masked = line.startswith("#") and "NOT USED IN THE FIT" in line
        if line.startswith("#") and not masked:
            continue
        if not profile:
            raise AssertionError("data row appears before its profile section")
        numeric = line.removeprefix("#").split("#", 1)[0]
        values = tuple(float(value) for value in numeric.split())
        if not values or any(not math.isfinite(value) for value in values):
            raise AssertionError(f"{profile}: data contains a non-finite value")
        rows.append(DataRow(profile, masked, values))
    return tuple(profiles), tuple(rows)


def _parameter_values(text: str) -> tuple[ParameterValue, ...]:
    section = ""
    result: list[ParameterValue] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line
            continue
        key, separator, remainder = line.partition("=")
        if not separator:
            continue
        value_text, _comment, comment = remainder.partition("#")
        match = UNCERTAINTY.search(comment)
        value = float(value_text.strip())
        uncertainty = None if match is None else float(match.group("value"))
        if not math.isfinite(value) or (
            uncertainty is not None and not math.isfinite(uncertainty)
        ):
            raise AssertionError(f"{section} {key.strip()}: non-finite parameter data")
        result.append(
            ParameterValue(
                (section, key.strip()),
                value,
                uncertainty,
            )
        )
    return tuple(result)


def _fit_schema(
    text: str,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[FitRow, ...]]:
    profiles: list[str] = []
    headers: list[str] = []
    rows: list[FitRow] = []
    profile = ""
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            profile = line
            profiles.append(line)
            continue
        if line.startswith("#"):
            headers.append(line)
            continue
        if not profile:
            raise AssertionError("fit row appears before its profile section")
        values = tuple(float(value) for value in line.split())
        if len(values) < 2 or any(not math.isfinite(value) for value in values):
            raise AssertionError(f"{profile}: fit curve has an invalid numeric row")
        rows.append(FitRow(profile, values[:-1]))
    return tuple(profiles), tuple(headers), tuple(rows)


def _statistics_schema(text: str) -> tuple[tuple[str, str], ...]:
    record = tomllib.loads(text)
    return tuple((key, type(value).__name__) for key, value in record.items())


def _compare_output_schema(path: str, legacy: str, native: str) -> None:
    if path.endswith(".fit"):
        if _fit_schema(native) != _fit_schema(legacy):
            raise AssertionError(f"{path}: fit curve schema or row order differs")
    elif Path(path).name == "statistics.toml":
        if _statistics_schema(native) != _statistics_schema(legacy):
            raise AssertionError(f"{path}: statistics schema differs")
    else:
        raise AssertionError(f"{path}: unsupported frozen output schema")


def _fallback_limit(identifier: tuple[str, str], legacy: float, case: Case) -> float:
    name = identifier[0].upper() + " " + identifier[1].upper()
    parameter_name = identifier[1].upper()
    policy = case.fallback
    if re.fullmatch(r"P[A-Z]", parameter_name) or "POPULATION" in name:
        return policy.population_abs
    if "DW_" in name or "CS_" in name:
        return policy.shift_abs
    if any(token in name for token in ("R1_", "R2_", "R1A_", "R2A_")):
        return policy.relaxation_abs
    if any(token in name for token in ("KEX", "KAB", "KBA", "KAC", "KCA")):
        return max(policy.other_abs, policy.rate_relative * abs(legacy))
    return max(policy.other_abs, policy.other_relative * abs(legacy))


def _compare_data(
    path: str, legacy_text: str, native_text: str
) -> tuple[float, str | None]:
    legacy_profiles, legacy_rows = _data_rows(legacy_text)
    native_profiles, native_rows = _data_rows(native_text)
    if native_profiles != legacy_profiles:
        raise AssertionError(f"{path}: selected profile order differs")
    if len(native_rows) != len(legacy_rows):
        raise AssertionError(f"{path}: selected data-point count differs")
    worst_fraction = 0.0
    worst_detail: str | None = None
    for index, (legacy, native) in enumerate(
        zip(legacy_rows, native_rows, strict=True)
    ):
        if (
            native.profile != legacy.profile
            or native.masked != legacy.masked
            or len(legacy.values) < 4
            or len(native.values) != len(legacy.values)
        ):
            raise AssertionError(f"{path}:{index}: data schema differs")
        if native.values[:-1] != legacy.values[:-1]:
            raise AssertionError(f"{path}:{index}: metadata/exp/error ordering differs")
        limit = 0.1 * abs(legacy.values[-2])
        difference = abs(native.values[-1] - legacy.values[-1])
        if limit:
            fraction = difference / limit
            if fraction > worst_fraction:
                worst_fraction = fraction
                worst_detail = (
                    f"{path}:{index}: calculated difference {difference} versus "
                    f"0.1 * experimental error = {limit}"
                )
    return worst_fraction, worst_detail


def _compare_parameters(
    path: str,
    legacy_text: str,
    native_text: str,
    case: Case,
) -> tuple[float, str | None]:
    legacy = _parameter_values(legacy_text)
    native = _parameter_values(native_text)
    if tuple(item.identifier for item in native) != tuple(
        item.identifier for item in legacy
    ):
        raise AssertionError(f"{path}: parameter identifiers/order or roles differ")
    worst_fraction = 0.0
    worst_detail: str | None = None
    for legacy_item, native_item in zip(legacy, native, strict=True):
        uncertainty = legacy_item.uncertainty
        limit = (
            0.1 * uncertainty
            if uncertainty is not None
            and math.isfinite(uncertainty)
            and uncertainty > 0.0
            else _fallback_limit(legacy_item.identifier, legacy_item.value, case)
        )
        difference = abs(native_item.value - legacy_item.value)
        if limit:
            fraction = difference / limit
            if fraction > worst_fraction:
                worst_fraction = fraction
                worst_detail = (
                    f"{path}:{legacy_item.identifier}: parameter difference "
                    f"{difference} versus limit {limit}"
                )
    return worst_fraction, worst_detail


def _expand(case: Case, patterns: tuple[str, ...]) -> list[str]:
    return [
        str(path) for pattern in patterns for path in sorted(case.example.glob(pattern))
    ]


def _run_case(case: Case, output: Path) -> None:
    if output.exists():
        shutil.rmtree(output)
    command = [
        sys.executable,
        "-m",
        "chemex",
        "fit",
        "-e",
        *_expand(case, case.experiments),
        "-p",
        *_expand(case, case.parameters),
        "-m",
        *_expand(case, case.methods),
        "-d",
        case.model,
        "-o",
        str(output),
        "--workers",
        "1",
        "--native-threads",
        "1",
    ]
    subprocess.run(command, cwd=case.example, check=True)  # noqa: S603


def compare_case(
    archive: tarfile.TarFile, case: Case, output: Path
) -> dict[str, object]:
    members = _members(archive, case)
    expected = _comparable_paths(set(members))
    actual = _comparable_paths(
        {
            path.relative_to(output).as_posix()
            for path in output.rglob("*")
            if path.is_file()
        }
    )
    if actual != expected:
        raise AssertionError(
            f"{case.slug}: structured output paths differ; "
            f"missing={sorted(expected - actual)!r}, unexpected={sorted(actual - expected)!r}"
        )

    worst_data_fraction = 0.0
    worst_parameter_fraction = 0.0
    worst_data_detail: str | None = None
    worst_parameter_detail: str | None = None
    for path in sorted(expected):
        legacy = _read_member(archive, members[path])
        native = (output / path).read_text(encoding="utf-8")
        if path.endswith(".dat"):
            fraction, detail = _compare_data(path, legacy, native)
            if fraction > worst_data_fraction:
                worst_data_fraction = fraction
                worst_data_detail = detail
        elif "Parameters" in Path(path).parts:
            fraction, detail = _compare_parameters(path, legacy, native, case)
            if fraction > worst_parameter_fraction:
                worst_parameter_fraction = fraction
                worst_parameter_detail = detail
        else:
            _compare_output_schema(path, legacy, native)
    passed = worst_data_fraction <= 1.0 and worst_parameter_fraction <= 1.0
    return {
        "status": "PASS" if passed else "FAIL",
        "comparable_output_path_count": len(expected),
        "worst_calculated_fraction_of_limit": worst_data_fraction,
        "worst_parameter_fraction_of_limit": worst_parameter_fraction,
        "worst_calculated_comparison": worst_data_detail,
        "worst_parameter_comparison": worst_parameter_detail,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--run", action="store_true")
    parser.add_argument(
        "--case", choices=[case.slug for case in CASES], action="append"
    )
    args = parser.parse_args()
    selected = tuple(case for case in CASES if not args.case or case.slug in args.case)
    args.output_root.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[str, object]] = {}
    with tarfile.open(args.archive, mode="r:xz") as archive:
        for case in selected:
            output = args.output_root / case.slug
            if args.run:
                _run_case(case, output)
            try:
                results[case.slug] = compare_case(archive, case, output)
            except AssertionError as error:
                results[case.slug] = {
                    "status": "FAIL",
                    "structural_failure": str(error),
                }
    print(json.dumps(results, indent=2, sort_keys=True))
    if any(result["status"] != "PASS" for result in results.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

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

# This is the one case-specific #657 disposition established by same-vector
# cross-evaluation on Asterix.  It is deliberately not a tolerance waiver:
# legacy-parity remains failed, while the independently diagnosed same objective
# has a reproducible lower native basin.  The published statistics values below
# are exact at ChemEx's output precision and prevent this disposition from
# accepting an unrelated future numerical change.
_CEST_ACCEPTED_DIFFERENCE = {
    "statistics_path": "STEP2/All/statistics.toml",
    "published_legacy_chi_square": 3140.90,
    "published_native_chi_square": 2849.45,
    "diagnosed_legacy_chi_square": 3140.903510372318,
    "diagnosed_native_chi_square": 2849.4551945314374,
    "diagnosed_legacy_l18cd1_chi_square": 1093.424072389286,
    "diagnosed_native_l18cd1_chi_square": 801.9806529406368,
}
_CEST_ACCEPTED_FAILURE_SCOPE = frozenset(
    {
        "data:STEP2/All/Data/23hz.dat:[L18CD1]",
        "data:STEP2/Groups/3_L18CD1/Data/23hz.dat:[L18CD1]",
        *(
            f"parameter:{root}/constrained.toml:{identifier!r}"
            for root in (
                "STEP2/All/Parameters",
                "STEP2/Groups/3_L18CD1/Parameters",
            )
            for identifier in (
                ('["R1_B, B0->598.8MHZ"]', "L18CD1"),
                ('["R2_B, B0->598.8MHZ"]', "L18CD1"),
                ("[CS_B]", "L18CD1"),
            )
        ),
        *(
            f"parameter:{root}/fitted.toml:{identifier!r}"
            for root in (
                "STEP2/All/Parameters",
                "STEP2/Groups/3_L18CD1/Parameters",
            )
            for identifier in (
                ('["R1_A, B0->598.8MHZ"]', "L18CD1"),
                ('["R2_A, B0->598.8MHZ"]', "L18CD1"),
                ("[CS_A]", "L18CD1"),
                ("[DW_AB]", "L18CD1"),
            )
        ),
    }
)
_CEST_ACCEPTED_NUMERICAL_SIGNATURE = {
    "legacy_cs_a": 24.9724,
    "native_cs_a": 24.9288,
    "legacy_dw_ab": -0.142414,
    "native_dw_ab": 0.193203,
    "legacy_cs_b": 24.8300,
    "native_cs_b": 25.1220,
    "legacy_l18cd1_chi_square": 1093.424158056536,
    "native_l18cd1_chi_square": 801.9787254040632,
}
# These identify the diagnosed basin at product-output precision; they do not
# relax or replace any legacy-parity tolerance declared in ``CASES``.
_CEST_SIGNATURE_ABSOLUTE_LIMITS = {
    "legacy_cs_a": 0.001,
    "native_cs_a": 0.001,
    "legacy_dw_ab": 0.001,
    "native_dw_ab": 0.001,
    "legacy_cs_b": 0.001,
    "native_cs_b": 0.001,
    "legacy_l18cd1_chi_square": 0.01,
    "native_l18cd1_chi_square": 0.01,
}


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


def _comparison_disposition(
    case: Case,
    parity_passed: bool,
    published_chi_squares: dict[str, tuple[float, float]],
    parity_failure_scope: frozenset[str],
    numerical_signature: dict[str, float],
) -> dict[str, object]:
    if parity_passed:
        return {"status": "PASS_PARITY"}
    if case.slug != "cest-13c-label-cn":
        return {"status": "FAIL"}

    path = str(_CEST_ACCEPTED_DIFFERENCE["statistics_path"])
    observed = published_chi_squares.get(path)
    expected = (
        float(_CEST_ACCEPTED_DIFFERENCE["published_legacy_chi_square"]),
        float(_CEST_ACCEPTED_DIFFERENCE["published_native_chi_square"]),
    )
    signature_matches = numerical_signature.keys() == (
        _CEST_ACCEPTED_NUMERICAL_SIGNATURE.keys()
    ) and all(
        math.isclose(
            numerical_signature[name],
            expected_value,
            rel_tol=0.0,
            abs_tol=_CEST_SIGNATURE_ABSOLUTE_LIMITS[name],
        )
        for name, expected_value in _CEST_ACCEPTED_NUMERICAL_SIGNATURE.items()
    )
    if (
        observed != expected
        or parity_failure_scope != _CEST_ACCEPTED_FAILURE_SCOPE
        or not signature_matches
    ):
        return {"status": "FAIL"}
    return {
        "status": "ACCEPTED_DIFFERENCE",
        "parity_status": "FAIL",
        "objective_semantics_parity": "PASS",
        "observed_numerical_signature": numerical_signature,
        "disposition": (
            "Frozen legacy observation is not the normative optimum; native and "
            "legacy evaluate the same objective, and native TRF reproducibly reaches "
            "the lower accepted basin."
        ),
        "diagnosed_legacy_chi_square": _CEST_ACCEPTED_DIFFERENCE[
            "diagnosed_legacy_chi_square"
        ],
        "diagnosed_native_chi_square": _CEST_ACCEPTED_DIFFERENCE[
            "diagnosed_native_chi_square"
        ],
        "diagnosed_legacy_l18cd1_chi_square": _CEST_ACCEPTED_DIFFERENCE[
            "diagnosed_legacy_l18cd1_chi_square"
        ],
        "diagnosed_native_l18cd1_chi_square": _CEST_ACCEPTED_DIFFERENCE[
            "diagnosed_native_l18cd1_chi_square"
        ],
    }


def _profile_chi_square(text: str, profile: str) -> float:
    _profiles, rows = _data_rows(text)
    return sum(
        ((row.values[-1] - row.values[-3]) / row.values[-2]) ** 2
        for row in rows
        if row.profile == profile and not row.masked
    )


def _cest_numerical_signature(
    path: str,
    legacy: str,
    native: str,
) -> dict[str, float]:
    if path == "STEP2/All/Data/23hz.dat":
        return {
            "legacy_l18cd1_chi_square": _profile_chi_square(legacy, "[L18CD1]"),
            "native_l18cd1_chi_square": _profile_chi_square(native, "[L18CD1]"),
        }
    if path not in {
        "STEP2/All/Parameters/fitted.toml",
        "STEP2/All/Parameters/constrained.toml",
    }:
        return {}
    legacy_values = {item.identifier: item.value for item in _parameter_values(legacy)}
    native_values = {item.identifier: item.value for item in _parameter_values(native)}
    identifiers = (
        {
            "cs_a": ("[CS_A]", "L18CD1"),
            "dw_ab": ("[DW_AB]", "L18CD1"),
        }
        if path.endswith("fitted.toml")
        else {"cs_b": ("[CS_B]", "L18CD1")}
    )
    return {
        f"{implementation}_{name}": values[identifier]
        for name, identifier in identifiers.items()
        for implementation, values in (
            ("legacy", legacy_values),
            ("native", native_values),
        )
    }


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
) -> tuple[float, str | None, frozenset[str]]:
    legacy_profiles, legacy_rows = _data_rows(legacy_text)
    native_profiles, native_rows = _data_rows(native_text)
    if native_profiles != legacy_profiles:
        raise AssertionError(f"{path}: selected profile order differs")
    if len(native_rows) != len(legacy_rows):
        raise AssertionError(f"{path}: selected data-point count differs")
    worst_fraction = 0.0
    worst_detail: str | None = None
    failures: set[str] = set()
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
            if fraction > 1.0:
                failures.add(f"data:{path}:{legacy.profile}")
            if fraction > worst_fraction:
                worst_fraction = fraction
                worst_detail = (
                    f"{path}:{index}: calculated difference {difference} versus "
                    f"0.1 * experimental error = {limit}"
                )
    return worst_fraction, worst_detail, frozenset(failures)


def _compare_parameters(
    path: str,
    legacy_text: str,
    native_text: str,
    case: Case,
) -> tuple[float, str | None, frozenset[str]]:
    legacy = _parameter_values(legacy_text)
    native = _parameter_values(native_text)
    if tuple(item.identifier for item in native) != tuple(
        item.identifier for item in legacy
    ):
        raise AssertionError(f"{path}: parameter identifiers/order or roles differ")
    worst_fraction = 0.0
    worst_detail: str | None = None
    failures: set[str] = set()
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
            if fraction > 1.0:
                failures.add(f"parameter:{path}:{legacy_item.identifier!r}")
            if fraction > worst_fraction:
                worst_fraction = fraction
                worst_detail = (
                    f"{path}:{legacy_item.identifier}: parameter difference "
                    f"{difference} versus limit {limit}"
                )
    return worst_fraction, worst_detail, frozenset(failures)


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
    published_chi_squares: dict[str, tuple[float, float]] = {}
    parity_failure_scope: set[str] = set()
    numerical_signature: dict[str, float] = {}
    for path in sorted(expected):
        legacy = _read_member(archive, members[path])
        native = (output / path).read_text(encoding="utf-8")
        if path.endswith(".dat"):
            fraction, detail, failures = _compare_data(path, legacy, native)
            parity_failure_scope.update(failures)
            if case.slug == "cest-13c-label-cn":
                numerical_signature.update(
                    _cest_numerical_signature(path, legacy, native)
                )
            if fraction > worst_data_fraction:
                worst_data_fraction = fraction
                worst_data_detail = detail
        elif "Parameters" in Path(path).parts:
            fraction, detail, failures = _compare_parameters(path, legacy, native, case)
            parity_failure_scope.update(failures)
            if case.slug == "cest-13c-label-cn":
                numerical_signature.update(
                    _cest_numerical_signature(path, legacy, native)
                )
            if fraction > worst_parameter_fraction:
                worst_parameter_fraction = fraction
                worst_parameter_detail = detail
        else:
            _compare_output_schema(path, legacy, native)
            if Path(path).name == "statistics.toml":
                published_chi_squares[path] = (
                    float(tomllib.loads(legacy)["chi-square"]),
                    float(tomllib.loads(native)["chi-square"]),
                )
    passed = worst_data_fraction <= 1.0 and worst_parameter_fraction <= 1.0
    return {
        **_comparison_disposition(
            case,
            passed,
            published_chi_squares,
            frozenset(parity_failure_scope),
            numerical_signature,
        ),
        "comparable_output_path_count": len(expected),
        "parity_failure_scope": sorted(parity_failure_scope),
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
    successful = {"PASS_PARITY", "ACCEPTED_DIFFERENCE"}
    if any(result["status"] not in successful for result in results.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

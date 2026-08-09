"""Attested, content-addressed numerical execution lanes.

Numerical lanes are baseline evidence infrastructure.  They deliberately do
not select a solver or alter the product's execution settings.  A lane gains
authority only after its already-imported numerical stack matches every claim
frozen in its definition.
"""

from __future__ import annotations

import ctypes
import hashlib
import json
import math
import os
import platform as platform_module
import struct
import sys
from collections.abc import Sequence
from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path
from typing import Any, Literal

from chemex.runtime.execution import NATIVE_THREAD_ENV_VARS

_SCHEMA_VERSION = 1
_SEMANTIC_VERSION = "chemex-numerical-lane-v1"
_HASH_LENGTH = 64
_RECIPE_NAME = "numerical-lane.Dockerfile"
_PROVENANCE_PATH = Path("/opt/chemex-numerical-lane/provenance.json")
_PLATFORM_MANIFEST = (
    "debian:bookworm-slim@sha256:"
    "362e64223cc0da95422b3b13c045186fc0a81250e765d31c025fbddf257f6143"
)
_LOCKFILE_HASH = "f0fb2ffc7b1a5ecd1bf7ac43956fc4861b96c058d158948b68b4e97027a6086a"
_NUMERICAL_LIBRARY = "scipy-openblas-0.3.33.112.0"
_ISA = "x86-64-v3"
_FLOATING_POINT_MODE = "ieee754-binary64-round-to-nearest"
type ComparisonScopeKind = Literal[
    "WITHIN_LANE_BITWISE", "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON"
]


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {
            "kind": kind,
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "record": record,
        },
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return _sha256(encoded)


def _digest(value: str, name: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != _HASH_LENGTH
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{name} must be a SHA-256 digest")
    return value


def _non_empty(value: str, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a non-empty string")
    return value


def _pinned_platform_manifest(value: str) -> str:
    repository, separator, digest = value.rpartition("@sha256:")
    if not separator or not repository:
        raise ValueError("platform manifest must use an immutable SHA-256 digest")
    _digest(digest, "platform manifest digest")
    return value


def _recipe_hash() -> str:
    return _sha256(files(__package__).joinpath(_RECIPE_NAME).read_bytes())


def _canonical_architecture(value: str) -> str:
    return {"x86_64": "amd64", "aarch64": "arm64"}.get(value, value)


def _unavailable_provenance() -> tuple[str, str, str, str]:
    unavailable = "0" * _HASH_LENGTH
    return (
        unavailable,
        unavailable,
        unavailable,
        f"unavailable@sha256:{unavailable}",
    )


def _current_provenance(path: Path) -> tuple[str, str, str, str]:
    try:
        record = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return _unavailable_provenance()
    if not isinstance(record, dict) or set(record) != {
        "build_recipe_hash",
        "dependency_lock_hash",
        "platform_manifest",
        "python_source_hash",
    }:
        return _unavailable_provenance()
    build_recipe_hash = record["build_recipe_hash"]
    dependency_lock_hash = record["dependency_lock_hash"]
    platform_manifest = record["platform_manifest"]
    python_source_hash = record["python_source_hash"]
    if not all(
        isinstance(value, str)
        for value in (
            build_recipe_hash,
            dependency_lock_hash,
            platform_manifest,
            python_source_hash,
        )
    ):
        return _unavailable_provenance()
    try:
        return (
            _digest(build_recipe_hash, "build recipe hash"),
            _digest(dependency_lock_hash, "dependency lock hash"),
            _digest(python_source_hash, "Python source hash"),
            _pinned_platform_manifest(platform_manifest),
        )
    except ValueError:
        return _unavailable_provenance()


class LaneAuthorityError(RuntimeError):
    """Raised when an execution cannot claim a lane's authority."""


@dataclass(frozen=True, slots=True)
class RuntimeEnvironment:
    """The post-import facts that determine numerical-lane authority."""

    python_implementation: str
    python_version: str
    python_abi: str
    python_source_hash: str
    platform: str
    platform_manifest: str
    dependency_lock_hash: str
    build_recipe_hash: str
    numerical_library: str
    isa: str
    workers: int | None
    native_threads: int | None
    floating_point_mode: str
    imported_packages: tuple[str, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        for value, name in (
            (self.python_implementation, "Python implementation"),
            (self.python_version, "Python version"),
            (self.python_abi, "Python ABI"),
            (self.platform, "platform"),
            (self.numerical_library, "numerical library"),
            (self.isa, "ISA"),
            (self.floating_point_mode, "floating-point mode"),
        ):
            _non_empty(value, name)
        _digest(self.dependency_lock_hash, "dependency lock hash")
        _digest(self.build_recipe_hash, "build recipe hash")
        _digest(self.python_source_hash, "Python source hash")
        _pinned_platform_manifest(self.platform_manifest)
        for value, name in (
            (self.workers, "workers"),
            (self.native_threads, "native numerical threads"),
        ):
            if value is not None and value < 1:
                raise ValueError(f"{name} must be positive or unavailable")
        packages = self.imported_packages
        if (
            not packages
            or tuple(sorted(packages)) != packages
            or len(set(packages)) != len(packages)
            or any(not isinstance(package, str) or not package for package in packages)
        ):
            raise ValueError(
                "imported packages must be unique, ordered, non-empty strings"
            )
        object.__setattr__(
            self,
            "identity",
            _identity(
                "post-import-runtime-environment",
                (
                    self.python_implementation,
                    self.python_version,
                    self.python_abi,
                    self.python_source_hash,
                    self.platform,
                    self.platform_manifest,
                    self.dependency_lock_hash,
                    self.build_recipe_hash,
                    self.numerical_library,
                    self.isa,
                    self.workers,
                    self.native_threads,
                    self.floating_point_mode,
                    self.imported_packages,
                ),
            ),
        )

    @classmethod
    def from_current_process(
        cls, provenance_path: Path = _PROVENANCE_PATH
    ) -> RuntimeEnvironment:
        """Capture facts after importing NumPy and SciPy, without repairing them."""
        import numpy
        import scipy

        (
            build_recipe_hash,
            dependency_lock_hash,
            python_source_hash,
            platform_manifest,
        ) = _current_provenance(provenance_path)
        return cls(
            python_implementation=platform_module.python_implementation(),
            python_version=".".join(map(str, sys.version_info[:3])),
            python_abi=f"cp{sys.version_info.major}{sys.version_info.minor}",
            python_source_hash=python_source_hash,
            platform=(
                f"{sys.platform}/"
                f"{_canonical_architecture(platform_module.machine().lower())}"
            ),
            platform_manifest=platform_manifest,
            dependency_lock_hash=dependency_lock_hash,
            build_recipe_hash=build_recipe_hash,
            numerical_library=_current_numerical_library(numpy),
            isa=_current_isa(numpy),
            workers=_current_workers(),
            native_threads=_current_native_threads(),
            floating_point_mode=_current_floating_point_mode(),
            imported_packages=(
                f"numpy=={numpy.__version__}",
                f"scipy=={scipy.__version__}",
            ),
        )


def _current_numerical_library(numpy: Any) -> str:
    config = numpy.__config__.CONFIG
    dependencies = config.get("Build Dependencies", {})
    blas = dependencies.get("blas", {})
    name = blas.get("name")
    version = blas.get("version")
    if not isinstance(name, str) or not isinstance(version, str):
        return "unavailable"
    return f"{name}-{version}"


def _current_isa(numpy: Any) -> str:
    features = numpy._core._multiarray_umath.__cpu_features__
    if features.get("X86_V3") and not features.get("X86_V4"):
        return _ISA
    return "unavailable"


def _current_native_threads() -> int | None:
    values = tuple(os.environ.get(name) for name in NATIVE_THREAD_ENV_VARS)
    if any(value is None for value in values):
        return None
    try:
        counts = {int(value) for value in values if value is not None}
    except ValueError:
        return None
    if len(counts) != 1:
        return None
    count = counts.pop()
    return count if count > 0 else None


def _current_workers() -> int | None:
    value = os.environ.get("CHEMEX_NUMERICAL_LANE_WORKERS")
    try:
        workers = int(value) if value is not None else 0
    except ValueError:
        return None
    return workers if workers > 0 else None


def _current_floating_point_mode() -> str:
    try:
        fegetround = ctypes.CDLL(None).fegetround
        fegetround.restype = ctypes.c_int
        rounding_mode = fegetround()
    except (AttributeError, OSError):
        return "unavailable"
    if (
        rounding_mode == 0
        and sys.float_info.radix == 2
        and struct.pack(">d", 1.0) == b"?\xf0\x00\x00\x00\x00\x00\x00"
    ):
        return _FLOATING_POINT_MODE
    return "unavailable"


@dataclass(frozen=True, slots=True)
class LaneAttestation:
    """A content-addressed proof that one process satisfies one lane."""

    lane_identity: str
    environment_identity: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _digest(self.lane_identity, "lane identity")
        _digest(self.environment_identity, "environment identity")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "numerical-lane-attestation",
                (self.lane_identity, self.environment_identity),
            ),
        )


@dataclass(frozen=True, slots=True)
class NumericalLane:
    """One closed numerical environment and its replay promise."""

    name: str
    python_version: str
    python_abi: str
    python_source_hash: str
    platform: str
    platform_manifest: str
    dependency_lock_hash: str
    build_recipe_hash: str
    numerical_library: str
    isa: str
    workers: int
    native_threads: int
    floating_point_mode: str
    required_packages: tuple[str, ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        for value, name in (
            (self.name, "lane name"),
            (self.python_version, "Python version"),
            (self.python_abi, "Python ABI"),
            (self.platform, "platform"),
            (self.platform_manifest, "platform manifest"),
            (self.numerical_library, "numerical library"),
            (self.isa, "ISA"),
            (self.floating_point_mode, "floating-point mode"),
        ):
            _non_empty(value, name)
        _digest(self.python_source_hash, "Python source hash")
        _pinned_platform_manifest(self.platform_manifest)
        _digest(self.dependency_lock_hash, "dependency lock hash")
        _digest(self.build_recipe_hash, "build recipe hash")
        if self.workers < 1 or self.native_threads < 1:
            raise ValueError(
                "worker and native numerical thread counts must be positive"
            )
        if (
            not self.required_packages
            or tuple(sorted(self.required_packages)) != self.required_packages
        ):
            raise ValueError("required packages must be non-empty and ordered")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "numerical-lane",
                (
                    self.name,
                    self.python_version,
                    self.python_abi,
                    self.python_source_hash,
                    self.platform,
                    self.platform_manifest,
                    self.dependency_lock_hash,
                    self.build_recipe_hash,
                    self.numerical_library,
                    self.isa,
                    self.workers,
                    self.native_threads,
                    self.floating_point_mode,
                    self.required_packages,
                ),
            ),
        )

    def compatibility_delta(self, other: NumericalLane) -> dict[str, tuple[str, str]]:
        """Return the only fields allowed to differ between compatibility lanes."""
        values = {
            "python_version": (self.python_version, other.python_version),
            "python_abi": (self.python_abi, other.python_abi),
            "python_source_hash": (self.python_source_hash, other.python_source_hash),
            "platform": (self.platform, other.platform),
            "platform_manifest": (self.platform_manifest, other.platform_manifest),
            "dependency_lock_hash": (
                self.dependency_lock_hash,
                other.dependency_lock_hash,
            ),
            "build_recipe_hash": (self.build_recipe_hash, other.build_recipe_hash),
            "numerical_library": (self.numerical_library, other.numerical_library),
            "isa": (self.isa, other.isa),
            "workers": (str(self.workers), str(other.workers)),
            "native_threads": (str(self.native_threads), str(other.native_threads)),
            "floating_point_mode": (
                self.floating_point_mode,
                other.floating_point_mode,
            ),
            "required_packages": (
                str(self.required_packages),
                str(other.required_packages),
            ),
        }
        return {field: value for field, value in values.items() if value[0] != value[1]}

    def attest(self, environment: RuntimeEnvironment) -> LaneAttestation:
        """Grant authority only when every required post-import claim matches."""
        mismatches: list[str] = []
        expected = (
            ("Python implementation", "CPython", environment.python_implementation),
            ("Python version", self.python_version, environment.python_version),
            ("Python ABI", self.python_abi, environment.python_abi),
            (
                "Python source",
                self.python_source_hash,
                environment.python_source_hash,
            ),
            ("platform", self.platform, environment.platform),
            (
                "platform manifest",
                self.platform_manifest,
                environment.platform_manifest,
            ),
            (
                "dependency lock",
                self.dependency_lock_hash,
                environment.dependency_lock_hash,
            ),
            ("build recipe", self.build_recipe_hash, environment.build_recipe_hash),
            (
                "numerical library",
                self.numerical_library,
                environment.numerical_library,
            ),
            ("ISA", self.isa, environment.isa),
            ("workers", str(self.workers), str(environment.workers)),
            (
                "native numerical threads",
                str(self.native_threads),
                str(environment.native_threads),
            ),
            (
                "floating-point mode",
                self.floating_point_mode,
                environment.floating_point_mode,
            ),
        )
        mismatches.extend(
            f"{name}: expected {expected_value}, got {actual_value}"
            for name, expected_value, actual_value in expected
            if expected_value != actual_value
        )
        missing = sorted(
            set(self.required_packages) - set(environment.imported_packages)
        )
        if missing:
            mismatches.append(
                f"post-import package claims unavailable: {', '.join(missing)}"
            )
        if mismatches:
            raise LaneAuthorityError(
                "Numerical lane authority rejected: " + "; ".join(mismatches)
            )
        return LaneAttestation(self.identity, environment.identity)

    def attest_current_process(self) -> LaneAttestation:
        return self.attest(RuntimeEnvironment.from_current_process())


@dataclass(frozen=True, slots=True)
class ComparisonScope:
    """The distinct contract for replaying or comparing baseline evidence."""

    kind: ComparisonScopeKind
    left_lane_identity: str
    right_lane_identity: str

    def __post_init__(self) -> None:
        _digest(self.left_lane_identity, "left lane identity")
        _digest(self.right_lane_identity, "right lane identity")
        if self.kind == "WITHIN_LANE_BITWISE" and (
            self.left_lane_identity != self.right_lane_identity
        ):
            raise ValueError("Within-lane replay requires one lane identity")
        if self.kind == "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON" and (
            self.left_lane_identity == self.right_lane_identity
        ):
            raise ValueError("Cross-lane comparison requires distinct lane identities")


@dataclass(frozen=True, slots=True)
class CrossLaneNumericalPolicy:
    """Explicit tolerance policy required for cross-lane numerical comparison."""

    rtol: float
    atol: float

    def __post_init__(self) -> None:
        if not (
            math.isfinite(self.rtol)
            and math.isfinite(self.atol)
            and self.rtol >= 0.0
            and self.atol >= 0.0
        ):
            raise ValueError("cross-lane tolerances must be finite and non-negative")


@dataclass(frozen=True, slots=True)
class ComparisonOutcome:
    """The observed result of one scoped replay or numerical comparison."""

    scope: ComparisonScope
    equivalent: bool


def compare_values(
    scope: ComparisonScope,
    left: Sequence[float],
    right: Sequence[float],
    *,
    policy: CrossLaneNumericalPolicy | None = None,
) -> ComparisonOutcome:
    """Compare binary64 values without allowing tolerances into exact replay."""
    import numpy

    left_values = numpy.asarray(left, dtype=numpy.float64)
    right_values = numpy.asarray(right, dtype=numpy.float64)
    if scope.kind == "WITHIN_LANE_BITWISE":
        if policy is not None:
            raise ValueError("Within-lane replay does not permit a tolerance policy")
        equivalent = bool(
            left_values.shape == right_values.shape
            and numpy.array_equal(
                left_values.view(numpy.uint64), right_values.view(numpy.uint64)
            )
        )
    else:
        if policy is None:
            raise ValueError(
                "Cross-lane comparison requires an explicit numerical policy"
            )
        equivalent = bool(
            numpy.allclose(
                left_values,
                right_values,
                rtol=policy.rtol,
                atol=policy.atol,
                equal_nan=True,
            )
        )
    return ComparisonOutcome(scope, equivalent)


def comparison_scope(left: NumericalLane, right: NumericalLane) -> ComparisonScope:
    """Keep exact replay within a lane and numerical comparison across lanes."""
    if left.identity == right.identity:
        return ComparisonScope("WITHIN_LANE_BITWISE", left.identity, right.identity)
    return ComparisonScope(
        "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON", left.identity, right.identity
    )


def canonical_lanes() -> tuple[NumericalLane, NumericalLane]:
    """Return the 3.13 authority lane and the 3.14 compatibility replay lane."""
    shared = {
        "platform": "linux/amd64",
        "platform_manifest": _PLATFORM_MANIFEST,
        "dependency_lock_hash": _LOCKFILE_HASH,
        "build_recipe_hash": _recipe_hash(),
        "numerical_library": _NUMERICAL_LIBRARY,
        "isa": _ISA,
        "workers": 1,
        "native_threads": 1,
        "floating_point_mode": _FLOATING_POINT_MODE,
        "required_packages": ("numpy==2.5.1", "scipy==1.18.0"),
    }
    return (
        NumericalLane(
            "canonical-linux-amd64-python-3.13-v1",
            "3.13.5",
            "cp313",
            "e6190f52699b534ee203d9f417bdbca05a92f23e35c19c691a50ed2942835385",
            **shared,
        ),
        NumericalLane(
            "compatibility-linux-amd64-python-3.14-v1",
            "3.14.5",
            "cp314",
            "9c22bfe9939a6c5418fc74b289a5f1cc41859ae82ac6b163016b5844bd0a86bc",
            **shared,
        ),
    )

"""Attested, content-addressed numerical execution lanes.

The module has one authority-minting interface: :meth:`NumericalLane.attest_current_process`.
It combines an externally observed immutable OCI digest with facts measured from the
already-imported Python/NumPy/SciPy/OpenBLAS process.  Expected recipe metadata alone
never grants lane authority.
"""

from __future__ import annotations

import ctypes
import hashlib
import importlib.metadata
import json
import math
import os
import platform as platform_module
import struct
import subprocess
import sys
import sysconfig
import weakref
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, fields
from importlib.resources import files
from pathlib import Path, PurePosixPath
from typing import Any, Literal, NamedTuple, SupportsIndex, cast

_SCHEMA_VERSION = 2
_SEMANTIC_VERSION = "chemex-numerical-lane-v2"
_HASH_LENGTH = 64
_RECIPE_NAMES = ("build-image.sh", "numerical-lane.Dockerfile")
_PROVENANCE_PATH = Path("/opt/chemex-numerical-lane/provenance.json")
_BUILD_CONTEXT_MANIFEST_PATH = Path(
    "/opt/chemex-numerical-lane/build-context-manifest.txt"
)
_OS_PACKAGE_MANIFEST_PATH = Path("/opt/chemex-numerical-lane/os-package-manifest.txt")
_WHEEL_MANIFEST_PATH = Path("/opt/chemex-numerical-lane/wheel-manifest.txt")
_WHEELHOUSE_PATH = Path("/opt/chemex-numerical-lane/wheels")
_BUILD_ROOT = Path("/opt/chemex-numerical-lane/build-context")
_FPSTATE_PATH = Path("/opt/chemex-numerical-lane/libchemex_fpstate.so")
_HISTORICAL_RECIPE_HASH = (
    "af852fa2b39b40c39f9b84a8f8ae50246c66ba315656b22c5807bf59fb85e71a"
)
_HISTORICAL_LANE_NAMES = {
    "CANONICAL_NUMERICAL": "canonical-linux-amd64-python-3.13-v1",
    "PYTHON_COMPATIBILITY": "compatibility-linux-amd64-python-3.14-v1",
}
_PROSPECTIVE_LANE_NAMES = {
    "CANONICAL_NUMERICAL": "canonical-linux-amd64-python-3.13-v2",
    "PYTHON_COMPATIBILITY": "compatibility-linux-amd64-python-3.14-v2",
}
_PLATFORM_MANIFEST = (
    "debian:bookworm-slim@sha256:"
    "362e64223cc0da95422b3b13c045186fc0a81250e765d31c025fbddf257f6143"
)
_LOCKFILE_HASH = "cc7a8e08d8fb8f1ea4255b63452598f6dbe041a8b4024de0f3af065020088004"
_PROSPECTIVE_LOCKFILE_HASH = (
    "a632a77d8ef4d074b6c2bc26e1e05a826251d475468b7f5dc5c8efe30ee46735"
)
_UV_VERSION = "0.11.15"
_UV_WHEEL_HASH = "98edf1bdaf82447014852051d93e3ee95012509c567bf057fd117e6bdbd9a807"
_OPENBLAS_CORE = "Haswell"
_NUMPY_BLAS = (
    '{"name":"scipy-openblas","openblas configuration":'
    '"OpenBLAS 0.3.33.112.0  USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell '
    'MAX_THREADS=64","version":"0.3.33.112.0"}'
)
_SCIPY_BLAS_313 = (
    '{"cython blas ilp64":false,"has ilp64":false,"name":"scipy-openblas",'
    '"openblas configuration":"OpenBLAS 0.3.31.dev DYNAMIC_ARCH NO_AFFINITY '
    'SkylakeX MAX_THREADS=64","version":"0.3.31.dev"}'
)
_SCIPY_LAPACK_313 = (
    '{"has ilp64":false,"name":"scipy-openblas","openblas configuration":'
    '"OpenBLAS 0.3.31.dev DYNAMIC_ARCH NO_AFFINITY SkylakeX MAX_THREADS=64",'
    '"version":"0.3.31.dev"}'
)
_SCIPY_BLAS_314 = _SCIPY_BLAS_313.replace("SkylakeX", "Haswell")
_SCIPY_LAPACK_314 = _SCIPY_LAPACK_313.replace("SkylakeX", "Haswell")
_OPENBLAS_RUNTIME = (
    '[["numpy.libs/libscipy_openblas64_-017048f4.so",'
    '"OpenBLAS 0.3.33.112.0  USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell '
    'MAX_THREADS=64","Haswell",1],'
    '["scipy.libs/libscipy_openblas-5f890258.so",'
    '"OpenBLAS 0.3.31.dev DYNAMIC_ARCH NO_AFFINITY Haswell MAX_THREADS=64",'
    '"Haswell",1]]'
)
_ISA = "x86-64-v3"
_FLOATING_POINT_MODE = "binary64-round-nearest-gradual-underflow"
_NUMPY_DISPATCH_RESTRICTIONS = tuple(
    sorted(
        {
            "AVX512_ICL",
            "AVX512_SPR",
            "X86_V4",
        }
    )
)
NATIVE_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
_COMPATIBILITY_ALLOWED_FIELDS = {
    "image_digest",
    "loaded_library_manifest_hash",
    "numpy_installation_hash",
    "python_abi",
    "python_executable_hash",
    "python_source_hash",
    "python_version",
    "scipy_blas",
    "scipy_installation_hash",
    "scipy_lapack",
    "wheel_manifest_hash",
}
_PROVENANCE_KEYS = {
    "build_context_hash",
    "build_recipe_hash",
    "dependency_lock_hash",
    "os_package_manifest_hash",
    "platform_manifest",
    "python_source_hash",
    "uv_version",
    "uv_wheel_hash",
    "wheel_manifest_hash",
}
type LaneRole = Literal["CANONICAL_NUMERICAL", "PYTHON_COMPATIBILITY"]
type ComparisonScopeKind = Literal[
    "WITHIN_LANE_BITWISE", "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON"
]


class LaneAuthorityError(RuntimeError):
    """Raised when an execution cannot claim a lane's authority."""


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {
            "kind": kind,
            "record": record,
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
        },
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return _sha256(encoded)


def _strict_json_loads(content: str) -> object:
    """Decode JSON without silently normalizing duplicate or non-finite values."""

    def object_from_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
        record: dict[str, object] = {}
        for key, value in pairs:
            if key in record:
                raise ValueError(f"duplicate JSON member: {key}")
            record[key] = value
        return record

    def reject_constant(value: str) -> object:
        raise ValueError(f"non-finite JSON number: {value}")

    return json.loads(
        content,
        object_pairs_hook=object_from_pairs,
        parse_constant=reject_constant,
    )


def _digest(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != _HASH_LENGTH
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{name} must be a SHA-256 digest")
    return value


def _image_digest(value: object) -> str:
    if not isinstance(value, str) or not value.startswith("sha256:"):
        raise ValueError("image digest must use the sha256:<digest> form")
    _digest(value.removeprefix("sha256:"), "image digest")
    return value


def _non_empty(value: object, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a non-empty string")
    return value


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _ordered_strings(value: object, name: str) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        raise TypeError(f"{name} must be a sequence")
    result = tuple(_non_empty(item, name) for item in value)
    if not result or tuple(sorted(result)) != result or len(set(result)) != len(result):
        raise ValueError(f"{name} must be non-empty, unique, and canonically ordered")
    return result


def _exact_keys(record: Mapping[str, object], expected: set[str], name: str) -> None:
    if set(record) != expected:
        raise ValueError(f"{name} record keys are not canonical")


def _record_mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a record")
    return cast("Mapping[str, object]", value)


def _pinned_platform_manifest(value: object) -> str:
    text = _non_empty(value, "platform manifest")
    repository, separator, digest = text.rpartition("@sha256:")
    if not separator or not repository:
        raise ValueError("platform manifest must use an immutable SHA-256 digest")
    _digest(digest, "platform manifest digest")
    return text


def _recipe_hash() -> str:
    recipe = hashlib.sha256()
    package_files = files(__package__)
    for name in _RECIPE_NAMES:
        member_hash = _sha256(package_files.joinpath(name).read_bytes())
        recipe.update(
            f"{member_hash}  src/chemex/numerical_lanes/{name}\n".encode("ascii")
        )
    return recipe.hexdigest()


def _canonical_architecture(value: str) -> str:
    return {"x86_64": "amd64", "aarch64": "arm64"}.get(value, value)


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def _manifest_hash(entries: Sequence[str], kind: str) -> str:
    return _identity(kind, tuple(entries))


def _parse_file_manifest(content: bytes, name: str) -> tuple[tuple[str, str], ...]:
    try:
        lines = content.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise LaneAuthorityError(f"{name} is not ASCII") from error
    entries: list[tuple[str, str]] = []
    for line in lines:
        digest, separator, relative_name = line.partition("  ")
        path = PurePosixPath(relative_name)
        if (
            not separator
            or path.is_absolute()
            or not path.parts
            or any(part in {"", ".", ".."} for part in path.parts)
        ):
            raise LaneAuthorityError(f"{name} is malformed or path-dependent")
        try:
            entries.append((_digest(digest, f"{name} member hash"), path.as_posix()))
        except ValueError as error:
            raise LaneAuthorityError(f"{name} is malformed") from error
    if not entries or len({entry[1] for entry in entries}) != len(entries):
        raise LaneAuthorityError(f"{name} is empty or has duplicate members")
    return tuple(entries)


def _verify_file_manifest(
    manifest_path: Path,
    root: Path,
    actual_relative_names: Sequence[str],
    name: str,
) -> str:
    try:
        content = manifest_path.read_bytes()
    except OSError as error:
        raise LaneAuthorityError(f"{name} is unavailable") from error
    entries = _parse_file_manifest(content, name)
    frozen_names = tuple(sorted(relative_name for _, relative_name in entries))
    actual_names = tuple(sorted(actual_relative_names))
    if frozen_names != actual_names:
        raise LaneAuthorityError(f"{name} does not describe the running filesystem")
    for expected_hash, relative_name in entries:
        path = root.joinpath(*PurePosixPath(relative_name).parts)
        try:
            observed_hash = _sha256(path.read_bytes())
        except OSError as error:
            raise LaneAuthorityError(f"{name} member is unavailable") from error
        if observed_hash != expected_hash:
            raise LaneAuthorityError(f"{name} member content does not match")
    return _sha256(content)


def _current_build_context_hash() -> str:
    relative_names = ["pyproject.toml", "uv.lock"]
    try:
        relative_names.extend(
            path.relative_to(_BUILD_ROOT).as_posix()
            for path in (_BUILD_ROOT / "src/chemex/numerical_lanes").rglob("*")
            if path.is_file()
            and path.suffix != ".pyc"
            and not (
                path.parent.name == "manifests"
                and path.parent.parent.name == "numerical_lanes"
                and path.suffix == ".json"
            )
        )
    except OSError as error:
        raise LaneAuthorityError("build context is unavailable") from error
    return _verify_file_manifest(
        _BUILD_CONTEXT_MANIFEST_PATH,
        _BUILD_ROOT,
        relative_names,
        "build-context manifest",
    )


def _current_wheel_manifest_hash() -> str:
    try:
        names = tuple(path.name for path in _WHEELHOUSE_PATH.glob("*.whl"))
    except OSError as error:
        raise LaneAuthorityError("wheelhouse is unavailable") from error
    return _verify_file_manifest(
        _WHEEL_MANIFEST_PATH,
        _WHEELHOUSE_PATH,
        names,
        "wheel manifest",
    )


def _current_os_package_manifest_hash() -> str:
    environment = dict(os.environ)
    environment["LC_ALL"] = "C"
    try:
        result = subprocess.run(
            [
                "/usr/bin/dpkg-query",
                "--show",
                "--showformat=${Package}\\t${Version}\\t${Architecture}\\n",
            ],
            check=True,
            capture_output=True,
            env=environment,
        )
        frozen = _OS_PACKAGE_MANIFEST_PATH.read_bytes()
    except (OSError, subprocess.CalledProcessError) as error:
        raise LaneAuthorityError("OS package observation is unavailable") from error
    try:
        observed_lines = sorted(result.stdout.decode("utf-8").splitlines())
    except UnicodeDecodeError as error:
        raise LaneAuthorityError("OS package observation is not UTF-8") from error
    observed = ("\n".join(observed_lines) + "\n").encode("utf-8")
    if frozen != observed:
        raise LaneAuthorityError("OS package manifest does not match the running image")
    return _sha256(frozen)


def _read_build_provenance(path: Path) -> dict[str, str]:
    try:
        raw = _strict_json_loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise LaneAuthorityError("build provenance is unavailable") from error
    if not isinstance(raw, Mapping) or set(raw) != _PROVENANCE_KEYS:
        raise LaneAuthorityError("build provenance is malformed or noncanonical")
    record = cast("Mapping[str, object]", raw)
    try:
        provenance = {
            "build_context_hash": _digest(
                record["build_context_hash"], "build context hash"
            ),
            "build_recipe_hash": _digest(
                record["build_recipe_hash"], "build recipe hash"
            ),
            "dependency_lock_hash": _digest(
                record["dependency_lock_hash"], "dependency lock hash"
            ),
            "os_package_manifest_hash": _digest(
                record["os_package_manifest_hash"], "OS package manifest hash"
            ),
            "platform_manifest": _pinned_platform_manifest(record["platform_manifest"]),
            "python_source_hash": _digest(
                record["python_source_hash"], "Python source hash"
            ),
            "uv_version": _non_empty(record["uv_version"], "uv version"),
            "uv_wheel_hash": _digest(record["uv_wheel_hash"], "uv wheel hash"),
            "wheel_manifest_hash": _digest(
                record["wheel_manifest_hash"], "wheel manifest hash"
            ),
        }
    except (KeyError, TypeError, ValueError) as error:
        raise LaneAuthorityError(
            "build provenance is malformed or noncanonical"
        ) from error
    try:
        dependency_lock_hash = _sha256((_BUILD_ROOT / "uv.lock").read_bytes())
    except OSError as error:
        raise LaneAuthorityError("dependency lock is unavailable") from error
    observed = {
        "build_context_hash": _current_build_context_hash(),
        "build_recipe_hash": _recipe_hash(),
        "dependency_lock_hash": dependency_lock_hash,
        "os_package_manifest_hash": _current_os_package_manifest_hash(),
        "wheel_manifest_hash": _current_wheel_manifest_hash(),
    }
    mismatches = sorted(
        name for name, value in observed.items() if provenance[name] != value
    )
    if mismatches:
        raise LaneAuthorityError(
            "build provenance does not match the running filesystem: "
            + ", ".join(mismatches)
        )
    return provenance


def _installed_package_hash(package: str) -> str:
    distribution = importlib.metadata.distribution(package)
    package_files = distribution.files
    if not package_files:
        raise LaneAuthorityError(f"installed {package} file manifest is unavailable")
    members: list[tuple[str, str, int]] = []
    for package_file in sorted(package_files, key=str):
        location = Path(package_file.locate())
        if location.is_file():
            content = location.read_bytes()
            members.append(
                (
                    PurePosixPath(str(package_file)).as_posix(),
                    _sha256(content),
                    len(content),
                )
            )
    if not members:
        raise LaneAuthorityError(f"installed {package} content is unavailable")
    return _identity("installed-package-content", tuple(members))


def _dependency_identity(module: Any, dependency: str) -> str:
    try:
        config = module.show_config(mode="dicts")
        data = config["Build Dependencies"][dependency]
    except (AttributeError, KeyError, TypeError) as error:
        raise LaneAuthorityError(
            f"{module.__name__} {dependency} build identity is unavailable"
        ) from error
    if not isinstance(data, Mapping) or data.get("found") is not True:
        raise LaneAuthorityError(
            f"{module.__name__} {dependency} build identity is unavailable"
        )
    semantic_keys = (
        "cython blas ilp64",
        "has ilp64",
        "name",
        "openblas configuration",
        "version",
    )
    semantic = {key: data[key] for key in semantic_keys if key in data}
    if not isinstance(semantic.get("name"), str) or not isinstance(
        semantic.get("version"), str
    ):
        raise LaneAuthorityError(
            f"{module.__name__} {dependency} build identity is ambiguous"
        )
    return _canonical_json(semantic)


def _load_numerical_stack() -> tuple[Any, Any]:
    import numpy
    import scipy
    import scipy.linalg
    import scipy.optimize

    left = numpy.eye(2, dtype=numpy.float64)
    numpy.dot(left, left)
    scipy.linalg.solve(left, numpy.ones(2), assume_a="gen")
    return numpy, scipy


def _loaded_library_paths() -> tuple[Path, ...]:
    try:
        lines = Path("/proc/self/maps").read_text(encoding="ascii").splitlines()
    except (OSError, UnicodeDecodeError) as error:
        raise LaneAuthorityError("loaded-library observation is unavailable") from error
    paths: set[Path] = set()
    for line in lines:
        candidate = line.rsplit(maxsplit=1)[-1]
        if not candidate.startswith("/"):
            continue
        path = Path(candidate.removesuffix(" (deleted)"))
        path_text = path.as_posix().lower()
        if path.is_file() and (
            "/numpy/" in path_text
            or "/scipy/" in path_text
            or "/numpy.libs/" in path_text
            or "/scipy.libs/" in path_text
            or "openblas" in path_text
        ):
            paths.add(path)
    if not paths:
        raise LaneAuthorityError("loaded numerical libraries are unavailable")
    return tuple(sorted(paths))


def _library_key(path: Path) -> str:
    parts = path.parts
    if "site-packages" in parts:
        index = parts.index("site-packages")
        return PurePosixPath(*parts[index + 1 :]).as_posix()
    return path.name


def _loaded_library_manifest(paths: Sequence[Path]) -> str:
    entries = tuple(
        sorted(
            f"{_library_key(path)}@sha256:{_sha256(path.read_bytes())}"
            for path in paths
        )
    )
    if len({entry.split("@sha256:", 1)[0] for entry in entries}) != len(entries):
        raise LaneAuthorityError("loaded numerical library names are ambiguous")
    return _manifest_hash(entries, "loaded-numerical-libraries")


def _openblas_runtime(paths: Sequence[Path]) -> tuple[str, str, int]:
    candidates = [path for path in paths if "openblas" in path.name.lower()]
    if not candidates:
        raise LaneAuthorityError("loaded OpenBLAS libraries are unavailable")
    symbol_families = (
        (
            "openblas_get_config",
            "openblas_get_corename",
            "openblas_get_num_threads",
        ),
        (
            "scipy_openblas_get_config",
            "scipy_openblas_get_corename",
            "scipy_openblas_get_num_threads",
        ),
        (
            "scipy_openblas_get_config64_",
            "scipy_openblas_get_corename64_",
            "scipy_openblas_get_num_threads64_",
        ),
    )
    observed: list[tuple[str, str, str, int]] = []
    for path in candidates:
        try:
            library = ctypes.CDLL(str(path))
        except OSError as error:
            raise LaneAuthorityError(
                "loaded OpenBLAS library cannot be inspected"
            ) from error
        for config_symbol, core_symbol, threads_symbol in symbol_families:
            try:
                get_config = getattr(library, config_symbol)
                get_core = getattr(library, core_symbol)
                get_threads = getattr(library, threads_symbol)
            except AttributeError:
                continue
            get_config.restype = ctypes.c_char_p
            get_core.restype = ctypes.c_char_p
            get_threads.restype = ctypes.c_int
            config_value = get_config()
            core_value = get_core()
            if config_value and core_value:
                observed.append(
                    (
                        _library_key(path),
                        config_value.decode("ascii"),
                        core_value.decode("ascii"),
                        int(get_threads()),
                    )
                )
                break
        else:
            raise LaneAuthorityError(
                f"OpenBLAS runtime introspection unavailable for {_library_key(path)}"
            )
    cores = {item[2] for item in observed}
    thread_counts = {item[3] for item in observed}
    if len(cores) != 1 or len(thread_counts) != 1:
        raise LaneAuthorityError("loaded OpenBLAS runtimes disagree")
    return (
        _canonical_json(sorted(observed)),
        cores.pop(),
        thread_counts.pop(),
    )


def _current_numpy_dispatch(numpy: Any) -> tuple[str, str]:
    value = os.environ.get("NPY_DISABLE_CPU_FEATURES")
    if value is None:
        raise LaneAuthorityError("NumPy dispatch restrictions are unavailable")
    restrictions = tuple(
        sorted(part.strip() for part in value.split(",") if part.strip())
    )
    if len(set(restrictions)) != len(restrictions):
        raise LaneAuthorityError("NumPy dispatch restrictions are ambiguous")
    try:
        multiarray = numpy._core._multiarray_umath
        features = multiarray.__cpu_features__
        dispatch = multiarray.__cpu_dispatch__
    except AttributeError as error:
        raise LaneAuthorityError(
            "NumPy effective CPU dispatch is unavailable"
        ) from error
    if not isinstance(dispatch, list) or not all(
        isinstance(name, str) for name in dispatch
    ):
        raise LaneAuthorityError("NumPy CPU dispatch target list is ambiguous")
    enabled_dispatch = sorted(name for name in dispatch if features.get(name) is True)
    if enabled_dispatch != ["X86_V3"]:
        raise LaneAuthorityError("effective NumPy ISA is not x86-64-v3")
    enabled_restricted = sorted(
        name for name in restrictions if features.get(name) is True
    )
    if enabled_restricted:
        raise LaneAuthorityError(
            "effective NumPy dispatch still enables restricted features: "
            + ", ".join(enabled_restricted)
        )
    return ",".join(restrictions), _ISA


def _current_workers() -> int:
    value = os.environ.get("CHEMEX_NUMERICAL_LANE_WORKERS")
    try:
        workers = int(value) if value is not None else 0
    except ValueError as error:
        raise LaneAuthorityError("worker policy is unavailable") from error
    if workers < 1:
        raise LaneAuthorityError("worker policy is unavailable")
    return workers


def _validate_thread_environment(actual_threads: int) -> None:
    values = tuple(os.environ.get(name) for name in NATIVE_THREAD_ENV_VARS)
    if any(value is None for value in values):
        raise LaneAuthorityError("native-thread controls are unavailable")
    try:
        declared = {int(value) for value in values if value is not None}
    except ValueError as error:
        raise LaneAuthorityError("native-thread controls are unavailable") from error
    if declared != {actual_threads}:
        raise LaneAuthorityError("declared and actual native-thread counts differ")


def _current_floating_point_mode() -> str:
    try:
        process = ctypes.CDLL(None)
        fegetround = process.fegetround
        fegetround.restype = ctypes.c_int
        helper = ctypes.CDLL(str(_FPSTATE_PATH))
        get_mxcsr = helper.chemex_get_mxcsr
        get_x87_control = helper.chemex_get_x87_control_word
        get_mxcsr.restype = ctypes.c_uint
        get_x87_control.restype = ctypes.c_ushort
        rounding_mode = int(fegetround())
        mxcsr_controls = int(get_mxcsr()) & 0xFFC0
        x87_controls = int(get_x87_control()) & 0x1F3F
    except (AttributeError, OSError) as error:
        raise LaneAuthorityError(
            "hardware floating-point state is unavailable"
        ) from error
    if not (
        rounding_mode == 0
        and mxcsr_controls == 0x1F80
        and x87_controls == 0x033F
        and sys.float_info.radix == 2
        and struct.pack(">d", 1.0) == b"?\xf0\x00\x00\x00\x00\x00\x00"
    ):
        raise LaneAuthorityError("hardware floating-point state is not canonical")
    return _FLOATING_POINT_MODE


@dataclass(frozen=True, slots=True)
class LaneSemantics:
    """Complete numerical semantics shared by a lane and its observed process."""

    python_implementation: str
    python_version: str
    python_abi: str
    python_source_hash: str
    python_executable_hash: str
    platform: str
    platform_manifest: str
    dependency_lock_hash: str
    build_recipe_hash: str
    build_context_hash: str
    uv_version: str
    uv_wheel_hash: str
    wheel_manifest_hash: str
    os_package_manifest_hash: str
    image_digest: str
    numpy_version: str
    numpy_installation_hash: str
    scipy_version: str
    scipy_installation_hash: str
    numpy_blas: str
    numpy_lapack: str
    scipy_blas: str
    scipy_lapack: str
    openblas_configuration: str
    openblas_core: str
    loaded_library_manifest_hash: str
    numpy_dispatch_restrictions: str
    isa: str
    workers: int
    native_threads: int
    floating_point_mode: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        for attribute in fields(self):
            if attribute.name == "identity":
                continue
            value = getattr(self, attribute.name)
            if attribute.name in {"workers", "native_threads"}:
                _positive_int(value, attribute.name)
            elif attribute.name == "platform_manifest":
                _pinned_platform_manifest(value)
            elif attribute.name == "image_digest":
                _image_digest(value)
            elif attribute.name.endswith("_hash"):
                _digest(value, attribute.name)
            else:
                _non_empty(value, attribute.name)
        object.__setattr__(
            self,
            "identity",
            _identity("numerical-lane-semantics", self._identity_values()),
        )

    def _identity_values(self) -> tuple[object, ...]:
        return tuple(
            getattr(self, attribute.name)
            for attribute in fields(self)
            if attribute.name != "identity"
        )

    def to_record(self) -> dict[str, object]:
        return {
            attribute.name: getattr(self, attribute.name) for attribute in fields(self)
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> LaneSemantics:
        names = {attribute.name for attribute in fields(cls)}
        _exact_keys(record, names, "Lane semantics")
        semantics = cls(
            python_implementation=_non_empty(
                record["python_implementation"], "Python implementation"
            ),
            python_version=_non_empty(record["python_version"], "Python version"),
            python_abi=_non_empty(record["python_abi"], "Python ABI"),
            python_source_hash=_digest(
                record["python_source_hash"], "Python source hash"
            ),
            python_executable_hash=_digest(
                record["python_executable_hash"], "Python executable hash"
            ),
            platform=_non_empty(record["platform"], "platform"),
            platform_manifest=_pinned_platform_manifest(record["platform_manifest"]),
            dependency_lock_hash=_digest(
                record["dependency_lock_hash"], "dependency lock hash"
            ),
            build_recipe_hash=_digest(record["build_recipe_hash"], "build recipe hash"),
            build_context_hash=_digest(
                record["build_context_hash"], "build context hash"
            ),
            uv_version=_non_empty(record["uv_version"], "uv version"),
            uv_wheel_hash=_digest(record["uv_wheel_hash"], "uv wheel hash"),
            wheel_manifest_hash=_digest(
                record["wheel_manifest_hash"], "wheel manifest hash"
            ),
            os_package_manifest_hash=_digest(
                record["os_package_manifest_hash"], "OS package manifest hash"
            ),
            image_digest=_image_digest(record["image_digest"]),
            numpy_version=_non_empty(record["numpy_version"], "NumPy version"),
            numpy_installation_hash=_digest(
                record["numpy_installation_hash"], "NumPy installation hash"
            ),
            scipy_version=_non_empty(record["scipy_version"], "SciPy version"),
            scipy_installation_hash=_digest(
                record["scipy_installation_hash"], "SciPy installation hash"
            ),
            numpy_blas=_non_empty(record["numpy_blas"], "NumPy BLAS"),
            numpy_lapack=_non_empty(record["numpy_lapack"], "NumPy LAPACK"),
            scipy_blas=_non_empty(record["scipy_blas"], "SciPy BLAS"),
            scipy_lapack=_non_empty(record["scipy_lapack"], "SciPy LAPACK"),
            openblas_configuration=_non_empty(
                record["openblas_configuration"], "OpenBLAS configuration"
            ),
            openblas_core=_non_empty(record["openblas_core"], "OpenBLAS core"),
            loaded_library_manifest_hash=_digest(
                record["loaded_library_manifest_hash"],
                "loaded library manifest hash",
            ),
            numpy_dispatch_restrictions=_non_empty(
                record["numpy_dispatch_restrictions"],
                "NumPy dispatch restrictions",
            ),
            isa=_non_empty(record["isa"], "ISA"),
            workers=_positive_int(record["workers"], "workers"),
            native_threads=_positive_int(record["native_threads"], "native threads"),
            floating_point_mode=_non_empty(
                record["floating_point_mode"], "floating-point mode"
            ),
        )
        if record["identity"] != semantics.identity:
            raise ValueError("Lane semantics identity does not match payload")
        return semantics


@dataclass(frozen=True, slots=True)
class RuntimeEnvironment:
    """Facts measured from the actual already-imported numerical process."""

    semantics: LaneSemantics
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.semantics, LaneSemantics):
            raise TypeError("runtime environment requires canonical lane semantics")
        object.__setattr__(
            self,
            "identity",
            _identity("post-import-runtime-environment", self.semantics.identity),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "identity": self.identity,
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "semantics": self.semantics.to_record(),
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> RuntimeEnvironment:
        _exact_keys(
            record,
            {"identity", "schema_version", "semantic_version", "semantics"},
            "Runtime environment",
        )
        if (
            record["schema_version"] != _SCHEMA_VERSION
            or record["semantic_version"] != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported runtime-environment semantics")
        environment = cls(
            LaneSemantics.from_record(_record_mapping(record["semantics"], "semantics"))
        )
        if record["identity"] != environment.identity:
            raise ValueError("Runtime-environment identity does not match payload")
        return environment

    @classmethod
    def from_current_process(
        cls,
        image_digest: str,
        provenance_path: Path = _PROVENANCE_PATH,
    ) -> RuntimeEnvironment:
        """Observe required facts after loading the numerical stack; fail closed."""

        provenance = _read_build_provenance(provenance_path)
        numpy, scipy = _load_numerical_stack()
        library_paths = _loaded_library_paths()
        openblas_configuration, openblas_core, native_threads = _openblas_runtime(
            library_paths
        )
        _validate_thread_environment(native_threads)
        numpy_dispatch_restrictions, isa = _current_numpy_dispatch(numpy)
        abi = sysconfig.get_config_var("SOABI")
        if not isinstance(abi, str) or not abi:
            raise LaneAuthorityError("Python ABI is unavailable")
        semantics = LaneSemantics(
            python_implementation=platform_module.python_implementation(),
            python_version=".".join(map(str, sys.version_info[:3])),
            python_abi=abi,
            python_source_hash=provenance["python_source_hash"],
            python_executable_hash=_sha256(Path(sys.executable).read_bytes()),
            platform=(
                f"{sys.platform}/"
                f"{_canonical_architecture(platform_module.machine().lower())}"
            ),
            platform_manifest=provenance["platform_manifest"],
            dependency_lock_hash=provenance["dependency_lock_hash"],
            build_recipe_hash=provenance["build_recipe_hash"],
            build_context_hash=provenance["build_context_hash"],
            uv_version=provenance["uv_version"],
            uv_wheel_hash=provenance["uv_wheel_hash"],
            wheel_manifest_hash=provenance["wheel_manifest_hash"],
            os_package_manifest_hash=provenance["os_package_manifest_hash"],
            image_digest=_image_digest(image_digest),
            numpy_version=numpy.__version__,
            numpy_installation_hash=_installed_package_hash("numpy"),
            scipy_version=scipy.__version__,
            scipy_installation_hash=_installed_package_hash("scipy"),
            numpy_blas=_dependency_identity(numpy, "blas"),
            numpy_lapack=_dependency_identity(numpy, "lapack"),
            scipy_blas=_dependency_identity(scipy, "blas"),
            scipy_lapack=_dependency_identity(scipy, "lapack"),
            openblas_configuration=openblas_configuration,
            openblas_core=openblas_core,
            loaded_library_manifest_hash=_loaded_library_manifest(library_paths),
            numpy_dispatch_restrictions=numpy_dispatch_restrictions,
            isa=isa,
            workers=_current_workers(),
            native_threads=native_threads,
            floating_point_mode=_current_floating_point_mode(),
        )
        return cls(semantics)


@dataclass(frozen=True, slots=True)
class LaneAttestation:
    """Canonical evidence that one observed process matched one lane exactly.

    Parsing or constructing this record never grants live lane authority.  Only
    :meth:`NumericalLane.attest_current_process` can return that process-local
    capability.
    """

    lane_identity: str
    environment_identity: str
    workers: int
    native_threads: int
    method: Literal["POST_IMPORT_CURRENT_PROCESS"]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _digest(self.lane_identity, "lane identity")
        _digest(self.environment_identity, "environment identity")
        _positive_int(self.workers, "attested workers")
        _positive_int(self.native_threads, "attested native threads")
        if self.method != "POST_IMPORT_CURRENT_PROCESS":
            raise ValueError("Unknown lane-attestation method")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "numerical-lane-attestation",
                (
                    self.lane_identity,
                    self.environment_identity,
                    self.workers,
                    self.native_threads,
                    self.method,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "environment_identity": self.environment_identity,
            "identity": self.identity,
            "lane_identity": self.lane_identity,
            "method": self.method,
            "native_threads": self.native_threads,
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "workers": self.workers,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> LaneAttestation:
        _exact_keys(
            record,
            {
                "environment_identity",
                "identity",
                "lane_identity",
                "method",
                "native_threads",
                "schema_version",
                "semantic_version",
                "workers",
            },
            "Lane attestation",
        )
        if (
            record["schema_version"] != _SCHEMA_VERSION
            or record["semantic_version"] != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported lane-attestation semantics")
        method = _non_empty(record["method"], "lane-attestation method")
        if method != "POST_IMPORT_CURRENT_PROCESS":
            raise ValueError("Unknown lane-attestation method")
        record_identity = _digest(record["identity"], "lane-attestation identity")
        attestation = cls(
            _digest(record["lane_identity"], "lane identity"),
            _digest(record["environment_identity"], "environment identity"),
            _positive_int(record["workers"], "attested workers"),
            _positive_int(record["native_threads"], "attested native threads"),
            "POST_IMPORT_CURRENT_PROCESS",
        )
        if record_identity != attestation.identity:
            raise ValueError("Lane-attestation identity does not match payload")
        return attestation


class _LiveLaneBinding(NamedTuple):
    """Immutable process-owned facts bound to one opaque live capability."""

    lane_identity: str
    lane_role: LaneRole
    attestation_identity: str
    environment_identity: str
    workers: int
    native_threads: int


type _LiveLaneAuthorityFact = Literal[
    "lane_identity",
    "lane_role",
    "attestation_identity",
    "environment_identity",
    "workers",
    "native_threads",
]


class _LiveLaneAuthorityMismatch(LaneAuthorityError):
    def __init__(self, fact: _LiveLaneAuthorityFact) -> None:
        self.fact = fact
        super().__init__(f"Live lane authority does not match required {fact}")


class LiveLaneAuthority:
    """Non-serializable capability owned by one successfully probed process."""

    __slots__ = ("__weakref__",)

    def __new__(cls) -> LiveLaneAuthority:
        raise TypeError("Live lane authority is minted only by current-process probing")

    def to_record(self) -> dict[str, object]:
        """Serialize evidence only; deserialization cannot recreate this capability."""

        return _validate_live_lane_authority(self).to_record()

    def __copy__(self) -> LiveLaneAuthority:
        raise TypeError("Live lane authority cannot be copied")

    def __deepcopy__(self, memo: object) -> LiveLaneAuthority:
        _ = memo
        raise TypeError("Live lane authority cannot be copied")

    def __reduce_ex__(self, protocol: SupportsIndex) -> str | tuple[object, ...]:
        _ = protocol
        raise TypeError("Live lane authority cannot be pickled")


_LIVE_LANE_BINDINGS: weakref.WeakKeyDictionary[LiveLaneAuthority, _LiveLaneBinding] = (
    weakref.WeakKeyDictionary()
)


def _validate_live_lane_authority(
    authority: LiveLaneAuthority,
    *,
    required_lane_identity: str | None = None,
    required_lane_role: LaneRole | None = None,
    required_attestation_identity: str | None = None,
    required_environment_identity: str | None = None,
    required_workers: int | None = None,
    required_native_threads: int | None = None,
) -> LaneAttestation:
    """Validate one exact live token without exposing its registry binding."""

    if not isinstance(authority, LiveLaneAuthority):
        raise TypeError("Live lane validation requires current-process authority")
    try:
        binding = _LIVE_LANE_BINDINGS[authority]
    except KeyError as error:
        raise LaneAuthorityError("Live lane authority is not process-owned") from error
    required_facts = (
        ("lane_role", required_lane_role, binding.lane_role),
        ("lane_identity", required_lane_identity, binding.lane_identity),
        (
            "attestation_identity",
            required_attestation_identity,
            binding.attestation_identity,
        ),
        (
            "environment_identity",
            required_environment_identity,
            binding.environment_identity,
        ),
        ("workers", required_workers, binding.workers),
        ("native_threads", required_native_threads, binding.native_threads),
    )
    for fact, required, observed in required_facts:
        if required is not None and required != observed:
            raise _LiveLaneAuthorityMismatch(fact)
    evidence = LaneAttestation(
        binding.lane_identity,
        binding.environment_identity,
        binding.workers,
        binding.native_threads,
        "POST_IMPORT_CURRENT_PROCESS",
    )
    if evidence.identity != binding.attestation_identity:
        raise LaneAuthorityError(
            "Live lane authority binding is internally inconsistent"
        )
    return evidence


@dataclass(frozen=True, slots=True)
class NumericalLane:
    """One complete numerical environment and its replay promise."""

    name: str
    role: LaneRole
    semantics: LaneSemantics
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _non_empty(self.name, "lane name")
        if self.role not in {"CANONICAL_NUMERICAL", "PYTHON_COMPATIBILITY"}:
            raise ValueError("Unknown numerical-lane role")
        if not isinstance(self.semantics, LaneSemantics):
            raise TypeError("Numerical lane requires canonical semantics")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "numerical-lane",
                (self.name, self.role, self.semantics.identity),
            ),
        )

    @classmethod
    def capture_current_process(
        cls,
        name: str,
        role: LaneRole,
        image_digest: str,
        provenance_path: Path = _PROVENANCE_PATH,
    ) -> NumericalLane:
        """Capture a candidate matching an explicit historical or prospective lane."""

        environment = RuntimeEnvironment.from_current_process(
            image_digest, provenance_path
        )
        cls._validate_canonical_contract(name, role, environment.semantics)
        return cls(name, role, environment.semantics)

    @staticmethod
    def _validate_canonical_contract(
        name: str, role: LaneRole, semantics: LaneSemantics
    ) -> None:
        historical_name = _HISTORICAL_LANE_NAMES[role]
        prospective_name = _PROSPECTIVE_LANE_NAMES[role]
        if name == historical_name:
            expected_recipe_hash = _HISTORICAL_RECIPE_HASH
            expected_lockfile_hash = _LOCKFILE_HASH
        elif name == prospective_name:
            expected_recipe_hash = _recipe_hash()
            expected_lockfile_hash = _PROSPECTIVE_LOCKFILE_HASH
        else:
            raise LaneAuthorityError(
                f"Numerical lane contract rejected: unexpected {role} lane name {name}"
            )
        expected_version = "3.13.5" if role == "CANONICAL_NUMERICAL" else "3.14.5"
        expected_source = (
            "e6190f52699b534ee203d9f417bdbca05a92f23e35c19c691a50ed2942835385"
            if role == "CANONICAL_NUMERICAL"
            else "9c22bfe9939a6c5418fc74b289a5f1cc41859ae82ac6b163016b5844bd0a86bc"
        )
        expected = {
            "build_recipe_hash": expected_recipe_hash,
            "dependency_lock_hash": expected_lockfile_hash,
            "floating_point_mode": _FLOATING_POINT_MODE,
            "isa": _ISA,
            "native_threads": 1,
            "numpy_blas": _NUMPY_BLAS,
            "numpy_dispatch_restrictions": ",".join(_NUMPY_DISPATCH_RESTRICTIONS),
            "numpy_lapack": _NUMPY_BLAS,
            "numpy_version": "2.5.1",
            "openblas_configuration": _OPENBLAS_RUNTIME,
            "openblas_core": _OPENBLAS_CORE,
            "platform": "linux/amd64",
            "platform_manifest": _PLATFORM_MANIFEST,
            "python_implementation": "CPython",
            "python_abi": (
                "cpython-313-x86_64-linux-gnu"
                if role == "CANONICAL_NUMERICAL"
                else "cpython-314-x86_64-linux-gnu"
            ),
            "python_source_hash": expected_source,
            "python_version": expected_version,
            "scipy_version": "1.18.0",
            "scipy_blas": (
                _SCIPY_BLAS_313 if role == "CANONICAL_NUMERICAL" else _SCIPY_BLAS_314
            ),
            "scipy_lapack": (
                _SCIPY_LAPACK_313
                if role == "CANONICAL_NUMERICAL"
                else _SCIPY_LAPACK_314
            ),
            "uv_version": _UV_VERSION,
            "uv_wheel_hash": _UV_WHEEL_HASH,
            "workers": 1,
        }
        mismatches = [
            f"{name}: expected {expected_value}, got {getattr(semantics, name)}"
            for name, expected_value in expected.items()
            if getattr(semantics, name) != expected_value
        ]
        if mismatches:
            raise LaneAuthorityError(
                "Numerical lane contract rejected: " + "; ".join(mismatches)
            )

    def attest_current_process(
        self,
        image_digest: str,
        provenance_path: Path = _PROVENANCE_PATH,
    ) -> LiveLaneAuthority:
        """Mint authority only after an exact current-process observation."""

        environment = RuntimeEnvironment.from_current_process(
            image_digest, provenance_path
        )
        if environment.semantics != self.semantics:
            mismatches = [
                attribute.name
                for attribute in fields(LaneSemantics)
                if attribute.name != "identity"
                and getattr(environment.semantics, attribute.name)
                != getattr(self.semantics, attribute.name)
            ]
            raise LaneAuthorityError(
                "Numerical lane authority rejected; mismatched claims: "
                + ", ".join(mismatches)
            )
        evidence = LaneAttestation(
            self.identity,
            environment.identity,
            environment.semantics.workers,
            environment.semantics.native_threads,
            "POST_IMPORT_CURRENT_PROCESS",
        )
        authority = object.__new__(LiveLaneAuthority)
        _LIVE_LANE_BINDINGS[authority] = _LiveLaneBinding(
            lane_identity=self.identity,
            lane_role=self.role,
            attestation_identity=evidence.identity,
            environment_identity=environment.identity,
            workers=evidence.workers,
            native_threads=evidence.native_threads,
        )
        return authority

    def compatibility_delta(self, other: NumericalLane) -> dict[str, tuple[str, str]]:
        """Report every structural difference between the two lane records."""

        delta: dict[str, tuple[str, str]] = {}
        for attribute in fields(LaneSemantics):
            if attribute.name == "identity":
                continue
            left = str(getattr(self.semantics, attribute.name))
            right = str(getattr(other.semantics, attribute.name))
            if left != right:
                delta[attribute.name] = (left, right)
        return delta

    def validate_compatibility_lane(self, other: NumericalLane) -> None:
        """Reject non-Python structural drift between the 3.13 and 3.14 lanes."""

        if self.role != "CANONICAL_NUMERICAL" or other.role != "PYTHON_COMPATIBILITY":
            raise ValueError(
                "Compatibility validation requires canonical then compatibility"
            )
        unexpected = sorted(
            set(self.compatibility_delta(other)) - _COMPATIBILITY_ALLOWED_FIELDS
        )
        if unexpected:
            raise LaneAuthorityError(
                "Compatibility lane has non-Python structural differences: "
                + ", ".join(unexpected)
            )

    def to_record(self) -> dict[str, object]:
        return {
            "identity": self.identity,
            "name": self.name,
            "role": self.role,
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
            "semantics": self.semantics.to_record(),
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> NumericalLane:
        _exact_keys(
            record,
            {
                "identity",
                "name",
                "role",
                "schema_version",
                "semantic_version",
                "semantics",
            },
            "Numerical lane",
        )
        if (
            record["schema_version"] != _SCHEMA_VERSION
            or record["semantic_version"] != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported numerical-lane semantics")
        role = _non_empty(record["role"], "lane role")
        if role not in {"CANONICAL_NUMERICAL", "PYTHON_COMPATIBILITY"}:
            raise ValueError("Unknown numerical-lane role")
        lane = cls(
            _non_empty(record["name"], "lane name"),
            role,
            LaneSemantics.from_record(
                _record_mapping(record["semantics"], "lane semantics")
            ),
        )
        if record["identity"] != lane.identity:
            raise ValueError("Numerical-lane identity does not match payload")
        return lane


@dataclass(frozen=True, slots=True)
class ComparisonScope:
    """Qualified contract for replaying or comparing attested evidence."""

    kind: ComparisonScopeKind
    left_lane_identity: str
    right_lane_identity: str
    left_attestation_identity: str
    right_attestation_identity: str
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        if self.kind not in {
            "WITHIN_LANE_BITWISE",
            "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON",
        }:
            raise ValueError("Unknown comparison scope")
        for value, name in (
            (self.left_lane_identity, "left lane identity"),
            (self.right_lane_identity, "right lane identity"),
            (self.left_attestation_identity, "left attestation identity"),
            (self.right_attestation_identity, "right attestation identity"),
        ):
            _digest(value, name)
        same_lane = self.left_lane_identity == self.right_lane_identity
        if self.kind == "WITHIN_LANE_BITWISE" and not same_lane:
            raise ValueError("Within-lane replay requires one lane identity")
        if self.kind == "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON" and same_lane:
            raise ValueError("Cross-lane comparison requires distinct lane identities")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "numerical-comparison-scope",
                (
                    self.kind,
                    self.left_lane_identity,
                    self.right_lane_identity,
                    self.left_attestation_identity,
                    self.right_attestation_identity,
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "identity": self.identity,
            "kind": self.kind,
            "left_attestation_identity": self.left_attestation_identity,
            "left_lane_identity": self.left_lane_identity,
            "right_attestation_identity": self.right_attestation_identity,
            "right_lane_identity": self.right_lane_identity,
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> ComparisonScope:
        expected = {
            "identity",
            "kind",
            "left_attestation_identity",
            "left_lane_identity",
            "right_attestation_identity",
            "right_lane_identity",
            "schema_version",
            "semantic_version",
        }
        _exact_keys(record, expected, "Comparison scope")
        if (
            record["schema_version"] != _SCHEMA_VERSION
            or record["semantic_version"] != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported comparison-scope semantics")
        kind = _non_empty(record["kind"], "comparison scope")
        if kind not in {
            "WITHIN_LANE_BITWISE",
            "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON",
        }:
            raise ValueError("Unknown comparison scope")
        scope = cls(
            kind,
            _digest(record["left_lane_identity"], "left lane identity"),
            _digest(record["right_lane_identity"], "right lane identity"),
            _digest(record["left_attestation_identity"], "left attestation identity"),
            _digest(record["right_attestation_identity"], "right attestation identity"),
        )
        if record["identity"] != scope.identity:
            raise ValueError("Comparison-scope identity does not match payload")
        return scope


@dataclass(frozen=True, slots=True)
class CrossLaneNumericalPolicy:
    """Content-identified artifact-family policy for cross-lane comparisons."""

    artifact_family: str
    rtol: float
    atol: float
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _non_empty(self.artifact_family, "artifact family")
        if not (
            isinstance(self.rtol, float)
            and isinstance(self.atol, float)
            and math.isfinite(self.rtol)
            and math.isfinite(self.atol)
            and self.rtol >= 0.0
            and self.atol >= 0.0
        ):
            raise ValueError("cross-lane tolerances must be finite non-negative floats")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "cross-lane-numerical-policy",
                (self.artifact_family, self.rtol.hex(), self.atol.hex()),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "artifact_family": self.artifact_family,
            "atol": float(self.atol).hex(),
            "identity": self.identity,
            "rtol": float(self.rtol).hex(),
            "schema_version": _SCHEMA_VERSION,
            "semantic_version": _SEMANTIC_VERSION,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> CrossLaneNumericalPolicy:
        _exact_keys(
            record,
            {
                "artifact_family",
                "atol",
                "identity",
                "rtol",
                "schema_version",
                "semantic_version",
            },
            "Cross-lane policy",
        )
        if (
            record["schema_version"] != _SCHEMA_VERSION
            or record["semantic_version"] != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported cross-lane policy semantics")
        policy = cls(
            _non_empty(record["artifact_family"], "artifact family"),
            float.fromhex(_non_empty(record["rtol"], "relative tolerance")),
            float.fromhex(_non_empty(record["atol"], "absolute tolerance")),
        )
        if record["identity"] != policy.identity:
            raise ValueError("Cross-lane policy identity does not match payload")
        return policy


@dataclass(frozen=True, slots=True)
class ComparisonOutcome:
    """Content-identified result of one qualified comparison."""

    scope_identity: str
    policy_identity: str | None
    equivalent: bool
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        _digest(self.scope_identity, "comparison scope identity")
        if self.policy_identity is not None:
            _digest(self.policy_identity, "comparison policy identity")
        if not isinstance(self.equivalent, bool):
            raise TypeError("comparison outcome must be boolean")
        object.__setattr__(
            self,
            "identity",
            _identity(
                "numerical-comparison-outcome",
                (self.scope_identity, self.policy_identity, self.equivalent),
            ),
        )

    def to_record(self) -> dict[str, object]:
        return {
            "equivalent": self.equivalent,
            "identity": self.identity,
            "policy_identity": self.policy_identity,
            "schema_version": _SCHEMA_VERSION,
            "scope_identity": self.scope_identity,
            "semantic_version": _SEMANTIC_VERSION,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> ComparisonOutcome:
        _exact_keys(
            record,
            {
                "equivalent",
                "identity",
                "policy_identity",
                "schema_version",
                "scope_identity",
                "semantic_version",
            },
            "Comparison outcome",
        )
        if (
            record["schema_version"] != _SCHEMA_VERSION
            or record["semantic_version"] != _SEMANTIC_VERSION
        ):
            raise ValueError("Unsupported comparison-outcome semantics")
        policy_value = record["policy_identity"]
        if policy_value is not None:
            policy_value = _digest(policy_value, "comparison policy identity")
        equivalent = record["equivalent"]
        if not isinstance(equivalent, bool):
            raise TypeError("comparison outcome must be boolean")
        outcome = cls(
            _digest(record["scope_identity"], "comparison scope identity"),
            policy_value,
            equivalent,
        )
        if record["identity"] != outcome.identity:
            raise ValueError("Comparison-outcome identity does not match payload")
        return outcome


def compare_values(
    scope: ComparisonScope,
    left: Sequence[float],
    right: Sequence[float],
    *,
    policy: CrossLaneNumericalPolicy | None = None,
) -> ComparisonOutcome:
    """Compare equal-structure binary64 artifacts under one qualified scope."""

    import numpy

    left_values = numpy.asarray(left, dtype=numpy.float64)
    right_values = numpy.asarray(right, dtype=numpy.float64)
    if left_values.shape != right_values.shape:
        raise ValueError("Numerical comparisons require exactly matching structure")
    if not (numpy.isfinite(left_values).all() and numpy.isfinite(right_values).all()):
        raise ValueError("Numerical comparisons require finite values")
    if scope.kind == "WITHIN_LANE_BITWISE":
        if policy is not None:
            raise ValueError("Within-lane replay does not permit a tolerance policy")
        equivalent = bool(
            numpy.array_equal(
                left_values.view(numpy.uint64), right_values.view(numpy.uint64)
            )
        )
        return ComparisonOutcome(scope.identity, None, equivalent)
    if policy is None:
        raise ValueError("Cross-lane comparison requires an explicit numerical policy")
    difference = numpy.abs(left_values - right_values)
    scale = numpy.maximum(numpy.abs(left_values), numpy.abs(right_values))
    equivalent = bool(numpy.all(difference <= policy.atol + policy.rtol * scale))
    return ComparisonOutcome(scope.identity, policy.identity, equivalent)


def comparison_scope(
    left: LiveLaneAuthority, right: LiveLaneAuthority
) -> ComparisonScope:
    """Create a scope only from two qualified post-import attestations."""

    if not isinstance(left, LiveLaneAuthority) or not isinstance(
        right, LiveLaneAuthority
    ):
        raise TypeError("Comparison scope requires live current-process lane authority")
    left_evidence = _validate_live_lane_authority(left)
    right_evidence = _validate_live_lane_authority(right)
    kind: ComparisonScopeKind = (
        "WITHIN_LANE_BITWISE"
        if left_evidence.lane_identity == right_evidence.lane_identity
        else "CROSS_LANE_NUMERICAL_ARTIFACT_COMPARISON"
    )
    return ComparisonScope(
        kind,
        left_evidence.lane_identity,
        right_evidence.lane_identity,
        left_evidence.identity,
        right_evidence.identity,
    )


def _load_lane_pair(
    names: tuple[str, str], manifest_directory: Path | None
) -> tuple[NumericalLane, NumericalLane]:
    directory = (
        manifest_directory
        if manifest_directory is not None
        else Path(str(files(__package__).joinpath("manifests")))
    )
    lanes: list[NumericalLane] = []
    for lane_name in names:
        filename = f"{lane_name}.json"
        try:
            raw = _strict_json_loads((directory / filename).read_text(encoding="ascii"))
        except (OSError, UnicodeDecodeError, ValueError) as error:
            raise LaneAuthorityError(
                f"frozen numerical-lane manifest is unavailable: {filename}"
            ) from error
        lane = NumericalLane.from_record(_record_mapping(raw, "numerical lane"))
        if lane.name != lane_name:
            raise LaneAuthorityError("frozen numerical-lane name is invalid")
        NumericalLane._validate_canonical_contract(lane.name, lane.role, lane.semantics)
        lanes.append(lane)
    if (
        lanes[0].role != "CANONICAL_NUMERICAL"
        or lanes[1].role != "PYTHON_COMPATIBILITY"
    ):
        raise LaneAuthorityError("frozen numerical-lane roles are invalid")
    lanes[0].validate_compatibility_lane(lanes[1])
    return lanes[0], lanes[1]


def canonical_lanes(
    manifest_directory: Path | None = None,
) -> tuple[NumericalLane, NumericalLane]:
    """Load the immutable historical #588 v1 lane manifests."""

    return _load_lane_pair(
        (
            _HISTORICAL_LANE_NAMES["CANONICAL_NUMERICAL"],
            _HISTORICAL_LANE_NAMES["PYTHON_COMPATIBILITY"],
        ),
        manifest_directory,
    )


def prospective_lanes(
    manifest_directory: Path | None = None,
) -> tuple[NumericalLane, NumericalLane]:
    """Load prospective v2 manifests without reinterpreting historical authority."""

    return _load_lane_pair(
        (
            _PROSPECTIVE_LANE_NAMES["CANONICAL_NUMERICAL"],
            _PROSPECTIVE_LANE_NAMES["PYTHON_COMPATIBILITY"],
        ),
        manifest_directory,
    )

"""Small software-environment record used by current scientific evidence."""

from __future__ import annotations

import hashlib
import json
import platform
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version

import numpy as np

from chemex import __version__


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"kind": kind, "schema_version": 1, "record": record},
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "unavailable"


@dataclass(frozen=True, slots=True)
class RuntimeEnvironment:
    """Resolved versions and native libraries relevant to numerical evidence."""

    chemex_version: str
    python_version: str
    python_implementation: str
    platform: str
    numpy_version: str
    scipy_version: str
    emcee_version: str
    numerical_libraries: tuple[tuple[str, str, str], ...]
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        libraries = tuple(
            tuple(str(value) for value in item) for item in self.numerical_libraries
        )
        if any(len(item) != 3 for item in libraries):
            raise ValueError("Numerical libraries require kind, name, and version")
        object.__setattr__(self, "numerical_libraries", libraries)
        object.__setattr__(
            self,
            "identity",
            _identity(
                "native-provenance-environment",
                (
                    self.chemex_version,
                    self.python_version,
                    self.python_implementation,
                    self.platform,
                    self.numpy_version,
                    self.scipy_version,
                    self.emcee_version,
                    libraries,
                ),
            ),
        )

    @classmethod
    def from_current_process(cls) -> RuntimeEnvironment:
        """Capture the current product and native numerical runtime."""
        build_dependencies = getattr(np.__config__, "CONFIG", {}).get(
            "Build Dependencies", {}
        )
        libraries = tuple(
            (
                kind,
                str(details.get("name", "unknown")),
                str(details.get("version", "unknown")),
            )
            for kind, details in sorted(build_dependencies.items())
            if isinstance(details, dict) and details.get("found", False)
        )
        return cls(
            chemex_version=__version__,
            python_version=platform.python_version(),
            python_implementation=platform.python_implementation(),
            platform=platform.platform(),
            numpy_version=_package_version("numpy"),
            scipy_version=_package_version("scipy"),
            emcee_version=_package_version("emcee"),
            numerical_libraries=libraries,
        )

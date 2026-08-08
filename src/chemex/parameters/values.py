"""Revisioned ChemEx-native central values for one analysis session."""

from __future__ import annotations

import json
import math
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from numbers import Real
from threading import RLock
from types import MappingProxyType
from uuid import uuid4

from chemex.parameters.sealed import SealedConfiguration


class AnalysisValuesCommitError(ValueError):
    """Base error for a rejected atomic central-value commit."""


class IncompatibleAnalysisValuesError(AnalysisValuesCommitError):
    """Raised when a freshness snapshot belongs to another analysis state."""


class StaleAnalysisValuesError(AnalysisValuesCommitError):
    """Raised when the analysis-wide revision has advanced."""


class InvalidAnalysisValuesCommitError(AnalysisValuesCommitError):
    """Raised when commit scope or candidate central values are invalid."""


@dataclass(frozen=True, slots=True)
class AnalysisValuesSnapshot(Mapping[str, float]):
    """Immutable central values and the identities needed for freshness checks."""

    occurrence_identity: str
    model_identity: str
    definitions_identity: str
    configuration_identity: str
    revision: int
    _items: tuple[tuple[str, float], ...]
    _index: Mapping[str, float] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        items = tuple((str(param_id), float(value)) for param_id, value in self._items)
        index = MappingProxyType(dict(items))
        if len(index) != len(items):
            msg = "Analysis values snapshot contains duplicate parameter IDs"
            raise ValueError(msg)
        if isinstance(self.revision, bool) or self.revision < 0:
            msg = "Analysis values snapshot revision must be non-negative"
            raise ValueError(msg)
        if any(not math.isfinite(value) for _param_id, value in items):
            msg = "Analysis values snapshot contains a non-finite value"
            raise ValueError(msg)
        object.__setattr__(self, "_items", items)
        object.__setattr__(self, "_index", index)

    def __iter__(self) -> Iterator[str]:
        return (param_id for param_id, _value in self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, key: str) -> float:
        return self._index[key]

    def to_json(self) -> str:
        """Serialize with stable field order and exact binary64 value tokens."""
        payload = {
            "schema_version": 1,
            "occurrence_identity": self.occurrence_identity,
            "model_identity": self.model_identity,
            "definitions_identity": self.definitions_identity,
            "configuration_identity": self.configuration_identity,
            "revision": self.revision,
            "values": [
                {"id": param_id, "value": float(value).hex()}
                for param_id, value in self._items
            ],
        }
        return json.dumps(payload, ensure_ascii=True, separators=(",", ":"))

    @classmethod
    def from_json(cls, encoded: str) -> AnalysisValuesSnapshot:
        """Restore a snapshot from its deterministic JSON representation."""
        payload = json.loads(encoded)
        if not isinstance(payload, dict) or payload.get("schema_version") != 1:
            msg = "Unsupported analysis values snapshot schema"
            raise ValueError(msg)
        values = payload.get("values")
        if not isinstance(values, list):
            msg = "Analysis values snapshot values must be a list"
            raise TypeError(msg)
        try:
            items = tuple(
                (entry["id"], float.fromhex(entry["value"])) for entry in values
            )
            occurrence_identity = payload["occurrence_identity"]
            model_identity = payload["model_identity"]
            definitions_identity = payload["definitions_identity"]
            configuration_identity = payload["configuration_identity"]
            revision = payload["revision"]
        except (KeyError, TypeError, ValueError) as error:
            msg = "Invalid analysis values snapshot payload"
            raise ValueError(msg) from error
        if (
            not all(
                isinstance(identity, str)
                for identity in (
                    occurrence_identity,
                    model_identity,
                    definitions_identity,
                    configuration_identity,
                )
            )
            or isinstance(revision, bool)
            or not isinstance(revision, int)
        ):
            msg = "Invalid analysis values snapshot identity or revision"
            raise ValueError(msg)
        return cls(
            occurrence_identity=occurrence_identity,
            model_identity=model_identity,
            definitions_identity=definitions_identity,
            configuration_identity=configuration_identity,
            revision=revision,
            _items=items,
        )


class AnalysisValues:
    """Session-owned mutable central values with one analysis-wide revision."""

    def __init__(self) -> None:
        self._occurrence_identity = ""
        self._model_identity = ""
        self._definitions_identity = ""
        self._configuration_identity = ""
        self._items: tuple[tuple[str, float], ...] = ()
        self._bounds: Mapping[str, tuple[float, float]] = MappingProxyType({})
        self._revision: int | None = None
        self._lock = RLock()

    @property
    def is_initialized(self) -> bool:
        """Whether sealed configuration has initialized central values."""
        with self._lock:
            return self._revision is not None

    def initialize(
        self,
        model_identity: str,
        configuration: SealedConfiguration,
    ) -> None:
        """Initialize revision zero from immutable configured effective values."""
        with self._lock:
            if self._revision is not None:
                msg = "Analysis values are already initialized"
                raise RuntimeError(msg)

            items: list[tuple[str, float]] = []
            for config in configuration:
                value = config.effective_value
                if value is None or not math.isfinite(value):
                    msg = (
                        f"Invalid initial analysis value for {config.param_id!r}: "
                        f"{value!r}"
                    )
                    raise ValueError(msg)
                items.append((config.param_id, float(value)))

            self._occurrence_identity = uuid4().hex
            self._model_identity = str(model_identity)
            self._definitions_identity = configuration.definitions_identity
            self._configuration_identity = configuration.identity
            self._items = tuple(items)
            self._bounds = MappingProxyType(
                {
                    config.param_id: (config.lower_bound, config.upper_bound)
                    for config in configuration
                }
            )
            self._revision = 0

    def snapshot(self) -> AnalysisValuesSnapshot:
        """Return an immutable point-in-time copy of all configured values."""
        with self._lock:
            return self._snapshot_unlocked()

    def _snapshot_unlocked(self) -> AnalysisValuesSnapshot:
        if self._revision is None:
            msg = "Analysis values are not initialized"
            raise RuntimeError(msg)
        return AnalysisValuesSnapshot(
            occurrence_identity=self._occurrence_identity,
            model_identity=self._model_identity,
            definitions_identity=self._definitions_identity,
            configuration_identity=self._configuration_identity,
            revision=self._revision,
            _items=self._items,
        )

    def commit(
        self,
        candidate: Mapping[str, float],
        *,
        expected: AnalysisValuesSnapshot,
        scope: tuple[str, ...],
    ) -> AnalysisValuesSnapshot:
        """Atomically commit a complete finite candidate for an explicit scope."""
        with self._lock:
            return self._commit_unlocked(candidate, expected=expected, scope=scope)

    def _commit_unlocked(
        self,
        candidate: Mapping[str, float],
        *,
        expected: AnalysisValuesSnapshot,
        scope: tuple[str, ...],
    ) -> AnalysisValuesSnapshot:
        current = self._snapshot_unlocked()
        if expected.occurrence_identity != current.occurrence_identity:
            msg = "Analysis values snapshot belongs to another analysis occurrence"
            raise IncompatibleAnalysisValuesError(msg)
        current_identity = (
            current.model_identity,
            current.definitions_identity,
            current.configuration_identity,
        )
        expected_identity = (
            expected.model_identity,
            expected.definitions_identity,
            expected.configuration_identity,
        )
        if expected_identity != current_identity:
            msg = "Analysis values snapshot has incompatible model or configuration"
            raise IncompatibleAnalysisValuesError(msg)
        if expected.revision != current.revision:
            msg = (
                f"Analysis values revision is stale: expected {expected.revision}, "
                f"current {current.revision}"
            )
            raise StaleAnalysisValuesError(msg)
        if expected._items != current._items:
            msg = "Analysis values snapshot does not match the current revision"
            raise IncompatibleAnalysisValuesError(msg)

        scope_ids = tuple(scope)
        scope_set = set(scope_ids)
        candidate_ids = set(candidate)
        configured_ids = set(current)
        if (
            not scope_ids
            or len(scope_set) != len(scope_ids)
            or scope_set - configured_ids
            or candidate_ids != scope_set
        ):
            msg = (
                "Invalid analysis values commit scope: "
                f"scope={scope_ids!r}, candidate_ids={tuple(candidate)!r}"
            )
            raise InvalidAnalysisValuesCommitError(msg)

        replacements: dict[str, float] = {}
        for param_id in scope_ids:
            value = candidate[param_id]
            if isinstance(value, bool) or not isinstance(value, Real):
                msg = f"Invalid central value for {param_id!r}: {value!r}"
                raise InvalidAnalysisValuesCommitError(msg)
            numeric_value = float(value)
            lower_bound, upper_bound = self._bounds[param_id]
            if (
                not math.isfinite(numeric_value)
                or not lower_bound <= numeric_value <= upper_bound
            ):
                msg = f"Invalid central value for {param_id!r}: {value!r}"
                raise InvalidAnalysisValuesCommitError(msg)
            replacements[param_id] = numeric_value

        self._items = tuple(
            (param_id, replacements.get(param_id, value))
            for param_id, value in self._items
        )
        self._revision = current.revision + 1
        return self._snapshot_unlocked()

    def reset(self) -> None:
        """Return to the uninitialized state for a full analysis rebuild."""
        with self._lock:
            self._occurrence_identity = ""
            self._model_identity = ""
            self._definitions_identity = ""
            self._configuration_identity = ""
            self._items = ()
            self._bounds = MappingProxyType({})
            self._revision = None

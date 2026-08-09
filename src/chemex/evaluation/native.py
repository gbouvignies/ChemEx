"""Non-authoritative immutable native evaluation boundary (#585).

This module deliberately has no optimizer or workflow entry point.  It is a
qualification harness: callers build a frozen plan from already-selected
ChemEx experiments, bind trusted local profile implementations, then evaluate
one complete independent frame through a private single-owner workspace.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import OrderedDict
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass, field
from numbers import Real
from threading import get_ident
from typing import Any, Literal, cast

import numpy as np

from chemex.containers.experiments import Experiments
from chemex.containers.profile import Profile
from chemex.parameters.parameterization import (
    ActiveParameterization,
    IndependentValueFrame,
    ParameterizationError,
    ResolvedParameterValues,
)
from chemex.typing import Array

_SCHEMA_VERSION = 1
_NORMALIZATION_VERSION = "weighted-zero-intercept-v1"
_RESIDUAL_VERSION = "calculated-minus-experimental-weighted-masked-v1"
_COMPATIBILITY_VERSION = "chemex-profile-kernel-v1"


@dataclass(frozen=True, slots=True)
class EvaluationFrame:
    """The narrow #572 public independent-value input.

    Lifecycle occurrence/revision validation happens before this projection is
    made by ``from_lifecycle_frame``; evaluators never expose that metadata.
    """

    parameterization_identity: str
    _items: tuple[tuple[str, float], ...]

    def __post_init__(self) -> None:
        items: list[tuple[str, float]] = []
        for param_id, value in self._items:
            if not isinstance(param_id, str):
                raise TypeError("Evaluation-frame parameter IDs must be strings")
            items.append((param_id, _finite_evaluation_scalar(value)))
        if len({param_id for param_id, _value in items}) != len(items):
            raise ValueError("Evaluation-frame parameter IDs must be unique")
        object.__setattr__(self, "_items", tuple(items))

    @classmethod
    def from_lifecycle_frame(
        cls,
        parameterization: ActiveParameterization,
        frame: IndependentValueFrame,
    ) -> EvaluationFrame:
        """Project a lifecycle-checked #584 frame into the public #572 input."""
        if (
            frame.parameterization_identity != parameterization.identity
            or frame.program_fingerprint != parameterization.program.fingerprint
            or frame.occurrence_identity != parameterization.occurrence_identity
            or frame.revision != parameterization.source_revision
        ):
            raise ValueError("Lifecycle frame is incompatible with parameterization")
        if tuple(param_id for param_id, _value in frame._items) != (
            parameterization.independent_ids
        ):
            raise ValueError("Lifecycle frame has non-canonical independent ID order")
        return cls(
            parameterization.identity,
            frame._items,
        )


def _finite_evaluation_scalar(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError("Evaluation-frame values must be real binary64 scalars")
    try:
        scalar = float(value)
    except OverflowError as error:
        raise ValueError("Evaluation-frame values must be finite") from error
    if not math.isfinite(scalar):
        raise ValueError("Evaluation-frame values must be finite")
    return 0.0 if scalar == 0.0 else scalar


def _canonical_float(value: float) -> str:
    scalar = float(value)
    if not math.isfinite(scalar):
        raise ValueError("Native evaluation plans require finite scalar data")
    return (0.0 if scalar == 0.0 else scalar).hex()


def _identity(kind: str, record: object) -> str:
    encoded = json.dumps(
        {"schema": _SCHEMA_VERSION, "kind": kind, "record": record},
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _semantic_value(value: object) -> object:
    """Encode finite scientific values without repr or runtime identity."""
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, Real):
        if isinstance(value, bool):
            return value
        return {"binary64": _canonical_float(float(value))}
    if isinstance(value, np.ndarray):
        raw_values = cast("Any", value).tolist()
        return _semantic_value(raw_values)
    if isinstance(value, Mapping):
        return {
            str(key): _semantic_value(item)
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_semantic_value(item) for item in value]
    raise TypeError(
        f"Unsupported native scientific descriptor value: {type(value).__name__}"
    )


def _semantic_json(value: object) -> str:
    return json.dumps(
        _semantic_value(value), ensure_ascii=True, separators=(",", ":"), sort_keys=True
    )


def _kernel_configuration(profile: Profile) -> str:
    settings = getattr(profile.pulse_sequence, "settings", None)
    if not hasattr(settings, "model_dump"):
        raise ValueError("Native profile kernels must expose immutable settings")
    return _semantic_json(settings.model_dump(mode="json"))


def _observation_metadata(profile: Profile) -> str:
    return _semantic_json(profile.data.metadata.tolist())


def _spectrometer_configuration(profile: Profile) -> str:
    return _semantic_json(profile.spectrometer.native_kernel_descriptor())


def _readonly(values: Array) -> Array:
    result = np.array(values, dtype=np.float64, copy=True)
    result.flags.writeable = False
    return result


def _record_list(value: object, *, field_name: str) -> list[object]:
    if not isinstance(value, list):
        raise TypeError(f"Evaluation-plan {field_name} must be a list")
    return cast("list[object]", value)


def _record_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"Evaluation-plan {field_name} must be an integer")
    return value


def _record_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"Evaluation-plan {field_name} must be a scalar")
    return _finite_evaluation_scalar(value)


def _record_binary64(value: object, *, field_name: str) -> float:
    if not isinstance(value, str):
        raise TypeError(f"Evaluation-plan {field_name} must be a binary64 string")
    try:
        scalar = float.fromhex(value)
    except ValueError as error:
        raise ValueError(f"Evaluation-plan {field_name} is not binary64") from error
    if _canonical_float(scalar) != value:
        raise ValueError(f"Evaluation-plan {field_name} is not canonical binary64")
    return scalar


def _record_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"Evaluation-plan {field_name} must be a boolean")
    return value


def _record_local_inputs(value: object) -> tuple[tuple[str, str], ...]:
    items = _record_list(value, field_name="local_inputs")
    pairs: list[tuple[str, str]] = []
    for item in items:
        pair = _record_list(item, field_name="local_input")
        if len(pair) != 2:
            raise TypeError("Evaluation-plan local inputs must be string pairs")
        name, param_id = pair
        if not isinstance(name, str) or not isinstance(param_id, str):
            raise TypeError("Evaluation-plan local inputs must be string pairs")
        pairs.append((name, param_id))
    return tuple(pairs)


def _profile_source_identity(profile: Profile) -> str:
    """Fingerprint every scientific input that a trusted binding must retain."""
    data = profile.data
    return _identity(
        "profile-source",
        (
            tuple(
                (name, param_id) for name, param_id in sorted(profile.name_map.items())
            ),
            bool(profile.is_scaled),
            type(profile.pulse_sequence).__module__,
            type(profile.pulse_sequence).__qualname__,
            _kernel_configuration(profile),
            _spectrometer_configuration(profile),
            tuple(_canonical_float(value) for value in np.ravel(data.exp)),
            tuple(_canonical_float(value) for value in np.ravel(data.err)),
            tuple(bool(value) for value in np.ravel(data.mask)),
            _observation_metadata(profile),
        ),
    )


@dataclass(frozen=True, slots=True)
class ProfilePlan:
    """Frozen population and kernel contract for one ChemEx profile."""

    identity: str
    source_identity: str
    experiment_ordinal: int
    profile_ordinal: int
    observation_offset: int
    local_inputs: tuple[tuple[str, str], ...]
    is_scaled: bool
    experimental_values: tuple[float, ...]
    uncertainties: tuple[float, ...]
    mask: tuple[bool, ...]
    kernel_identity: str
    kernel_configuration: str
    spectrometer_configuration: str
    observation_metadata: str
    output_shape: tuple[int, ...]

    @property
    def param_ids(self) -> tuple[str, ...]:
        return tuple(sorted({param_id for _name, param_id in self.local_inputs}))

    @property
    def observation_count(self) -> int:
        return len(self.experimental_values)

    @property
    def retained_observation_indices(self) -> tuple[int, ...]:
        return tuple(index for index, included in enumerate(self.mask) if included)


@dataclass(frozen=True, slots=True)
class EvaluationPlan:
    """Serializable immutable specification of one exact flat population."""

    parameterization_identity: str
    constraint_program_identity: str
    profiles: tuple[ProfilePlan, ...]
    normalization_version: str = _NORMALIZATION_VERSION
    residual_version: str = _RESIDUAL_VERSION
    identity: str = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "identity",
            _identity(
                "evaluation-plan",
                (
                    self.parameterization_identity,
                    self.constraint_program_identity,
                    tuple(
                        (
                            item.identity,
                            item.source_identity,
                            item.experiment_ordinal,
                            item.profile_ordinal,
                            item.observation_offset,
                            item.local_inputs,
                            item.is_scaled,
                            tuple(
                                _canonical_float(value)
                                for value in item.experimental_values
                            ),
                            tuple(
                                _canonical_float(value) for value in item.uncertainties
                            ),
                            item.mask,
                            item.kernel_identity,
                            item.kernel_configuration,
                            item.spectrometer_configuration,
                            item.observation_metadata,
                            item.output_shape,
                        )
                        for item in self.profiles
                    ),
                    self.normalization_version,
                    self.residual_version,
                ),
            ),
        )

    @property
    def observation_count(self) -> int:
        return sum(item.observation_count for item in self.profiles)

    @property
    def retained_observation_count(self) -> int:
        return sum(len(item.retained_observation_indices) for item in self.profiles)

    def to_record(self) -> dict[str, object]:
        """Return an executable-object-free serialization payload."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "parameterization_identity": self.parameterization_identity,
            "constraint_program_identity": self.constraint_program_identity,
            "normalization_version": self.normalization_version,
            "residual_version": self.residual_version,
            "profiles": [
                {
                    "identity": item.identity,
                    "source_identity": item.source_identity,
                    "experiment_ordinal": item.experiment_ordinal,
                    "profile_ordinal": item.profile_ordinal,
                    "observation_offset": item.observation_offset,
                    "local_inputs": [list(value) for value in item.local_inputs],
                    "is_scaled": item.is_scaled,
                    "experimental_values": [
                        _canonical_float(value) for value in item.experimental_values
                    ],
                    "uncertainties": [
                        _canonical_float(value) for value in item.uncertainties
                    ],
                    "mask": list(item.mask),
                    "kernel_identity": item.kernel_identity,
                    "kernel_configuration": item.kernel_configuration,
                    "spectrometer_configuration": item.spectrometer_configuration,
                    "observation_metadata": item.observation_metadata,
                    "output_shape": list(item.output_shape),
                }
                for item in self.profiles
            ],
            "identity": self.identity,
        }

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> EvaluationPlan:
        """Restore and verify a plan before it is rebound to local kernels."""
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise ValueError("Unsupported native evaluation-plan schema")
        raw_profiles = record.get("profiles")
        if not isinstance(raw_profiles, list):
            raise TypeError("Evaluation-plan profiles must be a list")
        try:
            profiles = tuple(
                ProfilePlan(
                    identity=str(cast("Mapping[str, object]", item)["identity"]),
                    source_identity=str(
                        cast("Mapping[str, object]", item)["source_identity"]
                    ),
                    experiment_ordinal=_record_int(
                        cast("Mapping[str, object]", item)["experiment_ordinal"],
                        field_name="experiment_ordinal",
                    ),
                    profile_ordinal=_record_int(
                        cast("Mapping[str, object]", item)["profile_ordinal"],
                        field_name="profile_ordinal",
                    ),
                    observation_offset=_record_int(
                        cast("Mapping[str, object]", item)["observation_offset"],
                        field_name="observation_offset",
                    ),
                    local_inputs=_record_local_inputs(
                        cast("Mapping[str, object]", item)["local_inputs"]
                    ),
                    is_scaled=_record_bool(
                        cast("Mapping[str, object]", item)["is_scaled"],
                        field_name="is_scaled",
                    ),
                    experimental_values=tuple(
                        _record_binary64(value, field_name="experimental_values")
                        for value in _record_list(
                            cast("Mapping[str, object]", item)["experimental_values"],
                            field_name="experimental_values",
                        )
                    ),
                    uncertainties=tuple(
                        _record_binary64(value, field_name="uncertainties")
                        for value in _record_list(
                            cast("Mapping[str, object]", item)["uncertainties"],
                            field_name="uncertainties",
                        )
                    ),
                    mask=tuple(
                        _record_bool(value, field_name="mask")
                        for value in _record_list(
                            cast("Mapping[str, object]", item)["mask"],
                            field_name="mask",
                        )
                    ),
                    kernel_identity=str(
                        cast("Mapping[str, object]", item)["kernel_identity"]
                    ),
                    kernel_configuration=str(
                        cast("Mapping[str, object]", item)["kernel_configuration"]
                    ),
                    spectrometer_configuration=str(
                        cast("Mapping[str, object]", item)["spectrometer_configuration"]
                    ),
                    observation_metadata=str(
                        cast("Mapping[str, object]", item)["observation_metadata"]
                    ),
                    output_shape=tuple(
                        _record_int(value, field_name="output_shape")
                        for value in _record_list(
                            cast("Mapping[str, object]", item)["output_shape"],
                            field_name="output_shape",
                        )
                    ),
                )
                for item in raw_profiles
                if isinstance(item, Mapping)
            )
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Malformed native evaluation plan") from error
        if len(profiles) != len(raw_profiles):
            raise ValueError("Malformed native evaluation-profile record")
        plan = cls(
            parameterization_identity=str(record["parameterization_identity"]),
            constraint_program_identity=str(record["constraint_program_identity"]),
            profiles=profiles,
            normalization_version=str(record["normalization_version"]),
            residual_version=str(record["residual_version"]),
        )
        if record.get("identity") != plan.identity:
            raise ValueError("Evaluation-plan fingerprint does not match its payload")
        return plan


@dataclass(frozen=True, slots=True)
class ProfileEvaluation:
    """One complete immutable profile segment in canonical plan order."""

    profile_identity: str
    observation_offset: int
    observation_count: int
    residual_offset: int
    residual_count: int
    retained_observation_indices: tuple[int, ...]
    normalization_factor: float
    kernel_identity: str


def _profile_evaluation_from_record(
    record: object,
    descriptor: ProfilePlan,
) -> ProfileEvaluation:
    if not isinstance(record, Mapping):
        raise TypeError("Evaluation-result profile must be a mapping")
    values = cast("Mapping[str, object]", record)
    retained = tuple(
        _record_int(value, field_name="retained observation index")
        for value in _record_list(
            values.get("retained_observation_indices"),
            field_name="retained_observation_indices",
        )
    )
    if (
        values.get("profile_identity") != descriptor.identity
        or values.get("kernel_identity") != descriptor.kernel_identity
        or retained != descriptor.retained_observation_indices
    ):
        raise ValueError("Evaluation-result profile does not match frozen plan")
    return ProfileEvaluation(
        descriptor.identity,
        _record_int(values.get("observation_offset"), field_name="observation_offset"),
        _record_int(values.get("observation_count"), field_name="observation_count"),
        _record_int(values.get("residual_offset"), field_name="residual_offset"),
        _record_int(values.get("residual_count"), field_name="residual_count"),
        retained,
        _record_float(values.get("normalization_factor"), field_name="normalization"),
        descriptor.kernel_identity,
    )


@dataclass(frozen=True, slots=True)
class EvaluationResult:
    """The only successful native-evaluation outcome."""

    plan_identity: str
    parameterization_identity: str
    evaluator_compatibility_identity: str
    resolved_values: ResolvedParameterValues
    unscaled_calculations: Array
    normalized_calculations: Array
    residuals: Array
    profiles: tuple[ProfileEvaluation, ...]

    @property
    def identity(self) -> str:
        """Deterministic identity of this complete immutable scientific outcome."""
        return _identity(
            "evaluation-result",
            (
                self.plan_identity,
                self.parameterization_identity,
                self.evaluator_compatibility_identity,
                tuple(self.resolved_values.items()),
                tuple(_canonical_float(value) for value in self.unscaled_calculations),
                tuple(
                    _canonical_float(value) for value in self.normalized_calculations
                ),
                tuple(_canonical_float(value) for value in self.residuals),
                tuple(
                    (
                        item.profile_identity,
                        item.observation_offset,
                        item.observation_count,
                        item.residual_offset,
                        item.residual_count,
                        item.retained_observation_indices,
                        _canonical_float(item.normalization_factor),
                        item.kernel_identity,
                    )
                    for item in self.profiles
                ),
            ),
        )

    def to_record(self) -> dict[str, object]:
        """Serialize immutable scientific evidence without runtime machinery."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "plan_identity": self.plan_identity,
            "parameterization_identity": self.parameterization_identity,
            "evaluator_compatibility_identity": self.evaluator_compatibility_identity,
            "resolved": {
                "program_fingerprint": self.resolved_values.program_fingerprint,
                "occurrence_identity": self.resolved_values.occurrence_identity,
                "revision": self.resolved_values.revision,
                "items": list(self.resolved_values.items()),
            },
            "unscaled_calculations": self.unscaled_calculations.tolist(),
            "normalized_calculations": self.normalized_calculations.tolist(),
            "residuals": self.residuals.tolist(),
            "profiles": [
                {
                    "profile_identity": item.profile_identity,
                    "observation_offset": item.observation_offset,
                    "observation_count": item.observation_count,
                    "residual_offset": item.residual_offset,
                    "residual_count": item.residual_count,
                    "retained_observation_indices": list(
                        item.retained_observation_indices
                    ),
                    "normalization_factor": item.normalization_factor,
                    "kernel_identity": item.kernel_identity,
                }
                for item in self.profiles
            ],
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        plan: EvaluationPlan,
    ) -> EvaluationResult:
        """Restore a result only when its population matches a frozen plan."""
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise ValueError("Unsupported native evaluation-result schema")
        if record.get("plan_identity") != plan.identity:
            raise ValueError("Evaluation result belongs to another plan")
        resolved_record = record.get("resolved")
        if not isinstance(resolved_record, Mapping):
            raise TypeError("Evaluation-result resolved values must be a mapping")
        resolved_values = cast("Mapping[str, object]", resolved_record)
        items = _record_list(resolved_values.get("items"), field_name="resolved.items")
        resolved = ResolvedParameterValues(
            parameterization_identity=str(record["parameterization_identity"]),
            program_fingerprint=str(resolved_values["program_fingerprint"]),
            occurrence_identity=str(resolved_values["occurrence_identity"]),
            revision=_record_int(resolved_values["revision"], field_name="revision"),
            _items=tuple(
                (
                    str(_record_list(item, field_name="resolved item")[0]),
                    _record_float(
                        _record_list(item, field_name="resolved item")[1],
                        field_name="resolved value",
                    ),
                )
                for item in items
            ),
        )
        profiles_record = _record_list(record.get("profiles"), field_name="profiles")
        if len(profiles_record) != len(plan.profiles):
            raise ValueError("Evaluation result has the wrong profile population")
        profiles = tuple(
            _profile_evaluation_from_record(item, descriptor)
            for item, descriptor in zip(profiles_record, plan.profiles, strict=True)
        )
        unscaled = _readonly(
            np.asarray(
                _record_list(
                    record.get("unscaled_calculations"),
                    field_name="unscaled calculations",
                ),
                dtype=np.float64,
            )
        )
        normalized = _readonly(
            np.asarray(
                _record_list(
                    record.get("normalized_calculations"),
                    field_name="normalized calculations",
                ),
                dtype=np.float64,
            )
        )
        residuals = _readonly(
            np.asarray(
                _record_list(record.get("residuals"), field_name="residuals"),
                dtype=np.float64,
            )
        )
        if (
            unscaled.shape != (plan.observation_count,)
            or normalized.shape != (plan.observation_count,)
            or residuals.shape != (plan.retained_observation_count,)
            or not np.all(np.isfinite(unscaled))
            or not np.all(np.isfinite(normalized))
            or not np.all(np.isfinite(residuals))
        ):
            raise ValueError("Malformed native evaluation-result arrays")
        return cls(
            plan.identity,
            str(record["parameterization_identity"]),
            str(record["evaluator_compatibility_identity"]),
            resolved,
            unscaled,
            normalized,
            residuals,
            profiles,
        )


@dataclass(frozen=True, slots=True)
class EvaluationFailure:
    """Atomic typed non-success outcome; no usable vectors accompany it."""

    plan_identity: str
    parameterization_identity: str
    stage: Literal[
        "frame", "resolution", "kernel", "normalization", "residual", "binding"
    ]
    category: str
    validity: Literal[
        "INVALID_REQUEST",
        "INVALID_TRIAL",
        "INVALID_PLAN_OR_BINDING",
        "IMPLEMENTATION_FAILURE",
    ]
    profile_identity: str | None = None
    message: str = ""

    def to_record(self) -> dict[str, object]:
        """Serialize stable failure information without exception objects."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "plan_identity": self.plan_identity,
            "parameterization_identity": self.parameterization_identity,
            "stage": self.stage,
            "category": self.category,
            "validity": self.validity,
            "profile_identity": self.profile_identity,
            "message": self.message,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        plan: EvaluationPlan,
    ) -> EvaluationFailure:
        """Restore a sanitized failure only for the matching frozen plan."""
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise ValueError("Unsupported native evaluation-failure schema")
        if record.get("plan_identity") != plan.identity:
            raise ValueError("Evaluation failure belongs to another plan")
        stage = str(record.get("stage"))
        validity = str(record.get("validity"))
        if stage not in {
            "frame",
            "resolution",
            "kernel",
            "normalization",
            "residual",
            "binding",
        }:
            raise ValueError("Unknown native evaluation-failure stage")
        if validity not in {
            "INVALID_REQUEST",
            "INVALID_TRIAL",
            "INVALID_PLAN_OR_BINDING",
            "IMPLEMENTATION_FAILURE",
        }:
            raise ValueError("Unknown native evaluation-failure validity")
        profile_identity = record.get("profile_identity")
        if profile_identity is not None and profile_identity not in {
            item.identity for item in plan.profiles
        }:
            raise ValueError("Evaluation failure names a profile outside its plan")
        return cls(
            plan.identity,
            str(record["parameterization_identity"]),
            stage,
            str(record["category"]),
            validity,
            None if profile_identity is None else str(profile_identity),
            str(record.get("message", "")),
        )


@dataclass(frozen=True, slots=True)
class CacheStatistics:
    hits: int
    misses: int
    entries: int


@dataclass(frozen=True, slots=True)
class _CachedProfile:
    unscaled: Array
    normalized: Array
    residuals: Array
    normalization_factor: float


@dataclass(slots=True)
class _Workspace:
    templates: tuple[Profile, ...]
    profiles: tuple[Profile, ...]
    cache: OrderedDict[tuple[object, ...], _CachedProfile] = field(
        default_factory=OrderedDict
    )
    hits: int = 0
    misses: int = 0
    poisoned: bool = False

    def reset(self) -> None:
        self.profiles = tuple(deepcopy(profile) for profile in self.templates)
        self.cache.clear()


def _build_plan(
    experiments: Experiments,
    parameterization: ActiveParameterization,
) -> tuple[EvaluationPlan, tuple[tuple[int, int, Profile], ...]]:
    profile_plans: list[ProfilePlan] = []
    sources: list[tuple[int, int, Profile]] = []
    offset = 0
    for experiment_ordinal, experiment in enumerate(experiments):
        for profile_ordinal, profile in enumerate(experiment.profiles):
            exp = np.asarray(profile.data.exp, dtype=np.float64)
            err = np.asarray(profile.data.err, dtype=np.float64)
            mask = np.asarray(profile.data.mask, dtype=np.bool_)
            if exp.ndim != 1 or err.shape != exp.shape or mask.shape != exp.shape:
                raise ValueError(
                    "Native evaluation requires one-dimensional symmetric data"
                )
            retained = mask
            if not np.all(np.isfinite(exp[retained])) or not np.all(
                np.isfinite(err[retained]) & (err[retained] > 0.0)
            ):
                raise ValueError(
                    "Retained native observations require finite positive errors"
                )
            if profile.is_scaled and not np.any(retained):
                raise ValueError(
                    "A scaled native profile cannot retain zero observations"
                )
            source_identity = _profile_source_identity(profile)
            kernel_identity = (
                f"{type(profile.pulse_sequence).__module__}:"
                f"{type(profile.pulse_sequence).__qualname__}:{_COMPATIBILITY_VERSION}"
            )
            identity = _identity(
                "profile-plan",
                (source_identity, experiment_ordinal, profile_ordinal, offset),
            )
            profile_plans.append(
                ProfilePlan(
                    identity=identity,
                    source_identity=source_identity,
                    experiment_ordinal=experiment_ordinal,
                    profile_ordinal=profile_ordinal,
                    observation_offset=offset,
                    local_inputs=tuple(sorted(profile.name_map.items())),
                    is_scaled=profile.is_scaled,
                    experimental_values=tuple(float(value) for value in exp),
                    uncertainties=tuple(float(value) for value in err),
                    mask=tuple(bool(value) for value in mask),
                    kernel_identity=kernel_identity,
                    kernel_configuration=_kernel_configuration(profile),
                    spectrometer_configuration=_spectrometer_configuration(profile),
                    observation_metadata=_observation_metadata(profile),
                    output_shape=tuple(exp.shape),
                )
            )
            sources.append((experiment_ordinal, profile_ordinal, profile))
            offset += exp.size
    return (
        EvaluationPlan(
            parameterization_identity=parameterization.identity,
            constraint_program_identity=parameterization.program.fingerprint,
            profiles=tuple(profile_plans),
        ),
        tuple(sources),
    )


class EvaluationEngine:
    """Binds one immutable plan to trusted local ChemEx profile kernels."""

    def __init__(
        self,
        plan: EvaluationPlan,
        parameterization: ActiveParameterization,
        sources: Sequence[tuple[int, int, Profile]],
    ) -> None:
        if plan.parameterization_identity != parameterization.identity:
            raise ValueError("Evaluation plan belongs to another parameterization")
        if plan.constraint_program_identity != parameterization.program.fingerprint:
            raise ValueError("Evaluation plan belongs to another constraint program")
        if (
            plan.normalization_version != _NORMALIZATION_VERSION
            or plan.residual_version != _RESIDUAL_VERSION
        ):
            raise ValueError("Evaluation plan has incompatible execution semantics")
        if len(plan.profiles) != len(sources):
            raise ValueError("Trusted profile bindings do not match evaluation plan")
        expected_offset = 0
        for descriptor, (experiment_ordinal, profile_ordinal, source) in zip(
            plan.profiles, sources, strict=True
        ):
            if (
                descriptor.experiment_ordinal != experiment_ordinal
                or descriptor.profile_ordinal != profile_ordinal
                or descriptor.observation_offset != expected_offset
                or descriptor.observation_count != source.data.exp.size
                or descriptor.identity
                != _identity(
                    "profile-plan",
                    (
                        _profile_source_identity(source),
                        experiment_ordinal,
                        profile_ordinal,
                        expected_offset,
                    ),
                )
            ):
                raise ValueError(
                    "Trusted population structure does not match frozen plan"
                )
            if descriptor.source_identity != _profile_source_identity(source):
                raise ValueError("Trusted profile binding does not match frozen plan")
            kernel_identity = (
                f"{type(source.pulse_sequence).__module__}:"
                f"{type(source.pulse_sequence).__qualname__}:{_COMPATIBILITY_VERSION}"
            )
            if (
                descriptor.kernel_identity != kernel_identity
                or descriptor.kernel_configuration != _kernel_configuration(source)
                or descriptor.spectrometer_configuration
                != _spectrometer_configuration(source)
                or descriptor.local_inputs != tuple(sorted(source.name_map.items()))
                or descriptor.observation_metadata != _observation_metadata(source)
                or descriptor.output_shape != tuple(source.data.exp.shape)
            ):
                raise ValueError("Trusted kernel descriptor does not match frozen plan")
            expected_offset += descriptor.observation_count
        self.plan = plan
        self._parameterization = parameterization
        # Copy construction-owned mutable machinery once.  Plans retain no
        # source objects, and callers cannot alter a future workspace by
        # mutating the legacy occurrence after it was bound.
        self._sources = tuple(deepcopy(source) for _e, _p, source in sources)
        self.compatibility_identity = _identity(
            "evaluation-compatibility",
            tuple(item.kernel_identity for item in plan.profiles),
        )

    @classmethod
    def from_experiments(
        cls,
        experiments: Experiments,
        parameterization: ActiveParameterization,
    ) -> EvaluationEngine:
        plan, sources = _build_plan(experiments, parameterization)
        return cls(plan, parameterization, sources)

    @classmethod
    def bind(
        cls,
        plan: EvaluationPlan,
        parameterization: ActiveParameterization,
        experiments: Experiments,
    ) -> EvaluationEngine:
        """Trusted local rebinding path for a deserialized plan."""
        sources = tuple(
            (experiment_ordinal, profile_ordinal, profile)
            for experiment_ordinal, experiment in enumerate(experiments)
            for profile_ordinal, profile in enumerate(experiment.profiles)
        )
        return cls(plan, parameterization, sources)

    def new_evaluator(self) -> BoundEvaluator:
        """Create an empty single-owner workspace from trusted implementations."""
        return BoundEvaluator(
            self.plan,
            self._parameterization,
            self.compatibility_identity,
            _Workspace(
                self._sources,
                tuple(deepcopy(profile) for profile in self._sources),
            ),
        )


class BoundEvaluator:
    """One non-reentrant evaluator with private reusable ChemEx workspaces."""

    def __init__(
        self,
        plan: EvaluationPlan,
        parameterization: ActiveParameterization,
        compatibility_identity: str,
        workspace: _Workspace,
    ) -> None:
        self.plan = plan
        self._parameterization = parameterization
        self.compatibility_identity = compatibility_identity
        self._workspace = workspace
        self._owner = get_ident()

    @property
    def cache_statistics(self) -> CacheStatistics:
        return CacheStatistics(
            self._workspace.hits,
            self._workspace.misses,
            len(self._workspace.cache),
        )

    def _failure(
        self,
        stage: Literal[
            "frame", "resolution", "kernel", "normalization", "residual", "binding"
        ],
        category: str,
        validity: Literal[
            "INVALID_REQUEST",
            "INVALID_TRIAL",
            "INVALID_PLAN_OR_BINDING",
            "IMPLEMENTATION_FAILURE",
        ],
        *,
        profile_identity: str | None = None,
        message: str = "",
    ) -> EvaluationFailure:
        if validity in {"INVALID_PLAN_OR_BINDING", "IMPLEMENTATION_FAILURE"}:
            self._workspace.poisoned = True
            self._workspace.cache.clear()
        elif validity == "INVALID_TRIAL":
            self._workspace.reset()
        return EvaluationFailure(
            self.plan.identity,
            self._parameterization.identity,
            stage,
            category,
            validity,
            profile_identity,
            message,
        )

    def _calculate_unscaled(
        self,
        descriptor: ProfilePlan,
        profile: Profile,
        resolved: ResolvedParameterValues,
    ) -> Array | EvaluationFailure:
        """Run and validate the narrow unscaled profile kernel."""
        try:
            unscaled = np.asarray(profile.calculate_unscaled(resolved))
        except Exception as error:  # noqa: BLE001 - native kernel fault boundary
            return self._failure(
                "kernel",
                "kernel_exception",
                "IMPLEMENTATION_FAILURE",
                profile_identity=descriptor.identity,
                message=str(error),
            )
        if unscaled.shape != (descriptor.observation_count,):
            return self._failure(
                "kernel",
                "unexpected_shape",
                "IMPLEMENTATION_FAILURE",
                profile_identity=descriptor.identity,
            )
        if not np.issubdtype(unscaled.dtype, np.floating):
            return self._failure(
                "kernel",
                "unexpected_dtype",
                "IMPLEMENTATION_FAILURE",
                profile_identity=descriptor.identity,
            )
        unscaled_values = _readonly(unscaled)
        if not np.all(np.isfinite(unscaled_values)):
            return self._failure(
                "kernel",
                "non_finite_calculation",
                "INVALID_TRIAL",
                profile_identity=descriptor.identity,
            )
        return unscaled_values

    def _normalize_profile(
        self,
        descriptor: ProfilePlan,
        unscaled: Array,
    ) -> _CachedProfile | EvaluationFailure:
        """Apply the frozen #572 normalization and retained residual contract."""
        exp = np.asarray(descriptor.experimental_values, dtype=np.float64)
        err = np.asarray(descriptor.uncertainties, dtype=np.float64)
        retained = np.asarray(descriptor.mask, dtype=np.bool_)
        if descriptor.is_scaled:
            selected_exp = exp[retained]
            selected_calc = unscaled[retained]
            selected_err = err[retained]
            numerator = float(
                np.dot(selected_exp / selected_err, selected_calc / selected_err)
            )
            denominator = float(
                np.dot(selected_calc / selected_err, selected_calc / selected_err)
            )
            if (
                not math.isfinite(numerator)
                or not math.isfinite(denominator)
                or denominator == 0.0
            ):
                return self._failure(
                    "normalization",
                    "invalid_normalization",
                    "INVALID_TRIAL",
                    profile_identity=descriptor.identity,
                )
            scale = numerator / denominator
            if not math.isfinite(scale):
                return self._failure(
                    "normalization",
                    "non_finite_normalization",
                    "INVALID_TRIAL",
                    profile_identity=descriptor.identity,
                )
        else:
            scale = 1.0
        normalized = _readonly(scale * unscaled)
        residuals = _readonly((normalized[retained] - exp[retained]) / err[retained])
        if not np.all(np.isfinite(normalized)) or not np.all(np.isfinite(residuals)):
            return self._failure(
                "residual",
                "non_finite_residual",
                "INVALID_TRIAL",
                profile_identity=descriptor.identity,
            )
        return _CachedProfile(unscaled, normalized, residuals, scale)

    def _cached_profile(
        self,
        descriptor: ProfilePlan,
        profile: Profile,
        resolved: ResolvedParameterValues,
    ) -> _CachedProfile | EvaluationFailure:
        key = (
            descriptor.identity,
            *(resolved[param_id] for param_id in descriptor.param_ids),
        )
        cached = self._workspace.cache.get(key)
        if cached is not None:
            self._workspace.cache.move_to_end(key)
            self._workspace.hits += 1
            return cached
        self._workspace.misses += 1
        unscaled = self._calculate_unscaled(descriptor, profile, resolved)
        if isinstance(unscaled, EvaluationFailure):
            return unscaled
        cached = self._normalize_profile(descriptor, unscaled)
        if isinstance(cached, EvaluationFailure):
            return cached
        self._workspace.cache[key] = cached
        self._workspace.cache.move_to_end(key)
        while len(self._workspace.cache) > 5:
            self._workspace.cache.popitem(last=False)
        return cached

    def evaluate(
        self,
        frame: EvaluationFrame,
    ) -> EvaluationResult | EvaluationFailure:
        """Resolve once and atomically return a complete result or typed failure."""
        if get_ident() != self._owner:
            return self._failure(
                "binding", "workspace_owner_violation", "INVALID_REQUEST"
            )
        if self._workspace.poisoned:
            return self._failure(
                "binding", "workspace_poisoned", "INVALID_PLAN_OR_BINDING"
            )
        try:
            if frame.parameterization_identity != self._parameterization.identity:
                return self._failure(
                    "frame", "parameterization_mismatch", "INVALID_REQUEST"
                )
            if (
                tuple(param_id for param_id, _value in frame._items)
                != self._parameterization.independent_ids
            ):
                return self._failure(
                    "frame", "independent_value_order", "INVALID_REQUEST"
                )
            lifecycle_frame = IndependentValueFrame(
                self._parameterization.identity,
                self._parameterization.program.fingerprint,
                self._parameterization.occurrence_identity,
                self._parameterization.source_revision,
                frame._items,
            )
            resolved = self._parameterization.resolve(lifecycle_frame)
        except ParameterizationError as error:
            return self._failure(
                "resolution", error.code, "INVALID_REQUEST", message=str(error)
            )
        except Exception as error:  # noqa: BLE001 - trusted-program boundary
            return self._failure(
                "resolution",
                "unexpected_resolution_error",
                "IMPLEMENTATION_FAILURE",
                message=str(error),
            )
        unscaled_parts: list[Array] = []
        normalized_parts: list[Array] = []
        residual_parts: list[Array] = []
        profiles: list[ProfileEvaluation] = []
        residual_offset = 0
        for descriptor, profile in zip(
            self.plan.profiles, self._workspace.profiles, strict=True
        ):
            cached = self._cached_profile(descriptor, profile, resolved)
            if isinstance(cached, EvaluationFailure):
                return cached
            retained = descriptor.retained_observation_indices
            profiles.append(
                ProfileEvaluation(
                    descriptor.identity,
                    descriptor.observation_offset,
                    descriptor.observation_count,
                    residual_offset,
                    len(retained),
                    retained,
                    cached.normalization_factor,
                    descriptor.kernel_identity,
                )
            )
            residual_offset += len(retained)
            unscaled_parts.append(cached.unscaled)
            normalized_parts.append(cached.normalized)
            residual_parts.append(cached.residuals)
        return EvaluationResult(
            self.plan.identity,
            self._parameterization.identity,
            self.compatibility_identity,
            resolved,
            _readonly(
                np.concatenate(unscaled_parts) if unscaled_parts else np.empty(0)
            ),
            _readonly(
                np.concatenate(normalized_parts) if normalized_parts else np.empty(0)
            ),
            _readonly(
                np.concatenate(residual_parts) if residual_parts else np.empty(0)
            ),
            tuple(profiles),
        )

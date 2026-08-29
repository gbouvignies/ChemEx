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
import os
from collections import OrderedDict
from collections.abc import Iterator, Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass, field
from numbers import Real
from threading import get_ident
from typing import Any, Literal, cast

import numpy as np

from chemex.containers.data import Data
from chemex.containers.experiments import Experiments
from chemex.containers.profile import Profile, PulseSequence
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.parameterization import (
    ActiveParameterization,
    IndependentValueFrame,
    ParameterizationError,
)
from chemex.typing import Array

_SCHEMA_VERSION = 1
_SCALAR_VERSION = "canonical-binary64-v1"
_NORMALIZATION_VERSION = "weighted-zero-intercept-ordered-binary64-v2"
_RESIDUAL_VERSION = "calculated-minus-experimental-weighted-masked-v1"
_MASKING_VERSION = "retained-true-observations-v1"
_ORDERING_VERSION = "experiment-profile-row-v1"
_FAILURE_VERSION = "typed-failure-v1"
_DIAGNOSTICS_VERSION = "no-contractual-diagnostics-v1"
_COMPATIBILITY_VERSION = "chemex-profile-kernel-v1"
_PROFILE_CACHE_CAPACITY = 5
type FailureStage = Literal[
    "frame",
    "resolution",
    "projection",
    "cache",
    "kernel",
    "normalization",
    "residual",
    "result",
    "binding",
]
type FailureValidity = Literal[
    "INVALID_REQUEST",
    "INVALID_TRIAL",
    "INVALID_PLAN_OR_BINDING",
    "IMPLEMENTATION_FAILURE",
]
_FAILURE_TAXONOMY = {
    "frame": {
        "parameterization_mismatch": "INVALID_REQUEST",
        "independent_value_order": "INVALID_REQUEST",
    },
    "resolution": {
        "no_match": "INVALID_PLAN_OR_BINDING",
        "ambiguity": "INVALID_PLAN_OR_BINDING",
        "self_reference": "INVALID_PLAN_OR_BINDING",
        "cycle": "INVALID_PLAN_OR_BINDING",
        "domain_error": "INVALID_TRIAL",
        "evaluation_error": "INVALID_TRIAL",
        "non_finite": "INVALID_TRIAL",
        "incomplete_dependencies": "INVALID_PLAN_OR_BINDING",
        "incompatible_input": "INVALID_PLAN_OR_BINDING",
        "program_mismatch": "INVALID_PLAN_OR_BINDING",
        "unsupported_expression": "INVALID_PLAN_OR_BINDING",
        "model_derivation_override": "INVALID_PLAN_OR_BINDING",
        "unexpected_resolution_error": "IMPLEMENTATION_FAILURE",
    },
    "projection": {"missing_local_parameter": "INVALID_PLAN_OR_BINDING"},
    "cache": {"cache_key_exception": "IMPLEMENTATION_FAILURE"},
    "kernel": {
        "kernel_exception": "IMPLEMENTATION_FAILURE",
        "unexpected_container": "IMPLEMENTATION_FAILURE",
        "unexpected_shape": "IMPLEMENTATION_FAILURE",
        "unexpected_dtype": "IMPLEMENTATION_FAILURE",
        "non_finite_calculation": "INVALID_TRIAL",
    },
    "normalization": {"invalid_normalization": "INVALID_TRIAL"},
    "residual": {"non_finite_residual": "INVALID_TRIAL"},
    "result": {"result_assembly_exception": "IMPLEMENTATION_FAILURE"},
    "binding": {
        "workspace_owner_violation": "INVALID_REQUEST",
        "workspace_reentrant": "INVALID_REQUEST",
        "workspace_process_violation": "INVALID_PLAN_OR_BINDING",
        "workspace_poisoned": "INVALID_PLAN_OR_BINDING",
        "unexpected_evaluation_error": "IMPLEMENTATION_FAILURE",
    },
}


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
        return cls(parameterization.evaluator_identity, frame._items)

    def to_record(self) -> dict[str, object]:
        """Serialize only canonical evaluator-facing frame semantics."""
        return {
            "schema_version": _SCHEMA_VERSION,
            "parameterization_identity": self.parameterization_identity,
            "items": [
                [param_id, _canonical_float(value)] for param_id, value in self._items
            ],
            "identity": self.identity,
        }

    @property
    def identity(self) -> str:
        return _identity(
            "evaluation-frame",
            (
                self.parameterization_identity,
                tuple((key, _canonical_float(value)) for key, value in self._items),
            ),
        )

    @classmethod
    def from_record(cls, record: Mapping[str, object]) -> EvaluationFrame:
        _record_exact_keys(
            record,
            {"schema_version", "parameterization_identity", "items", "identity"},
            "Evaluation-frame",
        )
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise ValueError("Unsupported native evaluation-frame schema")
        identity = record.get("parameterization_identity")
        if not isinstance(identity, str):
            raise TypeError(
                "Evaluation-frame parameterization identity must be a string"
            )
        items = tuple(
            _record_frame_item(item)
            for item in _record_list(record.get("items"), field_name="frame items")
        )
        frame = cls(identity, items)
        if record.get("identity") != frame.identity:
            raise ValueError("Evaluation-frame fingerprint does not match its payload")
        return frame


def _finite_evaluation_scalar(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, np.float64)):
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


def _semantic_value(value: object, *, allow_infinity: bool = False) -> object:
    """Encode finite scientific values without repr or runtime identity."""
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, Real):
        if isinstance(value, bool):
            return value
        scalar = float(value)
        if allow_infinity and math.isinf(scalar):
            return {"infinity": "positive" if scalar > 0.0 else "negative"}
        return {"binary64": _canonical_float(scalar)}
    if isinstance(value, np.ndarray):
        raw_values = cast("Any", value).tolist()
        return _semantic_value(raw_values, allow_infinity=allow_infinity)
    if isinstance(value, Mapping):
        return {
            str(key): _semantic_value(item, allow_infinity=allow_infinity)
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [_semantic_value(item, allow_infinity=allow_infinity) for item in value]
    raise TypeError(
        f"Unsupported native scientific descriptor value: {type(value).__name__}"
    )


def _semantic_json(value: object, *, allow_infinity: bool = False) -> str:
    return json.dumps(
        _semantic_value(value, allow_infinity=allow_infinity),
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def _kernel_configuration(profile: Profile) -> str:
    settings = getattr(profile.pulse_sequence, "settings", None)
    if not hasattr(settings, "model_dump"):
        raise ValueError("Native profile kernels must expose immutable settings")
    # Pulse settings may use infinity as an explicit configuration sentinel
    # (for example, CEST's unbounded sweep width).  It is fingerprinted as a
    # stable token here; evaluated values and observations remain finite-only.
    return _semantic_json(
        settings.model_dump(mode="json"),
        allow_infinity=True,
    )


def _observation_metadata(profile: Profile) -> str:
    return _semantic_json(profile.data.metadata.tolist())


def _spectrometer_configuration(profile: Profile) -> str:
    return _semantic_json(profile.spectrometer.native_kernel_descriptor())


def _readonly(values: Array) -> Array:
    """Own immutable binary64 storage that cannot be made writable again."""
    copied = np.ascontiguousarray(values, dtype=np.float64)
    return np.frombuffer(copied.tobytes(), dtype=np.float64)


def _weighted_reduction(
    experimental: Array, calculation: Array, uncertainty: Array, *, square: bool
) -> float:
    """Left-to-right binary64 terms in canonical retained-observation order."""
    total = 0.0
    for exp, calc, err in zip(experimental, calculation, uncertainty, strict=True):
        left = float(calc) / float(err)
        right = left if square else float(exp) / float(err)
        term = left * right
        total = total + term
        if (
            not math.isfinite(left)
            or not math.isfinite(right)
            or not math.isfinite(total)
        ):
            return math.nan
    return total


def _normalization_factor(descriptor: ProfilePlan, unscaled: Array) -> float | None:
    if not descriptor.is_scaled:
        return 1.0
    retained = np.asarray(descriptor.mask, dtype=np.bool_)
    exp = np.asarray(descriptor.experimental_values, dtype=np.float64)[retained]
    err = np.asarray(descriptor.uncertainties, dtype=np.float64)[retained]
    calc = unscaled[retained]
    numerator = _weighted_reduction(exp, calc, err, square=False)
    denominator = _weighted_reduction(exp, calc, err, square=True)
    if (
        not math.isfinite(numerator)
        or not math.isfinite(denominator)
        or denominator == 0.0
    ):
        return None
    scale = numerator / denominator
    return scale if math.isfinite(scale) else None


def _parameterization_failure_validity(
    error: ParameterizationError,
) -> Literal["INVALID_REQUEST", "INVALID_TRIAL", "INVALID_PLAN_OR_BINDING"]:
    validity = _FAILURE_TAXONOMY["resolution"].get(error.code)
    if validity is None:
        raise RuntimeError("Undeclared parameterization failure")
    return cast(
        'Literal["INVALID_REQUEST", "INVALID_TRIAL", "INVALID_PLAN_OR_BINDING"]',
        validity,
    )


def _validate_plan_runtime_versions(plan: EvaluationPlan) -> None:
    if (
        plan.normalization_version != _NORMALIZATION_VERSION
        or plan.residual_version != _RESIDUAL_VERSION
        or plan.scalar_version != _SCALAR_VERSION
        or plan.masking_version != _MASKING_VERSION
        or plan.ordering_version != _ORDERING_VERSION
        or plan.failure_version != _FAILURE_VERSION
        or plan.diagnostics_version != _DIAGNOSTICS_VERSION
    ):
        raise ValueError("Evaluation plan has incompatible execution semantics")


def _record_failure_stage_validity(
    record: Mapping[str, object],
) -> tuple[FailureStage, FailureValidity]:
    stage = _record_string(record.get("stage"), "failure stage")
    validity = _record_string(record.get("validity"), "failure validity")
    if stage not in _FAILURE_TAXONOMY:
        raise ValueError("Unknown native evaluation-failure stage")
    if validity not in {
        "INVALID_REQUEST",
        "INVALID_TRIAL",
        "INVALID_PLAN_OR_BINDING",
        "IMPLEMENTATION_FAILURE",
    }:
        raise ValueError("Unknown native evaluation-failure validity")
    category = _record_string(record.get("category"), "failure category")
    if _FAILURE_TAXONOMY[stage].get(category) != validity:
        raise ValueError("Evaluation failure has incompatible stage and validity")
    return cast(FailureStage, stage), validity


def _record_exact_keys(
    record: Mapping[str, object], expected: set[str], name: str
) -> None:
    if set(record) != expected:
        raise ValueError(f"{name} record has unknown or missing fields")


def _record_frame_item(value: object) -> tuple[str, float]:
    pair = _record_list(value, field_name="frame item")
    if len(pair) != 2 or not isinstance(pair[0], str):
        raise TypeError("Evaluation-frame items must be [string, binary64] pairs")
    return pair[0], _record_binary64(pair[1], field_name="frame value")


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


def _record_string(value: object, field_name: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"Evaluation record {field_name} must be a string")
    return value


def _record_string_sequence(value: object, field_name: str) -> tuple[str, ...]:
    items = _record_list(value, field_name=field_name)
    result = tuple(_record_string(item, field_name) for item in items)
    if len(set(result)) != len(result):
        raise ValueError(f"Evaluation record {field_name} must be unique")
    return result


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
    resolved_ids: tuple[str, ...] = ()
    scalar_version: str = _SCALAR_VERSION
    normalization_version: str = _NORMALIZATION_VERSION
    residual_version: str = _RESIDUAL_VERSION
    masking_version: str = _MASKING_VERSION
    ordering_version: str = _ORDERING_VERSION
    failure_version: str = _FAILURE_VERSION
    diagnostics_version: str = _DIAGNOSTICS_VERSION
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
                    self.resolved_ids,
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
                    self.scalar_version,
                    self.masking_version,
                    self.ordering_version,
                    self.failure_version,
                    self.diagnostics_version,
                ),
            ),
        )

    @property
    def compatibility_identity(self) -> str:
        """Return the runtime compatibility identity sealed by this plan."""
        return _compatibility_identity(self)

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
            "resolved_ids": list(self.resolved_ids),
            "scalar_version": self.scalar_version,
            "normalization_version": self.normalization_version,
            "residual_version": self.residual_version,
            "masking_version": self.masking_version,
            "ordering_version": self.ordering_version,
            "failure_version": self.failure_version,
            "diagnostics_version": self.diagnostics_version,
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
        _record_exact_keys(
            record,
            {
                "schema_version",
                "parameterization_identity",
                "constraint_program_identity",
                "resolved_ids",
                "scalar_version",
                "normalization_version",
                "residual_version",
                "masking_version",
                "ordering_version",
                "failure_version",
                "diagnostics_version",
                "profiles",
                "identity",
            },
            "Evaluation-plan",
        )
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
        expected_offset = 0
        for descriptor in profiles:
            if (
                len(descriptor.uncertainties) != descriptor.observation_count
                or len(descriptor.mask) != descriptor.observation_count
                or descriptor.output_shape != (descriptor.observation_count,)
                or descriptor.observation_offset != expected_offset
                or descriptor.experiment_ordinal < 0
                or descriptor.profile_ordinal < 0
                or descriptor.identity
                != _identity(
                    "profile-plan",
                    (
                        descriptor.source_identity,
                        descriptor.experiment_ordinal,
                        descriptor.profile_ordinal,
                        descriptor.observation_offset,
                    ),
                )
            ):
                raise ValueError("Malformed native evaluation-profile semantics")
            expected_offset += descriptor.observation_count
        plan = cls(
            parameterization_identity=_record_string(
                record["parameterization_identity"], "parameterization identity"
            ),
            constraint_program_identity=_record_string(
                record["constraint_program_identity"], "constraint program identity"
            ),
            profiles=profiles,
            resolved_ids=_record_string_sequence(
                record["resolved_ids"], "resolved IDs"
            ),
            scalar_version=_record_string(record["scalar_version"], "scalar version"),
            normalization_version=_record_string(
                record["normalization_version"], "normalization version"
            ),
            residual_version=_record_string(
                record["residual_version"], "residual version"
            ),
            masking_version=_record_string(
                record["masking_version"], "masking version"
            ),
            ordering_version=_record_string(
                record["ordering_version"], "ordering version"
            ),
            failure_version=_record_string(
                record["failure_version"], "failure version"
            ),
            diagnostics_version=_record_string(
                record["diagnostics_version"], "diagnostics version"
            ),
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


@dataclass(frozen=True, slots=True)
class ResolvedEvaluationValues(Mapping[str, float]):
    """Evaluator-semantic resolved snapshot with no lifecycle authority."""

    parameterization_identity: str
    program_identity: str
    _items: tuple[tuple[str, float], ...]
    _index: Mapping[str, float] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        items = tuple(
            (key, _finite_evaluation_scalar(value)) for key, value in self._items
        )
        if len({key for key, _value in items}) != len(items):
            raise ValueError("Resolved evaluation values must have unique IDs")
        object.__setattr__(self, "_items", items)
        object.__setattr__(self, "_index", dict(items))

    def __iter__(self) -> Iterator[str]:
        return iter(self._index)

    def __len__(self) -> int:
        return len(self._index)

    def __getitem__(self, key: str) -> float:
        return self._index[key]

    def ordered_items(self) -> tuple[tuple[str, float], ...]:
        return self._items


def _profile_evaluation_from_record(
    record: object,
    descriptor: ProfilePlan,
) -> ProfileEvaluation:
    if not isinstance(record, Mapping):
        raise TypeError("Evaluation-result profile must be a mapping")
    values = cast("Mapping[str, object]", record)
    _record_exact_keys(
        values,
        {
            "profile_identity",
            "observation_offset",
            "observation_count",
            "residual_offset",
            "residual_count",
            "retained_observation_indices",
            "normalization_factor",
            "kernel_identity",
        },
        "Evaluation-result profile",
    )
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
        _record_binary64(
            values.get("normalization_factor"), field_name="normalization"
        ),
        descriptor.kernel_identity,
    )


def _compatibility_identity(plan: EvaluationPlan) -> str:
    return _identity(
        "evaluation-compatibility",
        tuple(item.kernel_identity for item in plan.profiles),
    )


def _validate_result_relationships(
    plan: EvaluationPlan,
    unscaled: Array,
    normalized: Array,
    residuals: Array,
    profiles: tuple[ProfileEvaluation, ...],
) -> None:
    residual_offset = 0
    for descriptor, item in zip(plan.profiles, profiles, strict=True):
        if (
            item.observation_offset != descriptor.observation_offset
            or item.observation_count != descriptor.observation_count
            or item.residual_offset != residual_offset
            or item.residual_count != len(descriptor.retained_observation_indices)
        ):
            raise ValueError("Evaluation result has inconsistent profile offsets")
        start = item.observation_offset
        stop = start + item.observation_count
        retained = np.asarray(item.retained_observation_indices, dtype=np.intp)
        expected_factor = _normalization_factor(descriptor, unscaled[start:stop])
        if expected_factor is None or _canonical_float(
            item.normalization_factor
        ) != _canonical_float(expected_factor):
            raise ValueError("Evaluation result has inconsistent normalization")
        with np.errstate(
            divide="raise",
            over="raise",
            under="ignore",
            invalid="raise",
        ):
            expected_normalized = item.normalization_factor * unscaled[start:stop]
            expected_residuals = (
                normalized[start:stop][retained]
                - np.asarray(descriptor.experimental_values)[retained]
            ) / np.asarray(descriptor.uncertainties)[retained]
        if not np.array_equal(
            normalized[start:stop], expected_normalized
        ) or not np.array_equal(
            residuals[residual_offset : residual_offset + item.residual_count],
            expected_residuals,
        ):
            raise ValueError("Evaluation result arrays do not match profile semantics")
        residual_offset += item.residual_count


@dataclass(frozen=True, slots=True)
class EvaluationResult:
    """The only successful native-evaluation outcome."""

    plan_identity: str
    parameterization_identity: str
    evaluator_compatibility_identity: str
    resolved_values: ResolvedEvaluationValues
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
                self.resolved_values.ordered_items(),
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
        record: dict[str, object] = {
            "schema_version": _SCHEMA_VERSION,
            "plan_identity": self.plan_identity,
            "parameterization_identity": self.parameterization_identity,
            "evaluator_compatibility_identity": self.evaluator_compatibility_identity,
            "resolved": {
                "program_identity": self.resolved_values.program_identity,
                "items": [
                    [param_id, _canonical_float(value)]
                    for param_id, value in self.resolved_values.ordered_items()
                ],
            },
            "unscaled_calculations": [
                _canonical_float(value) for value in self.unscaled_calculations
            ],
            "normalized_calculations": [
                _canonical_float(value) for value in self.normalized_calculations
            ],
            "residuals": [_canonical_float(value) for value in self.residuals],
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
                    "normalization_factor": _canonical_float(item.normalization_factor),
                    "kernel_identity": item.kernel_identity,
                }
                for item in self.profiles
            ],
        }
        record["identity"] = self.identity
        return record

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        plan: EvaluationPlan,
    ) -> EvaluationResult:
        """Restore a result only when its population matches a frozen plan."""
        _record_exact_keys(
            record,
            {
                "schema_version",
                "plan_identity",
                "parameterization_identity",
                "evaluator_compatibility_identity",
                "resolved",
                "unscaled_calculations",
                "normalized_calculations",
                "residuals",
                "profiles",
                "identity",
            },
            "Evaluation-result",
        )
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise ValueError("Unsupported native evaluation-result schema")
        if record.get("plan_identity") != plan.identity:
            raise ValueError("Evaluation result belongs to another plan")
        resolved_record = record.get("resolved")
        if not isinstance(resolved_record, Mapping):
            raise TypeError("Evaluation-result resolved values must be a mapping")
        resolved_values = cast("Mapping[str, object]", resolved_record)
        _record_exact_keys(
            resolved_values, {"program_identity", "items"}, "Evaluation-result resolved"
        )
        items = _record_list(resolved_values.get("items"), field_name="resolved.items")
        parameterization_identity = _record_string(
            record["parameterization_identity"], "parameterization identity"
        )
        resolved = ResolvedEvaluationValues(
            parameterization_identity,
            _record_string(
                resolved_values["program_identity"], "resolved program identity"
            ),
            tuple(_record_frame_item(item) for item in items),
        )
        if (
            parameterization_identity != plan.parameterization_identity
            or resolved.program_identity != plan.constraint_program_identity
            or tuple(key for key, _value in resolved.ordered_items())
            != plan.resolved_ids
        ):
            raise ValueError("Evaluation result has incompatible resolved values")
        profiles_record = _record_list(record.get("profiles"), field_name="profiles")
        if len(profiles_record) != len(plan.profiles):
            raise ValueError("Evaluation result has the wrong profile population")
        profiles = tuple(
            _profile_evaluation_from_record(item, descriptor)
            for item, descriptor in zip(profiles_record, plan.profiles, strict=True)
        )
        unscaled = _readonly(
            np.asarray(
                [
                    _record_binary64(value, field_name="unscaled calculation")
                    for value in _record_list(
                        record.get("unscaled_calculations"),
                        field_name="unscaled calculations",
                    )
                ]
            )
        )
        normalized = _readonly(
            np.asarray(
                [
                    _record_binary64(value, field_name="normalized calculation")
                    for value in _record_list(
                        record.get("normalized_calculations"),
                        field_name="normalized calculations",
                    )
                ]
            )
        )
        residuals = _readonly(
            np.asarray(
                [
                    _record_binary64(value, field_name="residual")
                    for value in _record_list(
                        record.get("residuals"), field_name="residuals"
                    )
                ]
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
        expected_compatibility = _compatibility_identity(plan)
        if record.get("evaluator_compatibility_identity") != expected_compatibility:
            raise ValueError("Evaluation result has incompatible evaluator binding")
        _validate_result_relationships(plan, unscaled, normalized, residuals, profiles)
        result = cls(
            plan.identity,
            parameterization_identity,
            expected_compatibility,
            resolved,
            unscaled,
            normalized,
            residuals,
            profiles,
        )
        if record.get("identity") != result.identity:
            raise ValueError("Evaluation-result fingerprint does not match its payload")
        return result


@dataclass(frozen=True, slots=True)
class EvaluationFailure:
    """Atomic typed non-success outcome; no usable vectors accompany it."""

    plan_identity: str
    parameterization_identity: str
    stage: Literal[
        "frame",
        "resolution",
        "projection",
        "cache",
        "kernel",
        "normalization",
        "residual",
        "result",
        "binding",
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

    @property
    def identity(self) -> str:
        return _identity(
            "evaluation-failure",
            (
                self.plan_identity,
                self.parameterization_identity,
                self.stage,
                self.category,
                self.validity,
                self.profile_identity,
                self.message,
            ),
        )

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
            "identity": self.identity,
        }

    @classmethod
    def from_record(
        cls,
        record: Mapping[str, object],
        plan: EvaluationPlan,
    ) -> EvaluationFailure:
        """Restore a sanitized failure only for the matching frozen plan."""
        _record_exact_keys(
            record,
            {
                "schema_version",
                "plan_identity",
                "parameterization_identity",
                "stage",
                "category",
                "validity",
                "profile_identity",
                "message",
                "identity",
            },
            "Evaluation-failure",
        )
        if record.get("schema_version") != _SCHEMA_VERSION:
            raise ValueError("Unsupported native evaluation-failure schema")
        if record.get("plan_identity") != plan.identity:
            raise ValueError("Evaluation failure belongs to another plan")
        if record.get("parameterization_identity") != plan.parameterization_identity:
            raise ValueError("Evaluation failure belongs to another parameterization")
        stage, validity = _record_failure_stage_validity(record)
        profile_identity = record.get("profile_identity")
        if profile_identity is not None and profile_identity not in {
            item.identity for item in plan.profiles
        }:
            raise ValueError("Evaluation failure names a profile outside its plan")
        category = _record_string(record.get("category"), "failure category")
        message = _record_string(record.get("message"), "failure message")
        if not category:
            raise ValueError("Evaluation failure category cannot be empty")
        if (stage in {"frame", "resolution", "result", "binding"}) != (
            profile_identity is None
        ):
            raise ValueError("Evaluation failure has invalid profile location")
        failure = cls(
            plan.identity,
            plan.parameterization_identity,
            stage,
            category,
            validity,
            None if profile_identity is None else str(profile_identity),
            message,
        )
        if record.get("identity") != failure.identity:
            raise ValueError(
                "Evaluation-failure fingerprint does not match its payload"
            )
        return failure


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


@dataclass(frozen=True, slots=True)
class _LocalParameterProjection(Mapping[str, float]):
    """The only resolved values a native profile calculation can observe."""

    _items: tuple[tuple[str, float], ...]
    _values: Mapping[str, float] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "_values", dict(self._items))

    def __iter__(self) -> Iterator[str]:
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)

    def __getitem__(self, key: str) -> float:
        return self._values[key]


@dataclass(slots=True)
class _NativeKernelCapability:
    """Private adapter that exposes only local values and copied metadata."""

    spectrometer: Spectrometer
    pulse_sequence: PulseSequence
    metadata: Array
    output_size: int

    @classmethod
    def from_profile(cls, profile: Profile) -> _NativeKernelCapability:
        return cls(
            profile.spectrometer.new_native_workspace(),
            deepcopy(profile.pulse_sequence),
            np.array(profile.data.metadata, copy=True),
            profile.data.exp.size,
        )

    def calculate(self, local_values: Mapping[str, float]) -> Array:
        self.spectrometer.update(dict(local_values))
        return self.pulse_sequence.calculate(
            self.spectrometer,
            Data(
                exp=np.zeros(self.output_size, dtype=np.float64),
                err=np.ones(self.output_size, dtype=np.float64),
                metadata=np.array(self.metadata, copy=True),
            ),
        )


@dataclass(slots=True)
class _ProfileWorkspace:
    template: _NativeKernelCapability
    profile: _NativeKernelCapability = field(init=False)
    cache: OrderedDict[tuple[float, ...], _CachedProfile] = field(
        default_factory=OrderedDict
    )

    def __post_init__(self) -> None:
        self.profile = deepcopy(self.template)

    def reset(self) -> None:
        self.profile = deepcopy(self.template)
        self.cache.clear()


@dataclass(slots=True)
class _Workspace:
    profile_workspaces: tuple[_ProfileWorkspace, ...]
    resolved_values: dict[str, float] | None = None
    hits: int = 0
    misses: int = 0
    poisoned: bool = False

    def clear_caches(self) -> None:
        for workspace in self.profile_workspaces:
            workspace.cache.clear()

    def clear_resolution(self) -> None:
        self.resolved_values = None

    def reset(self) -> None:
        for workspace in self.profile_workspaces:
            workspace.reset()
        self.clear_resolution()


@dataclass(frozen=True, slots=True)
class ResampledProfileBinding:
    """Exact root-kernel rows used by one transformed resampling profile."""

    root_profile_index: int
    root_observation_indices: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.root_profile_index < 0 or not self.root_observation_indices:
            raise ValueError("Resampled profile binding requires root observations")
        if any(index < 0 for index in self.root_observation_indices):
            raise ValueError("Resampled observation indices must be non-negative")


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
            parameterization_identity=parameterization.evaluator_identity,
            constraint_program_identity=parameterization.program.fingerprint,
            profiles=tuple(profile_plans),
            resolved_ids=parameterization.scope_ids,
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
        if plan.parameterization_identity != parameterization.evaluator_identity:
            raise ValueError("Evaluation plan belongs to another parameterization")
        if plan.constraint_program_identity != parameterization.program.fingerprint:
            raise ValueError("Evaluation plan belongs to another constraint program")
        _validate_plan_runtime_versions(plan)
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
                or descriptor.is_scaled != source.is_scaled
                or descriptor.experimental_values
                != tuple(float(value) for value in np.asarray(source.data.exp))
                or descriptor.uncertainties
                != tuple(float(value) for value in np.asarray(source.data.err))
                or descriptor.mask
                != tuple(bool(value) for value in np.asarray(source.data.mask))
            ):
                raise ValueError("Trusted kernel descriptor does not match frozen plan")
            if not set(descriptor.param_ids).issubset(parameterization.scope_ids):
                raise ValueError(
                    "Evaluation plan projects parameters outside its closure"
                )
            expected_offset += descriptor.observation_count
        if plan.resolved_ids != parameterization.scope_ids:
            raise ValueError("Evaluation plan has incompatible resolved-value ordering")
        self.plan = plan
        self._parameterization = parameterization
        # Copy construction-owned mutable machinery once.  Plans retain no
        # source objects, and callers cannot alter a future workspace by
        # mutating the legacy occurrence after it was bound.
        self._sources = tuple(deepcopy(source) for _e, _p, source in sources)
        self._templates = tuple(
            _NativeKernelCapability.from_profile(source) for source in self._sources
        )
        self.compatibility_identity = _compatibility_identity(plan)

    @classmethod
    def from_experiments(
        cls,
        experiments: Experiments,
        parameterization: ActiveParameterization,
    ) -> EvaluationEngine:
        plan, sources = _build_plan(experiments, parameterization)
        return cls(plan, parameterization, sources)

    def rebind_parameterization(
        self,
        parameterization: ActiveParameterization,
    ) -> EvaluationEngine:
        """Bind the frozen plan and copied kernels to a fresh state occurrence."""
        sources = tuple(
            (descriptor.experiment_ordinal, descriptor.profile_ordinal, source)
            for descriptor, source in zip(
                self.plan.profiles,
                self._sources,
                strict=True,
            )
        )
        return type(self)(self.plan, parameterization, sources)

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

    def _projected_population(
        self,
        profile_indices: Sequence[int],
    ) -> tuple[EvaluationPlan, tuple[tuple[int, int, Profile], ...]]:
        """Build the immutable plan and trusted sources for a profile projection."""
        indices = tuple(profile_indices)
        if (
            not indices
            or tuple(sorted(set(indices))) != indices
            or indices[-1] >= len(self.plan.profiles)
            or indices[0] < 0
        ):
            raise ValueError("Profile projection requires unique root-ordered indices")
        profiles: list[ProfilePlan] = []
        sources: list[tuple[int, int, Profile]] = []
        offset = 0
        for index in indices:
            root = self.plan.profiles[index]
            source = self._sources[index]
            profiles.append(
                ProfilePlan(
                    identity=_identity(
                        "profile-plan",
                        (
                            root.source_identity,
                            root.experiment_ordinal,
                            root.profile_ordinal,
                            offset,
                        ),
                    ),
                    source_identity=root.source_identity,
                    experiment_ordinal=root.experiment_ordinal,
                    profile_ordinal=root.profile_ordinal,
                    observation_offset=offset,
                    local_inputs=root.local_inputs,
                    is_scaled=root.is_scaled,
                    experimental_values=root.experimental_values,
                    uncertainties=root.uncertainties,
                    mask=root.mask,
                    kernel_identity=root.kernel_identity,
                    kernel_configuration=root.kernel_configuration,
                    spectrometer_configuration=root.spectrometer_configuration,
                    observation_metadata=root.observation_metadata,
                    output_shape=root.output_shape,
                )
            )
            sources.append((root.experiment_ordinal, root.profile_ordinal, source))
            offset += root.observation_count
        plan = EvaluationPlan(
            parameterization_identity=self.plan.parameterization_identity,
            constraint_program_identity=self.plan.constraint_program_identity,
            profiles=tuple(profiles),
            resolved_ids=self.plan.resolved_ids,
            scalar_version=self.plan.scalar_version,
            normalization_version=self.plan.normalization_version,
            residual_version=self.plan.residual_version,
            masking_version=self.plan.masking_version,
            ordering_version=self.plan.ordering_version,
            failure_version=self.plan.failure_version,
            diagnostics_version=self.plan.diagnostics_version,
        )
        return plan, tuple(sources)

    def project_profiles(self, profile_indices: Sequence[int]) -> EvaluationEngine:
        """Project complete Profiles in root order into an isolated child engine."""
        plan, sources = self._projected_population(profile_indices)
        return EvaluationEngine(plan, self._parameterization, sources)

    def project_plan(self, profile_indices: Sequence[int]) -> EvaluationPlan:
        """Project an immutable child plan without compiling its runtime engine."""
        plan, _sources = self._projected_population(profile_indices)
        return plan

    def resampled_observation_metadata(self, binding: ResampledProfileBinding) -> str:
        """Return the canonical metadata descriptor for exact selected root rows."""
        if binding.root_profile_index >= len(self._templates):
            raise ValueError("Resampled profile names a foreign root profile")
        metadata = self._templates[binding.root_profile_index].metadata
        if any(index >= len(metadata) for index in binding.root_observation_indices):
            raise ValueError("Resampled profile names a foreign root observation")
        return _semantic_json(metadata[list(binding.root_observation_indices)].tolist())

    def resampled(
        self,
        plan: EvaluationPlan,
        bindings: Sequence[ResampledProfileBinding],
    ) -> EvaluationEngine:
        """Bind a transformed plan to fresh copies of exact trusted root kernels."""
        if plan.parameterization_identity != self._parameterization.evaluator_identity:
            raise ValueError("Resampled plan belongs to another parameterization")
        if (
            plan.constraint_program_identity
            != self._parameterization.program.fingerprint
            or plan.resolved_ids != self._parameterization.scope_ids
            or len(plan.profiles) != len(bindings)
        ):
            raise ValueError("Resampled plan has incompatible native lineage")
        _validate_plan_runtime_versions(plan)
        templates: list[_NativeKernelCapability] = []
        expected_offset = 0
        for descriptor, binding in zip(plan.profiles, bindings, strict=True):
            if binding.root_profile_index >= len(self.plan.profiles):
                raise ValueError("Resampled binding names a foreign root profile")
            root = self.plan.profiles[binding.root_profile_index]
            if (
                descriptor.observation_offset != expected_offset
                or descriptor.observation_count != len(binding.root_observation_indices)
                or descriptor.local_inputs != root.local_inputs
                or descriptor.is_scaled != root.is_scaled
                or descriptor.kernel_identity != root.kernel_identity
                or descriptor.kernel_configuration != root.kernel_configuration
                or descriptor.spectrometer_configuration
                != root.spectrometer_configuration
                or descriptor.observation_metadata
                != self.resampled_observation_metadata(binding)
                or descriptor.output_shape != (descriptor.observation_count,)
            ):
                raise ValueError(
                    "Resampled plan differs from its trusted root-kernel binding"
                )
            template = deepcopy(self._templates[binding.root_profile_index])
            template.metadata = np.array(
                template.metadata[list(binding.root_observation_indices)], copy=True
            )
            template.output_size = descriptor.observation_count
            templates.append(template)
            expected_offset += descriptor.observation_count
        engine = object.__new__(EvaluationEngine)
        engine.plan = plan
        engine._parameterization = self._parameterization
        engine._sources = ()
        engine._templates = tuple(templates)
        engine.compatibility_identity = _compatibility_identity(plan)
        return engine

    def new_evaluator(self) -> BoundEvaluator:
        """Create an empty single-owner workspace from trusted implementations."""
        return BoundEvaluator(
            self.plan,
            self._parameterization,
            self.compatibility_identity,
            _Workspace(
                tuple(
                    _ProfileWorkspace(deepcopy(template))
                    for template in self._templates
                ),
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
        self._owner_pid = os.getpid()
        self._in_flight = False

    @property
    def cache_statistics(self) -> CacheStatistics:
        return CacheStatistics(
            self._workspace.hits,
            self._workspace.misses,
            sum(
                len(workspace.cache) for workspace in self._workspace.profile_workspaces
            ),
        )

    def _failure(
        self,
        stage: Literal[
            "frame",
            "resolution",
            "projection",
            "cache",
            "kernel",
            "normalization",
            "residual",
            "result",
            "binding",
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
        if _FAILURE_TAXONOMY[stage].get(category) != validity:
            raise RuntimeError("Undeclared native evaluation failure")
        if validity in {"INVALID_PLAN_OR_BINDING", "IMPLEMENTATION_FAILURE"}:
            self._workspace.poisoned = True
            self._workspace.clear_caches()
            self._workspace.clear_resolution()
        elif validity == "INVALID_TRIAL":
            self._workspace.reset()
        return EvaluationFailure(
            self.plan.identity,
            self._parameterization.evaluator_identity,
            stage,
            category,
            validity,
            profile_identity,
            message,
        )

    def _calculate_unscaled(
        self,
        descriptor: ProfilePlan,
        profile: _NativeKernelCapability,
        resolved: Mapping[str, float],
    ) -> Array | EvaluationFailure:
        """Run and validate the narrow unscaled profile kernel."""
        try:
            local = _LocalParameterProjection(
                tuple(
                    (local_name, resolved[param_id])
                    for local_name, param_id in descriptor.local_inputs
                )
            )
        except KeyError as error:
            return self._failure(
                "projection",
                "missing_local_parameter",
                "INVALID_PLAN_OR_BINDING",
                profile_identity=descriptor.identity,
                message=str(error),
            )
        except Exception as error:  # noqa: BLE001 - projection integrity fence
            return self._failure(
                "projection",
                "missing_local_parameter",
                "INVALID_PLAN_OR_BINDING",
                profile_identity=descriptor.identity,
                message=str(error),
            )
        try:
            unscaled = profile.calculate(local)
        except Exception as error:  # noqa: BLE001 - native kernel fault boundary
            return self._failure(
                "kernel",
                "kernel_exception",
                "IMPLEMENTATION_FAILURE",
                profile_identity=descriptor.identity,
                message=str(error),
            )
        if not isinstance(unscaled, np.ndarray):
            return self._failure(
                "kernel",
                "unexpected_container",
                "IMPLEMENTATION_FAILURE",
                profile_identity=descriptor.identity,
            )
        if unscaled.shape != (descriptor.observation_count,):
            return self._failure(
                "kernel",
                "unexpected_shape",
                "IMPLEMENTATION_FAILURE",
                profile_identity=descriptor.identity,
            )
        if unscaled.dtype != np.dtype(np.float64):
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
            scale = _normalization_factor(descriptor, unscaled)
            if scale is None:
                return self._failure(
                    "normalization",
                    "invalid_normalization",
                    "INVALID_TRIAL",
                    profile_identity=descriptor.identity,
                )
        else:
            scale = 1.0
        try:
            with np.errstate(
                divide="raise",
                over="raise",
                under="ignore",
                invalid="raise",
            ):
                normalized = _readonly(scale * unscaled)
                residuals = _readonly(
                    (normalized[retained] - exp[retained]) / err[retained]
                )
        except FloatingPointError:
            return self._failure(
                "residual",
                "non_finite_residual",
                "INVALID_TRIAL",
                profile_identity=descriptor.identity,
            )
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
        profile: _NativeKernelCapability,
        cache: OrderedDict[tuple[float, ...], _CachedProfile],
        resolved: Mapping[str, float],
    ) -> _CachedProfile | EvaluationFailure:
        try:
            key = tuple(resolved[param_id] for param_id in descriptor.param_ids)
        except Exception as error:  # noqa: BLE001 - cache integrity fence
            return self._failure(
                "cache",
                "cache_key_exception",
                "IMPLEMENTATION_FAILURE",
                profile_identity=descriptor.identity,
                message=str(error),
            )
        cached = cache.get(key)
        if cached is not None:
            cache.move_to_end(key)
            self._workspace.hits += 1
            return cached
        self._workspace.misses += 1
        unscaled = self._calculate_unscaled(descriptor, profile, resolved)
        if isinstance(unscaled, EvaluationFailure):
            return unscaled
        cached = self._normalize_profile(descriptor, unscaled)
        if isinstance(cached, EvaluationFailure):
            return cached
        cache[key] = cached
        cache.move_to_end(key)
        while len(cache) > _PROFILE_CACHE_CAPACITY:
            cache.popitem(last=False)
        return cached

    def _resolve_frame(
        self,
        frame: EvaluationFrame,
    ) -> dict[str, float] | EvaluationFailure:
        if get_ident() != self._owner:
            return self._failure(
                "binding", "workspace_owner_violation", "INVALID_REQUEST"
            )
        if self._workspace.poisoned:
            return self._failure(
                "binding", "workspace_poisoned", "INVALID_PLAN_OR_BINDING"
            )
        try:
            if (
                frame.parameterization_identity
                != self._parameterization.evaluator_identity
            ):
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
            resolved = self._parameterization._resolve_values(
                lifecycle_frame,
                self._workspace.resolved_values,
            )
        except ParameterizationError as error:
            return self._failure(
                "resolution",
                error.code,
                _parameterization_failure_validity(error),
                message=str(error),
            )
        except Exception as error:  # noqa: BLE001 - trusted-program boundary
            return self._failure(
                "resolution",
                "unexpected_resolution_error",
                "IMPLEMENTATION_FAILURE",
                message=str(error),
            )
        self._workspace.resolved_values = resolved
        return resolved

    def _evaluate_profiles(
        self,
        resolved: Mapping[str, float],
    ) -> list[_CachedProfile] | EvaluationFailure:
        results: list[_CachedProfile] = []
        for descriptor, workspace in zip(
            self.plan.profiles,
            self._workspace.profile_workspaces,
            strict=True,
        ):
            cached = self._cached_profile(
                descriptor,
                workspace.profile,
                workspace.cache,
                resolved,
            )
            if isinstance(cached, EvaluationFailure):
                return cached
            results.append(cached)
        return results

    def _complete_result(
        self,
        resolved: Mapping[str, float],
        cached_profiles: Sequence[_CachedProfile],
        residuals: Array,
    ) -> EvaluationResult:
        residual_offset = 0
        profiles: list[ProfileEvaluation] = []
        for descriptor, cached in zip(self.plan.profiles, cached_profiles, strict=True):
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
        return EvaluationResult(
            self.plan.identity,
            self._parameterization.evaluator_identity,
            self.compatibility_identity,
            ResolvedEvaluationValues(
                self._parameterization.evaluator_identity,
                self._parameterization.program.fingerprint,
                tuple(
                    (param_id, resolved[param_id])
                    for param_id in self._parameterization.scope_ids
                ),
            ),
            _readonly(
                np.concatenate([item.unscaled for item in cached_profiles])
                if cached_profiles
                else np.empty(0)
            ),
            _readonly(
                np.concatenate([item.normalized for item in cached_profiles])
                if cached_profiles
                else np.empty(0)
            ),
            residuals,
            tuple(profiles),
        )

    def _evaluate_impl(
        self,
        frame: EvaluationFrame,
        *,
        materialize: bool,
    ) -> EvaluationResult | Array | EvaluationFailure:
        """Run the shared scientific path, optionally materializing full evidence."""
        resolved = self._resolve_frame(frame)
        if isinstance(resolved, EvaluationFailure):
            return resolved
        cached_profiles = self._evaluate_profiles(resolved)
        if isinstance(cached_profiles, EvaluationFailure):
            return cached_profiles
        try:
            residuals = _readonly(
                np.concatenate([item.residuals for item in cached_profiles])
                if cached_profiles
                else np.empty(0)
            )
            if not materialize:
                return residuals
            return self._complete_result(resolved, cached_profiles, residuals)
        except Exception as error:  # noqa: BLE001 - result integrity fence
            return self._failure(
                "result",
                "result_assembly_exception",
                "IMPLEMENTATION_FAILURE",
                message=str(error),
            )

    def _fenced_evaluate(
        self,
        frame: EvaluationFrame,
        *,
        materialize: bool,
    ) -> EvaluationResult | Array | EvaluationFailure:
        if os.getpid() != self._owner_pid:
            return self._failure(
                "binding", "workspace_process_violation", "INVALID_PLAN_OR_BINDING"
            )
        if get_ident() != self._owner:
            return self._failure(
                "binding", "workspace_owner_violation", "INVALID_REQUEST"
            )
        if self._in_flight:
            return self._failure("binding", "workspace_reentrant", "INVALID_REQUEST")
        if self._workspace.poisoned:
            return self._failure(
                "binding", "workspace_poisoned", "INVALID_PLAN_OR_BINDING"
            )
        self._in_flight = True
        try:
            with np.errstate(
                divide="raise",
                over="raise",
                under="ignore",
                invalid="raise",
            ):
                return self._evaluate_impl(frame, materialize=materialize)
        except FloatingPointError as error:
            return self._failure(
                "residual", "non_finite_residual", "INVALID_TRIAL", message=str(error)
            )
        except Exception as error:  # noqa: BLE001 - complete native boundary
            return self._failure(
                "binding",
                "unexpected_evaluation_error",
                "IMPLEMENTATION_FAILURE",
                message=str(error),
            )
        finally:
            self._in_flight = False

    def evaluate(self, frame: EvaluationFrame) -> EvaluationResult | EvaluationFailure:
        """Fence one complete call under the declared single-owner contract."""
        return cast(
            "EvaluationResult | EvaluationFailure",
            self._fenced_evaluate(frame, materialize=True),
        )

    def evaluate_residuals(self, frame: EvaluationFrame) -> Array | EvaluationFailure:
        """Evaluate one ephemeral solver request without publication evidence."""
        return cast(
            "Array | EvaluationFailure",
            self._fenced_evaluate(frame, materialize=False),
        )

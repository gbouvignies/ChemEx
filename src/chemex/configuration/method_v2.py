from __future__ import annotations

from pathlib import Path
from typing import Any

from chemex.configuration.method_expressions import (
    parse_constraint,
    parse_strict_de_coordinate,
    parse_strict_grid_axis,
    parse_strict_selector,
)
from chemex.configuration.method_plan import (
    ConstrainAction,
    Constraint,
    DeSearch,
    FitAction,
    FixAction,
    FormatOrigin,
    GridSearch,
    McmcRequest,
    MethodFormatError,
    MethodPlan,
    ParameterSelector,
    ProfileSelection,
    ResamplingRequest,
    RoleAction,
    SourceRef,
    StatisticsPlan,
    StepPlan,
    profile_selection,
)


def _selector(text: str, source: SourceRef) -> ParameterSelector:
    return parse_strict_selector(text, source)


def _constraint(text: str, source: SourceRef) -> Constraint:
    return parse_constraint(text, source, _selector)


def _seed(value: object, source: SourceRef) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 0 <= value < 1 << 64
    ):
        raise MethodFormatError("SEED must be an unsigned 64-bit integer", source)
    return value


def _positive_int(value: object, name: str, source: SourceRef) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise MethodFormatError(f"{name} must be a positive integer", source)
    return value


def _string_list(value: object, source: SourceRef) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise MethodFormatError(f"{source.field} must be a list of strings", source)
    return [item for item in value if isinstance(item, str)]


def _selection(
    include: object,
    exclude: object,
    filename: Path,
    name: str,
) -> ProfileSelection:
    for field, value in (("INCLUDE", include), ("EXCLUDE", exclude)):
        valid_all = isinstance(value, str) and value.lower() in {"*", "all"}
        valid_list = isinstance(value, list) and all(
            isinstance(item, (int, str)) and not isinstance(item, bool)
            for item in value
        )
        if value is not None and not valid_all and not valid_list:
            raise MethodFormatError(
                f"{field} must be ALL or a list of profile names/residue numbers",
                SourceRef(filename, name, field),
            )
    return profile_selection(include, exclude)


def _mapping(
    value: object,
    allowed: set[str],
    source: SourceRef,
) -> dict[str, object]:
    if not isinstance(value, dict):
        raise MethodFormatError(f"{source.field} must be a TOML table", source)
    normalized: dict[str, object] = {}
    for raw_key, item in value.items():
        key = str(raw_key).lower()
        if key in normalized:
            raise MethodFormatError(f"Duplicate key {raw_key}", source)
        if key not in allowed:
            raise MethodFormatError(f"Unsupported v2 field {raw_key}", source)
        normalized[key] = item
    return normalized


def _resampling(
    settings: dict[str, object], name: str, source: SourceRef
) -> ResamplingRequest | None:
    request = settings.get(name)
    if request is None:
        return None
    request_source = SourceRef(
        source.filename, source.step, f"STATISTICS.{name.upper()}"
    )
    if isinstance(request, int) and not isinstance(request, bool):
        return ResamplingRequest(_positive_int(request, "REPLICATES", request_source))
    expanded = _mapping(request, {"replicates", "seed"}, request_source)
    if "replicates" not in expanded:
        raise MethodFormatError("Expanded request requires REPLICATES", request_source)
    seed = expanded.get("seed")
    return ResamplingRequest(
        _positive_int(expanded["replicates"], "REPLICATES", request_source),
        _seed(seed, request_source) if seed is not None else None,
    )


def _burn(value: object, steps: int, source: SourceRef) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise MethodFormatError("BURN must be a nonnegative integer", source)
    if value >= steps:
        raise MethodFormatError("BURN must be smaller than STEPS", source)
    return value


def _mcmc(value: object, source: SourceRef) -> McmcRequest | None:
    if value is None:
        return None
    mcmc_source = SourceRef(source.filename, source.step, "STATISTICS.MCMC")
    if isinstance(value, int) and not isinstance(value, bool):
        return McmcRequest(_positive_int(value, "STEPS", mcmc_source))
    expanded = _mapping(value, {"steps", "burn", "seed"}, mcmc_source)
    if "steps" not in expanded:
        raise MethodFormatError("Expanded MCMC requires STEPS", mcmc_source)
    steps = _positive_int(expanded["steps"], "STEPS", mcmc_source)
    seed = expanded.get("seed")
    return McmcRequest(
        steps,
        _burn(expanded.get("burn"), steps, mcmc_source),
        _seed(seed, mcmc_source) if seed is not None else None,
    )


def _statistics(value: object, source: SourceRef) -> StatisticsPlan:
    settings = _mapping(value, {"mc", "bs", "bsn", "mcmc"}, source)
    return StatisticsPlan(
        _resampling(settings, "mc", source),
        _resampling(settings, "bs", source),
        _resampling(settings, "bsn", source),
        _mcmc(settings.get("mcmc"), source),
    )


def _grid_search(value: object, filename: Path, name: str) -> GridSearch:
    source = SourceRef(filename, name, "SEARCH.GRID")
    settings = _mapping(value, {"axes"}, source)
    if "axes" not in settings:
        raise MethodFormatError("SEARCH.GRID requires AXES", source)
    axes = _string_list(settings["axes"], SourceRef(filename, name, "SEARCH.GRID.AXES"))
    if not axes:
        raise MethodFormatError(
            "SEARCH.GRID.AXES must contain at least one entry",
            SourceRef(filename, name, "SEARCH.GRID.AXES"),
        )
    return GridSearch(
        tuple(
            parse_strict_grid_axis(
                text, SourceRef(filename, name, "SEARCH.GRID.AXES", index)
            )
            for index, text in enumerate(axes)
        )
    )


def _de_search(value: object, filename: Path, name: str) -> DeSearch:
    source = SourceRef(filename, name, "SEARCH.DE")
    settings = _mapping(value, {"seed", "coordinates"}, source)
    if "seed" not in settings or "coordinates" not in settings:
        raise MethodFormatError("SEARCH.DE requires SEED and COORDINATES", source)
    coordinates = _string_list(
        settings["coordinates"], SourceRef(filename, name, "SEARCH.DE.COORDINATES")
    )
    if not coordinates:
        raise MethodFormatError(
            "SEARCH.DE.COORDINATES must contain at least one entry",
            SourceRef(filename, name, "SEARCH.DE.COORDINATES"),
        )
    return DeSearch(
        _seed(settings["seed"], SourceRef(filename, name, "SEARCH.DE.SEED")),
        tuple(
            parse_strict_de_coordinate(
                text, SourceRef(filename, name, "SEARCH.DE.COORDINATES", index)
            )
            for index, text in enumerate(coordinates)
        ),
    )


def _search(value: object, filename: Path, name: str) -> GridSearch | DeSearch:
    source = SourceRef(filename, name, "SEARCH")
    settings = _mapping(value, {"grid", "de"}, source)
    if len(settings) != 1:
        raise MethodFormatError("SEARCH must contain exactly one of GRID or DE", source)
    if "grid" in settings:
        return _grid_search(settings["grid"], filename, name)
    return _de_search(settings["de"], filename, name)


def _role_action(
    value: object, filename: Path, name: str, action_index: int
) -> RoleAction:
    action_source = SourceRef(filename, name, "ROLES", action_index)
    settings = _mapping(value, {"fix", "fit", "constrain"}, action_source)
    if len(settings) != 1:
        keys = ", ".join(str(key).upper() for key in settings)
        raise MethodFormatError(
            f"Each ROLES action must contain exactly one key; found {keys}",
            action_source,
        )
    key, value = next(iter(settings.items()))
    texts = _string_list(value, action_source)
    if not texts:
        raise MethodFormatError("ROLES actions cannot be empty", action_source)
    source = SourceRef(filename, name, f"ROLES.{key}", action_index)
    if key == "constrain":
        return ConstrainAction(
            tuple(
                _constraint(
                    text,
                    SourceRef(
                        filename, name, f"ROLES[{action_index}].CONSTRAIN[{index}]"
                    ),
                )
                for index, text in enumerate(texts)
            ),
            source,
        )
    selectors = tuple(
        _selector(
            text,
            SourceRef(filename, name, f"ROLES[{action_index}].{key.upper()}[{index}]"),
        )
        for index, text in enumerate(texts)
    )
    return (
        FixAction(selectors, source) if key == "fix" else FitAction(selectors, source)
    )


def _role_actions(value: object, filename: Path, name: str) -> tuple[RoleAction, ...]:
    if not isinstance(value, list):
        raise MethodFormatError(
            "ROLES must be an ordered array of single-key actions",
            SourceRef(filename, name, "ROLES"),
        )
    if not value:
        raise MethodFormatError(
            "ROLES must contain at least one action",
            SourceRef(filename, name, "ROLES"),
        )
    return tuple(
        _role_action(action, filename, name, index)
        for index, action in enumerate(value)
    )


def _roles_from(
    value: object, filename: Path, name: str, earlier: set[str]
) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or value not in earlier:
        raise MethodFormatError(
            "ROLES_FROM must name one unique earlier step",
            SourceRef(filename, name, "ROLES_FROM"),
        )
    return value


def _step(
    filename: Path,
    name: str,
    settings: dict[str, Any],
    earlier: set[str],
) -> StepPlan:
    normalized = _mapping(
        settings,
        {"include", "exclude", "roles", "roles_from", "search", "statistics"},
        SourceRef(filename, name, "<step>"),
    )
    return StepPlan(
        name=name,
        selection=_selection(
            normalized.get("include"), normalized.get("exclude"), filename, name
        ),
        roles_from=_roles_from(normalized.get("roles_from"), filename, name, earlier),
        role_actions=(
            _role_actions(normalized["roles"], filename, name)
            if "roles" in normalized
            else ()
        ),
        search=(
            _search(normalized["search"], filename, name)
            if "search" in normalized
            else None
        ),
        statistics=(
            _statistics(
                normalized["statistics"], SourceRef(filename, name, "STATISTICS")
            )
            if "statistics" in normalized
            else None
        ),
    )


def adapt_v2(raw_steps: list[tuple[Path, str, dict[str, Any]]]) -> MethodPlan:
    steps: list[StepPlan] = []
    earlier: set[str] = set()
    for filename, name, settings in raw_steps:
        steps.append(_step(filename, name, settings, earlier))
        earlier.add(name)
    return MethodPlan(FormatOrigin.V2, tuple(steps))

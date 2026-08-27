"""Publish raw factors and exact numerical profiles for profiled GRID."""

from __future__ import annotations

import re
import shutil
from collections.abc import Callable, Mapping
from functools import partial
from pathlib import Path

import numpy as np

from chemex.messages import GridOutputProgressReporter, console
from chemex.optimize.profiled_grid import ProfiledGridOutcome, ProfiledGridSurface
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import SealedParameterModel
from chemex.parameters.sealed import parameter_name_from_definition
from chemex.plotters.grid import GridResult, plot_grid_1d, plot_grid_2d

type _ParameterToken = Callable[[str], str]


def _parameter_token(
    param_id: str,
    parameter_names: Mapping[str, ParamName],
) -> str:
    folder = parameter_names[param_id].folder.lower()
    return re.sub(r"[^a-z0-9_.-]+", "_", folder).strip("_") or param_id.strip("_")


def _write_factor_tables(
    outcome: ProfiledGridOutcome,
    factors_path: Path,
    token: _ParameterToken,
) -> None:
    aggregate = outcome.aggregate
    if aggregate is None:
        raise RuntimeError("Accepted profiled GRID output lacks its aggregate")
    for result, selected_ordinal in zip(
        outcome.factors,
        aggregate.selection.factor_point_ordinals,
        strict=True,
    ):
        axis_tokens = tuple(token(param_id) for param_id in result.factor.grid_ids)
        axis_suffix = "__".join(axis_tokens) if axis_tokens else "constant"
        profile_keys = result.factor.profile_keys
        profile_label = ";".join(
            f"experiment_{experiment + 1}:profile_{profile + 1}"
            for experiment, profile in profile_keys
        )
        profile_suffix = (
            "-".join(
                f"e{experiment + 1:02d}p{profile + 1:03d}"
                for experiment, profile in profile_keys
            )
            if len(profile_keys) <= 4
            else f"{len(profile_keys)}_profiles"
        )
        suffix = "__".join(filter(None, (axis_suffix, profile_suffix)))
        filename = (
            factors_path / f"factor_{result.factor.ordinal + 1:02d}__{suffix}.tsv"
        )
        headers = (
            "point",
            "profiles",
            *(token(param_id) for param_id in result.factor.grid_ids),
            "chi_square",
            "status",
            "selected",
            "objective_evaluations",
            "failure",
        )
        with filename.open("w", encoding="utf-8") as output:
            output.write("\t".join(headers) + "\n")
            for point in result.points:
                failure = (
                    ""
                    if point.failure is None
                    else re.sub(r"[\t\r\n]+", " ", point.failure)
                )
                output.write(
                    "\t".join(
                        (
                            str(point.ordinal),
                            profile_label,
                            *(repr(value) for _param_id, value in point.axis_items),
                            "" if point.chi_square is None else repr(point.chi_square),
                            point.status.value,
                            str(point.ordinal == selected_ordinal).lower(),
                            str(point.objective_evaluations),
                            failure,
                        )
                    )
                    + "\n"
                )


def _write_surface(
    destination: Path,
    surface: ProfiledGridSurface,
    token: _ParameterToken,
    selected_grid: Mapping[str, float],
) -> None:
    axis_ids = surface.axis_ids
    axis_values = surface.axis_values
    chisqr = surface.chi_square
    with destination.open("w", encoding="utf-8") as output:
        output.write(
            "\t".join((*map(token, axis_ids), "chi_square", "selected")) + "\n"
        )
        for indices in np.ndindex(chisqr.shape):
            values = tuple(
                axis_values[index][coordinate]
                for index, coordinate in enumerate(indices)
            )
            is_selected = all(
                selected_grid[param_id] == value
                for param_id, value in zip(axis_ids, values, strict=True)
            )
            output.write(
                "\t".join(
                    (
                        *(repr(value) for value in values),
                        repr(float(chisqr[indices])),
                        str(is_selected).lower(),
                    )
                )
                + "\n"
            )


def write_grid_output(
    outcome: ProfiledGridOutcome,
    path: Path,
    *,
    parameter_model: SealedParameterModel,
    accepted_values: Mapping[str, float],
    progress: GridOutputProgressReporter | None = None,
) -> None:
    """Publish one complete, atomically replaced GRID product tree."""
    aggregate = outcome.aggregate
    accepted = outcome.accepted_result
    if aggregate is None or accepted is None:
        raise RuntimeError("Accepted profiled GRID output lacks its aggregate")
    parameter_names = {
        definition.param_id: parameter_name_from_definition(definition)
        for definition in parameter_model.definitions
    }
    token = partial(_parameter_token, parameter_names=parameter_names)
    grid_path = path / "Grid"
    staging = path / ".Grid.tmp"
    shutil.rmtree(staging, ignore_errors=True)
    factors_path = staging / "Factors"
    profiles_1d_path = staging / "Profiles" / "1D"
    profiles_2d_path = staging / "Profiles" / "2D"
    factors_path.mkdir(parents=True, exist_ok=True)
    profiles_1d_path.mkdir(parents=True, exist_ok=True)
    profiles_2d_path.mkdir(parents=True, exist_ok=True)

    progress = progress or GridOutputProgressReporter(console)
    progress.start_writing()
    _write_factor_tables(outcome, factors_path, token)
    selected_grid = dict(aggregate.selection.grid_items)
    grids_1d: list[GridResult] = []
    for param_id, surface in aggregate.profiles_1d.items():
        _write_surface(
            profiles_1d_path / f"{token(param_id)}.tsv",
            surface,
            token,
            selected_grid,
        )
        grids_1d.append(
            GridResult(
                {param_id: np.asarray(surface.axis_values[0], dtype=np.float64)},
                surface.chi_square,
            )
        )
    grids_2d: list[GridResult] = []
    for param_ids, surface in aggregate.profiles_2d.items():
        _write_surface(
            profiles_2d_path / f"{token(param_ids[0])}__{token(param_ids[1])}.tsv",
            surface,
            token,
            selected_grid,
        )
        grids_2d.append(
            GridResult(
                {
                    param_id: np.asarray(values, dtype=np.float64)
                    for param_id, values in zip(
                        surface.axis_ids, surface.axis_values, strict=True
                    )
                },
                surface.chi_square,
            )
        )
    with (staging / "summary.toml").open("w", encoding="utf-8") as output:
        output.write('status = "complete"\n')
        output.write(f"selected_chi_square = {accepted.chi_square!r}\n")
        output.write(f"factor_count = {len(outcome.factors)}\n")
        for param_id, value in aggregate.selection.grid_items:
            output.write("\n[[selected_axes]]\n")
            output.write(f"parameter_id = {param_id!r}\n")
            output.write(f"name = {str(parameter_names[param_id])!r}\n")
            output.write(f"value = {value!r}\n")
    progress.finish_writing()

    progress.start_plotting()
    plot_grid_1d(
        grids_1d,
        staging,
        parameter_names=parameter_names,
        accepted_values=accepted_values,
    )
    plot_grid_2d(
        grids_2d,
        staging,
        parameter_names=parameter_names,
        accepted_values=accepted_values,
    )
    progress.finish_plotting()

    shutil.rmtree(grid_path, ignore_errors=True)
    staging.replace(grid_path)

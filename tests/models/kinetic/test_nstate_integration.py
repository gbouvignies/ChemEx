"""Small production-stack tests for generic N-state model authority."""

from __future__ import annotations

import math
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pytest

from chemex.configuration.methods import Method, Selection
from chemex.configuration.parameters import (
    DefaultListType,
    DefaultSetting,
    read_defaults,
)
from chemex.experiments.builder import Experiments, build_experiments
from chemex.parameters.feasible_coordinates import compile_feasible_coordinates
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import (
    ConstraintDomainError,
    NoParameterMatchError,
    ParameterRole,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.printers.parameters import classify_parameters
from chemex.run_info import serialize_parameter_file
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parents[3]
EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
PROFILE = SpinSystem.from_name("G2N-HN")


def _build_session(
    model_name: str,
    defaults: DefaultListType | None = None,
) -> tuple[AnalysisSession, Experiments]:
    session = AnalysisSession.create()
    session.set_model(model_name)
    experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[PROFILE], exclude=None),
        session=session,
    )
    session.parameters.set_defaults([] if defaults is None else defaults)
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    return session, experiments


def _names_by_id(session: AnalysisSession) -> Mapping[str, str]:
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    return {
        definition.param_id: definition.name
        for definition in parameter_model.definitions
    }


@pytest.mark.parametrize(
    ("model_name", "expected_kex"),
    (
        ("4st", {"KEX_AB", "KEX_AC", "KEX_AD", "KEX_BC", "KEX_BD", "KEX_CD"}),
        ("4st_linear", {"KEX_AB", "KEX_BC", "KEX_CD"}),
        ("4st_fork", {"KEX_AB", "KEX_AC", "KEX_AD"}),
    ),
)
def test_named_topologies_reach_the_sealed_parameter_model(
    model_name: str,
    expected_kex: set[str],
) -> None:
    session, experiments = _build_session(model_name)
    names = set(_names_by_id(session).values())

    assert {name for name in names if name.startswith("KEX_")} == expected_kex
    resolved = session.resolve_current_values(experiments.param_ids)
    for name in expected_kex:
        edge = name.removeprefix("KEX_")
        left, right = edge.lower()
        p_left = resolved[f"__P{left.upper()}"]
        p_right = resolved[f"__P{right.upper()}"]
        kex = resolved[f"__{name}"]
        forward = resolved[f"__K{edge}"]
        reverse = resolved[f"__K{edge[::-1]}"]
        expected_forward = kex * p_right / (p_left + p_right)
        expected_reverse = kex * p_left / (p_left + p_right)

        assert forward == pytest.approx(expected_forward)
        assert reverse == pytest.approx(expected_reverse)
        assert forward + reverse == pytest.approx(kex)
        assert p_left * forward == pytest.approx(p_right * reverse)


def test_invalid_simplex_fails_during_resolution_before_profile_calculation() -> None:
    session, experiments = _build_session("4st")
    snapshot = session.analysis_values.snapshot()
    session.analysis_values.commit(
        {"__PB": 0.7, "__PC": 0.7},
        expected=snapshot,
        scope=("__PB", "__PC"),
    )

    with pytest.raises(ConstraintDomainError) as error:
        session.resolve_current_values(experiments.param_ids)

    assert error.value.context["function_id"] == "population_complement"


def test_negative_constrained_kex_fails_before_negative_directional_rates() -> None:
    session, experiments = _build_session("4st")
    parameterization = session.compile_parameterization(
        Method(constraints=["[KEX_AB] = -1.0"]),
        experiments.param_ids,
    )

    with pytest.raises(ConstraintDomainError) as error:
        parameterization.resolve(
            parameterization.frame_from_snapshot(session.analysis_values.snapshot())
        )

    assert error.value.context["function_id"] == "pair_rates"


def test_absent_edge_method_selector_uses_normal_no_match_diagnostic() -> None:
    session, experiments = _build_session("4st_linear")

    with pytest.raises(NoParameterMatchError):
        session.compile_parameterization(Method(fit=["KEX_AD"]), experiments.param_ids)


def test_resolved_rates_reach_exchange_matrices_with_conservative_orientation() -> None:
    session, experiments = _build_session("4st")
    snapshot = session.analysis_values.snapshot()
    session.analysis_values.commit(
        {"__PB": 0.2, "__PC": 0.3, "__PD": 0.1, "__KEX_AB": 360.0},
        expected=snapshot,
        scope=("__PB", "__PC", "__PD", "__KEX_AB"),
    )
    resolved = session.resolve_current_values(experiments.param_ids)
    profile = next(iter(experiments)).profiles[0]

    profile.calculate_from_values(resolved)

    basis = profile.spectrometer.basis
    exchange = sum(
        resolved[param_id] * basis.matrices[name.lower()]
        for param_id, name in _names_by_id(session).items()
        if name.startswith("K") and not name.startswith("KEX_")
    )
    block_size = len(basis)
    assert resolved["__KAB"] == pytest.approx(120.0)
    assert resolved["__KBA"] == pytest.approx(240.0)
    assert exchange[block_size, 0] == pytest.approx(120.0)
    assert exchange[0, block_size] == pytest.approx(240.0)
    np.testing.assert_allclose(exchange.sum(axis=0), 0.0, atol=1.0e-13)


def test_complete_model_cycles_are_consistent_with_population_authority() -> None:
    session, experiments = _build_session("4st")
    snapshot = session.analysis_values.snapshot()
    session.analysis_values.commit(
        {"__PB": 0.2, "__PC": 0.3, "__PD": 0.1},
        expected=snapshot,
        scope=("__PB", "__PC", "__PD"),
    )
    values = session.resolve_current_values(experiments.param_ids)

    triangle = (
        values["__KAB"]
        / values["__KBA"]
        * values["__KBC"]
        / values["__KCB"]
        * values["__KCA"]
        / values["__KAC"]
    )
    four_state_cycle = (
        values["__KAB"]
        / values["__KBA"]
        * values["__KBD"]
        / values["__KDB"]
        * values["__KDC"]
        / values["__KCD"]
        * values["__KCA"]
        / values["__KAC"]
    )
    assert triangle == pytest.approx(1.0)
    assert four_state_cycle == pytest.approx(1.0)


def test_disconnected_five_state_mixture_preserves_populations_and_is_finite() -> None:
    session, experiments = _build_session("5st_linear")
    snapshot = session.analysis_values.snapshot()
    session.analysis_values.commit(
        {
            "__PB": 0.2,
            "__PC": 0.1,
            "__PD": 0.25,
            "__PE": 0.25,
            "__KEX_BC": 0.0,
        },
        expected=snapshot,
        scope=("__PB", "__PC", "__PD", "__PE", "__KEX_BC"),
    )
    resolved = session.resolve_current_values(experiments.param_ids)
    profile = next(iter(experiments)).profiles[0]
    calculation = profile.calculate_from_values(resolved)

    assert resolved["__KBC"] == resolved["__KCB"] == 0.0
    assert math.fsum(resolved[f"__P{state}"] for state in "ABCDE") == pytest.approx(1.0)
    assert resolved["__PA"] + resolved["__PB"] == pytest.approx(0.4)
    assert math.fsum(resolved[f"__P{state}"] for state in "CDE") == pytest.approx(0.6)
    assert np.all(np.isfinite(calculation))

    basis = profile.spectrometer.basis
    exchange = sum(
        resolved[param_id] * basis.matrices[name.lower()]
        for param_id, name in _names_by_id(session).items()
        if name.startswith("K") and not name.startswith("KEX_")
    )
    np.testing.assert_allclose(exchange.sum(axis=0), 0.0, atol=1.0e-13)

    profile.update_spectrometer_from_values(resolved)
    start = profile.spectrometer.get_start_magnetization(("iz",))
    initialized = {
        state: float(
            (profile.spectrometer.basis.vectors[f"iz_{state}"].T @ start).item()
        )
        for state in "abcde"
    }
    assert initialized == pytest.approx(
        {state: resolved[f"__P{state.upper()}"] for state in "abcde"}
    )


def test_simplex_feasible_coordinates_respect_fixed_populations_and_boundaries() -> (
    None
):
    session, experiments = _build_session("4st")
    snapshot = session.analysis_values.snapshot()
    committed = session.analysis_values.commit(
        {"__PB": 0.2, "__PC": 0.1, "__PD": 0.3},
        expected=snapshot,
        scope=("__PB", "__PC", "__PD"),
    )
    parameterization = session.compile_parameterization(
        Method(fix=["PB", "PC", "PD"], fit=["PC"]),
        experiments.param_ids,
    )
    frame = parameterization.frame_from_snapshot(committed)
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    controlled_ids = tuple(
        param_id
        for param_id in parameterization.independent_ids
        if parameterization.role(param_id) is ParameterRole.FIT
    )
    lower = tuple(
        parameter_model.configuration[param_id].lower_bound
        for param_id in controlled_ids
    )
    upper = tuple(
        parameter_model.configuration[param_id].upper_bound
        for param_id in controlled_ids
    )

    chart = compile_feasible_coordinates(
        parameterization,
        frame,
        controlled_ids,
        lower,
        upper,
    )

    assert chart is not None
    assert chart.has_coordinate_transform
    assert not chart.uses_private_relaxation_coordinates
    pc_index = controlled_ids.index("__PC")
    lower_vector = list(chart.solver_start)
    upper_vector = list(chart.solver_start)
    lower_vector[pc_index] = 0.0
    upper_vector[pc_index] = 1.0
    lower_point = chart.decode(lower_vector)
    upper_point = chart.decode(upper_vector)
    lower_values = parameterization.resolve(lower_point.frame)
    upper_values = parameterization.resolve(upper_point.frame)
    assert (lower_values["__PB"], lower_values["__PC"], lower_values["__PD"]) == (
        0.2,
        0.0,
        0.3,
    )
    assert lower_values["__PA"] == 0.5
    assert (upper_values["__PB"], upper_values["__PC"], upper_values["__PD"]) == (
        0.2,
        0.5,
        0.3,
    )
    assert upper_values["__PA"] == 0.0


def test_simplex_coordinates_respect_joint_fit_bounds_and_held_mass() -> None:
    session, experiments = _build_session(
        "4st",
        [
            (ParamName.from_section("PB"), DefaultSetting(0.2, 0.1, 0.3)),
            (ParamName.from_section("PC"), DefaultSetting(0.3, 0.2, 0.4)),
            (ParamName.from_section("PD"), DefaultSetting(0.25)),
        ],
    )
    parameterization = session.compile_parameterization(
        Method(fix=["PB", "PC", "PD"], fit=["PB", "PC"]),
        experiments.param_ids,
    )
    frame = parameterization.frame_from_snapshot(session.analysis_values.snapshot())
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    controlled_ids = tuple(
        param_id
        for param_id in parameterization.independent_ids
        if parameterization.role(param_id) is ParameterRole.FIT
    )
    bounds = parameter_model.configuration
    chart = compile_feasible_coordinates(
        parameterization,
        frame,
        controlled_ids,
        tuple(bounds[param_id].lower_bound for param_id in controlled_ids),
        tuple(bounds[param_id].upper_bound for param_id in controlled_ids),
    )

    assert chart is not None
    lower_vector = list(chart.solver_start)
    upper_vector = list(chart.solver_start)
    for name in ("__PB", "__PC"):
        index = controlled_ids.index(name)
        lower_vector[index] = 0.0
        upper_vector[index] = 1.0
    lower = parameterization.resolve(chart.decode(lower_vector).frame)
    upper = parameterization.resolve(chart.decode(upper_vector).frame)

    assert (
        lower["__PB"],
        lower["__PC"],
        lower["__PD"],
        lower["__PA"],
    ) == pytest.approx((0.1, 0.2, 0.25, 0.45))
    assert (
        upper["__PB"],
        upper["__PC"],
        upper["__PD"],
        upper["__PA"],
    ) == pytest.approx((0.3, 0.4, 0.25, 0.05))


@pytest.mark.parametrize("model_name", ("4st", "4st_linear", "4st_fork"))
def test_restart_and_output_expose_only_structural_public_coordinates(
    model_name: str,
) -> None:
    session, experiments = _build_session(model_name)
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    parameterization = session.compile_parameterization(Method(), experiments.param_ids)
    snapshot = session.analysis_values.snapshot()
    resolved = parameterization.resolve(parameterization.frame_from_snapshot(snapshot))
    restart = serialize_parameter_file(parameter_model, snapshot, state_kind="restart")
    classified = classify_parameters(
        parameter_model,
        resolved,
        parameterization,
        tuple(
            param_id
            for param_id in parameterization.independent_ids
            if parameterization.role(param_id) is ParameterRole.FIT
        ),
    )
    restart_keys = {
        line.split("=", 1)[0].strip().strip('"')
        for line in restart.splitlines()
        if "=" in line
    }
    constrained_ids = {
        parameter.param_id
        for parameter in (
            *classified.constrained.global_.values(),
            *classified.constrained.local.values(),
        )
    }
    names = _names_by_id(session)

    assert {"PB", "PC", "PD"} <= restart_keys
    assert {name for name in restart_keys if name.startswith("KEX_")} == {
        name for name in names.values() if name.startswith("KEX_")
    }
    assert "PA" not in restart_keys
    assert not {
        name
        for name in restart_keys
        if name.startswith("K") and not name.startswith("KEX_")
    }
    assert names.keys() >= constrained_ids
    assert {names[param_id] for param_id in constrained_ids} >= {"PA"}
    kinetic_constrained_names = {
        names[param_id]
        for param_id in constrained_ids
        if names[param_id] == "PA" or names[param_id].startswith("K")
    }
    assert all(
        name == "PA" or (name.startswith("K") and not name.startswith("KEX_"))
        for name in kinetic_constrained_names
    )


def test_legacy_restart_preserves_zero_values_but_not_historical_fixed_roles(
    tmp_path: Path,
) -> None:
    session, experiments = _build_session("4st")
    snapshot = session.analysis_values.snapshot()
    zeroed = session.analysis_values.commit(
        {"__KEX_AD": 0.0, "__KEX_BD": 0.0},
        expected=snapshot,
        scope=("__KEX_AD", "__KEX_BD"),
    )
    restart_path = tmp_path / "restart.toml"
    restart_path.write_text(
        serialize_parameter_file(
            session.parameter_factory.sealed_parameter_model,
            zeroed,
            state_kind="restart",
        )
    )

    continued = AnalysisSession.create()
    continued.set_model("4st")
    continued_experiments = build_experiments(
        [EXPERIMENT],
        Selection(include=[PROFILE], exclude=None),
        session=continued,
    )
    continued.parameters.set_defaults(read_defaults([restart_path]))
    assert continued.try_build_analysis_values()
    ordinary = continued.compile_parameterization(
        Method(),
        continued_experiments.param_ids,
    )
    migrated = continued.compile_parameterization(
        Method(fix=["KEX_AD", "KEX_BD"]),
        continued_experiments.param_ids,
    )
    values = ordinary.resolve(
        ordinary.frame_from_snapshot(continued.analysis_values.snapshot())
    )

    assert values["__KEX_AD"] == values["__KEX_BD"] == 0.0
    assert ordinary.role("__KEX_AD") is ParameterRole.FIT
    assert ordinary.role("__KEX_BD") is ParameterRole.FIT
    assert migrated.role("__KEX_AD") is ParameterRole.FIX
    assert migrated.role("__KEX_BD") is ParameterRole.FIX

"""Qualification of bounded, optional starts and authoritative final refinement."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest
from benchmarks.basin_discovery import build_context, refine

import chemex.optimize.profile_preparation as preparation
from chemex.configuration.method_plan import MethodFormatError, ProfilePreparation
from chemex.configuration.methods import read_method_plan
from chemex.optimize.direct_trf import DirectTrfConstructionError
from chemex.optimize.profile_preparation import (
    mirrored_endpoints,
    prepare_profile_start,
)
from chemex.optimize.profiled_grid import ProfiledGridPoint, ProfiledGridPointStatus
from tests.test_native_direct_trf import _qualification_fit


def _coordinates(context):
    return (
        tuple(key for key, *_ in context.coordinates if not key.startswith("__DW_")),
        tuple(key for key in context.problem.controlled_ids if key.startswith("__DW_")),
    )


def test_profile_method_roundtrip_and_legacy_de_diagnostic(tmp_path: Path) -> None:
    method = tmp_path / "method.toml"
    method.write_text(
        'FORMAT_VERSION = 2\n[STEP.SEARCH.PREPARE]\nHOLD = ["PB", "KEX_AB"]\nTRY_DW_SIGNS = ["DW_AB"]\n'
    )
    plan = read_method_plan([method])
    assert isinstance(plan.steps[0].search, ProfilePreparation)
    method.write_text(plan.render())
    assert read_method_plan([method]).steps == plan.steps
    method.write_text("FORMAT_VERSION = 2\n[STEP.SEARCH.DE]\nSEED = 1\n")
    with pytest.raises(
        MethodFormatError, match="SEARCH.DE has been retired.*SEARCH.PREPARE"
    ):
        read_method_plan([method])


@pytest.mark.parametrize(
    "settings", ("", 'HOLD = "PB"', "TRY_DW_SIGNS = [1]", "RANGES = []")
)
def test_invalid_profile_structure_is_rejected(tmp_path: Path, settings: str) -> None:
    method = tmp_path / "method.toml"
    method.write_text(f"FORMAT_VERSION = 2\n[STEP.SEARCH.PREPARE]\n{settings}\n")
    with pytest.raises(MethodFormatError):
        read_method_plan([method])


def test_endpoint_magnitudes_zero_and_connected_sign_combinations() -> None:
    assert mirrored_endpoints((("dw", 0.0),), ("dw",)) == ()
    assert mirrored_endpoints((("dw", 1.25), ("r2", 17.0)), ("dw",)) == (
        (("dw", -1.25), ("r2", 17.0)),
    )
    starts = mirrored_endpoints((("ab", 2.0), ("ac", -3.0), ("r2", 17.0)), ("ab", "ac"))
    assert len(starts) == 3
    assert {tuple(dict(start)[key] for key in ("ab", "ac")) for start in starts} == {
        (-2.0, -3.0),
        (-2.0, 3.0),
        (2.0, 3.0),
    }
    assert all(dict(start)["r2"] == 17.0 for start in starts)
    assert (
        sum(len(mirrored_endpoints(((f"dw{i}", 2.0),), (f"dw{i}",))) for i in range(20))
        == 20
    )


@pytest.mark.parametrize("zero_start", (False, True))
def test_real_cpmg_endpoint_branches_reach_lower_complete_basin(
    zero_start: bool,
) -> None:
    context = build_context("cpmg")
    hold, mirror = _coordinates(context)
    root = context.problem
    if zero_start:
        root = root.restart_from(
            tuple(
                0.0 if key in mirror else value
                for key, value in zip(root.controlled_ids, root.start, strict=True)
            )
        )
    result = prepare_profile_start(
        root, hold, mirror, context.parameterization, context.engine
    )
    assert not result.incomplete
    assert len(result.attempts) == 10  # five independent normal + mirrored fits
    assert all(
        point.objective_evaluations <= preparation.preparation_request_budget(3)
        for point in result.attempts
    )
    assert all(
        dict(point.nuisance_items)[key] != 0
        for point in result.attempts[::2]
        for key in mirror
        if key in dict(point.nuisance_items)
    )
    final = refine(context, result.vector)
    assert final["terminal"] == "accepted"
    assert final["chi2"] == pytest.approx(429.6700333, abs=1e-3)


def test_real_cest_singleton_preparation_reaches_lower_complete_basin() -> None:
    context = build_context("cest")
    hold, _ = _coordinates(context)
    result = prepare_profile_start(
        context.problem, hold, (), context.parameterization, context.engine
    )
    assert not result.incomplete
    assert len(result.attempts) == 3
    assert sum(point.chi_square for point in result.attempts) == pytest.approx(
        1516.1486504, abs=1e-3
    )
    final = refine(context, result.vector)
    assert final["terminal"] == "accepted"
    assert final["chi2"] == pytest.approx(1460.0390471, abs=1e-3)


def test_dcest_exhaustion_falls_back_and_complete_trf_still_converges() -> None:
    context = build_context("dcest")
    hold, mirror = _coordinates(context)
    with patch.object(preparation, "preparation_request_budget", return_value=2):
        result = prepare_profile_start(
            context.problem, hold, mirror, context.parameterization, context.engine
        )
    assert result.incomplete
    assert result.fallback_factors == (0,)
    assert len(result.attempts) == 1  # do not branch a failed normal fit
    assert result.attempts[0].objective_evaluations <= 2
    assert result.vector == context.problem.start
    final = refine(context, result.vector)
    assert final["terminal"] == "accepted"
    assert final["chi2"] == pytest.approx(1008.6594653, abs=1e-3)


def test_failed_alternate_retains_normal_endpoint_and_reports_incomplete() -> None:
    context = build_context("cpmg")
    hold, mirror = _coordinates(context)
    original = preparation._fit_factor_point
    calls = 0

    def fail_alternates(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls % 2 == 0:
            return ProfiledGridPoint(
                args[4], args[5], ProfiledGridPointStatus.FAILED, failure="exhausted"
            )
        return original(*args, **kwargs)

    with patch.object(preparation, "_fit_factor_point", side_effect=fail_alternates):
        result = prepare_profile_start(
            context.problem, hold, mirror, context.parameterization, context.engine
        )
    assert result.incomplete and not result.fallback_factors
    assert len(result.attempts) == 10
    for point in result.attempts[::2]:
        assert (
            dict(
                zip(context.problem.controlled_ids, result.vector, strict=True)
            ).items()
            >= dict(point.nuisance_items).items()
        )


def test_invalid_reconstruction_restores_original_root_start() -> None:
    context = build_context("cpmg")
    hold, mirror = _coordinates(context)
    with patch.object(preparation, "materialize_root_candidate", return_value=None):
        result = prepare_profile_start(
            context.problem, hold, mirror, context.parameterization, context.engine
        )
    assert result.incomplete and result.root_fallback
    assert result.vector == context.problem.start


def test_explicit_child_starts_preserve_root_context_and_reject_held_updates() -> None:
    _, _, parameterization, engine, root, _ = _qualification_fit()
    vector = tuple(value + 0.2 for value in root.start)
    restarted = root.restart_from(vector)
    child = restarted.derive_profiled_grid_point(
        factor_identity="factor",
        point_ordinal=0,
        projected_plan_identity=engine.plan.identity,
        grid_items=(),
        controlled_ids=root.controlled_ids,
        nuisance_start_items=tuple(zip(root.controlled_ids, vector, strict=True)),
    )
    assert child.start == vector
    assert child.source_snapshot is root.source_snapshot
    assert child.commit_scope == root.commit_scope
    assert not child.acceptance_authority
    restarted.validate_derived_problem(child)
    with pytest.raises(DirectTrfConstructionError, match="canonical factor FIT"):
        restarted.derive_profiled_grid_point(
            factor_identity="factor",
            point_ordinal=0,
            projected_plan_identity=engine.plan.identity,
            grid_items=(),
            controlled_ids=root.controlled_ids,
            nuisance_start_items=(("held", 3.0),),
        )
    # Captured items cannot be forged to alter a root-owned held value.
    with pytest.raises(DirectTrfConstructionError):
        forged = replace(
            child,
            independent_items=tuple(
                (key, value + 1) for key, value in child.independent_items
            ),
            derivation=replace(
                child.derivation,
                captured_independent_items=tuple(
                    (key, value + 1) for key, value in child.independent_items
                ),
            ),
        )
        restarted.validate_derived_problem(forged)


def test_out_of_bounds_mirrored_child_start_is_rejected_safely() -> None:
    context = build_context("cpmg")
    root = context.problem
    key = next(key for key in root.controlled_ids if key.startswith("__DW_"))
    with pytest.raises(DirectTrfConstructionError):
        root.derive_profiled_grid_point(
            factor_identity="factor",
            point_ordinal=0,
            projected_plan_identity=context.engine.plan.identity,
            grid_items=(),
            controlled_ids=(key,),
            nuisance_start_items=(
                (key, root.upper_bounds[root.controlled_ids.index(key)] + 1),
            ),
        )


@pytest.mark.parametrize(
    "settings,diagnostic",
    (
        ('HOLD = ["PB"]\nTRY_DW_SIGNS = ["PB"]', "overlapping"),
        ('TRY_DW_SIGNS = ["R1A_A"]', "only constant"),
        ('HOLD = ["UNKNOWN"]', "No parameter matches"),
    ),
)
def test_profile_coordinate_validation(
    tmp_path: Path, settings: str, diagnostic: str
) -> None:
    from tests.test_native_production_fitting import _programmatic_fit_context

    _, session, _ = _programmatic_fit_context(tmp_path / "Output")
    method = tmp_path / "method.toml"
    method.write_text(f"FORMAT_VERSION = 2\n[STEP.SEARCH.PREPARE]\n{settings}\n")
    model = session.parameter_factory.sealed_parameter_model
    assert model is not None
    with pytest.raises(MethodFormatError, match=diagnostic):
        read_method_plan([method]).validate(model)


def test_preparation_failure_cannot_skip_normal_final_fit_or_commit(
    tmp_path: Path, capsys
) -> None:
    from chemex.chemex import run
    from chemex.runtime import AnalysisSession
    from tests.test_native_production_fitting import _fit_arguments

    method = tmp_path / "method.toml"
    method.write_text(
        'FORMAT_VERSION = 2\n[STEP]\nROLES = [{FIX = ["KEX_AB"]}, {FIT = ["PB"]}]\n[STEP.SEARCH.PREPARE]\nHOLD = ["PB"]\n'
    )
    session = AnalysisSession.create()
    with patch.object(preparation, "preparation_request_budget", return_value=1):
        run(_fit_arguments(tmp_path / "Output", method), session=session)
    assert session.analysis_values.snapshot().revision == 1
    assert (tmp_path / "Output/Parameters/fitted.toml").exists()
    assert "continuing with the normal final fit" in capsys.readouterr().out


def test_preparation_cannot_commit_if_final_trf_fails(tmp_path: Path) -> None:
    from chemex.chemex import run
    from chemex.runtime import AnalysisSession
    from tests.test_native_production_fitting import _fit_arguments

    method = tmp_path / "method.toml"
    method.write_text(
        'FORMAT_VERSION = 2\n[STEP]\nROLES = [{FIX = ["KEX_AB"]}, {FIT = ["PB"]}]\n[STEP.SEARCH.PREPARE]\nHOLD = ["PB"]\n'
    )
    session = AnalysisSession.create()
    with (
        patch(
            "chemex.optimize.direct_trf.least_squares",
            side_effect=RuntimeError("solver failed"),
        ),
        pytest.raises(RuntimeError, match="did not commit"),
    ):
        run(_fit_arguments(tmp_path / "Output", method), session=session)
    assert session.analysis_values.snapshot().revision == 0
    assert not (tmp_path / "Output/Parameters/fitted.toml").exists()


@pytest.mark.parametrize(
    "factor_count,dw_count,expected", ((20, 1, 40), (1, 2, 4), (2, 0, 2), (1, 3, 2))
)
def test_factor_execution_is_additive_and_connected_combinations_are_local(
    factor_count: int, dw_count: int, expected: int
) -> None:
    from types import SimpleNamespace

    from chemex.optimize.profiled_grid import ProfiledGridFactor

    factors = tuple(
        ProfiledGridFactor(
            i, (i,), ("g",), tuple(f"dw{i}_{j}" for j in range(max(dw_count, 1)))
        )
        for i in range(factor_count)
    )
    ids = ("g", *(key for factor in factors for key in factor.nuisance_ids))
    root = SimpleNamespace(
        controlled_ids=ids,
        start=(100.0, *((0.0,) * (len(ids) - 1))),
        acceptance_authority=True,
        evaluation_plan_identity="plan",
        identity="root",
        validate_parameterization=lambda _: None,
    )
    engine = SimpleNamespace(
        plan=SimpleNamespace(identity="plan"), project_profiles=lambda _: None
    )
    calls = []
    budgets = []

    def endpoint(_root, factor, _engine, _par, ordinal, held, **kwargs):
        calls.append(kwargs["nuisance_start_items"])
        budgets.append(kwargs["objective_request_budget"])
        values = tuple(
            (key, 0.0 if dw_count == 0 else 2.2) for key in factor.nuisance_ids
        )
        return ProfiledGridPoint(
            ordinal, held, ProfiledGridPointStatus.SUCCESS, 1.0, values, 2
        )

    mirror = ids[1:]
    with (
        patch.object(preparation, "build_profiled_factors", return_value=factors),
        patch.object(preparation, "_fit_factor_point", side_effect=endpoint),
        patch.object(preparation, "materialize_root_candidate", return_value=None),
    ):
        result = prepare_profile_start(root, ("g",), mirror, None, engine)
    assert len(result.attempts) == expected
    assert len(calls) == (1 if dw_count > 2 else expected)
    if dw_count > 2:
        assert result.incomplete
        assert "at most two" in result.attempts[-1].failure
        assert result.attempts[0].status is ProfiledGridPointStatus.SUCCESS
    elif dw_count:
        assert budgets[0] == preparation.preparation_request_budget(max(dw_count, 1))
        assert budgets[1] == preparation.preparation_request_budget(
            max(dw_count, 1), alternate=True
        )
        assert any(any(value == -2.2 for _, value in start) for start in calls)


def test_budget_is_separate_from_full_trf_and_dimension_limited() -> None:
    from chemex.optimize.native_deterministic import _objective_request_budget

    _, _, _, _, root, _ = _qualification_fit()
    assert preparation.preparation_request_budget(4) == 3000
    assert preparation.preparation_request_budget(50) == 3000
    assert preparation.preparation_request_budget(4, alternate=True) == 500
    assert preparation.preparation_request_budget(50, alternate=True) == 500
    assert _objective_request_budget(root) == 2000 * (
        max(1, len(root.controlled_ids)) + 1
    )


def test_invalid_alternate_start_is_nonfatal_and_cannot_erase_normal_candidate() -> (
    None
):
    context = build_context("cpmg")
    hold, mirror = _coordinates(context)

    def invalid(endpoint, _mirror):
        return (
            tuple((key, 1e30 if key in mirror else value) for key, value in endpoint),
        )

    with patch.object(preparation, "mirrored_endpoints", side_effect=invalid):
        result = prepare_profile_start(
            context.problem, hold, mirror, context.parameterization, context.engine
        )
    assert result.incomplete
    assert not result.fallback_factors and not result.root_fallback
    assert all(
        point.status is ProfiledGridPointStatus.SUCCESS
        for point in result.attempts[::2]
    )
    assert all(
        point.status is ProfiledGridPointStatus.FAILED
        for point in result.attempts[1::2]
    )


def test_tc_dw_coefficients_are_not_silently_mirrored(tmp_path: Path) -> None:
    from chemex.parameters.sealed import ParamDefinition
    from tests.configuration.test_method_plan import _parameter_model

    model = _parameter_model(
        ParamDefinition("dw0", "DW0_AB", "15N", (), 0.0, -10.0, 10.0)
    )
    method = tmp_path / "method.toml"
    method.write_text(
        'FORMAT_VERSION = 2\n[STEP.SEARCH.PREPARE]\nTRY_DW_SIGNS = ["DW0_AB"]\n'
    )
    with pytest.raises(MethodFormatError, match="only constant"):
        read_method_plan([method]).validate(model)


def test_unavailable_factor_proof_restores_original_start_for_final_trf() -> None:
    from chemex.optimize.profiled_grid import ProfiledGridConstructionError

    context = build_context("cpmg")
    hold, mirror = _coordinates(context)
    with patch.object(
        preparation,
        "build_profiled_factors",
        side_effect=ProfiledGridConstructionError("cannot prove factors"),
    ):
        result = prepare_profile_start(
            context.problem, hold, mirror, context.parameterization, context.engine
        )
    assert result.incomplete and result.root_fallback
    assert result.vector == context.problem.start


@pytest.mark.parametrize("values", ((0.0, 0.0, 0.0), (1.0, -2.0, 3.0)))
def test_three_requested_mirrors_rejected_before_sign_enumeration(values) -> None:
    endpoint = tuple(zip(("ab1", "ab2", "ac"), values, strict=True))
    with (
        patch.object(
            preparation, "product", side_effect=AssertionError("unbounded enumeration")
        ),
        pytest.raises(
            DirectTrfConstructionError, match="at most two.*same coupled fit"
        ),
    ):
        mirrored_endpoints(endpoint, ("ab1", "ab2", "ac"))


def test_unqualified_factor_keeps_normal_endpoint_and_reports_reason(capsys) -> None:
    from chemex.messages import print_profile_preparation_result

    normal = ProfiledGridPoint(
        0,
        (),
        ProfiledGridPointStatus.SUCCESS,
        1.0,
        (("a", 2.0), ("b", -3.0), ("c", 4.0)),
        10,
    )
    attempts = [normal]
    assert preparation._factor_branch_starts(normal, ("a", "b", "c"), attempts) == ()
    assert attempts[0] is normal
    assert len(attempts) == 2
    failure = attempts[1].failure
    assert failure is not None and "at most two" in failure
    print_profile_preparation_result(True, 0, False, (failure,))
    assert "at most two selected DW parameters within the same coupled fit" in " ".join(
        capsys.readouterr().out.split()
    )


def test_alternate_budget_exhaustion_preserves_normal_endpoint_and_work_bound() -> None:
    context = build_context("cpmg")
    hold, mirror = _coordinates(context)
    budget = preparation.preparation_request_budget
    with patch.object(
        preparation,
        "preparation_request_budget",
        side_effect=lambda n, *, alternate=False: 1 if alternate else budget(n),
    ):
        result = prepare_profile_start(
            context.problem, hold, mirror, context.parameterization, context.engine
        )
    assert (
        result.incomplete and not result.fallback_factors and not result.root_fallback
    )
    assert len(result.attempts) == 10
    for normal, alternate in zip(
        result.attempts[::2], result.attempts[1::2], strict=True
    ):
        assert normal.status is ProfiledGridPointStatus.SUCCESS
        assert alternate.status is ProfiledGridPointStatus.FAILED
        assert alternate.objective_evaluations <= 1
        assert (
            dict(
                zip(context.problem.controlled_ids, result.vector, strict=True)
            ).items()
            >= dict(normal.nuisance_items).items()
        )


def test_real_connected_factor_over_limit_preserves_normal_and_final_trf() -> None:
    context = build_context("cpmg")
    _, mirror = _coordinates(context)
    assert len(mirror) == 5
    # Without holding coupling coordinates these five DWs share one factor.
    with patch.object(
        preparation, "product", side_effect=AssertionError("unqualified sign product")
    ):
        result = prepare_profile_start(
            context.problem, (), mirror, context.parameterization, context.engine
        )
    assert (
        result.incomplete and not result.fallback_factors and not result.root_fallback
    )
    assert len(result.attempts) == 2
    normal, rejected = result.attempts
    assert normal.status is ProfiledGridPointStatus.SUCCESS
    assert dict(normal.nuisance_items) == dict(
        zip(context.problem.controlled_ids, result.vector, strict=True)
    )
    assert rejected.status is ProfiledGridPointStatus.FAILED
    assert rejected.objective_evaluations == 0
    assert "at most two" in rejected.failure
    final = refine(context, result.vector)
    assert final["terminal"] == "accepted"
    assert final["chi2"] == pytest.approx(434.5551945, abs=1e-3)


@pytest.mark.parametrize(
    "table,settings", (("PROFILE", 'HOLD = ["PB"]'), ("PREPARE", 'MIRROR = ["DW_AB"]'))
)
def test_unshipped_preparation_spellings_have_no_compatibility_aliases(
    tmp_path: Path, table: str, settings: str
) -> None:
    method = tmp_path / "method.toml"
    method.write_text(f"FORMAT_VERSION = 2\n[STEP.SEARCH.{table}]\n{settings}\n")
    with pytest.raises(MethodFormatError, match="Unsupported v2 field"):
        read_method_plan([method])


def test_prepare_public_rendering_and_user_progress_words(
    tmp_path: Path, capsys
) -> None:
    from chemex.messages import (
        print_profile_preparation,
        print_profile_preparation_result,
    )

    method = tmp_path / "method.toml"
    method.write_text(
        'FORMAT_VERSION = 2\n[STEP.SEARCH.PREPARE]\nHOLD = ["PB", "KEX_AB"]\nTRY_DW_SIGNS = ["DW_AB"]\n'
    )
    rendered = read_method_plan([method]).render()
    assert "[STEP.SEARCH.PREPARE]" in rendered
    assert 'TRY_DW_SIGNS = ["DW_AB"]' in rendered
    assert "SEARCH.PROFILE" not in rendered and "MIRROR" not in rendered
    print_profile_preparation()
    print_profile_preparation_result(
        True, 1, False, ("TerminalFailure(category='objective_budget_exhausted')",)
    )
    output = capsys.readouterr().out
    assert "Improving starting values before fitting" in output
    assert "preliminary fit did not converge within its work limit" in output
    assert "normal final fit" in output
    assert "nuisance" not in output and "TerminalFailure" not in output

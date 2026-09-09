"""Product regressions for retained Direct-TRF Jacobians in binding fits."""

from __future__ import annotations

import math
import sys
import tomllib
from pathlib import Path
from unittest.mock import patch

import pytest

from chemex.chemex import run
from chemex.cli import build_parser
from chemex.configuration.methods import Method, Selection, read_method_plan
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import EvaluationEngine, EvaluationFrame, EvaluationResult
from chemex.experiments.builder import build_experiments
from chemex.optimize import direct_trf as direct_trf_module
from chemex.optimize import native_deterministic as native_deterministic_module
from chemex.optimize.deterministic_uncertainty import (
    AcceptedDeterministicFitFacts,
    ContinuousTrfBasis,
    DeterministicUncertainty,
    derive_deterministic_uncertainty,
)
from chemex.optimize.direct_trf import (
    AcceptedFitResult,
    DirectTrfInvocation,
    FinalResidualJacobianEvidence,
    OptimizationProblem,
    ResidualJacobianSource,
    execute_direct_trf,
)
from chemex.optimize.grouped_direct_trf import FitDecomposition
from chemex.optimize.native_mcmc import (
    McmcOperationTerminal,
    McmcPlan,
    execute_mcmc_evidence,
    resolve_product_mcmc_policy,
)
from chemex.optimize.uncertainty import ParameterUnit
from chemex.parameters.parameterization import ActiveParameterization
from chemex.parameters.sealed import parameter_name_from_definition
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
BINDING = ROOT / "examples/Combinations/2stBinding"
EXPERIMENTS = tuple(sorted((BINDING / "Experiments").glob("*.toml")))
PARAMETERS = BINDING / "Parameters/params.toml"

# The boundary/Jacobian regression needs two real binding modalities and a
# second condition, but not every shipped dataset.  These explicit profile
# identities retain the high-cost Step 1 stencil and spread Step 2 over the
# available CEST/CPMG profiles so the reduced witness keeps its binding
# constraints, rank, and boundary-warning population.
REPRESENTATIVE_EXPERIMENTS = (
    BINDING / "Experiments/cest_10hz_10p_1.toml",
    BINDING / "Experiments/cpmg_13p.toml",
)
BINDING_STEP1_PROFILES = (
    486,
    488,
    489,
    491,
    492,
    493,
    494,
    495,
    496,
    497,
    498,
    499,
)
REPRESENTATIVE_STEP2_PROFILES = (
    480,
    484,
    *BINDING_STEP1_PROFILES,
    500,
    504,
    511,
    515,
    522,
    527,
    533,
    540,
    546,
    552,
    557,
    566,
    571,
    575,
)


def _arguments(
    output: Path,
    method: Path,
    parameters: tuple[Path, ...] = (PARAMETERS,),
    *,
    experiments: tuple[Path, ...] = EXPERIMENTS,
):
    return build_parser().parse_args(
        [
            "fit",
            "-e",
            *(str(path) for path in experiments),
            "-p",
            *(str(path) for path in parameters),
            "-m",
            str(method),
            "-d",
            "2st_binding",
            "-o",
            str(output),
            "--plot",
            "nothing",
            "--workers",
            "1",
        ]
    )


def _capture_product_uncertainty(
    output: Path,
    method: Path,
    parameters: tuple[Path, ...] = (PARAMETERS,),
    *,
    experiments: tuple[Path, ...] = EXPERIMENTS,
) -> tuple[AnalysisSession, tuple[DeterministicUncertainty, ...]]:
    captured: list[DeterministicUncertainty] = []
    real_derive = native_deterministic_module.derive_deterministic_uncertainty

    def derive(facts):
        uncertainty = real_derive(facts)
        captured.append(uncertainty)
        return uncertainty

    session = AnalysisSession.create()
    with patch.object(
        native_deterministic_module,
        "derive_deterministic_uncertainty",
        derive,
    ):
        run(
            _arguments(output, method, parameters, experiments=experiments),
            session=session,
        )
    return session, tuple(captured)


def _write_representative_method(path: Path) -> Path:
    """Write the reduced two-step binding witness used by the boundary test."""

    step1_profiles = ", ".join(
        f'"{profile_id}"' for profile_id in BINDING_STEP1_PROFILES
    )
    step2_profiles = ", ".join(
        f'"{profile_id}"' for profile_id in REPRESENTATIVE_STEP2_PROFILES
    )
    path.write_text(
        f"""FORMAT_VERSION = 2

[STEP1]
INCLUDE = [{step1_profiles}]
ROLES = [
  {{ FIX = ["KD"] }},
]

[STEP2]
INCLUDE = [{step2_profiles}]
ROLES_FROM = "STEP1"
ROLES = [
  {{ FIX = ["KOFF"] }},
]
""",
        encoding="utf-8",
    )
    return path


def _write_extreme_binding_experiment(path: Path) -> Path:
    data_path = path.with_name("shift-data.txt")
    data_path.write_text(
        "486N-HN 121.34550 0.01\n487N-HN 120.00000 0.01\n",
        encoding="utf-8",
    )
    path.write_text(
        f"""[experiment]
name = "shift_15n_sq"

[conditions]
h_larmor_frq = 800.0
p_total = 1e-300
l_total = 1e-300

[data]
path = "{data_path}"
""",
        encoding="utf-8",
    )
    return path


def _qualification_anchor_with_unit_residual_jacobian(
    problem: OptimizationProblem,
    parameterization: ActiveParameterization,
    engine: EvaluationEngine,
) -> AcceptedFitResult:
    """Build exact accepted-point evidence while isolating constraint propagation."""
    vector = problem.start
    frame = EvaluationFrame.from_lifecycle_frame(
        parameterization,
        problem.lifecycle_frame(vector, parameterization),
    )
    evaluation = engine.new_evaluator().evaluate(frame)
    assert isinstance(evaluation, EvaluationResult)
    residuals = tuple(float(value) for value in evaluation.residuals)
    retained_jacobian = FinalResidualJacobianEvidence(
        ResidualJacobianSource.SCIPY_FINAL_2_POINT,
        problem.controlled_ids,
        vector,
        residuals,
        (len(residuals), len(problem.controlled_ids)),
        tuple(
            tuple(
                1.0 if row_index == column_index else 0.0
                for column_index in range(len(problem.controlled_ids))
            )
            for row_index, _residual in enumerate(residuals)
        ),
    )
    occurrence_identity = "exact-binding-boundary-qualification-occurrence"
    return AcceptedFitResult(
        occurrence_identity,
        problem.identity,
        "exact-binding-boundary-qualification-invocation",
        "exact-binding-boundary-qualification-execution",
        "exact-binding-boundary-qualification-materialization",
        problem.parameterization_identity,
        problem.evaluator_parameterization_identity,
        problem.source_snapshot.occurrence_identity,
        problem.source_snapshot.revision,
        problem.controlled_ids,
        vector,
        math.fsum(value * value for value in residuals),
        evaluation,
        problem.controlled_ids,
        tuple(zip(problem.controlled_ids, vector, strict=True)),
        "exact-binding-boundary-qualification-origin",
        final_residual_jacobian=retained_jacobian,
        occurrence_witness=direct_trf_module._mint_accepted_occurrence_witness(
            occurrence_identity
        ),
    )


def _assert_minimum_category_c_uncertainty_is_qualified_or_fail_closed(
    uncertainty: DeterministicUncertainty,
    output_ids: dict[str, str],
    fitted_names: tuple[str, ...],
    population_ids: dict[str, str],
) -> None:
    conclusions = {
        name: uncertainty.parameter(param_id) for name, param_id in output_ids.items()
    }
    assert all(
        conclusions[name] is not None
        and conclusions[name].reportable
        and conclusions[name].standard_error is not None
        and math.isfinite(conclusions[name].standard_error)
        and conclusions[name].standard_error > 0.0
        for name in fitted_names
    )
    for name in population_ids:
        conclusion = conclusions[name]
        assert conclusion is not None
        if conclusion.reportable:
            assert conclusion.standard_error is not None
            assert math.isfinite(conclusion.standard_error)
            assert conclusion.standard_error > 0.0
        else:
            assert conclusion.standard_error is None
            assert conclusion.unavailable_kind is not None
    assert uncertainty.root_evidence is not None
    constraint_jacobian = uncertainty.root_evidence.constraint_jacobian
    if constraint_jacobian is not None:
        assert all(
            math.isfinite(coefficient)
            for gradient in constraint_jacobian.gradient_by_target.values()
            for coefficient in gradient
        )


@pytest.mark.parametrize(
    (
        "model_name",
        "fitted_names",
        "output_names",
        "parameter_values",
        "kd_value",
    ),
    (
        (
            "2st_binding",
            ("KD",),
            ("PA", "PB"),
            "KOFF = 0.0\nKD = 5e-324",
            math.nextafter(0.0, 1.0),
        ),
        (
            "3st_double_binding",
            ("KD_AB",),
            ("PA", "PB", "PC"),
            "KOFF_AB = 0.0\nKOFF_AC = 0.0\nKD_AB = 5e-324\nKD_AC = 1e-6",
            math.nextafter(0.0, 1.0),
        ),
        (
            "3st_binding_cs",
            ("KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB"),
            (
                "KD_APP",
                "KOFF_BC",
                "KEQ_AB",
                "KEX_AB",
                "KAB",
                "KBA",
                "KD_BC",
                "KBC",
                "KCB",
                "PA",
                "PB",
                "PC",
            ),
            ("KOFF_BC = 0.0\nKEX_AB = 0.0\nKEQ_AB = 1.0\nKD_APP = 1e-6"),
            1.0e-6,
        ),
        (
            "3st_binding_cs",
            ("KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB"),
            (
                "KD_APP",
                "KOFF_BC",
                "KEQ_AB",
                "KEX_AB",
                "KAB",
                "KBA",
                "KD_BC",
                "KBC",
                "KCB",
                "PA",
                "PB",
                "PC",
            ),
            ("KOFF_BC = 100.0\nKEX_AB = 1e-10\nKEQ_AB = 1.0\nKD_APP = 1e-6"),
            1.0e-6,
        ),
        (
            "3st_binding_cs",
            ("KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB"),
            (
                "KD_APP",
                "KOFF_BC",
                "KEQ_AB",
                "KEX_AB",
                "KAB",
                "KBA",
                "KD_BC",
                "KBC",
                "KCB",
                "PA",
                "PB",
                "PC",
            ),
            ("KOFF_BC = 80.0\nKEX_AB = 250.0\nKEQ_AB = 3.0\nKD_APP = 2e-4"),
            2.0e-4,
        ),
        (
            "3st_binding_if",
            ("KD_APP", "KOFF_AB", "KEQ_BC", "KEX_BC"),
            (
                "KD_APP",
                "KOFF_AB",
                "KEQ_BC",
                "KEX_BC",
                "KBC",
                "KCB",
                "KD_AB",
                "KAB",
                "KBA",
                "PA",
                "PB",
                "PC",
            ),
            ("KOFF_AB = 0.0\nKEX_BC = 200.0\nKEQ_BC = 0.0\nKD_APP = 1e-6"),
            1.0e-6,
        ),
        (
            "3st_binding_if",
            ("KD_APP", "KOFF_AB", "KEQ_BC", "KEX_BC"),
            (
                "KD_APP",
                "KOFF_AB",
                "KEQ_BC",
                "KEX_BC",
                "KBC",
                "KCB",
                "KD_AB",
                "KAB",
                "KBA",
                "PA",
                "PB",
                "PC",
            ),
            ("KOFF_AB = 80.0\nKEX_BC = 250.0\nKEQ_BC = 3.0\nKD_APP = 2e-4"),
            2.0e-4,
        ),
        (
            "3st_binding_if",
            ("KD_APP", "KOFF_AB", "KEQ_BC", "KEX_BC"),
            (
                "KD_APP",
                "KOFF_AB",
                "KEQ_BC",
                "KEX_BC",
                "KBC",
                "KCB",
                "KD_AB",
                "KAB",
                "KBA",
                "PA",
                "PB",
                "PC",
            ),
            ("KOFF_AB = 100.0\nKEX_BC = 0.0\nKEQ_BC = 1.0\nKD_APP = 1e-6"),
            1.0e-6,
        ),
        (
            "3st_binding_cs",
            ("KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB"),
            ("KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB", "PA", "PB", "PC"),
            ("KOFF_BC = 0.0\nKEX_AB = 0.0\nKEQ_AB = 1.0\nKD_APP = 5e-324"),
            math.nextafter(0.0, 1.0),
        ),
        (
            "3st_binding_cs",
            ("KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB"),
            ("KD_APP", "KOFF_BC", "KEQ_AB", "KEX_AB", "PA", "PB", "PC"),
            ("KOFF_BC = 0.0\nKEX_AB = 0.0\nKEQ_AB = 5e-324\nKD_APP = 1e-3"),
            1.0e-3,
        ),
        (
            "4st_binding_3_bound_states",
            ("KD_APP", "KEQ_BC"),
            ("PA", "PB", "PC", "PD"),
            (
                "KOFF_AB = 0.0\nKEX_BC = 0.0\nKEX_CD = 0.0\n"
                "KD_APP = 2.2250738585072014e-308\nKEQ_BC = 0.0\nKEQ_CD = 1.0"
            ),
            sys.float_info.min,
        ),
    ),
    ids=(
        "simple-binding",
        "three-state-binding",
        "conformational-selection-zero-kex",
        "conformational-selection-near-zero-kex",
        "conformational-selection-interior",
        "induced-fit-zero-keq",
        "induced-fit-interior",
        "induced-fit-zero-kex",
        "conformational-selection-minimum-kd-app",
        "conformational-selection-minimum-keq",
        "four-state-binding",
    ),
)
def test_boundary_association_uncertainties_are_qualified(
    tmp_path: Path,
    model_name: str,
    fitted_names: tuple[str, ...],
    output_names: tuple[str, ...],
    parameter_values: str,
    kd_value: float,
) -> None:
    experiment = _write_extreme_binding_experiment(tmp_path / "experiment.toml")
    experiment.write_text(
        experiment.read_text(encoding="utf-8")
        .replace("p_total = 1e-300", "p_total = 1e-3")
        .replace("l_total = 1e-300", "l_total = 2e-3"),
        encoding="utf-8",
    )
    if model_name in {"3st_binding_cs", "3st_binding_if"}:
        data_path = experiment.with_name("shift-data.txt")
        data_path.write_text(
            data_path.read_text(encoding="utf-8")
            + "488N-HN 119.50000 0.01\n489N-HN 118.75000 0.01\n",
            encoding="utf-8",
        )
    parameters = tmp_path / "parameters.toml"
    parameters.write_text(
        f"[GLOBAL]\nR1_A = 1.5\nR2_A = 4.7\n{parameter_values}\n",
        encoding="utf-8",
    )
    session = AnalysisSession.create()
    session.set_model(model_name)
    experiments = build_experiments(
        [experiment],
        Selection(include="*", exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([parameters]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    parameter_model = session.parameter_factory.sealed_parameter_model
    configuration = session.parameter_factory.sealed_configuration
    assert parameter_model is not None
    assert configuration is not None
    output_ids = {
        definition.name: definition.param_id
        for definition in parameter_model.definitions
        if definition.name in output_names
    }
    assert set(output_ids) == set(output_names)
    population_ids = {
        name: param_id
        for name, param_id in output_ids.items()
        if name in {"PA", "PB", "PC", "PD"}
    }
    required_ids = experiments.param_ids | set(output_ids.values())
    baseline = session.compile_parameterization(Method(), required_ids)
    fitted_ids = {
        definition.param_id
        for definition in parameter_model.definitions
        if definition.name in fitted_names
    }
    assert len(fitted_ids) == len(fitted_names)
    fixed_names = sorted(
        {
            parameter_name_from_definition(parameter_model.definitions[param_id]).name
            for param_id in baseline.independent_ids
            if param_id not in fitted_ids
        }
    )
    parameterization = session.compile_parameterization(
        Method(fit=fitted_names, fix=fixed_names),
        required_ids,
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    start_by_name = {
        parameter_model.definitions[param_id].name: value
        for param_id, value in zip(
            problem.controlled_ids,
            problem.start,
            strict=True,
        )
    }
    assert start_by_name[fitted_names[0]] == kd_value
    if "KEQ_AB = 5e-324" in parameter_values:
        assert start_by_name["KEQ_AB"] == math.nextafter(0.0, 1.0)
    accepted = _qualification_anchor_with_unit_residual_jacobian(
        problem,
        parameterization,
        engine,
    )
    mcmc_kex_name = (
        "KEX_AB"
        if model_name == "3st_binding_cs"
        and "KEQ_AB = 1.0\nKD_APP = 1e-6" in parameter_values
        else "KEX_BC"
        if model_name == "3st_binding_if" and "KEX_BC = 0.0" in parameter_values
        else None
    )
    if mcmc_kex_name is not None:
        mcmc_accepted = AcceptedFitResult.for_qualification(
            occurrence_identity="category-c-boundary-mcmc-occurrence",
            problem_identity=problem.identity,
            invocation_identity="category-c-boundary-mcmc-invocation",
            execution_identity="category-c-boundary-mcmc-execution",
            materialization_identity="category-c-boundary-mcmc-materialization",
            parameterization_identity=parameterization.identity,
            evaluator_parameterization_identity=parameterization.evaluator_identity,
            source_occurrence_identity=problem.source_snapshot.occurrence_identity,
            source_revision=problem.source_snapshot.revision,
            controlled_ids=problem.controlled_ids,
            vector=problem.start,
            chi_square=accepted.chi_square,
            evaluation_result=accepted.evaluation_result,
            commit_scope=problem.commit_scope,
            commit_items=accepted.evaluation_result.resolved_values.ordered_items(),
            origin_context_identity="category-c-boundary-mcmc-origin",
        )
        mcmc_plan = McmcPlan.for_accepted(
            mcmc_accepted,
            source_problem=problem,
            parameterization=parameterization,
            source_engine=engine,
            policy=resolve_product_mcmc_policy(
                dimension=len(problem.controlled_ids),
                walkers=10,
                steps=2,
                root_seed=733,
            ),
            coordinate_units=tuple(
                (param_id, ParameterUnit.UNSPECIFIED)
                for param_id in problem.controlled_ids
            ),
        )
        [kex_index] = [
            index
            for index, param_id in enumerate(problem.controlled_ids)
            if parameter_model.definitions[param_id].name == mcmc_kex_name
        ]
        expected_kex = 1.0e-10 if "KEX_AB = 1e-10" in parameter_values else 0.0
        assert problem.start[kex_index] == expected_kex
        if expected_kex == 0.0:
            assert any(row[kex_index] == 0.0 for row in mcmc_plan.initial_ensemble)
        mcmc_operation = execute_mcmc_evidence(mcmc_accepted, mcmc_plan)
        assert mcmc_operation.terminal is McmcOperationTerminal.COMPLETED, (
            mcmc_operation.failure_message
        )
        assert mcmc_operation.evidence is not None
        initial_state = mcmc_operation.evidence.states[0]
        assert all(math.isfinite(value) for value in initial_state.log_densities)
        if expected_kex == 0.0:
            assert any(
                position[kex_index] == 0.0 for position in initial_state.positions
            )
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    uncertainty = derive_deterministic_uncertainty(
        AcceptedDeterministicFitFacts(
            accepted,
            problem,
            parameterization,
            engine,
            ContinuousTrfBasis(decomposition.partition_proof),
            "exact-binding-boundary-output-qualification",
        )
    )

    conclusions = {
        name: uncertainty.parameter(param_id) for name, param_id in output_ids.items()
    }
    minimum_category_c_boundary = (
        kd_value == math.nextafter(0.0, 1.0) or "KEQ_AB = 5e-324" in parameter_values
    ) and model_name == "3st_binding_cs"
    if minimum_category_c_boundary:
        _assert_minimum_category_c_uncertainty_is_qualified_or_fail_closed(
            uncertainty,
            output_ids,
            fitted_names,
            population_ids,
        )
        return
    assert all(
        conclusion is not None
        and conclusion.reportable
        and conclusion.standard_error is not None
        and math.isfinite(conclusion.standard_error)
        for conclusion in conclusions.values()
    ), (
        {
            name: (
                None
                if conclusion is None
                else (
                    conclusion.reportable,
                    conclusion.standard_error,
                    conclusion.unavailable_kind,
                )
            )
            for name, conclusion in conclusions.items()
        },
        ()
        if uncertainty.root_evidence is None
        else tuple(
            (failure.category, failure.message)
            for failure in uncertainty.root_evidence.failures
        ),
        None
        if uncertainty.root_evidence is None
        else (
            uncertainty.root_evidence.constraint_jacobian is not None,
            uncertainty.root_evidence.constrained_propagation is not None,
            uncertainty.root_evidence.constrained_marginal_errors is not None,
            uncertainty.root_evidence.constrained_correlations is not None,
        ),
        ()
        if uncertainty.root_evidence is None
        or uncertainty.root_evidence.constrained_marginal_errors is None
        else tuple(uncertainty.root_evidence.constrained_marginal_errors.entries),
    )
    assert uncertainty.root_evidence is not None
    constraint_jacobian = uncertainty.root_evidence.constraint_jacobian
    assert constraint_jacobian is not None
    gradients = {
        parameter_model.definitions[param_id].name: row
        for param_id, row in zip(
            constraint_jacobian.output_ids,
            constraint_jacobian.matrix,
            strict=True,
        )
        if param_id in population_ids.values()
    }
    kd_column = next(
        index
        for index, param_id in enumerate(problem.controlled_ids)
        if parameter_model.definitions[param_id].name == fitted_names[0]
    )
    if model_name in {"3st_binding_cs", "3st_binding_if"}:
        assert math.isfinite(gradients["PA"][kd_column])
        assert gradients["PA"][kd_column] > 0.0
    else:
        assert gradients["PA"][kd_column] == pytest.approx(
            1.0e3,
            rel=2.0e-7,
            abs=0.0,
        )
    assert gradients["PB"][kd_column] != 0.0
    for column in range(len(problem.controlled_ids)):
        column_gradients = tuple(gradient[column] for gradient in gradients.values())
        assert math.fsum(column_gradients) == pytest.approx(
            0.0,
            abs=4.0 * math.ulp(max(map(abs, column_gradients))),
        )
    complement_components = (
        {"pa", "pb", "pc"}
        if model_name in {"3st_binding_cs", "3st_binding_if"}
        else {"pb"}
    )
    assert any(
        diagnostic.component in complement_components
        and diagnostic.method == "normalized_population_complement"
        for diagnostic in constraint_jacobian.function_partial_diagnostics
    )


def test_unrepresentable_kon_does_not_block_simulation_or_fixed_fit_setup(
    tmp_path: Path,
) -> None:
    experiment = _write_extreme_binding_experiment(tmp_path / "experiment.toml")
    parameters = tmp_path / "parameters.toml"
    parameters.write_text(
        PARAMETERS.read_text(encoding="utf-8")
        .replace("KOFF = 50.0", "KOFF = 1.0")
        .replace("KD = 2.3E-6", "KD = 5e-324\nKON = 1.0"),
        encoding="utf-8",
    )
    method = tmp_path / "method.toml"
    method.write_text(
        """FORMAT_VERSION = 2

[STEP]
INCLUDE = [486]
ROLES = [
  { FIT = ["CS_A"] },
  { FIX = ["KD", "KOFF"] },
]
""",
        encoding="utf-8",
    )

    output = tmp_path / "Simulation"
    run(
        build_parser().parse_args(
            [
                "simulate",
                "-e",
                str(experiment),
                "-p",
                str(parameters),
                "-d",
                "2st_binding",
                "-o",
                str(output),
                "--plot",
                "nothing",
            ]
        ),
        session=AnalysisSession.create(),
    )
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    assert "KAB" in constrained
    assert "KON" not in constrained

    fit_session = AnalysisSession.create()
    fit_session.set_model("2st_binding")
    experiments = build_experiments(
        [experiment],
        Selection(include=[SpinSystem.from_name("486N-HN")], exclude=None),
        session=fit_session,
    )
    fit_session.parameters.set_defaults(read_defaults([parameters]))
    assert fit_session.try_build_analysis_values(), repr(
        fit_session.parameter_factory.native_construction_error
    )
    plan = read_method_plan([method])
    parameterization = fit_session.compile_parameterization_from_actions(
        plan.effective_role_actions()["STEP"],
        experiments.param_ids,
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = fit_session.parameter_factory.sealed_configuration
    assert configuration is not None
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        fit_session.analysis_values.snapshot(),
    )
    assert any(param_id.startswith("__KAB") for param_id in problem.commit_scope)
    assert all(not param_id.startswith("__KON") for param_id in problem.commit_scope)


@pytest.mark.parametrize(
    ("model_name", "fitted_name", "population_names"),
    (
        ("2st_monomer_dimer", "KD", ("PA", "PB")),
        ("3st_monomer_dimer_trimer", "KD1", ("PA", "PB", "PC")),
        ("3st_double_binding", "KD_AB", ("PA", "PB", "PC")),
        ("3st_binding_partner_2st", "KD_AB", ("PA", "PB", "PC")),
        (
            "4st_binding_3_bound_states",
            "KD_APP",
            ("PA", "PB", "PC", "PD", "KD_EFF"),
        ),
    ),
)
def test_association_population_uncertainties_are_reportable_in_structured_output(
    tmp_path: Path,
    model_name: str,
    fitted_name: str,
    population_names: tuple[str, ...],
) -> None:
    experiment = _write_extreme_binding_experiment(tmp_path / "experiment.toml")
    experiment.write_text(
        experiment.read_text(encoding="utf-8")
        .replace("p_total = 1e-300", "p_total = 1e-3")
        .replace("l_total = 1e-300", "l_total = 2e-3"),
        encoding="utf-8",
    )
    session = AnalysisSession.create()
    session.set_model(model_name)
    experiments = build_experiments(
        [experiment],
        Selection(include="*", exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    parameter_model = session.parameter_factory.sealed_parameter_model
    configuration = session.parameter_factory.sealed_configuration
    assert parameter_model is not None
    assert configuration is not None
    requested_ids = {
        definition.param_id
        for definition in parameter_model.definitions
        if definition.name in population_names
    }
    required_ids = experiments.param_ids | requested_ids
    baseline = session.compile_parameterization(Method(), required_ids)
    fitted_ids = {
        definition.param_id
        for definition in parameter_model.definitions
        if definition.name == fitted_name
    }
    assert len(fitted_ids) == 1
    fixed_names = sorted(
        {
            parameter_name_from_definition(parameter_model.definitions[param_id]).name
            for param_id in baseline.independent_ids
            if param_id not in fitted_ids
        }
    )
    parameterization = session.compile_parameterization(
        Method(fit=[fitted_name], fix=fixed_names),
        required_ids,
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )
    outcome = execute_direct_trf(
        problem,
        DirectTrfInvocation.for_problem(problem, objective_request_budget=100),
        parameterization,
        engine,
    )
    accepted = outcome.accepted_result
    assert accepted is not None
    decomposition = FitDecomposition.from_root(problem, parameterization, engine)
    uncertainty = derive_deterministic_uncertainty(
        AcceptedDeterministicFitFacts(
            accepted,
            problem,
            parameterization,
            engine,
            ContinuousTrfBasis(decomposition.partition_proof),
            "association-population-uncertainty-test",
        )
    )
    output_ids = tuple(
        param_id for param_id in parameterization.scope_ids if param_id in requested_ids
    )

    assert {
        parameter_model.definitions[param_id].name for param_id in output_ids
    } == set(population_names)
    assert all(
        (conclusion := uncertainty.parameter(param_id)) is not None
        and conclusion.reportable
        and conclusion.standard_error is not None
        and math.isfinite(conclusion.standard_error)
        for param_id in output_ids
    ), {
        parameter_model.definitions[param_id].name: uncertainty.parameter(param_id)
        for param_id in output_ids
    }


def _assert_representative_profiles_are_shipped() -> None:
    """Fail clearly if a selected residue's stable N profile disappears."""

    required_profiles = {
        f"{profile_id}N" for profile_id in REPRESENTATIVE_STEP2_PROFILES
    }
    for experiment in REPRESENTATIVE_EXPERIMENTS:
        with experiment.open("rb") as file:
            profiles = tomllib.load(file)["data"]["profiles"]
        missing = sorted(required_profiles.difference(profiles))
        assert not missing, (
            f"{experiment.name} is missing representative profiles: {missing}"
        )


def test_fitted_kd_uses_backend_jacobian_before_boundary_reporting(
    tmp_path: Path,
) -> None:
    method = tmp_path / "kd-method.toml"
    method.write_text(
        """[STEP]
INCLUDE = [486, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499]
FIT = ["KD"]
FIX = ["KOFF", "DW_AB", "R1_A", "R2_A", "CS_A"]
""",
        encoding="utf-8",
    )
    step1_method = tmp_path / "step1-method.toml"
    step1_method.write_text(
        """[STEP1]
INCLUDE = [486, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499]
FIX = ["KD"]
""",
        encoding="utf-8",
    )
    step1_output = tmp_path / "Step1"
    step1_session = AnalysisSession.create()
    run(_arguments(step1_output, step1_method), session=step1_session)
    assert step1_session.analysis_values.snapshot().revision == 1
    output = tmp_path / "Output"

    session, (uncertainty,) = _capture_product_uncertainty(
        output,
        method,
        (PARAMETERS, step1_output / "Parameters" / "fitted.toml"),
    )

    assert session.analysis_values.snapshot().revision == 1
    evidence = uncertainty.root_evidence
    assert evidence is not None
    jacobian = evidence.residual_jacobian
    assert jacobian is not None
    assert jacobian.method == "retained-scipy-final-2-point"
    assert jacobian.evaluation_count == 0
    assert evidence.rank_diagnostic is not None
    assert evidence.covariance is not None
    kd_index = jacobian.controlled_ids.index("__KD")
    kd = evidence.accepted_anchor.vector[kd_index]
    assert kd == pytest.approx(1.0156e-6, rel=2.0e-3)
    assert math.isfinite(evidence.covariance.covariance[kd_index][kd_index])
    fitted = (output / "Parameters" / "fitted.toml").read_text(encoding="utf-8")
    assert "KD" in fitted
    assert "Jacobian unavailable" not in fitted
    parameter_model = session.parameter_factory.sealed_parameter_model
    assert parameter_model is not None
    population_ids = tuple(
        definition.param_id
        for definition in parameter_model.definitions
        if definition.name in {"PA", "PB"}
        and uncertainty.parameter(definition.param_id) is not None
    )
    assert population_ids
    assert all(
        uncertainty.parameter(param_id).reportable for param_id in population_ids
    )
    constrained = (output / "Parameters" / "constrained.toml").read_text(
        encoding="utf-8"
    )
    population_lines = tuple(
        line for line in constrained.splitlines() if line.startswith(('"PA,', '"PB,'))
    )
    assert len(population_lines) == len(population_ids)
    assert all("# ±" in line for line in population_lines)


def test_representative_binding_step2_reports_boundary_warning_with_retained_jacobian(
    tmp_path: Path,
) -> None:
    output = tmp_path / "Output"
    method = _write_representative_method(tmp_path / "representative-method.toml")
    _assert_representative_profiles_are_shipped()

    session, uncertainty = _capture_product_uncertainty(
        output,
        method,
        experiments=REPRESENTATIVE_EXPERIMENTS,
    )

    assert session.analysis_values.snapshot().revision == 2
    assert len(uncertainty) == 2
    step2 = uncertainty[1].root_evidence
    assert step2 is not None
    assert step2.residual_jacobian is not None
    assert step2.residual_jacobian.method == "retained-scipy-final-2-point"
    assert step2.residual_jacobian.evaluation_count == 0
    assert step2.rank_diagnostic is not None
    assert all(failure.stage != "residual_linearization" for failure in step2.failures)
    fitted = (output / "STEP2" / "Parameters" / "fitted.toml").read_text(
        encoding="utf-8"
    )
    fitted_uncertainties = [line for line in fitted.splitlines() if "# ±" in line]
    warned_uncertainties = [
        line
        for line in fitted_uncertainties
        if "boundary may make uncertainty asymmetric" in line
    ]
    unavailable_uncertainties = [
        line for line in fitted.splitlines() if "error unavailable" in line
    ]
    assert step2.covariance is not None
    assert len(warned_uncertainties) == len(step2.covariance.simple_bound_warning_ids)
    assert 0 < len(warned_uncertainties) < len(fitted_uncertainties)
    assert unavailable_uncertainties == []
    assert "error unavailable: boundary limited" not in fitted
    assert "Jacobian unavailable" not in fitted

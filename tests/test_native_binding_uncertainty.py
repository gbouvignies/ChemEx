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


@pytest.mark.parametrize(
    (
        "model_name",
        "fitted_names",
        "population_names",
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
    ids=("simple-binding", "three-state-binding", "four-state-binding"),
)
def test_boundary_association_population_uncertainties_are_reportable(
    tmp_path: Path,
    model_name: str,
    fitted_names: tuple[str, ...],
    population_names: tuple[str, ...],
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
    population_ids = {
        definition.name: definition.param_id
        for definition in parameter_model.definitions
        if definition.name in population_names
    }
    assert set(population_ids) == set(population_names)
    required_ids = experiments.param_ids | set(population_ids.values())
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
    accepted = _qualification_anchor_with_unit_residual_jacobian(
        problem,
        parameterization,
        engine,
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
        name: uncertainty.parameter(param_id)
        for name, param_id in population_ids.items()
    }
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
    assert any(
        diagnostic.component == "pb"
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
    )


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

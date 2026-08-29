from __future__ import annotations

import dataclasses
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemex.configuration.methods import Selection, read_method_plan
from chemex.configuration.parameters import read_defaults
from chemex.evaluation.native import EvaluationEngine
from chemex.experiments.builder import build_experiments
from chemex.optimize.direct_trf import OptimizationProblem
from chemex.optimize.helper import execute_simulation
from chemex.parameters.feasible_coordinates import (
    FeasibleCoordinateConstructionError,
    ScientificFeasibilityError,
    _conditional_interval,
    _reference_coefficient,
    _validate_blocks,
    compile_feasible_coordinates,
    relaxation_state_is_on_boundary,
)
from chemex.parameters.parameterization import (
    BinaryExpression,
    CompiledConstraint,
    FunctionExpression,
    IndependentValueFrame,
    LiteralExpression,
    ReferenceExpression,
    UnaryExpression,
)
from chemex.parameters.relaxation import (
    RelaxationPsdBlock,
    SealedRelaxationDomains,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession


def test_two_state_shared_cross_rate_uses_the_limiting_psd_block() -> None:
    blocks = (
        RelaxationPsdBlock("transverse-a", "a", ("r2a", "r2aa"), ((0, 1, "eta"),)),
        RelaxationPsdBlock("transverse-b", "b", ("r2b", "r2ab"), ((0, 1, "eta"),)),
    )
    values = {"r2a": 34.0, "r2aa": 33.5, "r2b": 31.0, "r2ab": 30.5, "eta": 0.0}
    limits = [_conditional_interval(block, "eta", values)[1] for block in blocks]
    limit = min(limits)

    values["eta"] = limit
    _validate_blocks(blocks, values)
    values["eta"] = limit * (1.0 + 1.0e-9)
    with pytest.raises(ScientificFeasibilityError, match="positive semidefinite"):
        _validate_blocks(blocks, values)


def test_longitudinal_cross_interval_enforces_the_complete_three_mode_block() -> None:
    block = RelaxationPsdBlock(
        "longitudinal-a",
        "a",
        ("r1i", "r1s", "r1a"),
        ((0, 1, "sigma"), (0, 2, "etai"), (1, 2, "etas")),
    )
    values = {
        "r1i": 3.0,
        "r1s": 2.0,
        "r1a": 4.0,
        "sigma": 1.0,
        "etai": 0.0,
        "etas": 1.0,
    }
    lower, upper = _conditional_interval(block, "etai", values)

    values["etai"] = upper
    _validate_blocks((block,), values)
    values["etai"] = upper + 1.0e-6
    with pytest.raises(ScientificFeasibilityError, match="longitudinal-a"):
        _validate_blocks((block,), values)
    assert lower < upper < np.sqrt(values["r1i"] * values["r1a"])


@pytest.mark.parametrize("mu", [-12.0, 12.0])
def test_mu_boundary_is_absolute_r2mq(mu: float) -> None:
    block = RelaxationPsdBlock(
        "multiple-quantum-a",
        "a",
        ("r2mq", "r2mq"),
        ((0, 1, "mu"),),
    )
    values = {"r2mq": 12.0, "mu": mu}
    _validate_blocks((block,), values)
    values["mu"] = np.copysign(12.000001, mu)
    with pytest.raises(ScientificFeasibilityError, match="multiple-quantum-a"):
        _validate_blocks((block,), values)


def _affine_parameterization():
    block = RelaxationPsdBlock(
        "transverse-a",
        "a",
        ("r2", "r2a"),
        ((0, 1, "eta"),),
    )
    constraint = CompiledConstraint(
        "r2a",
        BinaryExpression(
            "subtract",
            ReferenceExpression("r2"),
            ReferenceExpression("r1"),
        ),
        ("r2", "r1"),
        "r2 - r1",
        "r2 - r1",
    )

    class Parameterization:
        scope_ids = ("r2", "r1", "eta", "r2a")
        evaluator_identity = "relaxation-test-evaluator"
        program = SimpleNamespace(
            constraints=(constraint,),
            relaxation_domains=SealedRelaxationDomains((block,)),
        )

        @staticmethod
        def resolve(frame):
            values = dict(frame.ordered_items())
            values.setdefault("r2", 0.0)
            values.setdefault("r1", 0.0)
            values["r2a"] = values["r2"] - values["r1"]
            return values

    return Parameterization()


def _dq_affine_parameterization():
    block = RelaxationPsdBlock(
        "transverse-a",
        "a",
        ("r2dq", "r2adq"),
        ((0, 1, "eta"),),
    )
    constraint = CompiledConstraint(
        "r2adq",
        BinaryExpression(
            "subtract",
            ReferenceExpression("r2dq"),
            BinaryExpression(
                "divide",
                ReferenceExpression("r1"),
                LiteralExpression(3.0),
            ),
        ),
        ("r2dq", "r1"),
        "method",
        "[R2DQ] - [R1] / 3.0",
    )

    class Parameterization:
        scope_ids = ("r2dq", "r1", "eta", "r2adq")
        evaluator_identity = "dq-relaxation-test-evaluator"
        program = SimpleNamespace(
            constraints=(constraint,),
            relaxation_domains=SealedRelaxationDomains((block,)),
        )

        @staticmethod
        def resolve(frame):
            values = dict(frame.ordered_items())
            values.setdefault("r2dq", 0.0)
            values.setdefault("r1", 0.0)
            values["r2adq"] = values["r2dq"] - values["r1"] / 3.0
            return values

    return Parameterization()


def test_compiler_recognizes_exact_dq_affine_nonnegative_floor() -> None:
    parameterization = _dq_affine_parameterization()
    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("r2dq", 20.0), ("r1", 1.5), ("eta", 0.0)),
    )
    maximum = float(np.finfo(np.float64).max)

    chart = compile_feasible_coordinates(
        parameterization,
        frame,
        ("r2dq",),
        (0.0,),
        (maximum,),
    )

    assert chart is not None
    assert chart.rate_excess_ids == (("r2dq", "r2adq"),)
    boundary = chart.decode((0.0,))
    values = parameterization.resolve(boundary.frame)
    assert values["r2dq"] == pytest.approx(0.5)
    assert values["r2adq"] == pytest.approx(0.0)
    _validate_blocks(chart.blocks, values)


def test_compiler_recognizes_exact_dq_affine_psd_floor() -> None:
    parameterization = _dq_affine_parameterization()
    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("r2dq", 20.0), ("r1", 1.5), ("eta", 3.0)),
    )
    maximum = float(np.finfo(np.float64).max)

    chart = compile_feasible_coordinates(
        parameterization,
        frame,
        ("r2dq",),
        (0.0,),
        (maximum,),
    )

    assert chart is not None
    assert chart.rate_excess_ids == (("r2dq", "r2adq"),)
    assert chart.static_rate_floors[0][0] == "r2dq"
    boundary = chart.decode((0.0,))
    values = parameterization.resolve(boundary.frame)
    assert values["r2dq"] == pytest.approx(0.25 + 0.5 * np.sqrt(36.25))
    assert values["r2dq"] * values["r2adq"] == pytest.approx(9.0)
    _validate_blocks(chart.blocks, values)


@pytest.mark.parametrize(
    ("expression", "expected"),
    [
        (
            BinaryExpression("add", ReferenceExpression("x"), LiteralExpression(4.0)),
            1.0,
        ),
        (
            BinaryExpression(
                "subtract", ReferenceExpression("x"), LiteralExpression(4.0)
            ),
            1.0,
        ),
        (
            BinaryExpression("add", LiteralExpression(4.0), ReferenceExpression("x")),
            1.0,
        ),
        (
            BinaryExpression(
                "multiply", LiteralExpression(2.5), ReferenceExpression("x")
            ),
            2.5,
        ),
        (
            BinaryExpression(
                "multiply", ReferenceExpression("x"), LiteralExpression(2.5)
            ),
            2.5,
        ),
        (
            BinaryExpression(
                "divide", ReferenceExpression("x"), LiteralExpression(4.0)
            ),
            0.25,
        ),
        (UnaryExpression("positive", ReferenceExpression("x")), 1.0),
        (UnaryExpression("negative", ReferenceExpression("x")), -1.0),
        (
            BinaryExpression(
                "subtract",
                BinaryExpression(
                    "multiply",
                    LiteralExpression(3.0),
                    BinaryExpression(
                        "add", ReferenceExpression("x"), LiteralExpression(2.0)
                    ),
                ),
                UnaryExpression("negative", ReferenceExpression("x")),
            ),
            4.0,
        ),
    ],
)
def test_reference_coefficient_supports_exact_affine_forms(
    expression: object,
    expected: float,
) -> None:
    assert _reference_coefficient(expression, "x") == expected


def test_reference_coefficient_expands_recursive_derived_references() -> None:
    constraints = {
        "y": CompiledConstraint(
            "y",
            BinaryExpression(
                "add",
                BinaryExpression(
                    "multiply", LiteralExpression(2.0), ReferenceExpression("x")
                ),
                LiteralExpression(1.0),
            ),
            ("x",),
            "method",
            "2.0 * [X] + 1.0",
        )
    }
    expression = BinaryExpression(
        "divide", ReferenceExpression("y"), LiteralExpression(2.0)
    )

    assert _reference_coefficient(expression, "x", constraints) == 1.0


def test_reference_coefficient_rejects_recursive_constraint_cycle() -> None:
    constraints = {
        "y": CompiledConstraint(
            "y",
            ReferenceExpression("z"),
            ("z",),
            "method",
            "[Z]",
        ),
        "z": CompiledConstraint(
            "z",
            ReferenceExpression("y"),
            ("y",),
            "method",
            "[Y]",
        ),
    }

    assert _reference_coefficient(ReferenceExpression("y"), "x", constraints) is None


@pytest.mark.parametrize(
    "expression",
    [
        BinaryExpression(
            "multiply", ReferenceExpression("x"), ReferenceExpression("y")
        ),
        BinaryExpression("divide", ReferenceExpression("x"), ReferenceExpression("y")),
        BinaryExpression(
            "multiply", ReferenceExpression("x"), ReferenceExpression("x")
        ),
        FunctionExpression("f", (ReferenceExpression("x"),)),
        BinaryExpression("divide", ReferenceExpression("x"), LiteralExpression(0.0)),
        BinaryExpression(
            "multiply", LiteralExpression(float("inf")), ReferenceExpression("x")
        ),
        BinaryExpression(
            "multiply",
            ReferenceExpression("x"),
            BinaryExpression(
                "subtract",
                BinaryExpression(
                    "add", ReferenceExpression("y"), LiteralExpression(1.0)
                ),
                ReferenceExpression("y"),
            ),
        ),
        BinaryExpression(
            "divide",
            ReferenceExpression("x"),
            BinaryExpression(
                "add",
                BinaryExpression(
                    "subtract", ReferenceExpression("y"), ReferenceExpression("y")
                ),
                LiteralExpression(1.0),
            ),
        ),
    ],
)
def test_reference_coefficient_rejects_non_affine_forms(expression: object) -> None:
    assert _reference_coefficient(expression, "x") is None


def test_compiler_rejects_nonlinear_relaxation_controller() -> None:
    parameterization = _dq_affine_parameterization()
    nonlinear = CompiledConstraint(
        "r2adq",
        BinaryExpression(
            "multiply", ReferenceExpression("r2dq"), ReferenceExpression("r1")
        ),
        ("r2dq", "r1"),
        "method",
        "[R2DQ] * [R1]",
    )
    parameterization.program.constraints = (nonlinear,)
    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("r2dq", 20.0), ("r1", 1.5), ("eta", 0.0)),
    )
    maximum = float(np.finfo(np.float64).max)

    with pytest.raises(
        FeasibleCoordinateConstructionError,
        match="unsupported affine controller",
    ):
        compile_feasible_coordinates(
            parameterization,
            frame,
            ("r2dq",),
            (0.0,),
            (maximum,),
        )


def _recursive_affine_parameterization(*, nonlinear: bool):
    block = RelaxationPsdBlock(
        "recursive-transverse-a",
        "a",
        ("held", "diagonal"),
        ((0, 1, "cross"),),
    )
    intermediate_expression = (
        BinaryExpression(
            "multiply",
            ReferenceExpression("controller"),
            ReferenceExpression("controller"),
        )
        if nonlinear
        else ReferenceExpression("controller")
    )
    constraints = (
        CompiledConstraint(
            "intermediate",
            intermediate_expression,
            ("controller",),
            "method",
            "[CONTROLLER] * [CONTROLLER]" if nonlinear else "[CONTROLLER]",
        ),
        CompiledConstraint(
            "diagonal",
            BinaryExpression(
                "add",
                ReferenceExpression("intermediate"),
                ReferenceExpression("offset"),
            ),
            ("intermediate", "offset"),
            "method",
            "[INTERMEDIATE] + [OFFSET]",
        ),
    )

    class Parameterization:
        scope_ids = (
            "held",
            "cross",
            "controller",
            "offset",
            "intermediate",
            "diagonal",
        )
        evaluator_identity = f"recursive-affine-{nonlinear}"
        program = SimpleNamespace(
            constraints=constraints,
            relaxation_domains=SealedRelaxationDomains((block,)),
        )

        @staticmethod
        def resolve(frame):
            values = dict(frame.ordered_items())
            controller = values["controller"]
            values["intermediate"] = (
                controller * controller if nonlinear else controller
            )
            values["diagonal"] = values["intermediate"] + values["offset"]
            return values

    return Parameterization()


def test_compiler_discovers_recursive_affine_controller() -> None:
    parameterization = _recursive_affine_parameterization(nonlinear=False)
    maximum = float(np.finfo(np.float64).max)
    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("held", 2.0), ("cross", 0.0), ("controller", 2.0), ("offset", -1.0)),
    )

    chart = compile_feasible_coordinates(
        parameterization,
        frame,
        ("controller",),
        (0.0,),
        (maximum,),
    )

    assert chart is not None
    assert chart.rate_excess_ids == (("controller", "diagonal"),)
    values = parameterization.resolve(chart.decode((0.0,)).frame)
    assert values["controller"] == pytest.approx(1.0)
    assert values["diagonal"] == pytest.approx(0.0)


def test_compiler_rejects_recursive_nonlinear_controller() -> None:
    parameterization = _recursive_affine_parameterization(nonlinear=True)
    maximum = float(np.finfo(np.float64).max)
    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("held", 2.0), ("cross", 0.0), ("controller", 2.0), ("offset", -1.0)),
    )

    with pytest.raises(
        FeasibleCoordinateConstructionError,
        match="unsupported affine controller",
    ):
        compile_feasible_coordinates(
            parameterization,
            frame,
            ("controller",),
            (0.0,),
            (maximum,),
        )


def test_shipped_v71hg2_compiles_complete_dq_tq_state_mapping() -> None:
    root = Path(__file__).parent.parent
    example = root / "examples/Combinations/CPMG_CH3_1H_DQ_TQ"
    session = AnalysisSession.create()
    session.set_model("2st")
    plan = read_method_plan([example / "Methods/method.toml"])
    experiments = build_experiments(
        sorted((example / "Experiments").glob("*.toml")),
        Selection(include=[SpinSystem.from_name("V71HG2")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(
        read_defaults([example / "Parameters/parameters.toml"])
    )
    assert session.try_build_analysis_values()
    parameterization = session.compile_parameterization_from_actions(
        plan.effective_role_actions()["STEP1"], experiments.param_ids
    )
    engine = EvaluationEngine.from_experiments(experiments, parameterization)
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None

    problem = OptimizationProblem.from_native(
        engine.plan,
        parameterization,
        configuration,
        session.analysis_values.snapshot(),
    )

    chart = problem.feasible_coordinates
    assert chart is not None
    dq_rate = "__R2DQ_A_V71HG2_700_2MHZ"
    tq_rate = "__R2TQ_A_V71HG2_700_2MHZ"
    assert chart.rate_excess_ids == (
        (dq_rate, "__R2ADQ_A_V71HG2_700_2MHZ"),
        (dq_rate, "__R2ADQ_B_V71HG2_700_2MHZ"),
        (tq_rate, "__R2ATQ_A_V71HG2_700_2MHZ"),
        (tq_rate, "__R2ATQ_B_V71HG2_700_2MHZ"),
    )
    private = list(chart.solver_start)
    private[chart.controlled_ids.index(dq_rate)] = 0.0
    private[chart.controlled_ids.index(tq_rate)] = 0.0
    resolved = parameterization.resolve(chart.decode(tuple(private)).frame)
    expected = {
        dq_rate: 0.5,
        "__R2DQ_B_V71HG2_700_2MHZ": 0.5,
        "__R2ADQ_A_V71HG2_700_2MHZ": 0.0,
        "__R2ADQ_B_V71HG2_700_2MHZ": 0.0,
        tq_rate: 1.5,
        "__R2TQ_B_V71HG2_700_2MHZ": 1.5,
        "__R2ATQ_A_V71HG2_700_2MHZ": 0.0,
        "__R2ATQ_B_V71HG2_700_2MHZ": 0.0,
    }
    assert {param_id: resolved[param_id] for param_id in expected} == expected


def test_compiler_returns_no_chart_when_feasibility_has_no_execution_effect() -> None:
    class Parameterization:
        evaluator_identity = "uncoupled-test-evaluator"
        scope_ids = ("rate",)
        program = SimpleNamespace(
            constraints=(),
            relaxation_domains=SealedRelaxationDomains(()),
        )

        @staticmethod
        def resolve(frame):
            return dict(frame.ordered_items())

    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("rate", 5.0),),
    )

    chart = compile_feasible_coordinates(
        Parameterization(),
        frame,
        ("rate",),
        (0.0,),
        (10.0,),
    )

    assert chart is None


def test_compiled_dynamic_affine_floor_spans_only_the_exact_psd_domain() -> None:
    parameterization = _affine_parameterization()
    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("r2", 5.0), ("r1", 2.0), ("eta", 3.0)),
    )
    maximum = float(np.finfo(np.float64).max)
    chart = compile_feasible_coordinates(
        parameterization,
        frame,
        ("r2", "r1"),
        (0.0, 0.0),
        (maximum, maximum),
    )

    assert chart is not None
    assert chart.has_coordinate_transform
    boundary = chart.decode((0.0, 10.0))
    values = parameterization.resolve(boundary.frame)
    assert values["r2"] == pytest.approx(0.5 * (10.0 + np.sqrt(136.0)))
    assert values["r2"] * values["r2a"] == pytest.approx(9.0)
    _validate_blocks(chart.blocks, values)


def test_compiled_static_affine_floor_supersedes_nonnegative_rate_excess() -> None:
    parameterization = _affine_parameterization()
    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("r2", 5.0), ("r1", 2.0), ("eta", 3.0)),
    )
    maximum = float(np.finfo(np.float64).max)
    chart = compile_feasible_coordinates(
        parameterization,
        frame,
        ("r2",),
        (0.0,),
        (maximum,),
    )

    assert chart is not None
    boundary = chart.decode((0.0,))
    values = parameterization.resolve(boundary.frame)
    assert values["r2"] == pytest.approx(0.5 * (2.0 + np.sqrt(40.0)))
    assert values["r2"] * values["r2a"] == pytest.approx(9.0)


def test_component_projection_rejects_a_non_closed_root_feasibility_domain() -> None:
    parameterization = _affine_parameterization()
    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("r2", 5.0), ("r1", 2.0), ("eta", 3.0)),
    )
    maximum = float(np.finfo(np.float64).max)
    root = compile_feasible_coordinates(
        parameterization,
        frame,
        ("r2", "eta"),
        (0.0, -maximum),
        (maximum, maximum),
    )

    assert root is not None
    with pytest.raises(TypeError):
        dataclasses.replace(
            root,
            controlled_domain_groups=(frozenset(),),
        )
    with pytest.raises(
        FeasibleCoordinateConstructionError,
        match="closed root-chart subset",
    ):
        root.project_component(
            frame,
            ("r2",),
            (0.0,),
            (maximum,),
        )


def test_collapsed_cross_rate_chart_has_singular_public_differential() -> None:
    parameterization = _affine_parameterization()
    frame = IndependentValueFrame(
        "parameterization",
        "program",
        "occurrence",
        0,
        (("eta", 0.0),),
    )
    maximum = float(np.finfo(np.float64).max)
    chart = compile_feasible_coordinates(
        parameterization,
        frame,
        ("eta",),
        (-maximum,),
        (maximum,),
    )

    assert chart is not None
    assert chart.has_coordinate_transform
    assert np.linalg.matrix_rank(chart.differential(chart.solver_start)) == 0
    values = parameterization.resolve(chart.decode((1.0,)).frame)
    assert values["eta"] == 0.0
    assert relaxation_state_is_on_boundary(
        parameterization,
        values,
    )


def test_invalid_simulation_state_fails_before_experiment_execution(
    tmp_path: Path,
) -> None:
    block = RelaxationPsdBlock(
        "transverse-a",
        "a",
        ("r2", "r2a"),
        ((0, 1, "eta"),),
    )
    parameterization = SimpleNamespace(
        scope_ids=("r2", "r2a", "eta"),
        program=SimpleNamespace(
            constraints=(),
            relaxation_domains=SealedRelaxationDomains((block,)),
        ),
    )

    class Experiments:
        @staticmethod
        def prepare_for_simulation(_values):
            raise AssertionError("scientific execution must not be reached")

    with pytest.raises(ScientificFeasibilityError, match="positive semidefinite"):
        execute_simulation(
            Experiments(),  # ty: ignore[invalid-argument-type]
            tmp_path,
            parameter_values={"r2": 1.0, "r2a": 1.0, "eta": 2.0},
            parameter_model=SimpleNamespace(),  # ty: ignore[invalid-argument-type]
            parameterization=parameterization,  # ty: ignore[invalid-argument-type]
        )


def test_unrelated_fixed_psd_boundary_does_not_block_fitted_uncertainty() -> None:
    block = RelaxationPsdBlock(
        "transverse-a",
        "a",
        ("r2", "r2a"),
        ((0, 1, "eta"),),
    )
    parameterization = SimpleNamespace(
        scope_ids=("r2", "r2a", "eta", "pb"),
        program=SimpleNamespace(
            constraints=(),
            relaxation_domains=SealedRelaxationDomains((block,)),
        ),
    )
    values = {"r2": 1.0, "r2a": 1.0, "eta": 1.0, "pb": 0.1}

    assert relaxation_state_is_on_boundary(parameterization, values)
    assert not relaxation_state_is_on_boundary(
        parameterization,
        values,
        controlled_ids=("pb",),
    )

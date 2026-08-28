from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemex.optimize.helper import execute_simulation
from chemex.parameters.feasible_coordinates import (
    ScientificFeasibilityError,
    _conditional_interval,
    _validate_blocks,
    compile_feasible_coordinates,
    relaxation_state_is_on_boundary,
)
from chemex.parameters.parameterization import (
    BinaryExpression,
    CompiledConstraint,
    IndependentValueFrame,
    ReferenceExpression,
)
from chemex.parameters.relaxation import (
    RelaxationPsdBlock,
    SealedRelaxationDomains,
)


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

    boundary = chart.decode((0.0,))
    values = parameterization.resolve(boundary.frame)
    assert values["r2"] == pytest.approx(0.5 * (2.0 + np.sqrt(40.0)))
    assert values["r2"] * values["r2a"] == pytest.approx(9.0)


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

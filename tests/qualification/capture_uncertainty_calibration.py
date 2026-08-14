"""Direct prospective acquisition runner for issue #601."""

# ruff: noqa: I001 -- scientific imports are intentionally delayed until attestation.

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections.abc import Iterator, Mapping
from dataclasses import replace
from itertools import product
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, cast

from chemex.numerical_lanes import (
    LaneAttestation,
    NumericalLane,
    RuntimeEnvironment,
    canonical_lanes,
)

# Frozen #581 inputs.  Keep literals compact and single-sourced.
# fmt: off
SPECIFICATION_ID = "chemex-uncertainty-calibration-v1"
SCALE_EXPONENTS = tuple(range(-8, 9))
TAU_REL = (0.0,) + tuple(2.0**-p for p in range(8, 49, 2))
KAPPA = tuple(2.0**p for p in range(13))
EXTENTS = (0, 2, 4, 8, 12, 16)
SVD_DRIVERS = ("gesvd", "gesdd")
RANK_GRID = (0.0,) + tuple(2.0**-p for p in range(16, 53, 2))
WEAK_GRID = (0.0,) + tuple(2.0**-p for p in range(2, 53, 2))
CLUSTER_GRID = WEAK_GRID
CONDITION_GRID = (1.0,) + tuple(2.0**p for p in range(2, 53, 2))
CORRELATION_GRID = (0.0,) + tuple(2.0**p for p in range(13))
COUNTS = (17, 10_296, 2, 400, 27, 27, 27, 14)
TOL = {"derivative": 2.0**-24, "covariance": 2.0e-6, "holdout_covariance": 5.0e-6, "projector": 2.0e-12, "holdout_projector": 5.0e-12, "conditioning": 2.0e-13, "holdout_conditioning": 5.0e-13, "correlation_interior_eps": 8.0, "boundary_eps": 128.0, "boundary_zeta": 3.0, "partial_eps": 64.0, "truth_uncertainty": 1.0 / 16.0}
CEILINGS = {"canonical": {"evaluation": 28_000, "scalar": 4_000, "svd": 32, "correlation": 24, "offline": 12_000_000}, "compatibility": {"evaluation": 751, "scalar": 266, "svd": 4, "correlation": 4, "offline": 0}}
PLANNED_MAXIMUM = {"evaluation": 26_360, "scalar": 3_586, "svd": 26, "correlation": 18, "offline": 12_000_000}
CASE_REQUESTS = {"centered": 273, "one_sided": 205, "three_coordinate": 409, "function": 266}
SUPPORTED = (
    ("solvent_fraction", "fraction", 1.0, "1", ("d2o",)),
    ("state_population", "fraction", 1.0, "1", ("pa", "pb", "pc", "pd", "pe", "pf", "pop_b")),
    ("equilibrium_ratio", "dimensionless", 1.0, "1", ("keq", "keq_bc", "keq_cd", "keq_l", "keq_pl")),
    ("phase_factor", "dimensionless", 1.0, "1", ("phi", "phi_a", "phi_b")),
    ("order_parameter", "dimensionless", 1.0, "1", ("s2", "s2_*")),
    ("directional_exchange_rate", "rate", 100.0, "s^-1", tuple(f"k{a}{b}" for a in "abcdef" for b in "abcdef" if a != b)),
    ("exchange_rate", "rate", 100.0, "s^-1", ("kex", "kex_ab", "kex_ac", "kex_ad", "kex_ae", "kex_af", "kex_bc", "kex_bd", "kex_be", "kex_bf", "kex_cd", "kex_ce", "kex_cf", "kex_de", "kex_df", "kex_ef")),
    ("dissociation_rate", "rate", 100.0, "s^-1", ("koff", "koff1", "koff2", "koff_ab", "koff_ac", "koff_bc")),
    ("hydrogen_exchange_rate", "rate", 1.0, "s^-1", ("kdh", "kdh_a", "kdh_b")),
    ("generated_hydrogen_exchange_rate", "rate", 1.0, "s^-1", ("khh", "khh_*")),
    ("longitudinal_relaxation", "rate", 1.0, "s^-1", ("r1", "r1_*")),
    ("transverse_relaxation", "rate", 10.0, "s^-1", ("r2", "r2_*", "r2a", "r2a_*", "r2mq", "r2mq_*", "r1a", "r1a_*")),
    ("cross_relaxation", "rate", 10.0, "s^-1", ("etaxy", "etaxy_*", "etaz", "etaz_*", "sigma", "sigma_*", "mu", "mu_*")),
    ("chemical_shift", "frequency", 1.0, "ppm", ("cs", "cs_*")),
    ("chemical_shift_difference", "frequency", 1.0, "ppm", ("dw", "dw_*", "dwp", "dwp_*")),
    ("coupling", "frequency", 100.0, "Hz", ("j", "j_*", "d", "d_*")),
)
UNSUPPORTED = frozenset(["c_dimer", "c_l", "c_monomer", "c_p", "c_pl", "c_pl1", "c_pl2", "c_pl3", "c_tetramer", "c_trimer", "dimer", "l1_free", "l_free", "monomer", "p_free", "pfree", "pl", "pl1", "pl2", "pl3", "kd", "kd1", "kd2", "kd_ab", "kd_ac", "kd_app", "kd_bc", "kd_eff", "kon", "kon1", "kon2", "kon_ab", "kon_ac", "kon_bc", "dh_ab", "dh_ac", "dh_ad", "dh_b", "dh_bc", "dh_bd", "dh_c", "dh_cd", "dh_d", "ds_ab", "ds_ac", "ds_ad", "ds_b", "ds_bc", "ds_bd", "ds_c", "ds_cd", "ds_d", "tauc", "dwm"])
UNSUPPORTED_PREFIXES = ("tauc_", "dwm_a", "dwm_i_a", "dwm_s_a")
# fmt: on
# fmt: off

# The remaining direct runner keeps long scientific tuples and typed constructor
# calls on one line so the qualification artifact stays reviewably small.

_MACHINERY_IMPORTED = False
np: Any = None
scipy_linalg: Any = None
migration_core_authority_selection: Any = None

if TYPE_CHECKING:
    from chemex.containers.data import Data
    from chemex.containers.profile import Profile
    from chemex.evaluation.native import EvaluationEngine, EvaluationFailure, EvaluationFrame, EvaluationResult
    from chemex.optimize.direct_trf import AcceptedFitResult, AffineHalfSpace, OptimizationProblem
    from chemex.optimize.uncertainty import ClaimState, FunctionAnalyticPartialCapability, FunctionFiniteDifferenceCapability, ParameterUnit, ResidualVarianceScaling, UncertaintyPolicy, _column_orientations, _expected_correlation_entries, _stencil_estimates, _stencil_vectors, compile_constraint_linearization_capabilities, derive_uncertainty_evidence
    from chemex.parameters.parameterization import ActiveParameterization, BinaryExpression, CompiledConstraint, ConstraintProgram, FunctionExpression, LiteralExpression, ParameterRole, ReferenceExpression, ScientificFunctionBinder
    from chemex.parameters.spin_system import SpinSystem
    from chemex.parameters.values import AnalysisValuesSnapshot
    from chemex.printers.data import Printer


def frozen_digest() -> str:
    record = (SCALE_EXPONENTS, TAU_REL, KAPPA, EXTENTS, SVD_DRIVERS, RANK_GRID, WEAK_GRID, CONDITION_GRID, CORRELATION_GRID, SUPPORTED, tuple(sorted(UNSUPPORTED)), UNSUPPORTED_PREFIXES, TOL, CEILINGS)
    return hashlib.sha256(repr(record).encode()).hexdigest()


def finite_difference_policies() -> Iterator[tuple[float, float, int, int]]:
    return product(TAU_REL, KAPPA, EXTENTS, EXTENTS)


def classify_scale_name(name: str) -> tuple[str, str, float, str]:
    for family, stratum, q0, unit, patterns in SUPPORTED:
        if name in patterns or any(item.endswith("*") and name.startswith(item[:-1]) for item in patterns):
            return family, stratum, q0, unit
    if name in UNSUPPORTED or name.startswith(UNSUPPORTED_PREFIXES):
        raise KeyError(f"catalogued but unsupported after calibration: {name}")
    raise KeyError(f"not catalogued: {name}")


def normalized_residuals(a: float, b: float, t: tuple[float, ...], y: tuple[float, ...], e: tuple[float, ...]) -> tuple[float, ...]:
    c, w = tuple(a + b * value for value in t), tuple(value**-2 for value in e)
    alpha = sum(obs * calc * weight for obs, calc, weight in zip(y, c, w, strict=True)) / sum(calc * calc * weight for calc, weight in zip(c, w, strict=True))
    return tuple((alpha * calc - obs) / error for calc, obs, error in zip(c, y, e, strict=True))


def normalization_jacobian(a: float, b: float, t: tuple[float, ...], y: tuple[float, ...], e: tuple[float, ...]) -> tuple[tuple[float, float], ...]:
    c, w = tuple(a + b * value for value in t), tuple(value**-2 for value in e)
    n = sum(obs * calc * weight for obs, calc, weight in zip(y, c, w, strict=True))
    d = sum(calc * calc * weight for calc, weight in zip(c, w, strict=True))
    alpha = n / d
    columns = []
    for dc in ((1.0,) * len(t), t):
        dn = sum(obs * value * weight for obs, value, weight in zip(y, dc, w, strict=True))
        dd = 2.0 * sum(calc * value * weight for calc, value, weight in zip(c, dc, w, strict=True))
        da = (dn * d - n * dd) / d**2
        columns.append(tuple((da * calc + alpha * value) / error for calc, value, error in zip(c, dc, e, strict=True)))
    return tuple(zip(*columns, strict=True))


def scientific_value(a: float, b: float) -> float:
    return math.exp(a / b) + math.log1p(b * b)


def scientific_partials(a: float, b: float) -> tuple[float, float]:
    exponential = math.exp(a / b)
    return exponential / b, -a * exponential / b**2 + 2.0 * b / (1.0 + b * b)


def scientific_partial_a(a: float, b: float) -> float:
    return scientific_partials(a, b)[0]


def scientific_partial_b(a: float, b: float) -> float:
    return scientific_partials(a, b)[1]


def holdout_value(a: float, b: float) -> float:
    return math.log1p(a * a) + math.sin(b)


def holdout_partials(a: float, b: float) -> tuple[float, float]:
    return 2.0 * a / (1.0 + a * a), math.cos(b)


def holdout_partial_a(a: float, b: float) -> float:
    return holdout_partials(a, b)[0]


def holdout_partial_b(a: float, b: float) -> float:
    return holdout_partials(a, b)[1]


def _matrix(name: str) -> tuple[tuple[float, ...], ...]:
    if name == "B1":
        return ((1.0e-9, 0.0), (0.0, 2.0e-21))
    if name == "B2":
        return ((1.0, 0.0), (0.0, 5.0e-13))
    if name == "B3":
        return ((1.0, 1.0), (2.0, 2.0), (3.0, 3.0))
    if name == "F2":
        return ((1.0, 0.0), (0.0, 2.0), (1.0, 1.0), (2.0, -1.0))
    if name == "H5":
        return ((1.0, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, -1.0))
    if name == "H4":
        r = np.asarray(((0.6, -0.8, 0.0), (0.8, 0.6, 0.0), (0.0, 0.0, 1.0)))
        s = np.asarray((3.0, 3.0 * (1.0 - 7.5e-11), 3.0e-7))
        return tuple(tuple(float(value) for value in row) for row in np.diag(s) @ r.T)
    raise KeyError(name)


def _setup(name: str, q0: float = 1.0) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[tuple[float, ...], ...]]:  # noqa: C901
    if name == "S_CENTER":
        return (0.25 * q0,), (q0,), (-10.0 * q0,), (10.0 * q0,), ()
    if name == "S_BOUND":
        return (0.0,), (q0,), (0.0,), (10.0 * q0,), ()
    if name == "A1":
        return (0.25, -20.0), (1.0, 100.0), (-10.0, -1000.0), (10.0, 1000.0), ()
    if name == "A2":
        return (0.0, 2.0), (1.0, 10.0), (0.0, -10.0), (10.0, 10.0), ((0.0, 1.0, 2.0),)
    if name == "A3":
        return (3.0, 1.0), (1.0, 1.0), (-10.0, -10.0), (10.0, 10.0), ()
    if name == "F2":
        return (0.5, 2.0), (1.0, 2.0), (-10.0, -10.0), (10.0, 10.0), ()
    if name == "H1":
        return (3.0, -0.002), (10.0, 0.01), (-100.0, -1.0), (100.0, 1.0), ()
    if name == "H2":
        return (1.0, -5.0), (1.0, 10.0), (-2.0, -10.0), (1.0, 10.0), ((0.0, -1.0, 5.0),)
    if name == "H3":
        return (2.0, 1.5), (1.0, 2.0), (-10.0, -10.0), (10.0, 10.0), ()
    if name == "H4":
        return (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (-10.0,) * 3, (10.0,) * 3, ()
    if name == "H5":
        return (0.25, -0.5), (0.5, 1.0), (-10.0, -10.0), (10.0, 10.0), ()
    raise KeyError(name)


def truth_jacobian(name: str, q0: float = 1.0) -> tuple[tuple[float, ...], ...]:
    if name == "S_CENTER":
        return ((math.exp(0.25) / q0,),)
    if name == "S_BOUND":
        return ((1.0 / q0,),)
    if name == "A1":
        return ((math.exp(0.25), -0.2 / 100.0), (math.cos(0.25), -2.0 / 100.0), (-0.2, 0.25 / 100.0))
    if name == "A2":
        return ((1.0, math.cos(0.2) / 10.0), (0.0, 0.3))
    if name == "A3":
        return normalization_jacobian(3.0, 1.0, (-2.0, -1.0, 1.0, 3.0), (1.4, 2.1, 3.8, 5.9), (1.0, 0.5, 2.0, 1.5))
    if name == "H1":
        u, v = 0.3, -0.2
        return (((1.0 - math.tanh(u) ** 2) / 10.0, math.exp(v) / 0.01), ((2.0 * u / (1.0 + u * u)) / 10.0, 3.0 * v * v / 0.01), (-math.sin(u - v) / 10.0, math.sin(u - v) / 0.01))
    if name == "H2":
        return ((math.e, math.cos(-0.5) / 10.0), (2.0, 0.3))
    if name == "H3":
        return normalization_jacobian(2.0, 1.5, (-3.0, 0.0, 2.0, 4.0, 5.0), (0.5, 2.2, 4.4, 7.1, 8.8), (0.75, 1.25, 0.5, 2.0, 1.0))
    return _matrix(name)


def _calculation(name: str, vector: tuple[float, ...], q0: float) -> tuple[float, ...]:
    if name in {"S_CENTER", "S_BOUND"}:
        return (math.exp(vector[0] / q0),)
    if name == "A1":
        u, v = vector[0], vector[1] / 100.0
        return (math.exp(u) + 0.5 * v * v, math.sin(u) - 2.0 * v, u * v + 0.25)
    if name in {"A2", "H2"}:
        return (math.exp(vector[0]) + math.sin(vector[1] / 10.0), vector[0] ** 2 + 0.3 * vector[1])
    if name in {"A3", "H3"}:
        return tuple(vector[0] + vector[1] * value for value in ((-2.0, -1.0, 1.0, 3.0) if name == "A3" else (-3.0, 0.0, 2.0, 4.0, 5.0)))
    if name == "H1":
        u, v = vector[0] / 10.0, vector[1] / 0.01
        return (math.tanh(u) + math.exp(v), math.log1p(u * u) + v**3, math.cos(u - v))
    matrix, anchor = _matrix(name), _setup(name)[0]
    return tuple(sum(coef * (value - base) for coef, value, base in zip(row, vector, anchor, strict=True)) for row in matrix)


def _import_machinery() -> None:
    global ActiveParameterization, AffineHalfSpace, AnalysisValuesSnapshot, AcceptedFitResult, BinaryExpression, ClaimState, CompiledConstraint, ConstraintProgram, Data, EvaluationEngine, EvaluationFailure, EvaluationFrame, EvaluationResult, FunctionAnalyticPartialCapability, FunctionExpression, FunctionFiniteDifferenceCapability, LiteralExpression, OptimizationProblem, ParameterRole, ParameterUnit, Printer, Profile, ReferenceExpression, ResidualVarianceScaling, ScientificFunctionBinder, SpinSystem, UncertaintyPolicy, compile_constraint_linearization_capabilities, derive_uncertainty_evidence, migration_core_authority_selection, np, scipy_linalg, _column_orientations, _expected_correlation_entries, _stencil_estimates, _stencil_vectors, _MACHINERY_IMPORTED
    if _MACHINERY_IMPORTED:
        return
    # isort: off -- imports must remain behind live lane attestation.
    import numpy as np
    from pydantic import BaseModel
    from scipy import linalg as scipy_linalg

    from chemex.containers.data import Data
    from chemex.containers.profile import Profile
    from chemex.evaluation.native import EvaluationEngine, EvaluationFailure, EvaluationFrame, EvaluationResult
    from chemex.migration_core import migration_core_authority_selection
    from chemex.optimize.direct_trf import AcceptedFitResult, AffineHalfSpace, OptimizationProblem
    from chemex.optimize.uncertainty import ClaimState, FunctionAnalyticPartialCapability, FunctionFiniteDifferenceCapability, ParameterUnit, ResidualVarianceScaling, UncertaintyPolicy, _column_orientations, _expected_correlation_entries, _stencil_estimates, _stencil_vectors, compile_constraint_linearization_capabilities, derive_uncertainty_evidence
    from chemex.parameters.parameterization import ActiveParameterization, BinaryExpression, CompiledConstraint, ConstraintProgram, FunctionExpression, LiteralExpression, ParameterRole, ReferenceExpression, ScientificFunctionBinder
    from chemex.parameters.spin_system import SpinSystem
    from chemex.parameters.values import AnalysisValuesSnapshot
    from chemex.printers.data import Printer
    # isort: on

    class _Settings(BaseModel):
        kind: str = "issue-601-truth"

    _Pulse.settings = _Settings()
    _MACHINERY_IMPORTED = True


class _Spectrometer:
    def __init__(self, name: str, q0: float, params: tuple[str, ...]) -> None:
        self.name, self.q0, self.params = name, q0, params
        self.spin_system = SpinSystem.from_name(name)
        self.values = {param.lower(): 0.0 for param in params}

    def update(self, values: dict[str, float]) -> None:
        self.values.update(values)

    def new_native_workspace(self) -> _Spectrometer:
        return _Spectrometer(self.name, self.q0, self.params)

    def native_kernel_descriptor(self) -> dict[str, object]:
        return {"kind": "issue-601", "case": self.name, "q0": self.q0}


class _Pulse:
    settings: Any = None

    def calculate(self, spectrometer: _Spectrometer, data: Any) -> Any:
        return np.asarray(_calculation(spectrometer.name, tuple(spectrometer.values[param.lower()] for param in spectrometer.params), spectrometer.q0))

    def is_reference(self, metadata: Any) -> Any:
        return np.zeros_like(metadata, dtype=bool)


def _derived_program(name: str, params: tuple[str, ...]) -> tuple[Any, tuple[str, ...], tuple[Any, ...]]:
    if name not in {"F2", "H5"}:
        binder = ScientificFunctionBinder(name, {})
        return binder, params, ()
    function = scientific_value if name == "F2" else holdout_value
    binder = ScientificFunctionBinder(name, {"f": function})
    a, b = ReferenceExpression("A"), ReferenceExpression("B")
    constraints = (
        CompiledConstraint("SUM", BinaryExpression("add", a, BinaryExpression("multiply", LiteralExpression(2.0), b)), ("A", "B"), "issue601", "A+2*B"),
        CompiledConstraint("SCI", FunctionExpression("f", (a, b)), ("A", "B"), "issue601", "f(A,B)"),
        CompiledConstraint("ZERO", BinaryExpression("add", LiteralExpression(1.0), BinaryExpression("multiply", LiteralExpression(0.0), a)), ("A",), "issue601", "1+0*A"),
    )
    return binder, (*params, "SUM", "SCI", "ZERO"), constraints


def _fixture(name: str, q0: float = 1.0) -> tuple[Any, Any, Any, Any]:
    anchor, _scales, lower, upper, halfspaces = _setup(name, q0)
    params = tuple(chr(65 + i) for i in range(len(anchor)))
    binder, scope, constraints = _derived_program(name, params)
    derived = tuple(item.target_id for item in constraints)
    program = ConstraintProgram("parameter-model", "model", "definitions", "configuration", binder.identity, scope, params, derived, constraints, derived)
    roles = tuple(
        (param, ParameterRole.FIT if param in params else ParameterRole.DERIVED)
        for param in scope
    )
    parameterization = ActiveParameterization(program, binder, f"issue-601-{name}", 601, roles)
    snapshot = AnalysisValuesSnapshot(parameterization.occurrence_identity, "model", "definitions", "configuration", 601, tuple(zip(params, anchor, strict=True)))
    if name == "A3":
        observed, errors, metadata, scaled = (1.4, 2.1, 3.8, 5.9), (1.0, 0.5, 2.0, 1.5), (-2.0, -1.0, 1.0, 3.0), True
    elif name == "H3":
        observed, errors, metadata, scaled = (0.5, 2.2, 4.4, 7.1, 8.8), (0.75, 1.25, 0.5, 2.0, 1.0), (-3.0, 0.0, 2.0, 4.0, 5.0), True
    else:
        size = len(_calculation(name, anchor, q0))
        observed, errors, metadata, scaled = (0.0,) * size, (1.0,) * size, tuple(float(i) for i in range(size)), False
    data = Data(exp=np.asarray(observed), err=np.asarray(errors), metadata=np.asarray(metadata))
    profile = Profile(data, cast("Any", _Spectrometer(name, q0, params)), cast("Any", _Pulse()), {param.lower(): param for param in params}, cast("Printer", None), is_scaled=scaled)
    engine = EvaluationEngine.from_experiments(cast("Any", (SimpleNamespace(profiles=(profile,)),)), parameterization)
    restrictions = tuple(AffineHalfSpace(f"{name}-{i}", tuple(row[:-1]), row[-1]) for i, row in enumerate(halfspaces))
    problem = OptimizationProblem(engine.plan.identity, parameterization.identity, parameterization.evaluator_identity, program.fingerprint, "configuration", snapshot, tuple(zip(params, anchor, strict=True)), params, (), anchor, lower, upper, scope, affine_half_spaces=restrictions)
    evaluation = engine.new_evaluator().evaluate(EvaluationFrame.from_lifecycle_frame(parameterization, problem.lifecycle_frame(anchor, parameterization)))
    assert isinstance(evaluation, EvaluationResult)
    accepted = AcceptedFitResult.for_qualification(occurrence_identity=f"issue-601-{name}-accepted", problem_identity=problem.identity, invocation_identity="issue-601", execution_identity="issue-601", materialization_identity="issue-601", parameterization_identity=parameterization.identity, evaluator_parameterization_identity=parameterization.evaluator_identity, source_occurrence_identity=snapshot.occurrence_identity, source_revision=snapshot.revision, controlled_ids=params, vector=anchor, chi_square=float(np.dot(evaluation.residuals, evaluation.residuals)), evaluation_result=evaluation, commit_scope=scope, commit_items=tuple(zip(params, anchor, strict=True)), origin_context_identity="issue-601")
    return parameterization, engine, problem, accepted


def _spend(counts: dict[str, int], role: str, kind: str, amount: int = 1) -> None:
    if counts[kind] + amount > CEILINGS[role][kind]:
        raise RuntimeError(f"{role} {kind} ceiling exceeded")
    counts[kind] += amount


def _accepted_copy(accepted: Any, problem_identity: str, occurrence_identity: str) -> Any:
    return AcceptedFitResult.for_qualification(occurrence_identity=occurrence_identity, problem_identity=problem_identity, invocation_identity=accepted.invocation_identity, execution_identity=accepted.execution_identity, materialization_identity=accepted.materialization_identity, parameterization_identity=accepted.parameterization_identity, evaluator_parameterization_identity=accepted.evaluator_parameterization_identity, source_occurrence_identity=accepted.source_occurrence_identity, source_revision=accepted.source_revision, controlled_ids=accepted.controlled_ids, vector=accepted.vector, chi_square=accepted.chi_square, evaluation_result=accepted.evaluation_result, commit_scope=accepted.commit_scope, commit_items=accepted.commit_items, origin_context_identity=accepted.origin_context_identity)


def _evaluate(
    context: tuple[Any, Any, Any, Any],
    evaluator: Any,
    vector: tuple[float, ...],
    counts: dict[str, int],
    role: str,
) -> tuple[float, ...]:
    parameterization, _engine, problem, _accepted = context
    _spend(counts, role, "evaluation")
    result = evaluator.evaluate(
        EvaluationFrame.from_lifecycle_frame(
            parameterization, problem.lifecycle_frame(vector, parameterization)
        )
    )
    if isinstance(result, EvaluationFailure):
        raise RuntimeError(f"truth evaluation failed: {result.category}")  # noqa: TRY004
    return tuple(float(value) for value in result.residuals)


def _union_trajectories(
    name: str, q0: float, counts: dict[str, int], role: str
) -> tuple[tuple[object, ...], ...]:
    context = _fixture(name, q0)
    _parameterization, engine, problem, accepted = context
    evaluator = engine.new_evaluator()
    base = _evaluate(context, evaluator, accepted.vector, counts, role)
    truth = truth_jacobian(name, q0)
    trajectories = []
    for column, scale in enumerate(_setup(name, q0)[1]):
        line = problem.coordinate_line_feasibility(accepted.vector, column)
        orientation, maximum = _column_orientations(
            accepted.vector[column],
            accepted.vector[column] + line.minimum_displacement,
            accepted.vector[column] + line.maximum_displacement,
        )[0]
        nominal = (2.0**-52) ** (1.0 / 3.0) * scale
        steps = tuple(
            dict.fromkeys(
                (
                    *tuple(
                        math.ldexp(nominal, exponent)
                        for exponent in range(-24, 25)
                        if math.ldexp(nominal, exponent) <= maximum
                    ),
                    maximum,
                )
            )
        )
        observations = {}
        for step in steps:
            vectors = _stencil_vectors(accepted.vector, column, step, orientation)
            assert vectors is not None
            values = tuple(
                _evaluate(context, evaluator, vector, counts, role)
                for vector in vectors
            )
            displacements = tuple(
                vector[column] - accepted.vector[column] for vector in vectors
            )
            observations[float(step).hex()] = (values, displacements)
        trajectories.append(
            (
                orientation,
                base,
                maximum,
                observations,
                tuple(row[column] for row in truth),
                1.0,
            )
        )
    return tuple(trajectories)


def _scalar_trajectories(counts: dict[str, int], role: str) -> tuple[tuple[object, ...], ...]:
    anchor, scales, output = (0.5, 2.0), (1.0, 2.0), 4.0
    base = scientific_value(*anchor)
    trajectories = []
    for column, scale in enumerate(scales):
        _spend(counts, role, "scalar")
        nominal = (2.0**-52) ** (1.0 / 3.0) * scale
        maximum = math.ldexp(nominal, 24)
        steps = tuple(math.ldexp(nominal, exponent) for exponent in range(-24, 24)) + (maximum,)
        observations = {}
        for step in steps:
            values = []
            for multiplier in (-2.0, -1.0, 1.0, 2.0):
                trial = list(anchor)
                trial[column] += multiplier * step
                _spend(counts, role, "scalar")
                values.append((scientific_value(*trial),))
            observations[float(step).hex()] = (tuple(values), (-2.0 * step, -step, step, 2.0 * step))
        trajectories.append(("centered", (base,), maximum, observations, (scientific_partials(*anchor)[column],), output))
    return tuple(trajectories)


def _policy_steps(q: float, maximum: float, policy: tuple[float, float, int, int]) -> tuple[float, ...]:
    nominal = (2.0**-52) ** (1.0 / 3.0) * q
    smaller, larger = policy[2:]
    candidates = [
        (abs(float(exponent)), math.ldexp(nominal, exponent))
        for exponent in range(-smaller, larger + 1)
        if math.ldexp(nominal, exponent) <= maximum
    ]
    candidates.append((abs(math.log2(maximum / nominal)), maximum))
    return tuple(value for _distance, value in sorted({float(value).hex(): (distance, value) for distance, value in candidates}.values()))


def _replay(trajectory: tuple[object, ...], q: float, output_scale: float, policy: tuple[float, float, int, int], counts: dict[str, int]) -> tuple[bool, bool, float, int]:
    _spend(counts, "canonical", "offline")
    orientation, base, maximum, observations, truth, _q0 = trajectory
    base = cast("tuple[float, ...]", base)
    truth = cast("tuple[float, ...]", truth)
    observations = cast("dict[str, tuple[tuple[tuple[float, ...], ...], tuple[float, ...]]]", observations)
    for attempt, step in enumerate(_policy_steps(q, cast("float", maximum), policy), 1):
        values, displacements = observations[float(step).hex()]
        fine, coarse = _stencil_estimates(cast("str", orientation), values, base, displacements)
        discrepancy = max(abs(left - right) for left, right in zip(fine, coarse, strict=True))
        derivative_scale = max(*(abs(value) for value in fine), *(abs(value) for value in coarse))
        function_scale = max(output_scale, *(abs(value) for value in base), *(abs(value) for row in values for value in row))
        roundoff = policy[1] * 2.0**-52 * function_scale / min(abs(value) for value in displacements)
        if discrepancy <= policy[0] * derivative_scale + roundoff:
            utilization = max(abs(value - exact) / (TOL["derivative"] * (output_scale / q + abs(exact))) for value, exact in zip(fine, truth, strict=True))
            return utilization <= 1.0, utilization > 1.0, utilization, attempt * len(values)
    return False, False, math.inf, len(_policy_steps(q, cast("float", maximum), policy)) * (4 if orientation == "centered" else 3)


def _score(trajectories: tuple[tuple[object, ...], ...], scales: tuple[float, ...], outputs: tuple[float, ...], policy: tuple[float, float, int, int], counts: dict[str, int]) -> tuple[int, float, float, int] | None:
    results = tuple(_replay(item, scale, output, policy, counts) for item, scale, output in zip(trajectories, scales, outputs, strict=True))
    false = sum(result[1] for result in results)
    if false or not all(result[0] for result in results):
        return None
    return false, max(result[2] for result in results), sum(result[2] for result in results), sum(result[3] for result in results)


def select_scales(counts: dict[str, int]) -> tuple[dict[str, float] | None, dict[str, object], tuple[tuple[object, ...], ...], tuple[tuple[object, ...], ...]]:  # noqa: C901
    selected, summaries, retained = {}, {}, []
    for family, _stratum, q0, _unit, _patterns in SUPPORTED:
        pair = _union_trajectories("S_CENTER", q0, counts, "canonical") + _union_trajectories("S_BOUND", q0, counts, "canonical")
        qualified = []
        for exponent in SCALE_EXPONENTS:
            q = math.ldexp(q0, exponent)
            best = None
            for policy in finite_difference_policies():
                score = _score(pair, (q, q), (1.0, 1.0), policy, counts)
                if score is not None:
                    key = (*score, policy[0], policy[1], policy[3], policy[2])
                    best = key if best is None or key < best else best
            if best is not None:
                qualified.append((best, abs(exponent), exponent, q))
        if not qualified:
            return None, {"status": "UNSUPPORTED", "family": family}, (), ()
        if all(item[2] in {-8, 8} for item in qualified):
            return None, {"status": "GRID_SATURATED", "family": family}, (), ()
        choice = min(qualified)
        selected[family] = choice[3]
        summaries[family] = {
            "q0": q0,
            "exponent": choice[2],
            "scale": choice[3],
            "neighbors": tuple((item[2], item[3], "not_selected") for item in qualified if abs(item[2] - choice[2]) == 1),
        }
        retained.extend(pair)
    scalar = _scalar_trajectories(counts, "canonical")
    function_scales = []
    for column, q0 in enumerate((1.0, 2.0)):
        choices = []
        for exponent in SCALE_EXPONENTS:
            q = math.ldexp(q0, exponent)
            best = min(((*score, policy, exponent, q) for policy in finite_difference_policies() if (score := _score((scalar[column],), (q,), (4.0,), policy, counts)) is not None), default=None)
            if best is not None:
                choices.append(best)
        if not choices:
            return None, {"status": "UNSUPPORTED", "family": f"A4-argument-{column}"}, (), ()
        choice = min(choices)
        function_scales.append(choice[-1])
        summaries[f"A4-argument-{column}"] = {"q0": q0, "exponent": choice[-2], "scale": choice[-1]}
    output_choices = []
    for exponent in SCALE_EXPONENTS:
        output = math.ldexp(4.0, exponent)
        best = min(((*score, policy, exponent, output) for policy in finite_difference_policies() if (score := _score(scalar, tuple(function_scales), (output, output), policy, counts)) is not None), default=None)
        if best is not None:
            output_choices.append(best)
    if not output_choices:
        return None, {"status": "UNSUPPORTED", "family": "A4-output"}, (), ()
    output_choice = min(output_choices)
    selected["A4-argument-0"], selected["A4-argument-1"], selected["A4-output"] = *function_scales, output_choice[-1]
    summaries["A4-output"] = {"Q0": 4.0, "exponent": output_choice[-2], "scale": output_choice[-1]}
    return selected, {"status": "SELECTED", "families": summaries}, tuple(retained), scalar


def calibrate_finite_differences(scales: Mapping[str, float], family_trajectories: tuple[tuple[object, ...], ...], scalar: tuple[tuple[object, ...], ...], counts: dict[str, int]) -> tuple[tuple[float, float, int, int] | None, dict[str, object]]:
    normalization = _union_trajectories("A3", 1.0, counts, "canonical")
    trajectories = family_trajectories + normalization + scalar
    family_values = tuple(scales[family] for family, *_rest in SUPPORTED)
    q = tuple(value for value in family_values for _ in (0, 1)) + (1.0, 1.0, scales["A4-argument-0"], scales["A4-argument-1"])
    outputs = (1.0,) * (len(trajectories) - 2) + (scales["A4-output"],) * 2
    best = None
    frontier = []
    non_edge = False
    for policy in finite_difference_policies():
        score = _score(trajectories, q, outputs, policy, counts)
        if score is None:
            continue
        edge = policy == (2.0**-8, 2.0**12, 16, 16)
        non_edge |= not edge
        key = (*score, policy[0], policy[1], policy[3], policy[2], policy)
        if best is None or key < best:
            best = key
        point = (score[1], score[3], policy)
        frontier = [item for item in frontier if not (point[0] <= item[0] and point[1] <= item[1])]
        if not any(item[0] <= point[0] and item[1] <= point[1] for item in frontier):
            frontier.append(point)
    if best is None:
        return None, {"status": "UNSUPPORTED"}
    if not non_edge:
        return None, {"status": "GRID_SATURATED"}
    policy = tuple(best[-1])
    neighbors = []
    axes = (TAU_REL, KAPPA, EXTENTS, EXTENTS)
    for axis, values in enumerate(axes):
        index = values.index(policy[axis])
        for neighbor_index in (index - 1, index + 1):
            if 0 <= neighbor_index < len(values):
                candidate = list(policy)
                candidate[axis] = values[neighbor_index]
                candidate_tuple = cast("tuple[float, float, int, int]", tuple(candidate))
                reason = "qualified" if _score(trajectories, q, outputs, candidate_tuple, counts) else "truth_or_reliability"
                neighbors.append((candidate_tuple, reason))
    return policy, {"status": "SELECTED", "policy": policy, "metrics": best[:4], "frontier": tuple(sorted(frontier)[:32]), "rejected_neighbors": tuple(neighbors)}


def select_svd_driver(counts: dict[str, int]) -> tuple[str | None, dict[str, object]]:
    records = []
    for driver in SVD_DRIVERS:
        spectra, error = [], 0.0
        for name, exact in (("B1", (1.0e-9, 2.0e-21)), ("B2", (1.0, 5.0e-13)), ("B3", (math.sqrt(28.0), 0.0))):
            _spend(counts, "canonical", "svd")
            values = tuple(float(value) for value in scipy_linalg.svd(np.asarray(_matrix(name)), full_matrices=False, compute_uv=False, lapack_driver=driver))
            spectra.append(values)
            error = max(error, max(abs(value - truth) / max(1.0, abs(truth)) for value, truth in zip(values, exact, strict=True)))
        records.append((error, SVD_DRIVERS.index(driver), driver, tuple(spectra)))
    selected = min(records)
    return selected[2], {"status": "SELECTED", "driver": selected[2], "metrics": tuple((item[2], item[0]) for item in records), "spectra": selected[3]}


def select_rank_policy(driver_record: Mapping[str, object], counts: dict[str, int]) -> tuple[tuple[float, float] | None, dict[str, object]]:
    spectra = cast("tuple[tuple[float, ...], ...]", driver_record["spectra"])
    truths = (2, 1, 1)
    qualified = []
    for absolute, relative in product(RANK_GRID, RANK_GRID):
        _spend(counts, "canonical", "offline")
        observed = tuple(sum(value > absolute + relative * spectrum[0] for value in spectrum) for spectrum in spectra)
        if observed == truths:
            qualified.append((absolute, relative))
    if not qualified:
        return None, {"status": "UNSUPPORTED"}
    if all(absolute == 2.0**-16 or relative == 2.0**-16 for absolute, relative in qualified):
        return None, {"status": "GRID_SATURATED"}
    selected = min(qualified)
    index = qualified.index(selected)
    return selected, {"status": "SELECTED", "policy": selected, "qualified_count": len(qualified), "neighbors": tuple(qualified[max(0, index - 1) : index] + qualified[index + 1 : index + 2])}


def _smallest_threshold(grid: tuple[float, ...], classify: Any, edge: float) -> tuple[float | None, dict[str, object]]:
    qualified = tuple(value for value in grid if classify(value))
    if not qualified:
        return None, {"status": "UNSUPPORTED"}
    if qualified == (edge,):
        return None, {"status": "GRID_SATURATED"}
    selected = qualified[0]
    index = grid.index(selected)
    neighbors = tuple((grid[i], "qualified" if grid[i] in qualified else "misclassified") for i in (index - 1, index + 1) if 0 <= i < len(grid))
    return selected, {"status": "SELECTED", "value": selected, "neighbors": neighbors}


def select_weak_policy() -> tuple[float | None, dict[str, object]]:
    ratios, truth = (1.0, 2.0e-6, 2.0e-7), (False, False, True)
    return _smallest_threshold(WEAK_GRID, lambda value: tuple(item <= value for item in ratios) == truth, 2.0**-2)


def select_cluster_policy() -> tuple[float | None, dict[str, object]]:
    spectrum = (4.0, 4.0 * (1.0 - 5.0e-11), 2.0, 2.0 * (1.0 - 5.0e-9), 1.0)
    gaps = tuple((left - right) / spectrum[0] for left, right in zip(spectrum, spectrum[1:], strict=False))
    truth = (True, False, False, False)
    return _smallest_threshold(CLUSTER_GRID, lambda value: tuple(item <= value for item in gaps) == truth, 2.0**-2)


def select_conditioning_policy() -> tuple[float | None, dict[str, object]]:
    return _smallest_threshold(CONDITION_GRID, lambda value: 5.0e5 <= value < 2.0e6, 2.0**52)


def select_correlation_policy(counts: dict[str, int]) -> tuple[float | None, dict[str, object]]:
    eps = 2.0**-52
    values = (1.0 + 64.0 * eps, -(1.0 + 64.0 * eps), 1.0 + 128.0 * eps)

    def classify(multiplier: float) -> bool:
        _spend(counts, "canonical", "correlation")
        return abs(values[0] - 1.0) <= multiplier * eps and abs(values[1] + 1.0) <= multiplier * eps and abs(values[2] - 1.0) > multiplier * eps

    return _smallest_threshold(CORRELATION_GRID, classify, 2.0**12)


def compose_uncertainty_policy(name: str, scales: tuple[float, ...], fd: tuple[float, float, int, int], driver: str, rank: tuple[float, float], weak: float, cluster: float, conditioning: float, correlation: float) -> Any:
    params = tuple(chr(65 + i) for i in range(len(scales)))
    units = tuple((param, ParameterUnit.DIMENSIONLESS) for param in params)
    return UncertaintyPolicy(calibration_identity=SPECIFICATION_ID, numerical_compatibility_requirement="canonical-linux-amd64-python-3.13;bounded-python-3.14-replay", coordinate_scales=tuple(zip(params, scales, strict=True)), coordinate_units=units, residual_variance_scaling=ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES, relative_step_tolerance=fd[0], roundoff_multiplier=fd[1], smaller_step_extent=fd[2], larger_step_extent=fd[3], svd_driver=cast("Any", driver), rank_absolute_tolerance=rank[0], rank_relative_tolerance=rank[1], weak_relative_tolerance=weak, singular_value_cluster_relative_tolerance=cluster, conditioning_limit=conditioning, correlation_roundoff_multiplier=correlation, affine_feasibility_policy="canonical-root-affine-halfspace-zeta-gt-3-v1")


def policy_record(policy: Any) -> dict[str, object]:
    return {"identity": policy.identity, "calibration_identity": policy.calibration_identity, "numerical_compatibility_requirement": policy.numerical_compatibility_requirement, "coordinate_scales": policy.coordinate_scales, "coordinate_units": tuple((name, unit.value) for name, unit in policy.coordinate_units), "residual_variance_scaling": policy.residual_variance_scaling.value, "relative_step_tolerance": policy.relative_step_tolerance, "roundoff_multiplier": policy.roundoff_multiplier, "smaller_step_extent": policy.smaller_step_extent, "larger_step_extent": policy.larger_step_extent, "svd_driver": policy.svd_driver, "rank_absolute_tolerance": policy.rank_absolute_tolerance, "rank_relative_tolerance": policy.rank_relative_tolerance, "weak_relative_tolerance": policy.weak_relative_tolerance, "singular_value_cluster_relative_tolerance": policy.singular_value_cluster_relative_tolerance, "conditioning_limit": policy.conditioning_limit, "correlation_roundoff_multiplier": policy.correlation_roundoff_multiplier, "affine_feasibility_policy": policy.affine_feasibility_policy}


def policy_from_record(record: Mapping[str, object]) -> Any:
    policy = UncertaintyPolicy(calibration_identity=cast("str", record["calibration_identity"]), numerical_compatibility_requirement=cast("str", record["numerical_compatibility_requirement"]), coordinate_scales=tuple((cast("str", item[0]), float(cast("Any", item[1]))) for item in cast("list[list[object]] | tuple[tuple[object, object], ...]", record["coordinate_scales"])), coordinate_units=tuple((cast("str", item[0]), ParameterUnit(cast("str", item[1]))) for item in cast("list[list[object]] | tuple[tuple[object, object], ...]", record["coordinate_units"])), residual_variance_scaling=ResidualVarianceScaling(cast("str", record["residual_variance_scaling"])), relative_step_tolerance=float(cast("float", record["relative_step_tolerance"])), roundoff_multiplier=float(cast("float", record["roundoff_multiplier"])), smaller_step_extent=int(cast("int", record["smaller_step_extent"])), larger_step_extent=int(cast("int", record["larger_step_extent"])), svd_driver=cast("Any", record["svd_driver"]), rank_absolute_tolerance=float(cast("float", record["rank_absolute_tolerance"])), rank_relative_tolerance=float(cast("float", record["rank_relative_tolerance"])), weak_relative_tolerance=float(cast("float", record["weak_relative_tolerance"])), singular_value_cluster_relative_tolerance=float(cast("float", record["singular_value_cluster_relative_tolerance"])), conditioning_limit=float(cast("float", record["conditioning_limit"])), correlation_roundoff_multiplier=float(cast("float", record["correlation_roundoff_multiplier"])), affine_feasibility_policy=cast("str", record["affine_feasibility_policy"]))
    if policy.identity != record.get("identity"):
        raise RuntimeError("selected UncertaintyPolicy identity mismatch")
    return policy


def _case_policy(name: str, selected: Any) -> Any:
    scales = _setup(name)[1]
    fd = (selected.relative_step_tolerance, selected.roundoff_multiplier, selected.smaller_step_extent, selected.larger_step_extent)
    return compose_uncertainty_policy(name, scales, fd, selected.svd_driver, (selected.rank_absolute_tolerance, selected.rank_relative_tolerance), selected.weak_relative_tolerance, selected.singular_value_cluster_relative_tolerance, selected.conditioning_limit, selected.correlation_roundoff_multiplier)


def _derive_case(name: str, selected: Any, environment: str, counts: dict[str, int], role: str, numerical_partial: bool = False, boundary_zeta: float | None = None) -> tuple[Any, dict[str, object]]:
    context = _fixture(name)
    parameterization, engine, problem, accepted = context
    policy = _case_policy(name, selected)
    _spend(counts, role, "evaluation", CASE_REQUESTS["one_sided" if name in {"A2", "H2"} else "centered"])
    scope: Any = ()
    units: Any = ()
    scales: Any = ()
    compiled: Any = None
    if name in {"F2", "H5"}:
        function = scientific_value if name == "F2" else holdout_value
        partials = (scientific_partial_a, scientific_partial_b) if name == "F2" else (holdout_partial_a, holdout_partial_b)
        capability = FunctionFiniteDifferenceCapability("f", None, _setup(name)[1], 4.0) if numerical_partial else FunctionAnalyticPartialCapability("f", None, hashlib.sha256(function.__name__.encode()).hexdigest(), partials)
        if numerical_partial:
            _spend(counts, role, "scalar", CASE_REQUESTS["function"])
        scope = parameterization.scope_ids
        compiled = compile_constraint_linearization_capabilities(parameterization, scope, (capability,))
        units = tuple((item, ParameterUnit.DIMENSIONLESS) for item in scope)
        scales = tuple((item, value) for item, value in zip(scope, (*_setup(name)[1], 3.0, 4.0, 1.0), strict=True))
    if boundary_zeta is not None:
        boundary = AffineHalfSpace("F2-boundary", (1.0, 0.0), accepted.vector[0] + boundary_zeta * math.sqrt(6.0 / 35.0))
        problem = replace(problem, affine_half_spaces=(boundary,))
        accepted = _accepted_copy(accepted, problem.identity, f"F2-boundary-{float(boundary_zeta).hex()}")
    _spend(counts, role, "svd")
    if role == "compatibility":
        _spend(counts, role, "correlation")
    evidence = derive_uncertainty_evidence(accepted, problem=problem, parameterization=parameterization, engine=engine, policy=policy, constrained_scope=scope, constrained_units=units, constrained_scales=scales, compiled_constraint_linearization=compiled, resolved_environment_identity=environment)
    if evidence.residual_jacobian is None:
        return evidence, {"case": name, "passed": False, "reason": "missing_jacobian"}
    observed, truth = np.asarray(evidence.residual_jacobian.matrix), np.asarray(truth_jacobian(name))
    envelope = TOL["derivative"] * (1.0 / np.asarray(_setup(name)[1])[None, :] + np.abs(truth))
    passed = bool(np.all(np.abs(observed - truth) <= envelope))
    if name == "F2" and evidence.covariance is not None:
        expected = np.asarray(((6.0 / 35.0, 1.0 / 35.0), (1.0 / 35.0, 6.0 / 35.0)))
        passed &= bool(np.max(np.abs(np.asarray(evidence.covariance.covariance) - expected)) / max(1.0, float(np.max(np.abs(expected)))) <= TOL["covariance"])
        passed &= evidence.constrained_marginal_errors is not None and evidence.constrained_marginal_errors.entries[-1].raw_value == 0.0
        passed &= evidence.constrained_propagation is not None and evidence.constrained_propagation.covariance[-1][-1] == 0.0
        if boundary_zeta is not None:
            expected_claim = ClaimState.VIOLATED if boundary_zeta == 2.0 else ClaimState.SATISFIED
            passed &= evidence.covariance.claim("AFFINE_FEASIBILITY") is expected_claim
    if name in {"A3", "H3"} and evidence.covariance is not None:
        passed &= evidence.covariance.profiled_normalization_count == 1
    return evidence, {"case": name if boundary_zeta is None else f"F2-boundary-{boundary_zeta:g}", "passed": passed, "policy_identity": policy.identity, "evidence_identity": evidence.identity}


def _correlation_case(raw: float, policy: Any, counts: dict[str, int], role: str) -> bool:
    if role == "compatibility":
        _spend(counts, role, "correlation")
    entries = _expected_correlation_entries(((1.0, raw), (raw, 1.0)), residual_variance_degenerate=False, policy=policy)
    return entries[0][1].outcome == "ENDPOINT_CANONICALIZED" and entries[0][1].value == math.copysign(1.0, raw)


def _h4_spectral_case(selected: Any, counts: dict[str, int], role: str) -> dict[str, object]:
    _spend(counts, role, "svd")
    _u, spectrum, vh = scipy_linalg.svd(np.asarray(_matrix("H4")), full_matrices=False, lapack_driver=selected.svd_driver)
    threshold = selected.rank_absolute_tolerance + selected.rank_relative_tolerance * spectrum[0]
    rank_ok = int(np.count_nonzero(spectrum > threshold)) == 3
    weak_ok = spectrum[-1] / spectrum[0] <= selected.weak_relative_tolerance
    cluster_ok = (spectrum[0] - spectrum[1]) / spectrum[0] <= selected.singular_value_cluster_relative_tolerance
    conditioning_ok = spectrum[0] / spectrum[-1] > selected.conditioning_limit
    r = np.asarray(((0.6, -0.8, 0.0), (0.8, 0.6, 0.0), (0.0, 0.0, 1.0)))
    projector_error = float(np.linalg.norm(vh[:2].T @ vh[:2] - r[:, :2] @ r[:, :2].T, ord=2))
    passed = rank_ok and weak_ok and cluster_ok and conditioning_ok and projector_error <= TOL["holdout_projector"]
    return {"case": "H4", "passed": passed, "projector_error": projector_error}


def validate_composed_cases(selected: Any, environment: str, counts: dict[str, int]) -> tuple[bool, tuple[dict[str, object], ...]]:
    records = []
    for name, numerical, zeta in (("A1", False, None), ("A2", False, None), ("A3", False, None), ("F2", True, None), ("F2", False, 2.0), ("F2", False, 4.0)):
        _evidence, record = _derive_case(name, selected, environment, counts, "canonical", numerical, zeta)
        records.append(record)
        if not record["passed"]:
            return False, tuple(records)
    return True, tuple(records)


def run_holdouts(selected: Any, environment: str, counts: dict[str, int]) -> tuple[bool, tuple[dict[str, object], ...]]:
    records: list[dict[str, object]] = []
    for name in ("H1", "H2", "H3"):
        _evidence, record = _derive_case(name, selected, environment, counts, "canonical")
        records.append(record)
        if not record["passed"]:
            return False, tuple(records)
    record = _h4_spectral_case(selected, counts, "canonical")
    records.append(record)
    if not record["passed"]:
        return False, tuple(records)
    _evidence, record = _derive_case("H5", selected, environment, counts, "canonical")
    records.append(record)
    if not record["passed"]:
        return False, tuple(records)
    raw = -(1.0 + 48.0 * 2.0**-52)
    passed = _correlation_case(raw, _case_policy("H5", selected), counts, "canonical")
    records.append({"case": "H6", "passed": passed})
    return passed, tuple(records)


def run_compatibility_replay(selected: Any, environment: str, counts: dict[str, int]) -> dict[str, object]:
    records = []
    for name in ("A1", "A2", "F2"):
        _evidence, record = _derive_case(name, selected, environment, counts, "compatibility", numerical_partial=name == "F2")
        records.append(record)
        if not record["passed"]:
            return {"status": "UNAVAILABLE", "results": tuple(records), "retuning": False}
    h4 = _h4_spectral_case(selected, counts, "compatibility")
    records.append(h4)
    if not h4["passed"]:
        return {"status": "UNAVAILABLE", "results": tuple(records), "retuning": False}
    passed = _correlation_case(1.0 + 64.0 * 2.0**-52, _case_policy("F2", selected), counts, "compatibility")
    records.append({"case": "C5", "passed": passed})
    return {"status": "COMPATIBLE" if passed else "UNAVAILABLE", "results": tuple(records), "retuning": False}


def _lane(role: str) -> NumericalLane:
    if role == "canonical":
        return canonical_lanes()[0]
    if role == "compatibility":
        return canonical_lanes()[1]
    raise ValueError("lane role must be canonical or compatibility")


def attest_lane(role: str, image_digest: str) -> dict[str, object]:
    lane = _lane(role)
    authority = lane.attest_current_process(image_digest)
    records = {"numerical_lane": lane.to_record(), "lane_attestation": authority.to_record(), "runtime_environment": RuntimeEnvironment(lane.semantics).to_record()}
    reconstruct_lane_records(records, role)
    return records


def reconstruct_lane_records(records: Mapping[str, object], role: str) -> tuple[Any, Any, Any]:
    if set(records) != {"numerical_lane", "lane_attestation", "runtime_environment"} or any(not isinstance(value, Mapping) for value in records.values()):
        raise RuntimeError("complete typed lane records required")
    lane = NumericalLane.from_record(cast("Mapping[str, object]", records["numerical_lane"]))
    attestation = LaneAttestation.from_record(cast("Mapping[str, object]", records["lane_attestation"]))
    environment = RuntimeEnvironment.from_record(cast("Mapping[str, object]", records["runtime_environment"]))
    if lane != _lane(role) or attestation.lane_identity != lane.identity or attestation.environment_identity != environment.identity or environment.semantics != lane.semantics:
        raise RuntimeError("lane records are inconsistent")
    return lane, attestation, environment


def validate_lane_records(records: Mapping[str, object], role: str) -> str:
    lane, attestation, environment = reconstruct_lane_records(records, role)
    if role == "canonical":
        selection = migration_core_authority_selection()
        if lane.identity != selection.lane_identity or attestation.identity != selection.attestation_identity or environment.identity != selection.environment_identity:
            raise RuntimeError("canonical lane does not match live #588 authority")
    return environment.identity


def _compact_result(status: str, counts: Mapping[str, int], **items: object) -> dict[str, object]:
    return {"specification_id": SPECIFICATION_ID, "specification_digest": frozen_digest(), "status": status, "counts": dict(counts), **items}


def acquire_canonical(environment: str) -> dict[str, object]:
    counts = dict.fromkeys(CEILINGS["canonical"], 0)
    scales, scale_record, trajectories, scalar = select_scales(counts)
    if scales is None:
        return _compact_result(cast("str", scale_record["status"]), counts, scale_selection=scale_record)
    fd, fd_record = calibrate_finite_differences(scales, trajectories, scalar, counts)
    if fd is None:
        return _compact_result(cast("str", fd_record["status"]), counts, scale_selection=scale_record, finite_difference=fd_record)
    driver, driver_record = select_svd_driver(counts)
    if driver is None:
        return _compact_result("UNSUPPORTED", counts, finite_difference=fd_record, svd_driver=driver_record)
    rank, rank_record = select_rank_policy(driver_record, counts)
    if rank is None:
        return _compact_result(cast("str", rank_record["status"]), counts, svd_driver=driver_record, rank=rank_record)
    weak, weak_record = select_weak_policy()
    cluster, cluster_record = select_cluster_policy()
    conditioning, conditioning_record = select_conditioning_policy()
    correlation, correlation_record = select_correlation_policy(counts)
    if None in {weak, cluster, conditioning, correlation}:
        return _compact_result("UNSUPPORTED_OR_SATURATED", counts, weak=weak_record, cluster=cluster_record, conditioning=conditioning_record, correlation=correlation_record)
    selected = compose_uncertainty_policy("A1", (scales["state_population"], scales["exchange_rate"]), fd, driver, rank, cast("float", weak), cast("float", cluster), cast("float", conditioning), cast("float", correlation))
    valid, composed = validate_composed_cases(selected, environment, counts)
    if not valid:
        return _compact_result("COMPOSED_VALIDATION_FAILED", counts, policy=policy_record(selected), composed=composed)
    passed, holdouts = run_holdouts(selected, environment, counts)
    phases = {"scales": scales, "finite_difference": fd, "svd_driver": driver, "rank": rank, "weak": weak, "cluster": cluster, "conditioning": conditioning, "correlation": correlation}
    metrics = {"scale": scale_record, "finite_difference": fd_record, "driver": driver_record, "rank": rank_record, "weak": weak_record, "cluster": cluster_record, "conditioning": conditioning_record, "correlation": correlation_record}
    return _compact_result("QUALIFIED" if passed else "HOLDOUT_FAILED_POLICY_UNAVAILABLE", counts, selected_phase_policies=phases, policy=policy_record(selected), decisive_metrics=metrics, composed=composed, holdouts=holdouts)


def acquire(role: str, image_digest: str, selected_policy: Mapping[str, object] | None = None) -> dict[str, object]:
    if role == "compatibility" and selected_policy is None:
        raise RuntimeError("compatibility replay requires selected canonical policy")
    lane_records = attest_lane(role, image_digest)
    _import_machinery()
    environment = validate_lane_records(lane_records, role)
    if role == "canonical":
        result = acquire_canonical(environment)
    else:
        counts = dict.fromkeys(CEILINGS["compatibility"], 0)
        selected = policy_from_record(cast("Mapping[str, object]", selected_policy))
        compatibility = run_compatibility_replay(selected, environment, counts)
        result = _compact_result(cast("str", compatibility["status"]), counts, policy_identity=selected.identity, compatibility=compatibility)
    return {**result, "canonical_provenance": lane_records}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--role", choices=("canonical", "compatibility"), required=True)
    parser.add_argument("--image-digest", required=True)
    parser.add_argument("--selected-policy", type=Path)
    arguments = parser.parse_args()
    selected = None if arguments.selected_policy is None else cast("Mapping[str, object]", json.loads(arguments.selected_policy.read_text()))
    arguments.output.write_text(json.dumps(acquire(arguments.role, arguments.image_digest, selected), allow_nan=False, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

"""Direct prospective acquisition runner for issue #601."""

# ruff: noqa: C901, E701, E702, I001 -- compact qualification artifact; imports await attestation.
from __future__ import annotations

# fmt: off
import argparse
import hashlib
import json
import math
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import replace
from itertools import product
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, cast
from chemex.numerical_lanes import LaneAttestation, NumericalLane, RuntimeEnvironment, canonical_lanes

# Frozen #581 inputs.  Keep literals compact and single-sourced.
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
BC_SCALE_FAMILY = "equilibrium_ratio"; COMPOSED_CASES = ("A1", "A2", "F1-absolute", "F1-estimated", "F2", "F2-boundary-2", "F2-boundary-4"); HOLDOUT_CASES = ("H1", "H2", "H3", "H4", "H5", "H6"); COMPATIBILITY_CASES = ("A1", "A2", "F2", "H4", "C5")
TOL = {"derivative": 2.0**-24, "covariance": 2.0e-6, "holdout_covariance": 5.0e-6, "projector": 2.0e-12, "holdout_projector": 5.0e-12, "conditioning": 2.0e-13, "holdout_conditioning": 5.0e-13, "correlation_interior_eps": 8.0, "boundary_eps": 128.0, "boundary_zeta": 3.0, "partial_eps": 64.0, "truth_uncertainty": 1.0 / 16.0}
CEILINGS = {"canonical": {"evaluation": 28_000, "scalar": 4_000, "svd": 32, "correlation": 24, "offline": 12_000_000}, "compatibility": {"evaluation": 751, "scalar": 266, "svd": 4, "correlation": 4, "offline": 0}}
PLANNED_BY_ROLE = {"canonical": {"evaluation": 25_609, "scalar": 3_320, "svd": 30, "correlation": 23}, "compatibility": {"evaluation": 751, "scalar": 266, "svd": 4, "correlation": 3}}
PLANNED_MAXIMUM = {"evaluation": 26_360, "scalar": 3_586, "svd": 34, "correlation": 26, "offline": 12_000_000}
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
UNSUPPORTED = frozenset(["c_dimer", "c_l", "c_monomer", "c_p", "c_pl", "c_pl1", "c_pl2", "c_pl3", "c_tetramer", "c_trimer", "dimer", "l1_free", "l2_free", "l_free", "monomer", "p_free", "pfree", "pl", "pl1", "pl2", "pl3", "kd", "kd1", "kd2", "kd_ab", "kd_ac", "kd_app", "kd_bc", "kd_eff", "kon", "kon1", "kon2", "kon_ab", "kon_ac", "kon_bc", "dh_ab", "dh_ac", "dh_ad", "dh_b", "dh_bc", "dh_bd", "dh_c", "dh_cd", "dh_d", "ds_ab", "ds_ac", "ds_ad", "ds_b", "ds_bc", "ds_bd", "ds_c", "ds_cd", "ds_d", "tauc", "dwm"])
UNSUPPORTED_PREFIXES = ("tauc_", "dwm_a", "dwm_i_a", "dwm_s_a")
# fmt: on
# fmt: off
# The remaining direct runner keeps long scientific tuples and typed constructor
# calls on one line so the qualification artifact stays reviewably small.
_MACHINERY_IMPORTED = False
np: Any = None
migration_core_authority_selection: Any = None
uncertainty_module: Any = None
if TYPE_CHECKING:
    from chemex.containers.data import Data; from chemex.containers.profile import Profile
    from chemex.evaluation.native import EvaluationEngine, EvaluationFailure, EvaluationFrame, EvaluationResult; from chemex.optimize.direct_trf import AcceptedFitResult, AffineHalfSpace, OptimizationProblem
    from chemex.optimize.uncertainty import ClaimState, FunctionAnalyticPartialCapability, FunctionFiniteDifferenceCapability, OperationTerminal, ParameterUnit, ResidualVarianceScaling, UncertaintyPolicy, _column_orientations, _correlations, _stencil_estimates, _stencil_vectors, compile_constraint_linearization_capabilities, derive_uncertainty_evidence
    from chemex.parameters.parameterization import ActiveParameterization, BinaryExpression, CompiledConstraint, ConstraintProgram, FunctionExpression, LiteralExpression, ParameterRole, ReferenceExpression, ScientificFunctionBinder; from chemex.parameters.spin_system import SpinSystem; from chemex.parameters.values import AnalysisValuesSnapshot; from chemex.printers.data import Printer
def frozen_digest() -> str:
    record = (SCALE_EXPONENTS, TAU_REL, KAPPA, EXTENTS, SVD_DRIVERS, RANK_GRID, WEAK_GRID, CONDITION_GRID, CORRELATION_GRID, SUPPORTED, tuple(sorted(UNSUPPORTED)), UNSUPPORTED_PREFIXES, BC_SCALE_FAMILY, COMPOSED_CASES, HOLDOUT_CASES, COMPATIBILITY_CASES, TOL, CEILINGS, PLANNED_BY_ROLE, PLANNED_MAXIMUM)
    return hashlib.sha256(repr(record).encode()).hexdigest()
def finite_difference_policies() -> Iterator[tuple[float, float, int, int]]:
    return product(TAU_REL, KAPPA, EXTENTS, EXTENTS)
def _scale_candidate_identity(family: str, q0: float, exponent: int, scale: float) -> str:
    return hashlib.sha256(repr((family, float(q0).hex(), exponent, float(scale).hex())).encode()).hexdigest()
def _scale_catalogue_digest(scales: Mapping[str, float]) -> str:
    return hashlib.sha256(repr(tuple((name, float(value).hex()) for name, value in sorted(scales.items()))).encode()).hexdigest()
def classify_scale_name(name: str) -> tuple[str, str, float, str]:
    for family, stratum, q0, unit, patterns in SUPPORTED:
        if name in patterns or any(item.endswith("*") and name.startswith(item[:-1]) for item in patterns):
            return family, stratum, q0, unit
    if name in UNSUPPORTED or name.startswith(UNSUPPORTED_PREFIXES):
        raise KeyError(f"catalogued but unsupported after calibration: {name}")
    raise KeyError(f"not catalogued: {name}")
def actual_catalogue_names() -> frozenset[str]:
    """Generate the current complete non-WIP CLI-reachable local-name catalogue."""
    from chemex.configuration.conditions import Conditions
    from chemex.models.factory import model_factory
    from chemex.models.loader import register_kinetic_settings
    from chemex.models.model import ModelSpec
    from chemex.nmr.basis import Basis
    from chemex.parameters.spins import create_base_param_settings
    register_kinetic_settings()
    conditions = Conditions(h_larmor_frq=600.0, temperature=25.0, p_total=1.0e-3, l_total=2.0e-3, d2o=0.1)
    names: set[str] = set()
    for model_name in model_factory.set:
        names.update(model_factory.create(model_name, conditions))
    model = ModelSpec(name="6st", states="abcdef", temp_coef=True)
    basis = Basis(type="ixyzsxyz", spin_system="nh", model=model)
    for state in model.states:
        names.update(create_base_param_settings(basis, state, conditions))
    return frozenset(names)
def catalogue_partition() -> tuple[frozenset[str], frozenset[str]]:
    supported, unsupported = set(), set()
    for name in actual_catalogue_names():
        try:
            classify_scale_name(name)
        except KeyError as error:
            if str(error).strip("'").startswith("catalogued but unsupported"):
                unsupported.add(name)
            else:
                raise
        else:
            supported.add(name)
    return frozenset(supported), frozenset(unsupported)
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
_ACTIVE_COUNTS: dict[str, int] | None = None
_ACTIVE_ROLE = "canonical"
def _charge_scalar() -> None:
    if _ACTIVE_COUNTS is not None:
        _spend(_ACTIVE_COUNTS, _ACTIVE_ROLE, "scalar")
def counted_scientific_value(a: float, b: float) -> float:
    _charge_scalar(); return scientific_value(a, b)
def counted_holdout_value(a: float, b: float) -> float:
    _charge_scalar(); return holdout_value(a, b)
def counted_scientific_partial_a(a: float, b: float) -> float:
    _charge_scalar(); return scientific_partial_a(a, b)
def counted_scientific_partial_b(a: float, b: float) -> float:
    _charge_scalar(); return scientific_partial_b(a, b)
def counted_holdout_partial_a(a: float, b: float) -> float:
    _charge_scalar(); return holdout_partial_a(a, b)
def counted_holdout_partial_b(a: float, b: float) -> float:
    _charge_scalar(); return holdout_partial_b(a, b)
def _matrix(name: str) -> tuple[tuple[float, ...], ...]:
    if name == "B1":
        return ((1.0e-9, 0.0), (0.0, 2.0e-21))
    if name == "B2":
        return ((1.0, 0.0), (0.0, 5.0e-13))
    if name == "B3":
        return ((1.0, 1.0), (2.0, 2.0), (3.0, 3.0))
    if name == "C1":
        return ((1.0, 0.0, 0.0), (0.0, 2.0e-6, 0.0), (0.0, 0.0, 2.0e-7))
    if name == "C2":
        return tuple(tuple(float(value) for value in row) for row in np.diag((4.0, 4.0 * (1.0 - 5.0e-11), 2.0, 2.0 * (1.0 - 5.0e-9), 1.0)))
    if name == "C3":
        return ((1.0, 0.0), (0.0, 2.0e-6))
    if name == "C4":
        return ((1.0, 0.0), (0.0, 5.0e-7))
    if name in {"C5", "H6"}:
        return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    if name == "F2":
        return ((1.0, 0.0), (0.0, 2.0), (1.0, 1.0), (2.0, -1.0))
    if name == "H5":
        return ((1.0, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, -1.0))
    if name == "H4":
        r = np.asarray(((0.6, -0.8, 0.0), (0.8, 0.6, 0.0), (0.0, 0.0, 1.0)))
        s = np.asarray((3.0, 3.0 * (1.0 - 7.5e-11), 3.0e-7))
        return tuple(tuple(float(value) for value in row) for row in np.diag(s) @ r.T)
    raise KeyError(name)
def _setup(name: str, q0: float = 1.0) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[tuple[float, ...], ...]]:
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
    if name in {"B1", "B2", "B3", "C1", "C2", "C3", "C4", "C5", "H6"}:
        size = len(_matrix(name)[0])
        return (0.0,) * size, (1.0,) * size, (-10.0,) * size, (10.0,) * size, ()
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
    global ActiveParameterization, AffineHalfSpace, AnalysisValuesSnapshot, AcceptedFitResult, BinaryExpression, ClaimState, CompiledConstraint, ConstraintProgram, Data, EvaluationEngine, EvaluationFailure, EvaluationFrame, EvaluationResult, FunctionAnalyticPartialCapability, FunctionExpression, FunctionFiniteDifferenceCapability, LiteralExpression, OperationTerminal, OptimizationProblem, ParameterRole, ParameterUnit, Printer, Profile, ReferenceExpression, ResidualVarianceScaling, ScientificFunctionBinder, SpinSystem, UncertaintyPolicy, compile_constraint_linearization_capabilities, derive_uncertainty_evidence, migration_core_authority_selection, np, uncertainty_module, _column_orientations, _correlations, _stencil_estimates, _stencil_vectors, _MACHINERY_IMPORTED
    if _MACHINERY_IMPORTED:
        return
    # isort: off -- imports must remain behind live lane attestation.
    import numpy as np; from pydantic import BaseModel
    from chemex.containers.data import Data; from chemex.containers.profile import Profile
    from chemex.evaluation.native import EvaluationEngine, EvaluationFailure, EvaluationFrame, EvaluationResult
    from chemex.migration_core import migration_core_authority_selection; from chemex.optimize import uncertainty as uncertainty_module; from chemex.optimize.direct_trf import AcceptedFitResult, AffineHalfSpace, OptimizationProblem
    from chemex.optimize.uncertainty import ClaimState, FunctionAnalyticPartialCapability, FunctionFiniteDifferenceCapability, OperationTerminal, ParameterUnit, ResidualVarianceScaling, UncertaintyPolicy, _column_orientations, _correlations, _stencil_estimates, _stencil_vectors, compile_constraint_linearization_capabilities, derive_uncertainty_evidence
    from chemex.parameters.parameterization import ActiveParameterization, BinaryExpression, CompiledConstraint, ConstraintProgram, FunctionExpression, LiteralExpression, ParameterRole, ReferenceExpression, ScientificFunctionBinder; from chemex.parameters.spin_system import SpinSystem; from chemex.parameters.values import AnalysisValuesSnapshot; from chemex.printers.data import Printer
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
def _linear_expression(coefficients: tuple[float, ...], params: tuple[str, ...]) -> Any:
    terms = tuple(BinaryExpression("multiply", LiteralExpression(value), ReferenceExpression(param)) for value, param in zip(coefficients, params, strict=True))
    expression = terms[0]
    for term in terms[1:]:
        expression = BinaryExpression("add", expression, term)
    return expression
def _derived_program(name: str, params: tuple[str, ...], counts: dict[str, int], role: str) -> tuple[Any, tuple[str, ...], tuple[Any, ...]]:
    if name not in {"F2", "H5", "C5", "H6"}:
        binder = ScientificFunctionBinder(name, {})
        return binder, params, ()
    if name in {"C5", "H6"}:
        x64 = (-7.538105255798958e-83, 5.347758357992345e-73, -4.026953966760176e-85)
        xbad = (-4.752313720786381e-76, 1.0091182916234204e-87, 3.426969796398343e-71)
        x48 = x64
        rows = ((x48, tuple(-2.1085807432963563e-11 * value for value in x48)) if name == "H6" else ((1.0, 0.0, 0.0), (0.75, math.sqrt(7.0) / 4.0, 0.0), x64, tuple(2.1085807432963563e-11 * value for value in x64), x64, tuple(-2.1085807432963563e-11 * value for value in x64), xbad, tuple(2.392941441612749e-15 * value for value in xbad)))
        constraints = tuple(CompiledConstraint(f"R{index}", _linear_expression(row, params), params, "issue601", f"C5-row-{index}") for index, row in enumerate(rows))
        binder = ScientificFunctionBinder(name, {})
        return binder, tuple(item.target_id for item in constraints), constraints
    function = counted_scientific_value if name == "F2" else counted_holdout_value
    binder = ScientificFunctionBinder(name, {"f": function})
    a, b = ReferenceExpression("A"), ReferenceExpression("B")
    constraints = (
        CompiledConstraint("SUM", BinaryExpression("add", a, BinaryExpression("multiply", LiteralExpression(2.0), b)), ("A", "B"), "issue601", "A+2*B"),
        CompiledConstraint("SCI", FunctionExpression("f", (a, b)), ("A", "B"), "issue601", "f(A,B)"),
        CompiledConstraint("ZERO", BinaryExpression("add", LiteralExpression(1.0), BinaryExpression("multiply", LiteralExpression(0.0), a)), ("A",), "issue601", "1+0*A"),
    )
    return binder, (*params, "SUM", "SCI", "ZERO"), constraints
def _fixture(name: str, counts: dict[str, int], role: str, q0: float = 1.0) -> tuple[Any, Any, Any, Any, tuple[str, ...]]:
    anchor, _scales, lower, upper, halfspaces = _setup(name, q0)
    params = tuple(chr(65 + i) for i in range(len(anchor)))
    binder, requested, constraints = _derived_program(name, params, counts, role)
    derived = tuple(item.target_id for item in constraints)
    scope = (*params, *derived)
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
    _spend(counts, role, "evaluation")
    evaluation = engine.new_evaluator().evaluate(EvaluationFrame.from_lifecycle_frame(parameterization, problem.lifecycle_frame(anchor, parameterization)))
    assert isinstance(evaluation, EvaluationResult)
    accepted = AcceptedFitResult.for_qualification(occurrence_identity=f"issue-601-{name}-accepted", problem_identity=problem.identity, invocation_identity="issue-601", execution_identity="issue-601", materialization_identity="issue-601", parameterization_identity=parameterization.identity, evaluator_parameterization_identity=parameterization.evaluator_identity, source_occurrence_identity=snapshot.occurrence_identity, source_revision=snapshot.revision, controlled_ids=params, vector=anchor, chi_square=float(np.dot(evaluation.residuals, evaluation.residuals)), evaluation_result=evaluation, commit_scope=scope, commit_items=tuple(zip(params, anchor, strict=True)), origin_context_identity="issue-601")
    return parameterization, engine, problem, accepted, requested
def _spend(counts: dict[str, int], role: str, kind: str, amount: int = 1) -> None:
    if counts[kind] + amount > CEILINGS[role][kind]:
        raise RuntimeError(f"{role} {kind} ceiling exceeded")
    counts[kind] += amount
class _CountingEvaluator:
    def __init__(self, evaluator: Any, owner: Any) -> None:
        self._evaluator, self._owner = evaluator, owner
    def evaluate(self, frame: Any) -> Any:
        if self._owner.charging:
            _spend(self._owner.counts, self._owner.role, "evaluation")
        return self._evaluator.evaluate(frame)
class _CountingEngine:
    def __init__(self, engine: Any, counts: dict[str, int], role: str) -> None:
        self._engine, self.plan, self.counts, self.role, self.charging = engine, engine.plan, counts, role, True
    def new_evaluator(self) -> Any:
        return _CountingEvaluator(self._engine.new_evaluator(), self)
@contextmanager
def _charged_svd(counts: dict[str, int], role: str, verification: Any = None) -> Iterator[None]:
    original = uncertainty_module.svd
    def wrapper(*args: object, **kwargs: object) -> Any:
        if verification is None or not verification():
            _spend(counts, role, "svd")
        return original(*args, **kwargs)
    uncertainty_module.svd = wrapper
    try:
        yield
    finally:
        uncertainty_module.svd = original
@contextmanager
def _charged_derivation(counts: dict[str, int], role: str, stop_after: str | None, engine: Any) -> Iterator[Any]:
    original_correlation = uncertainty_module._correlations
    original_subspaces = uncertainty_module._invariant_singular_subspaces
    original_validation = uncertainty_module.ResidualJacobianEvidence.__post_init__; original_constraint_validation = uncertainty_module.ConstraintJacobianEvidence.__post_init__
    original_marginal_validation = uncertainty_module.MarginalErrorEvidence.__post_init__; svd_classes = (uncertainty_module.RankDiagnostic, uncertainty_module.CovarianceEvidence, uncertainty_module.ConstrainedPropagationEvidence); original_svd_validations = tuple(item.__post_init__ for item in svd_classes); residual_complete = rank_complete = constraint_complete = constrained_marginal_complete = verifying = False
    def correlation(*args: object, **kwargs: object) -> Any:
        _spend(counts, role, "correlation"); return original_correlation(*args, **kwargs)
    def subspaces(*args: object, **kwargs: object) -> Any:
        nonlocal rank_complete
        result = original_subspaces(*args, **kwargs); rank_complete = True; return result
    def validate(artifact: Any) -> None:
        nonlocal residual_complete
        engine.charging = False
        try:
            original_validation(artifact)
        finally:
            engine.charging = True
        residual_complete = True
    def validate_constraint(artifact: Any) -> None:
        nonlocal constraint_complete; original_constraint_validation(artifact); constraint_complete = True
    def validate_marginal(artifact: Any) -> None:
        nonlocal constrained_marginal_complete; original_marginal_validation(artifact); constrained_marginal_complete |= artifact.source_family == "constrained_propagation"
    def validate_svd(index: int, artifact: Any) -> None:
        nonlocal verifying; verifying = True
        try: original_svd_validations[index](artifact)
        finally: verifying = False
    uncertainty_module._correlations = correlation
    uncertainty_module.ResidualJacobianEvidence.__post_init__ = validate; uncertainty_module.ConstraintJacobianEvidence.__post_init__ = validate_constraint; uncertainty_module.MarginalErrorEvidence.__post_init__ = validate_marginal
    for index, item in enumerate(svd_classes): item.__post_init__ = lambda artifact, index=index: validate_svd(index, artifact)
    if stop_after == "rank":
        uncertainty_module._invariant_singular_subspaces = subspaces
    try:
        with _charged_svd(counts, role, lambda: verifying):
            yield lambda: OperationTerminal.CANCELLED if (stop_after == "residual" and residual_complete) or (stop_after == "rank" and rank_complete) or (stop_after == "constraint" and constraint_complete) or (stop_after == "constrained_marginal" and constrained_marginal_complete) else None
    finally:
        uncertainty_module._correlations = original_correlation
        uncertainty_module._invariant_singular_subspaces = original_subspaces
        uncertainty_module.ResidualJacobianEvidence.__post_init__ = original_validation; uncertainty_module.ConstraintJacobianEvidence.__post_init__ = original_constraint_validation; uncertainty_module.MarginalErrorEvidence.__post_init__ = original_marginal_validation
        for item, validation in zip(svd_classes, original_svd_validations, strict=True): item.__post_init__ = validation
def _accepted_copy(accepted: Any, problem_identity: str, occurrence_identity: str) -> Any:
    return AcceptedFitResult.for_qualification(occurrence_identity=occurrence_identity, problem_identity=problem_identity, invocation_identity=accepted.invocation_identity, execution_identity=accepted.execution_identity, materialization_identity=accepted.materialization_identity, parameterization_identity=accepted.parameterization_identity, evaluator_parameterization_identity=accepted.evaluator_parameterization_identity, source_occurrence_identity=accepted.source_occurrence_identity, source_revision=accepted.source_revision, controlled_ids=accepted.controlled_ids, vector=accepted.vector, chi_square=accepted.chi_square, evaluation_result=accepted.evaluation_result, commit_scope=accepted.commit_scope, commit_items=accepted.commit_items, origin_context_identity=accepted.origin_context_identity)
def _evaluate(context: tuple[Any, ...], evaluator: Any, vector: tuple[float, ...], counts: dict[str, int], role: str) -> tuple[float, ...]:
    parameterization, _engine, problem, _accepted, _requested = context
    _spend(counts, role, "evaluation")
    result = evaluator.evaluate(EvaluationFrame.from_lifecycle_frame(parameterization, problem.lifecycle_frame(vector, parameterization)))
    if isinstance(result, EvaluationFailure):
        raise RuntimeError(f"truth evaluation failed: {result.category}")  # noqa: TRY004
    return tuple(float(value) for value in result.residuals)
def _union_trajectories(
    name: str, q0: float, counts: dict[str, int], role: str
) -> tuple[tuple[object, ...], ...]:
    context = _fixture(name, counts, role, q0)
    _parameterization, engine, problem, accepted, _requested = context
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
    _spend(counts, role, "scalar")
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
def _qualifier_failure(rows: Mapping[tuple[object, ...], tuple[bool, tuple[object, ...] | None, str | None]], axes: tuple[tuple[object, ...], ...]) -> str | None:
    identities = tuple(identity for passed, _key, identity in rows.values() if passed and identity is not None)
    if len(identities) != len(set(identities)):
        return "AMBIGUOUS_SELECTION"
    for axis, grid in enumerate(axes):
        groups: dict[tuple[object, ...], list[tuple[int, bool]]] = {}
        for candidate, (passed, _key, _identity) in rows.items():
            groups.setdefault(candidate[:axis] + candidate[axis + 1 :], []).append((grid.index(candidate[axis]), passed))
        for group in groups.values():
            flags = tuple(passed for _index, passed in sorted(group))
            if any(not flags[index] for index in range(flags.index(True), len(flags) - tuple(reversed(flags)).index(True))) if any(flags) else False:
                return "NON_MONOTONE_ADEQUACY"
    return None
def _select_one_scale(label: str, q0: float, trajectories: tuple[tuple[object, ...], ...], outputs: tuple[float, ...], counts: dict[str, int]) -> tuple[float | None, dict[str, object]]:
    candidates, outcomes = [], {}
    for exponent in SCALE_EXPONENTS:
        q = math.ldexp(q0, exponent)
        summaries = []
        for policy in finite_difference_policies():
            score = _score(trajectories, (q,) * len(trajectories), outputs, policy, counts)
            if score is not None:
                summaries.append((score[3], score[1], score[2], policy))
        if summaries:
            best = min(summaries)
            candidates.append((best[0], abs(exponent), exponent, q, best))
            outcomes[(exponent,)] = (True, (best[0], abs(exponent), exponent), hashlib.sha256(repr(("scale", label, float(q).hex())).encode()).hexdigest())
        else:
            outcomes[(exponent,)] = (False, None, None)
    if failure := _qualifier_failure(outcomes, (SCALE_EXPONENTS,)):
        return None, {"status": failure, "family": label}
    if not candidates:
        return None, {"status": "UNSUPPORTED", "family": label}
    if all(item[2] in {-8, 8} for item in candidates):
        return None, {"status": "GRID_SATURATED", "family": label}
    selected = min(candidates, key=lambda item: item[:3])
    by_exponent = {item[2]: item for item in candidates}
    neighbors = []
    for exponent in (selected[2] - 1, selected[2] + 1):
        item = by_exponent.get(exponent)
        neighbors.append({"exponent": exponent, "status": "QUALIFIED_REJECTED" if item else "DISQUALIFIED", "cost": None if item is None else item[0], "truth_error": None if item is None else item[4][1], "reasons": ("cost_abs_k_k",) if item else ("truth_or_reliability",)})
    return selected[3], {"status": "SELECTED", "q0": q0, "exponent": selected[2], "scale": selected[3], "candidate_identity": _scale_candidate_identity(label, q0, selected[2], selected[3]), "cost": selected[0], "truth_error": selected[4][1], "rejected_neighbors": tuple(neighbors)}
def select_scales(counts: dict[str, int]) -> tuple[dict[str, float] | None, dict[str, object], tuple[tuple[object, ...], ...], tuple[tuple[object, ...], ...]]:
    selected, summaries, retained = {}, {}, []
    trajectories_by_q0 = {q0: _union_trajectories("S_CENTER", q0, counts, "canonical") + _union_trajectories("S_BOUND", q0, counts, "canonical") for q0 in {item[2] for item in SUPPORTED}}
    selections = {}
    for q0, pair in trajectories_by_q0.items():
        value, record = _select_one_scale(f"q0={q0:g}", q0, pair, (1.0, 1.0), counts)
        if value is None:
            return None, record, (), ()
        selections[q0] = value, record
    for family, _stratum, q0, _unit, _patterns in SUPPORTED:
        selected[family] = selections[q0][0]; summaries[family] = {**selections[q0][1], "family": family, "candidate_identity": _scale_candidate_identity(family, q0, cast("int", selections[q0][1]["exponent"]), selected[family])}
        retained.extend(trajectories_by_q0[q0])
    scalar = _scalar_trajectories(counts, "canonical")
    function_scales = []
    for column, q0 in enumerate((1.0, 2.0)):
        value, record = _select_one_scale(f"A4-argument-{column}", q0, (scalar[column],), (4.0,), counts)
        if value is None:
            return None, record, (), ()
        function_scales.append(value)
        summaries[f"A4-argument-{column}"] = record
    # Output scale uses the same A0 rule, with Q replacing q in the truth envelope.
    output_candidates = []
    output_outcomes: dict[tuple[object, ...], tuple[bool, tuple[object, ...] | None, str | None]] = {}
    for exponent in SCALE_EXPONENTS:
        output = math.ldexp(4.0, exponent)
        best = min(((score[3], score[1], score[2], policy) for policy in finite_difference_policies() if (score := _score(scalar, tuple(function_scales), (output, output), policy, counts)) is not None), default=None)
        if best is not None:
            output_candidates.append((best[0], abs(exponent), exponent, output, best))
            output_outcomes[(exponent,)] = (True, (best[0], abs(exponent), exponent), hashlib.sha256(repr(("scale", "A4-output", float(output).hex())).encode()).hexdigest())
        else:
            output_outcomes[(exponent,)] = (False, None, None)
    if failure := _qualifier_failure(output_outcomes, (SCALE_EXPONENTS,)):
        return None, {"status": failure, "family": "A4-output"}, (), ()
    if not output_candidates:
        return None, {"status": "UNSUPPORTED", "family": "A4-output"}, (), ()
    if all(item[2] in {-8, 8} for item in output_candidates):
        return None, {"status": "GRID_SATURATED", "family": "A4-output"}, (), ()
    output_choice = min(output_candidates, key=lambda item: item[:3])
    selected["A4-argument-0"], selected["A4-argument-1"], selected["A4-output"] = *function_scales, output_choice[3]
    output_by_k = {item[2]: item for item in output_candidates}
    summaries["A4-output"] = {"status": "SELECTED", "Q0": 4.0, "exponent": output_choice[2], "scale": output_choice[3], "candidate_identity": _scale_candidate_identity("A4-output", 4.0, output_choice[2], output_choice[3]), "cost": output_choice[0], "truth_error": output_choice[4][1], "rejected_neighbors": tuple({"exponent": exponent, "status": "QUALIFIED_REJECTED" if exponent in output_by_k else "DISQUALIFIED", "cost": None if exponent not in output_by_k else output_by_k[exponent][0], "truth_error": None if exponent not in output_by_k else output_by_k[exponent][4][1], "reasons": ("cost_abs_k_k",) if exponent in output_by_k else ("truth_or_reliability",)} for exponent in (output_choice[2] - 1, output_choice[2] + 1))}
    return selected, {"status": "SELECTED", "families": summaries}, tuple(retained), scalar
def calibrate_finite_differences(scales: Mapping[str, float], family_trajectories: tuple[tuple[object, ...], ...], scalar: tuple[tuple[object, ...], ...], counts: dict[str, int]) -> tuple[tuple[float, float, int, int] | None, dict[str, object]]:
    normalization = _union_trajectories("A3", 1.0, counts, "canonical")
    trajectories = family_trajectories + normalization + scalar
    family_values = tuple(scales[family] for family, *_rest in SUPPORTED)
    q = tuple(value for value in family_values for _ in (0, 1)) + (1.0, 1.0, scales["A4-argument-0"], scales["A4-argument-1"])
    outputs = (1.0,) * (len(trajectories) - 2) + (scales["A4-output"],) * 2
    best = None
    frontier, outcomes = [], {}
    non_edge = False
    for policy in finite_difference_policies():
        score = _score(trajectories, q, outputs, policy, counts)
        if score is None:
            outcomes[policy] = (False, None, None)
            continue
        edge = policy == (2.0**-8, 2.0**12, 16, 16)
        non_edge |= not edge
        key = (*score, policy[0], policy[1], policy[3], policy[2], policy)
        outcomes[policy] = (True, key, hashlib.sha256(repr(("finite_difference", policy)).encode()).hexdigest())
        if best is None or key < best:
            best = key
        point = (score[1], score[3], policy)
        frontier = [item for item in frontier if not (point[0] <= item[0] and point[1] <= item[1])]
        if not any(item[0] <= point[0] and item[1] <= point[1] for item in frontier):
            frontier.append(point)
    if failure := _qualifier_failure(outcomes, (TAU_REL, KAPPA, EXTENTS, EXTENTS)):
        return None, {"status": failure}
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
                score = _score(trajectories, q, outputs, candidate_tuple, counts)
                neighbors.append({"candidate": candidate_tuple, "status": "QUALIFIED_REJECTED" if score else "DISQUALIFIED", "metrics": score, "reasons": ("pareto_or_tie_break",) if score else ("truth_or_reliability",)})
    return policy, {"status": "SELECTED", "policy": policy, "metrics": best[:4], "frontier": tuple(sorted(frontier)[:32]), "rejected_neighbors": tuple(neighbors)}
def _bc_scales(name: str, scales: Mapping[str, float]) -> tuple[float, ...]:
    return (scales[BC_SCALE_FAMILY],) * len(_setup(name)[1])
def _spectral_observation(name: str, scales: Mapping[str, float], phase_a: Any, driver: str, rank: tuple[float, float], environment: str, counts: dict[str, int]) -> tuple[Any | None, dict[str, object]]:
    parameterization, engine, problem, accepted, _requested = _fixture(name, counts, "canonical")
    policy = compose_uncertainty_policy(name, _bc_scales(name, scales), (phase_a.relative_step_tolerance, phase_a.roundoff_multiplier, phase_a.smaller_step_extent, phase_a.larger_step_extent), driver, rank, max(rank), 0.0, 2.0**52, 0.0)
    counted_engine = _CountingEngine(engine, counts, "canonical")
    with _charged_derivation(counts, "canonical", "rank", counted_engine) as cancellation:
        evidence = derive_uncertainty_evidence(accepted, problem=problem, parameterization=parameterization, engine=cast("EvaluationEngine", counted_engine), policy=policy, cancellation_probe=cancellation, resolved_environment_identity=environment)
    if evidence.residual_jacobian is None or evidence.rank_diagnostic is None:
        return None, {"case": name, "status": "DISQUALIFIED", "reasons": tuple(item.category for item in evidence.failures)}
    observed = np.asarray(evidence.residual_jacobian.matrix)
    truth = np.asarray(truth_jacobian(name))
    derivative_error = float(np.max(np.abs(observed - truth) / (TOL["derivative"] * (1.0 + np.abs(truth)))))
    diagnostic = evidence.rank_diagnostic
    return diagnostic, {"case": name, "status": "ACQUIRED", "policy_identity": policy.identity, "coordinate_scales": policy.coordinate_scales, "jacobian_identity": evidence.residual_jacobian.identity, "rank_identity": diagnostic.identity, "spectrum": diagnostic.singular_values, "rank": diagnostic.rank, "derivative_error": derivative_error}
def select_svd_driver(scales: Mapping[str, float], phase_a: Any, environment: str, counts: dict[str, int]) -> tuple[str | None, dict[str, object], tuple[Any, ...]]:
    records = []
    q = scales[BC_SCALE_FAMILY]
    expected = ((1.0e-9 * q, 2.0e-21 * q), (q, 5.0e-13 * q), (math.sqrt(28.0) * q, 0.0))
    for driver in SVD_DRIVERS:
        observations, summaries = [], []
        for name in ("B1", "B2", "B3"):
            diagnostic, summary = _spectral_observation(name, scales, phase_a, driver, (0.0, 0.0), environment, counts)
            if diagnostic is None or summary.get("status") != "ACQUIRED":
                observations = []
                summaries.append(summary)
                break
            observations.append(diagnostic)
            summaries.append(summary)
        if observations:
            case_errors = tuple(max(abs(value - truth) / max(1.0, abs(truth)) for value, truth in zip(diagnostic.singular_values, exact, strict=True)) for diagnostic, exact in zip(observations, expected, strict=True)); summaries = [{**summary, "status": "QUALIFIED" if summary["derivative_error"] <= 1.0 and error <= TOL["covariance"] else "DISQUALIFIED", "spectrum_error": error} for summary, error in zip(summaries, case_errors, strict=True)]
            if all(summary["status"] == "QUALIFIED" for summary in summaries):
                records.append((max(case_errors), 3, SVD_DRIVERS.index(driver), driver, tuple(observations), tuple(summaries)))
    if not records:
        return None, {"status": "UNSUPPORTED", "reasons": ("no_driver_qualified",)}, ()
    selected = min(records, key=lambda item: item[:4])
    return selected[3], {"status": "SELECTED", "driver": selected[3], "truth_error": selected[0], "cost": selected[1], "frontier": ((selected[3], selected[0], selected[1]),), "rejected_neighbors": tuple({"driver": item[3], "status": "QUALIFIED_REJECTED", "truth_error": item[0], "cost": item[1], "reasons": ("spectrum_error_cost_driver_order",)} for item in records if item is not selected), "cases": selected[5]}, selected[4]
def select_rank_policy(observations: tuple[Any, ...], counts: dict[str, int]) -> tuple[tuple[float, float] | None, dict[str, object]]:
    truths = (2, 1, 1)
    summaries: dict[tuple[float, float], tuple[str, tuple[str, ...], int, float]] = {}
    qualified, outcomes = [], {}
    for absolute, relative in product(RANK_GRID, RANK_GRID):
        _spend(counts, "canonical", "offline")
        ranks = tuple(sum(value > absolute + relative * item.singular_values[0] for value in item.singular_values) for item in observations)
        error = sum(abs(value - truth) for value, truth in zip(ranks, truths, strict=True))
        reasons = () if error == 0 else ("rank_misclassification",)
        status = "QUALIFIED" if not reasons else "DISQUALIFIED"
        summaries[(absolute, relative)] = status, reasons, 3, float(error)
        outcomes[(absolute, relative)] = (not reasons, (0.0, 0, absolute, relative, repr((absolute, relative))) if not reasons else None, hashlib.sha256(repr(("rank", absolute, relative)).encode()).hexdigest() if not reasons else None)
        if not reasons:
            qualified.append((absolute, relative))
    if failure := _qualifier_failure(outcomes, (RANK_GRID, RANK_GRID)):
        return None, {"status": failure}
    if not qualified:
        return None, {"status": "UNSUPPORTED"}
    if all(absolute == 2.0**-16 or relative == 2.0**-16 for absolute, relative in qualified):
        return None, {"status": "GRID_SATURATED"}
    selected = min(qualified, key=lambda item: (0.0, 0, item[0], item[1], repr(item)))
    neighbors = []
    for axis, grid in enumerate((RANK_GRID, RANK_GRID)):
        index = grid.index(selected[axis])
        for neighbor_index in (index - 1, index + 1):
            if 0 <= neighbor_index < len(grid):
                candidate = list(selected)
                candidate[axis] = grid[neighbor_index]
                key = cast("tuple[float, float]", tuple(candidate))
                status, reasons, cost, error = summaries[key]
                neighbors.append({"candidate": key, "status": status, "decisive_reasons": reasons or ("deterministic_tie_break",), "error": error, "cost": cost})
            else:
                neighbors.append({"candidate": None, "status": "GRID_EDGE", "decisive_reasons": ("no_immediate_neighbor",), "error": None, "cost": 0})
    return selected, {"status": "SELECTED", "policy": selected, "qualified_count": len(qualified), "truth_error": 0.0, "cost": 3, "frontier": ((selected, 0.0, 3),), "rejected_neighbors": tuple(neighbors)}
def _b2_perturbation_truth(diagnostic: Any, selected: Any, counts: dict[str, int]) -> tuple[bool, Any, Any, dict[str, object]]:
    scaled = np.asarray(diagnostic.source_jacobian.matrix) * np.asarray(diagnostic.source_jacobian.coordinate_scales)[None, :]; largest = diagnostic.singular_values[0]; delta = 1.0e-12 * largest
    matrices = tuple(scaled + np.asarray(((0.0, 0.0), (0.0, sign * delta))) for sign in (-1.0, 1.0)); spectra, rights = [], []
    with _charged_svd(counts, "canonical"):
        for matrix in matrices:
            _left, singular, right = uncertainty_module.svd(matrix, full_matrices=False, compute_uv=True, overwrite_a=False, check_finite=True, lapack_driver=selected.svd_driver); spectra.append(tuple(float(value) for value in singular)); rights.append(np.asarray(right))
    ranks = tuple(sum(value > diagnostic.threshold for value in spectrum) for spectrum in spectra); unstable = tuple(index for index in range(len(spectra[0])) if (spectra[0][index] > diagnostic.threshold) != (spectra[1][index] > diagnostic.threshold))
    null = np.outer(rights[0][unstable[0]], rights[0][unstable[0]]) if unstable else np.zeros((2, 2)); identifiable = np.eye(2) - null
    weak_entries = tuple(float(matrix[1, 1] / largest) for matrix in matrices); unperturbed_error = float(np.max(np.abs(scaled / largest - np.diag((1.0, 5.0e-13))))); truth_null = np.diag((0.0, 1.0)); projector_error = max(float(np.linalg.norm(null - truth_null, ord=2)), float(np.linalg.norm(np.outer(rights[1][1], rights[1][1]) - truth_null, ord=2)))
    passed = ranks == (1, 2) and unstable == (1,) and weak_entries[0] < 0.0 < weak_entries[1] and unperturbed_error <= TOL["derivative"] and projector_error <= TOL["projector"] and diagnostic.rank == 1
    return passed, null, identifiable, {"relative_envelope": 1.0e-12, "unperturbed_error": unperturbed_error, "weak_entries": weak_entries, "ranks": ranks, "unstable_modes": unstable}
def validate_rank_truth_cases(scales: Mapping[str, float], selected: Any, environment: str, counts: dict[str, int]) -> tuple[bool, tuple[dict[str, object], ...]]:
    records = []
    for name in ("B2", "B3"):
        evidence, record = _derive_case(name, selected, environment, counts, "canonical", coordinate_scales=_bc_scales(name, scales))
        diagnostic = evidence.rank_diagnostic
        perturbation_passed, perturbation_null, perturbation_identifiable, perturbation = _b2_perturbation_truth(diagnostic, selected, counts) if name == "B2" and diagnostic is not None else (True, None, None, None)
        expected_null = perturbation_null if name == "B2" else np.asarray(((0.5, -0.5), (-0.5, 0.5))); expected_identifiable = perturbation_identifiable if name == "B2" else np.eye(2) - expected_null
        projector_error = math.inf if diagnostic is None else max(float(np.linalg.norm(np.asarray(diagnostic.null_projector) - expected_null, ord=2)), float(np.linalg.norm(np.asarray(diagnostic.identifiable_projector) - expected_identifiable, ord=2)))
        typed = diagnostic is not None and diagnostic.rank == 1 and diagnostic.controlled_ids == ("A", "B") and any(item.classification.endswith("null") for item in diagnostic.subspaces)
        unavailable = evidence.covariance is None and evidence.marginal_errors is None and evidence.correlations is None and any(item.category == "rank_deficient" for item in evidence.failures)
        perturbation_error = math.inf if diagnostic is None else abs(diagnostic.singular_values[1] - (5.0e-13 * scales[BC_SCALE_FAMILY] if name == "B2" else 0.0))
        passed = bool(record["passed"]) and typed and unavailable and perturbation_passed and projector_error <= TOL["projector"] and perturbation_error <= 512.0 * 2.0**-52 * max(1.0, scales[BC_SCALE_FAMILY])
        records.append({**record, "passed": passed, "rank_identity": None if diagnostic is None else diagnostic.identity, "projector_error": projector_error, "perturbation_error": perturbation_error, "perturbation": perturbation, "covariance_available": evidence.covariance is not None})
    return all(record["passed"] for record in records), tuple(records)
def _threshold_outcome(grid: tuple[float, ...], qualified: list[float], edge: float) -> str:
    if not qualified:
        return "UNSUPPORTED"
    indices = tuple(grid.index(value) for value in qualified)
    if indices != tuple(range(indices[0], indices[-1] + 1)):
        return "AMBIGUOUS_SELECTION"
    if qualified == [edge]:
        return "GRID_SATURATED"
    return "SELECTED"
def select_weak_policy(observation: Any, counts: dict[str, int]) -> tuple[float | None, dict[str, object]]:
    qualified, summaries = [], {}
    for value in WEAK_GRID:
        _spend(counts, "canonical", "offline")
        classification = tuple(item <= value for item in observation.normalized_singular_values)
        error = sum(left != right for left, right in zip(classification, (False, False, True), strict=True))
        summaries[value] = error
        if error == 0:
            qualified.append(value)
    status = _threshold_outcome(WEAK_GRID, qualified, 2.0**-2)
    if status != "SELECTED":
        return None, {"status": status}
    selected = min(qualified)
    index = WEAK_GRID.index(selected)
    neighbors = tuple({"candidate": WEAK_GRID[i], "status": "QUALIFIED" if WEAK_GRID[i] in qualified else "DISQUALIFIED", "error": summaries[WEAK_GRID[i]], "cost": 1, "reasons": ("weak_misclassification",) if summaries[WEAK_GRID[i]] else ("smallest_threshold_tie_break",)} for i in (index - 1, index + 1) if 0 <= i < len(WEAK_GRID))
    return selected, {"status": status, "value": selected, "truth_error": 0, "cost": 1, "frontier": ((selected, 0, 1),), "rejected_neighbors": neighbors}
def select_cluster_policy(observation: Any, counts: dict[str, int]) -> tuple[float | None, dict[str, object]]:
    spectrum = observation.singular_values
    isolated = tuple(np.asarray(item.projector) for item in observation.subspaces)
    expected_groups = ((0, 1), (2,), (3,), (4,))
    qualified, summaries = [], {}
    for value in CLUSTER_GRID:
        _spend(counts, "canonical", "offline")
        groups, start = [], 0
        for index, (left, right) in enumerate(zip(spectrum, spectrum[1:], strict=False)):
            if (left - right) / spectrum[0] > value:
                groups.append(tuple(range(start, index + 1)))
                start = index + 1
        groups.append(tuple(range(start, len(spectrum))))
        projector_error = max(float(np.linalg.norm(sum((isolated[i] for i in group), np.zeros_like(isolated[0])) - np.diag(tuple(1.0 if i in group else 0.0 for i in range(len(spectrum)))), ord=2)) for group in groups)
        error = sum(tuple(groups) != expected_groups for _ in (0,)) + int(projector_error > TOL["projector"])
        summaries[value] = error, projector_error
        if error == 0:
            qualified.append(value)
    status = _threshold_outcome(CLUSTER_GRID, qualified, 2.0**-2)
    if status != "SELECTED":
        return None, {"status": status}
    selected = min(qualified)
    index = CLUSTER_GRID.index(selected)
    neighbors = tuple({"candidate": CLUSTER_GRID[i], "status": "QUALIFIED" if CLUSTER_GRID[i] in qualified else "DISQUALIFIED", "error": summaries[CLUSTER_GRID[i]][0], "projector_error": summaries[CLUSTER_GRID[i]][1], "cost": 1, "reasons": ("cluster_or_projector_misclassification",) if summaries[CLUSTER_GRID[i]][0] else ("smallest_threshold_tie_break",)} for i in (index - 1, index + 1) if 0 <= i < len(CLUSTER_GRID))
    return selected, {"status": status, "value": selected, "truth_error": summaries[selected][1], "cost": 1, "frontier": ((selected, summaries[selected][1], 1),), "rejected_neighbors": neighbors}
def select_conditioning_policy(adequate: Any, inadequate: Any, counts: dict[str, int]) -> tuple[float | None, dict[str, object]]:
    conditions = (adequate.singular_values[0] / adequate.singular_values[-1], inadequate.singular_values[0] / inadequate.singular_values[-1])
    relative_error = max(abs(conditions[0] - 5.0e5) / 5.0e5, abs(conditions[1] - 2.0e6) / 2.0e6)
    qualified, summaries = [], {}
    for value in CONDITION_GRID:
        _spend(counts, "canonical", "offline")
        error = int(not (conditions[0] <= value < conditions[1])) + int(relative_error > TOL["conditioning"])
        summaries[value] = error
        if error == 0:
            qualified.append(value)
    status = _threshold_outcome(CONDITION_GRID, qualified, 2.0**52)
    if status != "SELECTED":
        return None, {"status": status}
    selected = qualified[0]
    index = CONDITION_GRID.index(selected)
    neighbors = tuple({"candidate": CONDITION_GRID[i], "status": "QUALIFIED" if CONDITION_GRID[i] in qualified else "DISQUALIFIED", "error": summaries[CONDITION_GRID[i]], "cost": 2, "reasons": ("conditioning_misclassification",) if summaries[CONDITION_GRID[i]] else ("smallest_limit_tie_break",)} for i in (index - 1, index + 1) if 0 <= i < len(CONDITION_GRID))
    return selected, {"status": status, "value": selected, "truth_error": relative_error, "cost": 2, "frontier": ((selected, relative_error, 2),), "rejected_neighbors": neighbors}
def select_correlation_policy(source: Any, scales: Mapping[str, float], selected: Any, counts: dict[str, int]) -> tuple[float | None, dict[str, object]]:
    pairs = ((0, 1), (2, 3), (4, 5), (6, 7))
    truth_raw = (0.75, 1.0 + 64.0 * 2.0**-52, -1.0 - 64.0 * 2.0**-52, 1.0 + 256.0 * 2.0**-52)
    expected_outcomes = ("AVAILABLE", "ENDPOINT_CANONICALIZED", "ENDPOINT_CANONICALIZED", "CORRELATION_RANGE_VIOLATION"); expected_values = (0.75, 1.0, -1.0, None)
    qualified, summaries = [], {}
    for multiplier in CORRELATION_GRID:
        _spend(counts, "canonical", "offline")
        policy = replace(_case_policy("C5", selected, _bc_scales("C5", scales)), correlation_roundoff_multiplier=multiplier)
        evidence = _correlation_evidence(source, policy, counts, "canonical")
        entries = tuple(evidence.entries[left][right] for left, right in pairs)
        claims = {item.name: item.state for item in evidence.claims}
        error = sum(item.outcome != outcome or item.value != value for item, outcome, value in zip(entries, expected_outcomes, expected_values, strict=True)) + sum(item.raw_value is None or abs(cast("float", item.raw_value) - truth) > TOL["correlation_interior_eps"] * 2.0**-52 for item, truth in zip(entries, truth_raw, strict=True))
        error += int(claims.get("CORRELATION_ENTRY_VALIDITY") is not ClaimState.VIOLATED or claims.get("CORRELATION_SCOPE_REPORTABILITY") is not ClaimState.VIOLATED or evidence.scope_reportable)
        summaries[multiplier] = error, tuple(item.outcome for item in entries), evidence.identity
        if error == 0:
            qualified.append(multiplier)
    status = _threshold_outcome(CORRELATION_GRID, qualified, 2.0**12)
    if status != "SELECTED":
        return None, {"status": status}
    chosen = qualified[0]
    index = CORRELATION_GRID.index(chosen)
    neighbors = tuple({"candidate": CORRELATION_GRID[i], "status": "QUALIFIED" if CORRELATION_GRID[i] in qualified else "DISQUALIFIED", "error": summaries[CORRELATION_GRID[i]][0], "cost": 1, "reasons": ("endpoint_misclassification",) if summaries[CORRELATION_GRID[i]][0] else ("smallest_multiplier_tie_break",)} for i in (index - 1, index + 1) if 0 <= i < len(CORRELATION_GRID))
    return chosen, {"status": status, "value": chosen, "truth_error": 0, "cost": 1, "frontier": ((chosen, 0, 1),), "typed_evidence_identity": summaries[chosen][2], "truth_raw": truth_raw, "outcomes": summaries[chosen][1], "rejected_neighbors": neighbors}
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
def _case_policy(name: str, selected: Any, scales: tuple[float, ...] | None = None) -> Any:
    scales = _setup(name)[1] if scales is None else scales
    fd = (selected.relative_step_tolerance, selected.roundoff_multiplier, selected.smaller_step_extent, selected.larger_step_extent)
    return replace(compose_uncertainty_policy(name, scales, fd, selected.svd_driver, (selected.rank_absolute_tolerance, selected.rank_relative_tolerance), selected.weak_relative_tolerance, selected.singular_value_cluster_relative_tolerance, selected.conditioning_limit, selected.correlation_roundoff_multiplier), calibration_identity=selected.calibration_identity)
def _correlation_evidence(source: Any, policy: Any, counts: dict[str, int], role: str) -> Any:
    _spend(counts, role, "correlation")
    return _correlations(source_identity=source.identity, accepted_result_identity=source.accepted_result_identity, accepted_occurrence_identity=source.accepted_occurrence_identity, source_family="constrained_propagation", output_ids=source.output_ids, units=source.output_units, covariance=source.covariance, residual_variance_degenerate=False, source_reportable=source.source_covariance.usable and source.claim("LOCAL_FIRST_ORDER_DEGENERACY") is ClaimState.SATISFIED, policy=policy, inherited_claims=source.claims, source_artifact=source)
def _normalization_truth(evidence: Any, name: str, scaling: Any) -> tuple[bool, dict[str, object]]:
    covariance = evidence.covariance; truth = np.asarray(truth_jacobian(name)); m, n, g, nu = len(truth), truth.shape[1], 1, len(truth) - truth.shape[1] - 1
    t, y, errors = ((-2.0, -1.0, 1.0, 3.0), (1.4, 2.1, 3.8, 5.9), (1.0, 0.5, 2.0, 1.5)) if name == "A3" else ((-3.0, 0.0, 2.0, 4.0, 5.0), (0.5, 2.2, 4.4, 7.1, 8.8), (0.75, 1.25, 0.5, 2.0, 1.0)); anchor = _setup(name)[0]; expected_chi = sum(value**2 for value in normalized_residuals(anchor[0], anchor[1], t, y, errors))
    chi_square = evidence.accepted_anchor.chi_square; phi = 1.0 if scaling is ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES else chi_square / nu
    if np.linalg.matrix_rank(truth) < n:
        passed = abs(chi_square - expected_chi) <= TOL["covariance"] and covariance is None and evidence.rank_diagnostic is not None and evidence.rank_diagnostic.rank < n and evidence.marginal_errors is None and evidence.correlations is None and any(item.category == "rank_deficient" for item in evidence.failures)
        return passed, {"m": m, "n": n, "g": g, "nu": nu, "chi_square": chi_square, "phi": phi, "covariance_truth": None, "marginal_reportable": False, "correlation_reportable": False}
    if covariance is None or evidence.marginal_errors is None or evidence.correlations is None:
        return False, {"reason": "missing_typed_covariance_evidence"}
    expected = np.linalg.inv(truth.T @ truth) * phi
    error = float(np.max(np.abs(np.asarray(covariance.covariance) - expected)) / max(1.0, float(np.max(np.abs(expected)))))
    passed = abs(covariance.chi_square - expected_chi) <= TOL["covariance"] and (covariance.retained_residual_count, covariance.controlled_coordinate_count, covariance.profiled_normalization_count, covariance.nominal_residual_degrees_of_freedom) == (m, n, g, nu)
    passed &= covariance.residual_variance_scaling is scaling and abs(covariance.residual_variance_scale - phi) <= TOL["covariance"] * max(1.0, abs(phi))
    passed &= error <= (TOL["holdout_covariance"] if name == "H3" else TOL["covariance"])
    passed &= evidence.marginal_errors.scope_reportable and evidence.correlations.scope_reportable and all(item.value is not None for item in evidence.marginal_errors.entries)
    return passed, {"m": m, "n": n, "g": g, "nu": nu, "chi_square": covariance.chi_square, "phi": covariance.residual_variance_scale, "covariance_error": error, "marginal_identity": evidence.marginal_errors.identity, "correlation_identity": evidence.correlations.identity}
def _constraint_truth(evidence: Any, name: str, numerical_partial: bool) -> tuple[bool, dict[str, object]]:
    covariance, jacobian, propagation, marginals, correlations = evidence.covariance, evidence.constraint_jacobian, evidence.constrained_propagation, evidence.constrained_marginal_errors, evidence.constrained_correlations
    if covariance is None or jacobian is None or propagation is None or marginals is None or (name == "H5" and correlations is None):
        return False, {"reason": "missing_typed_constraint_evidence"}
    partials = scientific_partials(*_setup(name)[0]) if name == "F2" else holdout_partials(*_setup(name)[0])
    expected_g = np.asarray(((1.0, 0.0), (0.0, 1.0), (1.0, 2.0), partials, (0.0, 0.0)))
    partial_error = max((abs(item.estimate - partials[item.argument_index]) / max(1.0, abs(partials[item.argument_index])) for item in jacobian.function_partial_diagnostics), default=math.inf)
    truth = np.asarray(truth_jacobian(name)); expected_local = np.linalg.inv(truth.T @ truth); expected_covariance = expected_g @ expected_local @ expected_g.T
    local_error = float(np.max(np.abs(np.asarray(covariance.covariance) - expected_local)) / max(1.0, float(np.max(np.abs(expected_local)))))
    propagated_error = float(np.max(np.abs(np.asarray(propagation.covariance) - expected_covariance)) / max(1.0, float(np.max(np.abs(expected_covariance)))))
    zero = len(expected_g) - 1
    passed = bool(np.array_equal(np.asarray(jacobian.matrix)[zero], np.zeros(2))) and np.allclose(jacobian.matrix, expected_g, rtol=0.0, atol=64.0 * 2.0**-52)
    passed &= partial_error <= 64.0 * 2.0**-52 and bool(jacobian.function_partial_diagnostics) and all((item.method == "finite_difference") is numerical_partial for item in jacobian.function_partial_diagnostics)
    passed &= propagated_error <= (TOL["holdout_covariance"] if name == "H5" else TOL["covariance"])
    passed &= all(value == 0.0 for value in propagation.factor[zero]) and all(propagation.covariance[zero][i] == propagation.covariance[i][zero] == 0.0 for i in range(len(expected_g)))
    passed &= propagation.claim("LOCAL_FIRST_ORDER_DEGENERACY") is ClaimState.VIOLATED and propagation.claim(f"OUTPUT_FIRST_ORDER_NONDEGENERACY::{jacobian.output_ids[zero]}") is ClaimState.VIOLATED and all(propagation.claim(f"OUTPUT_FIRST_ORDER_NONDEGENERACY::{item}") is ClaimState.SATISFIED for item in jacobian.output_ids[:zero])
    passed &= marginals.entries[zero].value is None and marginals.entries[zero].raw_value == 0.0 and not marginals.scope_reportable
    passed &= all(item.value is not None for item in marginals.entries[:zero]) and (correlations is None or (not correlations.scope_reportable and all(correlations.entries[i][j].value is not None for i in range(zero) for j in range(zero)) and all(correlations.entries[zero][i].value is None for i in range(len(expected_g)))))
    if name == "H5":
        passed &= local_error <= TOL["holdout_covariance"] and evidence.marginal_errors is not None and evidence.correlations is not None and evidence.marginal_errors.scope_reportable and evidence.correlations.scope_reportable
    return passed, {"constraint_jacobian_identity": jacobian.identity, "local_covariance_error": local_error, "propagation_identity": propagation.identity, "marginal_identity": marginals.identity, "correlation_identity": None if correlations is None else correlations.identity, "partial_error": partial_error, "propagated_covariance_error": propagated_error, "zero_gradient_row": zero}
def _constraint_jacobian_truth(evidence: Any, name: str) -> tuple[bool, dict[str, object]]:
    jacobian = evidence.constraint_jacobian
    if jacobian is None:
        return False, {"reason": "missing_typed_constraint_jacobian"}
    partials = scientific_partials(*_setup(name)[0]); expected = np.asarray(((1.0, 0.0), (0.0, 1.0), (1.0, 2.0), partials, (0.0, 0.0)))
    error = float(np.max(np.abs(np.asarray(jacobian.matrix) - expected)))
    return error <= 64.0 * 2.0**-52 and bool(jacobian.function_partial_diagnostics), {"constraint_jacobian_identity": jacobian.identity, "constraint_jacobian_error": error}
def _derive_case(name: str, selected: Any, environment: str, counts: dict[str, int], role: str, numerical_partial: bool = False, boundary_zeta: float | None = None, scaling: Any = None, stop_after: str | None = None, coordinate_scales: tuple[float, ...] | None = None) -> tuple[Any, dict[str, object]]:
    global _ACTIVE_COUNTS, _ACTIVE_ROLE
    _ACTIVE_COUNTS, _ACTIVE_ROLE = counts, role
    parameterization, engine, problem, accepted, requested = _fixture(name, counts, role)
    policy = replace(_case_policy(name, selected, coordinate_scales), residual_variance_scaling=scaling or ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES)
    capabilities: tuple[Any, ...] = ()
    if name in {"F2", "H5"}:
        function = scientific_value if name == "F2" else holdout_value
        counted = (counted_scientific_partial_a, counted_scientific_partial_b) if name == "F2" else (counted_holdout_partial_a, counted_holdout_partial_b)
        capabilities = (FunctionFiniteDifferenceCapability("f", None, _setup(name)[1], 4.0),) if numerical_partial else (FunctionAnalyticPartialCapability("f", None, hashlib.sha256(function.__name__.encode()).hexdigest(), counted),)
    compiled = compile_constraint_linearization_capabilities(parameterization, requested, capabilities)
    units = tuple((item, ParameterUnit.DIMENSIONLESS) for item in requested)
    output_scales = tuple((item, 1.0) for item in requested)
    if boundary_zeta is not None:
        problem = replace(problem, affine_half_spaces=(AffineHalfSpace("F2-boundary", (1.0, 0.0), accepted.vector[0] + boundary_zeta * math.sqrt(6.0 / 35.0)),))
        accepted = _accepted_copy(accepted, problem.identity, f"F2-boundary-{float(boundary_zeta).hex()}")
    counted_engine = _CountingEngine(engine, counts, role)
    with _charged_derivation(counts, role, stop_after, counted_engine) as cancellation:
        evidence = derive_uncertainty_evidence(accepted, problem=problem, parameterization=parameterization, engine=cast("EvaluationEngine", counted_engine), policy=policy, constrained_scope=requested, constrained_units=units, constrained_scales=output_scales, compiled_constraint_linearization=compiled, cancellation_probe=cancellation, resolved_environment_identity=environment)
    case = ("F1-absolute" if policy.residual_variance_scaling is ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES else "F1-estimated") if name == "A3" else name if boundary_zeta is None else f"F2-boundary-{boundary_zeta:g}"
    if evidence.residual_jacobian is None:
        return evidence, {"case": case, "passed": False, "reason": "missing_jacobian"}
    truth = np.asarray(truth_jacobian(name)); observed = np.asarray(evidence.residual_jacobian.matrix)
    envelope = TOL["derivative"] * (1.0 / np.asarray(_setup(name)[1])[None, :] + np.abs(truth))
    passed, details = bool(np.all(np.abs(observed - truth) <= envelope)), {}
    if name in {"A3", "H3"}:
        extra, details = _normalization_truth(evidence, name, policy.residual_variance_scaling); passed &= extra
    if name in {"F2", "H5"}:
        extra, details = _constraint_jacobian_truth(evidence, name) if stop_after == "constraint" else _constraint_truth(evidence, name, numerical_partial); passed &= extra
    if boundary_zeta is not None and evidence.covariance is not None:
        passed &= evidence.covariance.claim("AFFINE_FEASIBILITY") is (ClaimState.VIOLATED if boundary_zeta == 2.0 else ClaimState.SATISFIED)
    return evidence, {"case": case, "passed": passed, "policy_identity": policy.identity, "evidence_identity": evidence.identity, **details}
def _h4_spectral_case(selected: Any, environment: str, counts: dict[str, int], role: str) -> dict[str, object]:
    evidence, record = _derive_case("H4", selected, environment, counts, role, stop_after="rank")
    diagnostic = evidence.rank_diagnostic
    if diagnostic is None:
        return {**record, "passed": False, "reason": "missing_rank_diagnostic"}
    expected = np.diag((1.0, 1.0, 0.0))
    cluster = next((item for item in diagnostic.subspaces if item.indices == (0, 1)), None)
    projector_error = math.inf if cluster is None else float(np.linalg.norm(np.asarray(cluster.projector) - expected, ord=2))
    passed = diagnostic.rank == 3 and diagnostic.weak_projector[2][2] > 1.0 - TOL["holdout_projector"] and diagnostic.singular_values[0] / diagnostic.singular_values[-1] > selected.conditioning_limit and projector_error <= TOL["holdout_projector"]
    return {**record, "passed": passed, "projector_error": projector_error, "rank_identity": diagnostic.identity}
def validate_composed_cases(selected: Any, environment: str, counts: dict[str, int]) -> tuple[bool, tuple[dict[str, object], ...]]:
    records = []
    cases = (("A1", False, None, None, "residual"), ("A2", False, None, None, "residual"), ("A3", False, None, ResidualVarianceScaling.ABSOLUTE_OBSERVATION_UNCERTAINTIES, None), ("A3", False, None, ResidualVarianceScaling.ESTIMATED_COMMON_RESIDUAL_VARIANCE, None), ("F2", True, None, None, "constrained_marginal"), ("F2", False, 2.0, None, "constrained_marginal"), ("F2", False, 4.0, None, "constrained_marginal"))
    for name, numerical, zeta, scaling, stop_after in cases:
        _evidence, record = _derive_case(name, selected, environment, counts, "canonical", numerical, zeta, scaling, stop_after)
        records.append(record)
        if not record["passed"]:
            return False, tuple(records)
    return True, tuple(records)
def run_holdouts(selected: Any, environment: str, counts: dict[str, int]) -> tuple[bool, tuple[dict[str, object], ...]]:
    records: list[dict[str, object]] = []
    for name in ("H1", "H2", "H3"):
        scaling = ResidualVarianceScaling.ESTIMATED_COMMON_RESIDUAL_VARIANCE if name == "H3" else None
        _evidence, record = _derive_case(name, selected, environment, counts, "canonical", scaling=scaling, stop_after="residual" if name in {"H1", "H2"} else None)
        records.append(record)
        if not record["passed"]:
            return False, tuple(records)
    record = _h4_spectral_case(selected, environment, counts, "canonical")
    records.append(record)
    if not record["passed"]:
        return False, tuple(records)
    _evidence, record = _derive_case("H5", selected, environment, counts, "canonical", numerical_partial=True)
    records.append(record)
    if not record["passed"]:
        return False, tuple(records)
    evidence, record = _derive_case("H6", selected, environment, counts, "canonical")
    correlation = None
    if evidence.constrained_propagation is not None:
        source = evidence.constrained_propagation; policy = _case_policy("H6", selected)
        correlation = _correlation_evidence(source, policy, counts, "canonical")
    passed = correlation is not None and correlation.entries[0][1].outcome == "ENDPOINT_CANONICALIZED" and correlation.entries[0][1].value == -1.0
    records.append({**record, "passed": passed, "correlation_identity": None if correlation is None else correlation.identity})
    return passed, tuple(records)
def run_compatibility_replay(selected: Any, scales: Mapping[str, float], environment: str, counts: dict[str, int]) -> dict[str, object]:
    records = []
    for name in ("A1", "A2"):
        _evidence, record = _derive_case(name, selected, environment, counts, "compatibility", stop_after="residual")
        records.append(record)
        if not record["passed"]:
            return {"status": "UNAVAILABLE", "results": tuple(records), "retuning": False}
    _evidence, f2 = _derive_case("F2", selected, environment, counts, "compatibility", numerical_partial=True, stop_after="constraint")
    records.append(f2)
    if not f2["passed"]:
        return {"status": "UNAVAILABLE", "results": tuple(records), "retuning": False}
    h4 = _h4_spectral_case(selected, environment, counts, "compatibility")
    records.append(h4)
    if not h4["passed"]:
        return {"status": "UNAVAILABLE", "results": tuple(records), "retuning": False}
    evidence, c5 = _derive_case("C5", selected, environment, counts, "compatibility", stop_after="constrained_marginal", coordinate_scales=_bc_scales("C5", scales))
    source = evidence.constrained_propagation
    passed = source is not None
    if source is not None:
        policy = _case_policy("C5", selected, _bc_scales("C5", scales))
        correlation = _correlation_evidence(source, policy, counts, "compatibility")
        passed = correlation.entries[2][3].outcome == "ENDPOINT_CANONICALIZED"
        c5["correlation_identity"] = correlation.identity
    records.append({**c5, "passed": passed})
    return {"status": "COMPATIBLE" if passed and tuple(item["case"] for item in records) == COMPATIBILITY_CASES else "UNAVAILABLE", "results": tuple(records), "retuning": False}
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
    return {"source_digest": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "specification_id": SPECIFICATION_ID, "specification_digest": frozen_digest(), "status": status, "counts": dict(counts), **items}
def acquire_canonical(environment: str) -> dict[str, object]:
    counts = dict.fromkeys(CEILINGS["canonical"], 0)
    scales, scale_record, trajectories, scalar = select_scales(counts)
    if scales is None:
        return _compact_result(cast("str", scale_record["status"]), counts, scale_selection=scale_record)
    fd, fd_record = calibrate_finite_differences(scales, trajectories, scalar, counts)
    if fd is None:
        return _compact_result(cast("str", fd_record["status"]), counts, scale_selection=scale_record, finite_difference=fd_record)
    phase_a = compose_uncertainty_policy("A1", (scales["state_population"], scales["exchange_rate"]), fd, "gesvd", (0.0, 0.0), 0.0, 0.0, 2.0**52, 0.0)
    driver, driver_record, rank_observations = select_svd_driver(scales, phase_a, environment, counts)
    if driver is None:
        return _compact_result("UNSUPPORTED", counts, finite_difference=fd_record, svd_driver=driver_record)
    rank, rank_record = select_rank_policy(rank_observations, counts)
    if rank is None:
        return _compact_result(cast("str", rank_record["status"]), counts, svd_driver=driver_record, rank=rank_record)
    phase_ab = replace(phase_a, svd_driver=driver, rank_absolute_tolerance=rank[0], rank_relative_tolerance=rank[1])
    rank_truth_passed, rank_truth = validate_rank_truth_cases(scales, phase_ab, environment, counts)
    rank_record = {**rank_record, "typed_truth_cases": rank_truth}
    if not rank_truth_passed:
        return _compact_result("UNSUPPORTED", counts, phase="B", svd_driver=driver_record, rank=rank_record)
    c_observations, c_records = [], []
    for name in ("C1", "C2", "C3", "C4"):
        diagnostic, record = _spectral_observation(name, scales, phase_ab, driver, rank, environment, counts)
        if diagnostic is None:
            return _compact_result("UNSUPPORTED", counts, phase="C", case=record)
        c_observations.append(diagnostic); c_records.append(record)
    c5_evidence, c5_record = _derive_case("C5", phase_ab, environment, counts, "canonical", stop_after="constrained_marginal", coordinate_scales=_bc_scales("C5", scales))
    source = c5_evidence.constrained_propagation
    if source is None:
        return _compact_result("UNSUPPORTED", counts, phase="C5", case=c5_record)
    weak, weak_record = select_weak_policy(c_observations[0], counts)
    cluster, cluster_record = select_cluster_policy(c_observations[1], counts)
    conditioning, conditioning_record = select_conditioning_policy(c_observations[2], c_observations[3], counts)
    correlation, correlation_record = select_correlation_policy(source, scales, phase_ab, counts)
    if None in {weak, cluster, conditioning, correlation}:
        return _compact_result("UNSUPPORTED_OR_SATURATED", counts, weak=weak_record, cluster=cluster_record, conditioning=conditioning_record, correlation=correlation_record)
    selected = replace(compose_uncertainty_policy("A1", (scales["state_population"], scales["exchange_rate"]), fd, driver, rank, cast("float", weak), cast("float", cluster), cast("float", conditioning), cast("float", correlation)), calibration_identity=f"{SPECIFICATION_ID}:{_scale_catalogue_digest(scales)}")
    valid, composed = validate_composed_cases(selected, environment, counts)
    if not valid:
        return _compact_result("COMPOSED_VALIDATION_FAILED", counts, policy=policy_record(selected), composed=composed)
    passed, holdouts = run_holdouts(selected, environment, counts)
    phases = {"scales": scales, "finite_difference": fd, "svd_driver": driver, "rank": rank, "weak": weak, "cluster": cluster, "conditioning": conditioning, "correlation": correlation}
    metrics = {"scale": scale_record, "finite_difference": fd_record, "driver": driver_record, "rank": rank_record, "phase_c_cases": (*c_records, c5_record), "weak": weak_record, "cluster": cluster_record, "conditioning": conditioning_record, "correlation": correlation_record}
    serialized = policy_record(selected)
    digest = hashlib.sha256(json.dumps(serialized, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return _compact_result("QUALIFIED" if passed else "HOLDOUT_FAILED_POLICY_UNAVAILABLE", counts, selected_phase_policies=phases, policy=serialized, policy_digest=digest, decisive_metrics=metrics, composed=composed, holdouts=holdouts)
def policy_from_canonical_record(record: Mapping[str, object]) -> Any:
    required = {"source_digest", "specification_id", "specification_digest", "status", "canonical_provenance", "selected_phase_policies", "policy", "policy_digest", "decisive_metrics", "composed", "holdouts"}
    if not required <= set(record) or record["source_digest"] != hashlib.sha256(Path(__file__).read_bytes()).hexdigest() or record["specification_id"] != SPECIFICATION_ID or record["specification_digest"] != frozen_digest() or record["status"] != "QUALIFIED":
        raise RuntimeError("compatibility requires the successful canonical acquisition record")
    validate_lane_records(cast("Mapping[str, object]", record["canonical_provenance"]), "canonical")
    metrics = cast("Mapping[str, Mapping[str, object]]", record["decisive_metrics"])
    if any(metrics[name]["status"] != "SELECTED" for name in ("scale", "finite_difference", "driver", "rank", "weak", "cluster", "conditioning", "correlation")):
        raise RuntimeError("canonical phase selection is incomplete")
    driver_cases = cast("tuple[Mapping[str, object], ...]", metrics["driver"].get("cases", ())); rank_cases = cast("tuple[Mapping[str, object], ...]", metrics["rank"].get("typed_truth_cases", ())); phase_c = cast("tuple[Mapping[str, object], ...]", metrics.get("phase_c_cases", ()))
    if tuple(item.get("case") for item in driver_cases) != ("B1", "B2", "B3") or not all(item.get("status") == "QUALIFIED" for item in driver_cases) or tuple(item.get("case") for item in rank_cases) != ("B2", "B3") or tuple(item.get("case") for item in phase_c) != ("C1", "C2", "C3", "C4", "C5") or not all(item.get("passed", item.get("status") == "ACQUIRED") is True for item in (*rank_cases, *phase_c)):
        raise RuntimeError("canonical typed phase evidence is incomplete")
    composed = cast("tuple[Mapping[str, object], ...]", record["composed"]); holdouts = cast("tuple[Mapping[str, object], ...]", record["holdouts"])
    if tuple(item.get("case") for item in composed) != COMPOSED_CASES or tuple(item.get("case") for item in holdouts) != HOLDOUT_CASES or not all(item.get("passed") is True for item in (*composed, *holdouts)):
        raise RuntimeError("canonical composed validation or holdouts did not pass")
    serialized = cast("Mapping[str, object]", record["policy"])
    digest = hashlib.sha256(json.dumps(serialized, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    policy = policy_from_record(serialized)
    phases = cast("Mapping[str, object]", record["selected_phase_policies"])
    fd = (policy.relative_step_tolerance, policy.roundoff_multiplier, policy.smaller_step_extent, policy.larger_step_extent)
    exact = set(phases) == {"scales", "finite_difference", "svd_driver", "rank", "weak", "cluster", "conditioning", "correlation"} and fd in finite_difference_policies() and fd == cast("tuple[float, float, int, int]", tuple(cast("tuple[object, ...]", metrics["finite_difference"]["policy"]))) == cast("tuple[float, float, int, int]", tuple(cast("tuple[object, ...]", phases["finite_difference"]))) and policy.svd_driver in SVD_DRIVERS and policy.svd_driver == metrics["driver"]["driver"] == phases["svd_driver"] and (policy.rank_absolute_tolerance, policy.rank_relative_tolerance) in product(RANK_GRID, RANK_GRID) and (policy.rank_absolute_tolerance, policy.rank_relative_tolerance) == cast("tuple[float, float]", tuple(cast("tuple[object, ...]", metrics["rank"]["policy"]))) == cast("tuple[float, float]", tuple(cast("tuple[object, ...]", phases["rank"])))
    exact &= policy.weak_relative_tolerance in WEAK_GRID and policy.weak_relative_tolerance == metrics["weak"]["value"] == phases["weak"] and policy.singular_value_cluster_relative_tolerance in CLUSTER_GRID and policy.singular_value_cluster_relative_tolerance == metrics["cluster"]["value"] == phases["cluster"] and policy.conditioning_limit in CONDITION_GRID and policy.conditioning_limit == metrics["conditioning"]["value"] == phases["conditioning"] and policy.correlation_roundoff_multiplier in CORRELATION_GRID and policy.correlation_roundoff_multiplier == metrics["correlation"]["value"] == phases["correlation"]
    scales = cast("Mapping[str, float]", phases["scales"])
    q0 = {family: base for family, _stratum, base, _unit, _patterns in SUPPORTED} | {"A4-argument-0": 1.0, "A4-argument-1": 2.0, "A4-output": 4.0}
    exact &= set(scales) == set(q0) and all(scales[name] in tuple(math.ldexp(base, exponent) for exponent in SCALE_EXPONENTS) for name, base in q0.items())
    scale_records = cast("Mapping[str, Mapping[str, object]]", metrics["scale"].get("families", {})); supported_q0 = {family: base for family, _stratum, base, _unit, _patterns in SUPPORTED}
    exact &= supported_q0.keys() <= scale_records.keys() and all((item := scale_records[family]).get("status") == "SELECTED" and item.get("family") == family and item.get("q0") == base and item.get("exponent") in SCALE_EXPONENTS and item.get("scale") == scales[family] == math.ldexp(base, cast("int", item["exponent"])) and item.get("candidate_identity") == _scale_candidate_identity(family, base, cast("int", item["exponent"]), scales[family]) for family, base in supported_q0.items())
    exact &= policy.calibration_identity == f"{SPECIFICATION_ID}:{_scale_catalogue_digest(scales)}" and dict(policy.coordinate_scales) == {"A": scales["state_population"], "B": scales["exchange_rate"]}
    if digest != record["policy_digest"] or not exact:
        raise RuntimeError("serialized policy is not the canonical selected policy")
    return policy
def acquire(role: str, image_digest: str, canonical_record: Mapping[str, object] | None = None) -> dict[str, object]:
    if role == "compatibility" and canonical_record is None:
        raise RuntimeError("compatibility replay requires the canonical acquisition record")
    lane_records = attest_lane(role, image_digest)
    _import_machinery()
    environment = validate_lane_records(lane_records, role)
    if role == "canonical":
        result = acquire_canonical(environment)
    else:
        counts = dict.fromkeys(CEILINGS["compatibility"], 0)
        selected = policy_from_canonical_record(cast("Mapping[str, object]", canonical_record))
        phases = cast("Mapping[str, object]", cast("Mapping[str, object]", canonical_record)["selected_phase_policies"])
        compatibility = run_compatibility_replay(selected, cast("Mapping[str, float]", phases["scales"]), environment, counts)
        result = _compact_result(cast("str", compatibility["status"]), counts, policy_identity=selected.identity, compatibility=compatibility)
    return {**result, f"{role}_provenance": lane_records}
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--role", choices=("canonical", "compatibility"), required=True)
    parser.add_argument("--image-digest", required=True)
    parser.add_argument("--canonical-record", type=Path)
    arguments = parser.parse_args()
    record = None if arguments.canonical_record is None else cast("Mapping[str, object]", json.loads(arguments.canonical_record.read_text()))
    arguments.output.write_text(json.dumps(acquire(arguments.role, arguments.image_digest, record), allow_nan=False, indent=2, sort_keys=True) + "\n")
if __name__ == "__main__":
    main()

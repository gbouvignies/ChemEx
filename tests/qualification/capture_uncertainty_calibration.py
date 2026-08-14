"""Frozen prospective acquisition specification for issue #601.

Commit A intentionally freezes inputs, staging, accounting, and compact reducers.
Scientific acquisition is not performed by importing or testing this module.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterable, Iterator, Mapping
from itertools import product
from pathlib import Path
from typing import Any, cast

from chemex.numerical_lanes import (
    LaneAttestation,
    NumericalLane,
    RuntimeEnvironment,
    canonical_lanes,
)

SPECIFICATION_ID = "chemex-uncertainty-calibration-v1"
SCALE_EXPONENTS = tuple(range(-8, 9))
RELATIVE_STEP_TOLERANCES = (0.0,) + tuple(2.0**-power for power in range(8, 49, 2))
ROUNDOFF_MULTIPLIERS = tuple(2.0**power for power in range(13))
STEP_EXTENTS = (0, 2, 4, 8, 12, 16)
SVD_DRIVERS = ("gesvd", "gesdd")
RANK_TOLERANCES = (0.0,) + tuple(2.0**-power for power in range(16, 53, 2))
WEAK_MODE_THRESHOLDS = (0.0,) + tuple(2.0**-power for power in range(2, 53, 2))
CLUSTER_THRESHOLDS = WEAK_MODE_THRESHOLDS
CONDITIONING_LIMITS = (1.0,) + tuple(2.0**power for power in range(2, 53, 2))
CORRELATION_MULTIPLIERS = (0.0,) + tuple(2.0**power for power in range(13))
CANDIDATE_COUNTS = {
    "characteristic_scale_per_family": 17,
    "finite_difference": 10_296,
    "svd_driver": 2,
    "rank": 400,
    "weak_mode": 27,
    "cluster": 27,
    "conditioning": 27,
    "correlation": 14,
}
SCALE_SELECTION_SCOPE = {
    "coordinate": "q=q0*2^k, k=-8..8, selected per supported semantic family/input",
    "output": "Q=Q0*2^k, k=-8..8, selected per required output family",
    "edge_only": "only k=-8 or k=8 qualifies => GRID_SATURATED and stop",
    "fallback": None,
}
FINITE_DIFFERENCE_SEMANTICS = {
    "tuple": "(tau_rel,kappa,E_minus,E_plus)",
    "exponents": "range(-E_minus,E_plus+1)",
    "independent_extents": True,
    "ordering": "increasing abs(exponent), then represented step; feasible bound anchor included",
    "centered_multipliers": (-2.0, -1.0, 1.0, 2.0),
    "one_sided_multipliers": (1.0, 2.0, 4.0),
    "negative_one_sided_multipliers": (-1.0, -2.0, -4.0),
}

# family, stratum, q0, unit, complete non-WIP CLI-reachable members/patterns
SUPPORTED_SCALE_FAMILIES = (
    ("solvent_fraction", "fraction", 1.0, "1", ("d2o",)),
    (
        "state_population",
        "fraction",
        1.0,
        "1",
        ("pa", "pb", "pc", "pd", "pe", "pf", "pop_b"),
    ),
    (
        "equilibrium_ratio",
        "dimensionless",
        1.0,
        "1",
        ("keq", "keq_bc", "keq_cd", "keq_l", "keq_pl"),
    ),
    ("phase_factor", "dimensionless", 1.0, "1", ("phi", "phi_a", "phi_b")),
    ("order_parameter", "dimensionless", 1.0, "1", ("s2", "s2_*")),
    (
        "directional_exchange_rate",
        "rate",
        100.0,
        "s^-1",
        tuple(
            f"k{left}{right}"
            for left in "abcdef"
            for right in "abcdef"
            if left != right
        ),
    ),
    (
        "exchange_rate",
        "rate",
        100.0,
        "s^-1",
        (
            "kex",
            "kex_ab",
            "kex_ac",
            "kex_ad",
            "kex_ae",
            "kex_af",
            "kex_bc",
            "kex_bd",
            "kex_be",
            "kex_bf",
            "kex_cd",
            "kex_ce",
            "kex_cf",
            "kex_de",
            "kex_df",
            "kex_ef",
        ),
    ),
    (
        "dissociation_rate",
        "rate",
        100.0,
        "s^-1",
        ("koff", "koff1", "koff2", "koff_ab", "koff_ac", "koff_bc"),
    ),
    ("hydrogen_exchange_rate", "rate", 1.0, "s^-1", ("kdh", "kdh_a", "kdh_b")),
    ("generated_hydrogen_exchange_rate", "rate", 1.0, "s^-1", ("khh", "khh_*")),
    ("longitudinal_relaxation", "rate", 1.0, "s^-1", ("r1", "r1_*")),
    (
        "transverse_relaxation",
        "rate",
        10.0,
        "s^-1",
        ("r2", "r2_*", "r2a", "r2a_*", "r2mq", "r2mq_*", "r1a", "r1a_*"),
    ),
    (
        "cross_relaxation",
        "rate",
        10.0,
        "s^-1",
        (
            "etaxy",
            "etaxy_*",
            "etaz",
            "etaz_*",
            "sigma",
            "sigma_*",
            "mu",
            "mu_*",
        ),
    ),
    ("chemical_shift", "frequency", 1.0, "ppm", ("cs", "cs_*")),
    (
        "chemical_shift_difference",
        "frequency",
        1.0,
        "ppm",
        ("dw", "dw_*", "dwp", "dwp_*"),
    ),
    ("coupling", "frequency", 100.0, "Hz", ("j", "j_*", "d", "d_*")),
)

UNSUPPORTED_SCALE_FAMILIES = frozenset(
    [
        "c_dimer",
        "c_l",
        "c_monomer",
        "c_p",
        "c_pl",
        "c_pl1",
        "c_pl2",
        "c_pl3",
        "c_tetramer",
        "c_trimer",
        "dimer",
        "l1_free",
        "l_free",
        "monomer",
        "p_free",
        "pfree",
        "pl",
        "pl1",
        "pl2",
        "pl3",
        "kd",
        "kd1",
        "kd2",
        "kd_ab",
        "kd_ac",
        "kd_app",
        "kd_bc",
        "kd_eff",
        "kon",
        "kon1",
        "kon2",
        "kon_ab",
        "kon_ac",
        "kon_bc",
        "dh_ab",
        "dh_ac",
        "dh_ad",
        "dh_b",
        "dh_bc",
        "dh_bd",
        "dh_c",
        "dh_cd",
        "dh_d",
        "ds_ab",
        "ds_ac",
        "ds_ad",
        "ds_b",
        "ds_bc",
        "ds_bd",
        "ds_c",
        "ds_cd",
        "ds_d",
        "tauc",
        "dwm",
    ]
)
UNSUPPORTED_JUSTIFICATIONS = {
    "concentration": "condition-dependent input/output concentration unit and q0",
    "dissociation_constant": "concentration unit and normalization unresolved",
    "association_rate": "compound concentration/time unit unresolved",
    "thermodynamic": "energy and entropy unit conventions unresolved",
    "tauc": "nanosecond input versus second propagation unresolved",
    "dwm": "ppm-per-degree-Celsius is not a frozen frequency family",
}

PHASE_OWNERSHIP = (
    ("A0", "characteristic_scale"),
    ("A1", "finite_difference"),
    ("B1", "svd_driver"),
    ("B2", "rank_threshold"),
    ("C1", "weak_mode"),
    ("C2", "cluster"),
    ("C3", "conditioning"),
    ("C4", "correlation_endpoint"),
    ("F", "composed_policy"),
    ("H", "holdout_validation_only"),
)

# Exact functions and coefficients are strings/tuples so the committed record is
# directly inspectable without importing NumPy, SciPy, or uncertainty machinery.
TRUTH_CASES = {
    "A1": {
        "anchor": (0.25, -20.0),
        "q": (1.0, 100.0),
        "uv": (0.25, -0.2),
        "f": "(exp(u)+v^2/2, sin(u)-2v, uv+1/4)",
        "j": "((exp(u)/qx,v/qy),(cos(u)/qx,-2/qy),(v/qx,u/qy))",
        "claims": ("centered", "mixed_scale", "discrepancy"),
    },
    "A2": {
        "anchor": (0.0, 2.0),
        "q": (1.0, 10.0),
        "bounds": ((0.0, 10.0), (-10.0, 10.0)),
        "affine": "y<=2",
        "f": "(exp(x)+sin(y/10),x^2+3y/10)",
        "j": "((1,cos(.2)/10),(0,.3))",
        "claims": ("positive_one_sided", "negative_one_sided", "boundary"),
    },
    "A3": {
        "t": (-2.0, -1.0, 1.0, 3.0),
        "y": (1.4, 2.1, 3.8, 5.9),
        "e": (1.0, 0.5, 2.0, 1.5),
        "anchor": (3.0, 1.0),
        "q": (1.0, 1.0),
        "formula": "c=a+bt;w=e^-2;alpha=sum(ycw)/sum(c^2w);r=(alpha*c-y)/e",
        "dof": (4, 2, 1, 1),
    },
    "A4": {
        "anchor": (0.5, 2.0),
        "q0": (1.0, 2.0),
        "Q0": 4.0,
        "candidates": "each q0*2^k and Q0*2^k for k=-8..8",
        "f": "exp(a/b)+log1p(b^2)",
        "partials": ("exp(a/b)/b", "-a*exp(a/b)/b^2+2b/(1+b^2)"),
    },
    "B1": {"matrix": ((1.0e-9, 0.0), (0.0, 2.0e-21)), "rank": 2},
    "B2": {
        "matrix": ((1.0, 0.0), (0.0, 5.0e-13)),
        "perturbation": 1.0e-12,
        "rank": 1,
        "null": (0.0, 1.0),
    },
    "B3": {
        "matrix": ((1.0, 1.0), (2.0, 2.0), (3.0, 3.0)),
        "rank": 1,
        "null": (1.0, -1.0),
        "covariance": None,
    },
    "C1": {"spectrum": (1.0, 2.0e-6, 2.0e-7), "weak": (False, False, True)},
    "C2": {
        "spectrum": (4.0, 4.0 * (1.0 - 5.0e-11), 2.0, 2.0 * (1.0 - 5.0e-9), 1.0),
        "clusters": ((0, 1), (2,), (3,), (4,)),
    },
    "C3": {"spectrum": (1.0, 2.0e-6), "condition": 5.0e5, "adequate": True},
    "C4": {"spectrum": (1.0, 5.0e-7), "condition": 2.0e6, "adequate": False},
    "C5": {
        "rho": (
            0.75,
            1.0 + 64.0 * 2.0**-52,
            -(1.0 + 64.0 * 2.0**-52),
            1.0 + 128.0 * 2.0**-52,
        ),
        "truth": ("interior", "+1", "-1", "invalid"),
    },
    "F1": {
        "reuses": "A3",
        "variance_scaling": (
            "absolute_observation_uncertainties",
            "estimated_common_residual_variance",
        ),
    },
    "F2": {
        "J": ((1.0, 0.0), (0.0, 2.0), (1.0, 1.0), (2.0, -1.0)),
        "anchor": (0.5, 2.0),
        "C": ((6.0 / 35.0, 1.0 / 35.0), (1.0 / 35.0, 6.0 / 35.0)),
        "outputs": ("x", "y", "x+2y", "A4", "zero_gradient"),
        "boundary": (
            "n=(1,0)",
            "b=.5+2sqrt(6/35) violated",
            "b=.5+4sqrt(6/35) satisfied",
        ),
    },
}

HOLDOUT_CASES = {
    "H1": {
        "anchor": (3.0, -0.002),
        "q": (10.0, 0.01),
        "uv": (0.3, -0.2),
        "f": "(tanh(u)+exp(v),log1p(u^2)+v^3,cos(u-v))",
    },
    "H2": {
        "anchor": (1.0, -5.0),
        "q": (1.0, 10.0),
        "bounds": ((-2.0, 1.0), (-10.0, 10.0)),
        "affine": "-y<=5",
        "f": "(exp(x)+sin(y/10),x^2+3y/10)",
    },
    "H3": {
        "t": (-3.0, 0.0, 2.0, 4.0, 5.0),
        "y": (0.5, 2.2, 4.4, 7.1, 8.8),
        "e": (0.75, 1.25, 0.5, 2.0, 1.0),
        "anchor": (2.0, 1.5),
        "q": (1.0, 2.0),
        "dof": (5, 2, 1, 2),
    },
    "H4": {
        "R": ((0.6, -0.8, 0.0), (0.8, 0.6, 0.0), (0.0, 0.0, 1.0)),
        "s": (3.0, 3.0 * (1.0 - 7.5e-11), 3.0e-7),
        "truth": ("first_pair_clustered", "last_weak", "condition_1e7_inadequate"),
    },
    "H5": {
        "J": ((1.0, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, -1.0)),
        "anchor": (0.25, -0.5),
        "q": (0.5, 1.0),
        "f": "log1p(a^2)+sin(b)",
        "degenerate": "zero_gradient",
    },
    "H6": {"rho": -(1.0 + 48.0 * 2.0**-52), "truth": "-1_endpoint"},
}
COMPATIBILITY_REPLAY_CASES = ("A1", "A2", "H4", "C5", "F2")

METRIC_TOLERANCES = {
    "derivative_relative_envelope": 2.0**-24,
    "truth_uncertainty_fraction": 1.0 / 16.0,
    "covariance_calibration_relative": 2.0e-6,
    "covariance_holdout_relative": 5.0e-6,
    "projector_calibration_norm": 2.0e-12,
    "projector_holdout_norm": 5.0e-12,
    "conditioning_calibration_relative": 2.0e-13,
    "conditioning_holdout_relative": 5.0e-13,
    "correlation_interior_epsilon": 8.0,
    "boundary_epsilon_multiplier": 128.0,
    "boundary_zeta": 3.0,
    "constraint_partial_epsilon_multiplier": 64.0,
    "zero_gradient_absolute": 0.0,
    "false_acceptance_count": 0,
}
METRIC_FORMULAS = {
    "derivative": "abs(D-D*)<=2^-24*(Q/q+abs(D*)); exact zero<=2^-24*Q/q",
    "reliability": "delta<=tau_rel*maxnorm(fine,coarse)+kappa*eps*F/min|dx|",
    "covariance": "maxabs(C-C*)/max(1,maxabs(C*))",
    "rank": "count(s_i>a_rank+rho_rank*s_1), strict",
    "weak": "s_i/s_1<=omega",
    "cluster": "(s_i-s_(i+1))/s_1<=gamma; chain; compare projector norm",
    "conditioning": "adequate iff s_1/s_n<=limit",
    "correlation": "clip iff abs(abs(rho)-1)<=multiplier*2^-52; otherwise reject",
    "boundary": "violated iff zeta<=3; geometry tolerance 128*eps*max(1,norm)",
    "normalization_dof": "nu=m-n-g",
    "constraint_partials": "abs(error)<=64*eps*max(1,abs(truth))",
    "propagation": "G*C*G.T; error=sqrt(diagonal)",
    "zero_gradient": "exact zero gradient => exact zero variance and error",
}

RESOURCE_LEDGER = {
    "canonical_python_313": {
        "evaluation_engine_requests": 25_609,
        "scalar_function_requests": 3_320,
        "svd_kernels": 22,
        "correlation_kernels": 14,
        "candidate_selection": 10_840,
        "holdouts": 6,
    },
    "compatibility_python_314": {
        "evaluation_engine_requests": 751,
        "scalar_function_requests": 266,
        "svd_kernels": 4,
        "correlation_kernels": 4,
        "candidate_selection": 0,
        "holdouts": 0,
    },
    "planned_maximum": {
        "evaluation_engine_requests": 26_360,
        "scalar_function_requests": 3_586,
        "svd_kernels": 26,
        "correlation_kernels": 18,
        "offline_candidate_comparisons": 12_000_000,
    },
    "ceiling": {
        "evaluation_engine_requests": 28_000,
        "scalar_function_requests": 4_000,
        "svd_kernels": 32,
        "correlation_kernels": 24,
        "offline_candidate_comparisons": 12_000_000,
    },
}
RESOURCE_BREAKDOWN = {
    "canonical_python_313": {
        "union_fd_evaluation_engine": 6_033,
        "selected_and_neighbor_evaluation_engine": 18_279,
        "composed_and_holdout_evaluation_engine": 1_297,
        "union_scalar_function": 394,
        "selected_and_neighbor_scalar_function": 2_394,
        "composed_and_holdout_scalar_function": 532,
    },
    "compatibility_python_314": {
        "selected_policy_evaluation_engine": 751,
        "selected_policy_scalar_function": 266,
    },
}
REQUEST_CEILINGS = {
    "union_centered_column": {"evaluation_engine": 201},
    "union_one_sided_column": {"evaluation_engine": 151},
    "union_normalization_case_A3": {"evaluation_engine": 401},
    "union_function_case_A4": {"scalar_function": 394},
    "direct_centered_column": {"evaluation_engine": 137},
    "direct_one_sided_column": {"evaluation_engine": 103},
    "direct_two_centered_columns": {"evaluation_engine": 273},
    "direct_two_one_sided_columns": {"evaluation_engine": 205},
    "direct_two_function_partials": {"scalar_function": 266},
    "offline_fd_candidate": {"evaluation_engine": 0, "comparisons": 36},
    "offline_rank_candidate": {"svd_kernels": 0, "comparisons": 3},
    "offline_threshold_candidate": {"numerical_kernels": 0},
    "svd_matrix_driver_pair": {"svd_kernels": 2},
}
ACQUISITION_PLAN = {
    "canonical": (
        ("A0", "acquire each family union trajectory once; replay 17 q0*2^k scales"),
        ("A1", "stream all 10,296 policies over selected-scale fixed observations"),
        ("B1", "run gesvd and gesdd once per fixed matrix; no fallback"),
        ("B2", "stream 400 rank pairs over retained selected-driver spectra"),
        ("C", "stream 27/27/27/14 thresholds over fixed spectra/covariance"),
        (
            "F",
            "confirm selected policy and decisive neighbors with typed #598 evidence",
        ),
        ("H", "run six holdouts in order; stop at first failure; no runner-up"),
    ),
    "compatibility": (
        ("P", "reuse #588/#591 claims, then replay only A1,A2,H4,C5,F2"),
        ("D", "failure marks policy unavailable on 3.14; never retune or run holdouts"),
    ),
}
SELECTION_RULES = {
    "hard_disqualifiers": (
        "non_finite",
        "false_acceptance",
        "truth_envelope",
        "incomplete_trajectory",
        "request_ceiling",
        "kernel_failure",
    ),
    "objectives": ("false_acceptance", "maximum_truth_error", "evaluation_requests"),
    "tie_break": (
        "false_acceptance",
        "truth_error",
        "cost",
        "tau_rel",
        "kappa",
        "larger_extent",
        "smaller_extent",
        "literal_identity",
    ),
    "saturation": {
        "scale": (-8, 8),
        "finite_difference": (2.0**-8, 2.0**12, 16, 16),
        "rank": 2.0**-16,
        "weak": 2.0**-2,
        "cluster": 2.0**-2,
        "conditioning": 2.0**52,
        "correlation": 2.0**12,
    },
    "rejected_neighbor": "only decisive Manhattan/Pareto neighbors with metric and reason",
}
PHASE_SELECTION = {
    "A0": {
        "truth_error": "max derivative-envelope utilization",
        "cost": "EvaluationEngine requests",
        "tie": "error,cost,abs(k),k",
        "edge": "all qualifiers at k=-8 or +8",
        "neighbors": "k-1,k+1",
    },
    "A1": {
        "truth_error": "max then sum derivative-envelope utilization",
        "cost": "stencil EvaluationEngine requests",
        "tie": "tau,kappa,E_plus,E_minus,identity",
        "edge": "tau=2^-8;kappa=2^12;E_minus=E_plus=16",
        "neighbors": "Manhattan grid neighbors",
    },
    "B1": {
        "truth_error": "max covariance/spectrum truth error",
        "cost": "SVD kernels",
        "tie": "gesvd before gesdd",
        "edge": None,
        "neighbors": "other driver",
    },
    "B2": {
        "truth_error": "rank misclassification count then margin error",
        "cost": "offline=0",
        "tie": "a_rank,rho_rank,identity",
        "edge": "a_rank or rho_rank=2^-16 only",
        "neighbors": "four Manhattan neighbors",
    },
    "C1": {
        "truth_error": "weak-mode misclassification count",
        "cost": "offline=0",
        "tie": "smallest omega",
        "edge": "omega=2^-2 only",
        "neighbors": "adjacent omega",
    },
    "C2": {
        "truth_error": "cluster/projector error",
        "cost": "offline=0",
        "tie": "smallest gamma",
        "edge": "gamma=2^-2 only",
        "neighbors": "adjacent gamma",
    },
    "C3": {
        "truth_error": "conditioning misclassification and relative error",
        "cost": "offline=0",
        "tie": "smallest adequate limit",
        "edge": "limit=2^52 only",
        "neighbors": "adjacent limit",
    },
    "C4": {
        "truth_error": "endpoint classification error",
        "cost": "offline=0",
        "tie": "smallest multiplier",
        "edge": "multiplier=2^12 only",
        "neighbors": "adjacent multiplier",
    },
}
PHASE_FAILURE_RULE = (
    "no qualifier, non-monotone/ambiguous qualifier set, or edge-only qualifier set "
    "stops the round; no later phase may alter an earlier selection"
)
EVIDENCE_FIELDS = (
    "selected_phase_policies",
    "composed_uncertainty_policy_identity",
    "qualification_scope",
    "decisive_candidate_metrics",
    "rejected_neighbors",
    "holdout_results",
    "canonical_provenance",
)

_SCIENTIFIC_MACHINERY_IMPORTED = False
UncertaintyPolicy: Any = None
migration_core_authority_selection: Any = None


def scale_candidates(q0: float) -> tuple[float, ...]:
    return tuple(q0 * 2.0**exponent for exponent in SCALE_EXPONENTS)


def finite_difference_policies() -> Iterator[tuple[float, float, int, int]]:
    return product(
        RELATIVE_STEP_TOLERANCES, ROUNDOFF_MULTIPLIERS, STEP_EXTENTS, STEP_EXTENTS
    )


def rank_candidates() -> Iterator[tuple[float, float]]:
    return product(RANK_TOLERANCES, RANK_TOLERANCES)


def characteristic_scale_family(name: str) -> tuple[object, ...]:
    for entry in SUPPORTED_SCALE_FAMILIES:
        members = cast("tuple[str, ...]", entry[4])
        if name in members or any(
            item.endswith("*") and name.startswith(item[:-1]) for item in members
        ):
            return entry
    if name in UNSUPPORTED_SCALE_FAMILIES:
        raise KeyError(f"unsupported characteristic-scale family: {name}")
    raise KeyError(f"characteristic-scale member not catalogued: {name}")


def advance_phase(
    selected: Mapping[str, object], phase: str, candidate: object
) -> dict[str, object]:
    order = tuple(item[0] for item in PHASE_OWNERSHIP)
    if (
        phase not in order
        or phase in selected
        or tuple(selected) != order[: len(selected)]
        or phase != order[len(selected)]
    ):
        raise RuntimeError(f"phase {phase} cannot retune or skip staged ownership")
    return {**selected, phase: candidate}


def phase_decision(
    qualified: list[Mapping[str, object]], *, edge_only: bool, ambiguous: bool
) -> str:
    if not qualified:
        return "UNSUPPORTED"
    if edge_only:
        return "GRID_SATURATED"
    if ambiguous:
        return "AMBIGUOUS_NON_MONOTONE"
    return "SELECTABLE"


def holdout_decision(passed: tuple[bool, ...]) -> str:
    if len(passed) != 6:
        raise ValueError("holdout validation requires exactly six frozen cases")
    return "QUALIFIED" if all(passed) else "HOLDOUT_FAILED_POLICY_UNAVAILABLE"


def _number(record: Mapping[str, object], name: str) -> float:
    value = record[name]
    if not isinstance(value, int | float):
        raise TypeError(f"candidate {name} must be numeric")
    return float(value)


def _dominates(left: Mapping[str, object], right: Mapping[str, object]) -> bool:
    left_pair = (_number(left, "truth_error"), _number(left, "cost"))
    right_pair = (_number(right, "truth_error"), _number(right, "cost"))
    return (
        left_pair[0] <= right_pair[0]
        and left_pair[1] <= right_pair[1]
        and left_pair != right_pair
    )


def reduce_candidates(candidates: Iterable[Mapping[str, object]]) -> dict[str, object]:
    """Keep a compact frontier and bounded decisive rejection provenance."""
    frontier: list[Mapping[str, object]] = []
    rejected: list[tuple[object, tuple[str, ...]]] = []
    count = 0
    for candidate in candidates:
        count += 1
        identity = candidate["candidate"]
        reasons = tuple(cast("Iterable[str]", candidate["reasons"]))
        if not candidate["qualified"]:
            if len(rejected) < 16:
                rejected.append((identity, reasons))
            continue
        if any(_dominates(item, candidate) for item in frontier):
            if len(rejected) < 16:
                rejected.append((identity, ("pareto_dominated",)))
            continue
        equivalent = next(
            (
                item
                for item in frontier
                if _number(item, "truth_error") == _number(candidate, "truth_error")
                and _number(item, "cost") == _number(candidate, "cost")
            ),
            None,
        )
        if equivalent is not None:
            existing_key = repr(equivalent.get("tie_break", equivalent["candidate"]))
            candidate_key = repr(candidate.get("tie_break", candidate["candidate"]))
            kept, dropped = (
                (candidate, equivalent)
                if candidate_key < existing_key
                else (equivalent, candidate)
            )
            frontier[frontier.index(equivalent)] = kept
            if len(rejected) < 16:
                rejected.append((dropped["candidate"], ("deterministic_tie_break",)))
            continue
        displaced = [item for item in frontier if _dominates(candidate, item)]
        rejected.extend(
            (item["candidate"], ("pareto_dominated",))
            for item in displaced[: max(0, 16 - len(rejected))]
        )
        frontier = [item for item in frontier if item not in displaced]
        frontier.append(candidate)
    frontier.sort(
        key=lambda item: (
            _number(item, "truth_error"),
            _number(item, "cost"),
            repr(item["candidate"]),
        )
    )
    return {
        "selected": None if not frontier else frontier[0]["candidate"],
        "considered_count": count,
        "frontier": tuple(
            (
                item["candidate"],
                _number(item, "truth_error"),
                int(_number(item, "cost")),
            )
            for item in frontier
        ),
        "rejected_neighbors": tuple(rejected),
    }


def _lane(role: str) -> NumericalLane:
    if role == "canonical":
        return canonical_lanes()[0]
    if role == "compatibility":
        return canonical_lanes()[1]
    raise ValueError("lane role must be canonical or compatibility")


def reconstruct_lane_records(
    records: Mapping[str, object], role: str
) -> tuple[NumericalLane, LaneAttestation, RuntimeEnvironment]:
    names = {"numerical_lane", "lane_attestation", "runtime_environment"}
    if set(records) != names or any(
        not isinstance(records[name], Mapping) for name in names
    ):
        raise RuntimeError("calibration requires complete typed lane records")
    lane = NumericalLane.from_record(
        cast("Mapping[str, object]", records["numerical_lane"])
    )
    attestation = LaneAttestation.from_record(
        cast("Mapping[str, object]", records["lane_attestation"])
    )
    environment = RuntimeEnvironment.from_record(
        cast("Mapping[str, object]", records["runtime_environment"])
    )
    if (
        lane != _lane(role)
        or attestation.lane_identity != lane.identity
        or attestation.environment_identity != environment.identity
        or environment.semantics != lane.semantics
    ):
        raise RuntimeError("calibration lane records are internally inconsistent")
    return lane, attestation, environment


def attest_lane(role: str, image_digest: str) -> dict[str, object]:
    lane = _lane(role)
    authority = lane.attest_current_process(image_digest)
    records = {
        "numerical_lane": lane.to_record(),
        "lane_attestation": authority.to_record(),
        "runtime_environment": RuntimeEnvironment(lane.semantics).to_record(),
    }
    reconstruct_lane_records(records, role)
    return records


def _import_scientific_machinery() -> None:
    global UncertaintyPolicy, migration_core_authority_selection
    global _SCIENTIFIC_MACHINERY_IMPORTED
    if _SCIENTIFIC_MACHINERY_IMPORTED:
        return
    from chemex.migration_core import migration_core_authority_selection
    from chemex.optimize.uncertainty import UncertaintyPolicy

    _SCIENTIFIC_MACHINERY_IMPORTED = True


def validate_lane_records(records: Mapping[str, object], role: str) -> None:
    lane, attestation, environment = reconstruct_lane_records(records, role)
    selection = migration_core_authority_selection()
    if role == "canonical" and (
        lane.identity != selection.lane_identity
        or attestation.identity != selection.attestation_identity
        or environment.identity != selection.environment_identity
    ):
        raise RuntimeError(
            "canonical calibration lane does not match live #588 authority"
        )


def prospective_record(
    role: str, selected_policy: Mapping[str, object] | None = None
) -> dict[str, object]:
    if role == "compatibility" and selected_policy is None:
        raise RuntimeError("compatibility replay requires a selected canonical policy")
    return {
        "specification_id": SPECIFICATION_ID,
        "status": "FROZEN_AWAITING_INDEPENDENT_REVIEW",
        "lane_role": role,
        "candidate_counts": CANDIDATE_COUNTS,
        "phase_ownership": PHASE_OWNERSHIP,
        "acquisition_plan": ACQUISITION_PLAN[role],
        "resource_ledger": RESOURCE_LEDGER,
        "selected_policy_input": None
        if selected_policy is None
        else dict(selected_policy),
    }


def acquire(
    role: str, image_digest: str, selected_policy: Mapping[str, object] | None = None
) -> dict[str, object]:
    if role not in {"canonical", "compatibility"}:
        raise ValueError("lane role must be canonical or compatibility")
    if role == "compatibility" and selected_policy is None:
        raise RuntimeError("compatibility replay requires a selected canonical policy")
    lane_records = attest_lane(role, image_digest)
    _import_scientific_machinery()
    validate_lane_records(lane_records, role)
    return {
        **prospective_record(role, selected_policy),
        "canonical_provenance": lane_records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--role", choices=("canonical", "compatibility"), required=True)
    parser.add_argument("--image-digest", required=True)
    parser.add_argument("--selected-policy", type=Path)
    arguments = parser.parse_args()
    selected = (
        None
        if arguments.selected_policy is None
        else cast(
            "Mapping[str, object]",
            json.loads(arguments.selected_policy.read_text(encoding="utf-8")),
        )
    )
    record = acquire(arguments.role, arguments.image_digest, selected)
    arguments.output.write_text(
        json.dumps(record, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

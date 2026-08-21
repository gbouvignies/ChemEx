from __future__ import annotations

from collections.abc import Callable

import pytest

from tests.qualification.compare_lmfit_removal_pr1 import (
    _BINDING_ACCEPTED_FAILURE_SCOPE,
    _BINDING_ACCEPTED_NUMERICAL_SIGNATURE,
    _CEST_ACCEPTED_FAILURE_SCOPE,
    _CEST_ACCEPTED_NUMERICAL_SIGNATURE,
    CASES,
    _comparison_disposition,
    _data_rows,
    _fallback_limit,
    _fit_schema,
    _parameter_values,
    _statistics_schema,
)


def test_population_fallback_uses_declared_absolute_tolerance() -> None:
    case = CASES[0]

    assert _fallback_limit(("[GLOBAL]", "PA"), 0.93, case) == pytest.approx(
        case.fallback.population_abs
    )


def test_data_reader_preserves_profile_order_and_masked_rows() -> None:
    profiles, rows = _data_rows(
        """[G2N-HN]
  1.0  2.0  0.5  2.1
# 3.0  4.0  0.6  4.1 # NOT USED IN THE FIT
"""
    )

    assert profiles == ("[G2N-HN]",)
    assert tuple((row.profile, row.masked, row.values) for row in rows) == (
        ("[G2N-HN]", False, (1.0, 2.0, 0.5, 2.1)),
        ("[G2N-HN]", True, (3.0, 4.0, 0.6, 4.1)),
    )


def test_fit_and_statistics_schema_preserve_structure_without_numerical_parity() -> (
    None
):
    legacy_fit = "[G2N-HN]\n# TIME CALC\n1.0 2.0\n3.0 4.0\n"
    native_fit = "[G2N-HN]\n# TIME CALC\n1.0 20.0\n3.0 40.0\n"

    assert _fit_schema(native_fit) == _fit_schema(legacy_fit)
    assert _statistics_schema('"number" = 1\n"score" = 2.0\n') == (
        ("number", "int"),
        ("score", "float"),
    )


def test_cest_accepted_difference_is_explicit_and_not_a_parity_pass() -> None:
    case = next(case for case in CASES if case.slug == "cest-13c-label-cn")

    disposition = _comparison_disposition(
        case,
        False,
        {"STEP2/All/statistics.toml": (3140.90, 2849.45)},
        _CEST_ACCEPTED_FAILURE_SCOPE,
        dict(_CEST_ACCEPTED_NUMERICAL_SIGNATURE),
    )

    assert disposition["status"] == "ACCEPTED_DIFFERENCE"
    assert disposition["parity_status"] == "FAIL"
    assert disposition["objective_semantics_parity"] == "PASS"
    assert disposition["diagnosed_legacy_l18cd1_chi_square"] == pytest.approx(
        1093.424072389286
    )
    assert disposition["diagnosed_native_l18cd1_chi_square"] == pytest.approx(
        801.9806529406368
    )


def test_cest_accepted_difference_rejects_an_unrecorded_objective() -> None:
    case = next(case for case in CASES if case.slug == "cest-13c-label-cn")

    disposition = _comparison_disposition(
        case,
        False,
        {"STEP2/All/statistics.toml": (3140.90, 2850.00)},
        _CEST_ACCEPTED_FAILURE_SCOPE,
        dict(_CEST_ACCEPTED_NUMERICAL_SIGNATURE),
    )

    assert disposition == {"status": "FAIL"}


def test_cest_accepted_difference_rejects_an_unrelated_parity_failure() -> None:
    case = next(case for case in CASES if case.slug == "cest-13c-label-cn")

    disposition = _comparison_disposition(
        case,
        False,
        {"STEP2/All/statistics.toml": (3140.90, 2849.45)},
        _CEST_ACCEPTED_FAILURE_SCOPE | {"data:STEP1/Data/23hz.dat:[I28CG2]"},
        dict(_CEST_ACCEPTED_NUMERICAL_SIGNATURE),
    )

    assert disposition == {"status": "FAIL"}


def test_cest_accepted_difference_rejects_an_unrelated_numerical_state() -> None:
    case = next(case for case in CASES if case.slug == "cest-13c-label-cn")
    signature = dict(_CEST_ACCEPTED_NUMERICAL_SIGNATURE)
    signature["native_dw_ab"] = 0.3

    disposition = _comparison_disposition(
        case,
        False,
        {"STEP2/All/statistics.toml": (3140.90, 2849.45)},
        _CEST_ACCEPTED_FAILURE_SCOPE,
        signature,
    )

    assert disposition == {"status": "FAIL"}


def test_binding_accepted_difference_is_explicit_and_not_a_parity_pass() -> None:
    case = next(case for case in CASES if case.slug == "2st-binding")

    disposition = _comparison_disposition(
        case,
        False,
        {"STEP2/All/statistics.toml": (37771.2, 35802.6)},
        _BINDING_ACCEPTED_FAILURE_SCOPE,
        dict(_BINDING_ACCEPTED_NUMERICAL_SIGNATURE),
    )

    assert disposition["status"] == "ACCEPTED_DIFFERENCE"
    assert disposition["parity_status"] == "FAIL"
    assert disposition["objective_semantics_parity"] == "PASS"
    assert disposition["diagnosed_legacy_chi_square"] == pytest.approx(
        37771.20405101328
    )
    assert disposition["diagnosed_native_chi_square"] == pytest.approx(
        35802.60988668793
    )


def test_binding_accepted_difference_rejects_a_step1_failure() -> None:
    case = next(case for case in CASES if case.slug == "2st-binding")

    disposition = _comparison_disposition(
        case,
        False,
        {"STEP2/All/statistics.toml": (37771.2, 35802.6)},
        _BINDING_ACCEPTED_FAILURE_SCOPE | {"data:STEP1/All/Data/cpmg_10p.dat:[513N]"},
        dict(_BINDING_ACCEPTED_NUMERICAL_SIGNATURE),
    )

    assert disposition == {"status": "FAIL"}


@pytest.mark.parametrize(
    "unexpected_failure",
    (
        "data:STEP2/All/Data/cest_20hz_10p.dat:[575N]",
        ("parameter:STEP2/All/Parameters/fitted.toml:('[DW_AB]', '575N')"),
    ),
)
def test_binding_accepted_difference_rejects_an_unrelated_step2_failure(
    unexpected_failure: str,
) -> None:
    case = next(case for case in CASES if case.slug == "2st-binding")

    disposition = _comparison_disposition(
        case,
        False,
        {"STEP2/All/statistics.toml": (37771.2, 35802.6)},
        _BINDING_ACCEPTED_FAILURE_SCOPE | {unexpected_failure},
        dict(_BINDING_ACCEPTED_NUMERICAL_SIGNATURE),
    )

    assert disposition == {"status": "FAIL"}


def test_binding_accepted_difference_rejects_a_changed_objective() -> None:
    case = next(case for case in CASES if case.slug == "2st-binding")

    disposition = _comparison_disposition(
        case,
        False,
        {"STEP2/All/statistics.toml": (37771.2, 35803.0)},
        _BINDING_ACCEPTED_FAILURE_SCOPE,
        dict(_BINDING_ACCEPTED_NUMERICAL_SIGNATURE),
    )

    assert disposition == {"status": "FAIL"}


def test_binding_accepted_difference_rejects_an_unrelated_numerical_state() -> None:
    case = next(case for case in CASES if case.slug == "2st-binding")
    signature = dict(_BINDING_ACCEPTED_NUMERICAL_SIGNATURE)
    signature["native_560n_dw_ab"] = -0.08

    disposition = _comparison_disposition(
        case,
        False,
        {"STEP2/All/statistics.toml": (37771.2, 35802.6)},
        _BINDING_ACCEPTED_FAILURE_SCOPE,
        signature,
    )

    assert disposition == {"status": "FAIL"}


def test_successful_parity_remains_pass_parity() -> None:
    case = next(case for case in CASES if case.slug == "2st-binding")

    disposition = _comparison_disposition(case, True, {}, frozenset(), {})

    assert disposition == {"status": "PASS_PARITY"}


@pytest.mark.parametrize(
    "reader,text",
    (
        (_data_rows, "[G2N-HN]\n1.0 2.0 0.5 nan\n"),
        (_parameter_values, "[GLOBAL]\nPB = nan\n"),
    ),
)
def test_comparison_readers_reject_non_finite_values(
    reader: Callable[[str], object],
    text: str,
) -> None:
    with pytest.raises(AssertionError, match="non-finite"):
        reader(text)

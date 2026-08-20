from __future__ import annotations

from collections.abc import Callable

import pytest

from tests.qualification.compare_lmfit_removal_pr1 import (
    CASES,
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

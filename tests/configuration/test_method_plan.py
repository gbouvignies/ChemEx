from __future__ import annotations

from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from chemex.configuration.method_input import prepare_method_plan
from chemex.configuration.method_plan import (
    ConstrainAction,
    DeSearch,
    FitAction,
    FixAction,
    FormatOrigin,
    MethodFormatError,
    ProfileSelection,
)
from chemex.configuration.method_validation import resolve_grid_axes
from chemex.configuration.methods import Method, Statistics, read_method_plan
from chemex.parameters.parameterization import (
    ParameterDeclaration,
    SealedParameterDeclarations,
    SealedParameterModel,
)
from chemex.parameters.sealed import (
    ParamConfig,
    ParamDefinition,
    SealedConfiguration,
    SealedDefinitions,
)

ROOT = Path(__file__).parents[2]


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _parameter_model(
    *definitions: ParamDefinition,
    supports_estimation: bool = True,
    model_expression: str = "",
    requires_independent: bool | None = None,
    fits_by_default: bool | None = None,
) -> SealedParameterModel:
    sealed_definitions = SealedDefinitions(tuple(definitions), {})
    configuration = SealedConfiguration(
        tuple(
            ParamConfig(
                definition.param_id,
                definition.default_value,
                definition.lower_bound,
                definition.upper_bound,
            )
            for definition in definitions
        ),
        {},
        sealed_definitions.identity,
    )
    declarations = SealedParameterDeclarations(
        tuple(
            ParameterDeclaration(
                definition.param_id,
                supports_estimation,
                model_expression,
                requires_independent=(
                    not bool(model_expression)
                    if requires_independent is None
                    else requires_independent
                ),
                fits_by_default=(
                    supports_estimation and not bool(model_expression)
                    if fits_by_default is None
                    else fits_by_default
                ),
            )
            for definition in definitions
        )
    )
    return SealedParameterModel(
        "2st",
        "model",
        sealed_definitions,
        configuration,
        declarations,
    )


def test_v1_and_v2_normalize_to_the_same_ordered_role_semantics(
    tmp_path: Path,
) -> None:
    v1 = _write(
        tmp_path / "v1.toml",
        """
[STEP]
CONSTRAINTS = ["[R2_B] = [R2_A]"]
FIX = ["PB"]
FIT = ["KEX_AB"]
""",
    )
    v2 = _write(
        tmp_path / "v2.toml",
        """
FORMAT_VERSION = 2

[STEP]
ROLES = [
  { CONSTRAIN = ["[R2_B] = [R2_A]"] },
  { FIX = ["PB"] },
  { FIT = ["KEX_AB"] },
]
""",
    )

    v1_plan = read_method_plan([v1])
    v2_plan = read_method_plan([v2])

    assert v1_plan.steps == v2_plan.steps


def test_v1_v2_and_programmatic_methods_prepare_equivalent_execution_semantics(
    tmp_path: Path,
) -> None:
    v1 = _write(
        tmp_path / "v1.toml",
        """
[STEP]
INCLUDE = ["1H"]

[STEP.STATISTICS]
MC = 2
""",
    )
    v2 = _write(
        tmp_path / "v2.toml",
        """
FORMAT_VERSION = 2

[STEP]
INCLUDE = ["1H"]

[STEP.STATISTICS.MC]
REPLICATES = 2
SEED = 0
""",
    )
    programmatic = {"STEP": Method(include=["1H"], statistics=Statistics(mc=2))}

    v1_plan = read_method_plan([v1])
    v2_plan = read_method_plan([v2])
    methods_plan = prepare_method_plan(programmatic, _parameter_model())

    assert v1_plan.steps == v2_plan.steps == methods_plan.steps


def test_prepare_method_plan_validates_canonical_input(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    plan = read_method_plan([])
    parameter_model = _parameter_model()
    validated: list[SealedParameterModel] = []
    monkeypatch.setattr(
        type(plan),
        "validate",
        lambda _plan, model: validated.append(model),
    )

    prepared = prepare_method_plan(plan, parameter_model)

    assert prepared is plan
    assert validated == [parameter_model]


@pytest.mark.parametrize("format_version", (1, 2))
def test_v1_and_v2_validate_same_group_companion_constraints(
    tmp_path: Path,
    format_version: int,
) -> None:
    model = _parameter_model(
        ParamDefinition(
            "r2adq-a", "R2ADQ_A", "I3HD1", (("h_larmor_frq", 700.2),), 0.0, 0.0, 100.0
        ),
        ParamDefinition(
            "r2dq-a", "R2DQ_A", "I3HD1", (("h_larmor_frq", 700.2),), 20.0, 0.0, 100.0
        ),
        ParamDefinition(
            "r1-a", "R1_A", "I3CD1", (("h_larmor_frq", 700.2),), 1.5, 0.0, 10.0
        ),
    )
    roles = (
        'CONSTRAINTS = ["[R2ADQ_A] = [R2DQ_A] - [R1_A] / 3.0"]'
        if format_version == 1
        else 'ROLES = [{ CONSTRAIN = ["[R2ADQ_A] = [R2DQ_A] - [R1_A] / 3.0"] }]'
    )
    version = "" if format_version == 1 else "FORMAT_VERSION = 2\n"
    method = _write(tmp_path / "method.toml", f"{version}[STEP]\n{roles}\n")

    read_method_plan([method]).validate(model)


@pytest.mark.parametrize("format_version", (1, 2))
def test_v1_and_v2_allow_supported_default_derivation_role_overrides(
    tmp_path: Path,
    format_version: int,
) -> None:
    model = _parameter_model(
        ParamDefinition("j-a", "J_A", "3N-H", (), -90.0, -120.0, -60.0),
        ParamDefinition("j-b", "J_B", "3N-H", (), -90.0, -120.0, -60.0),
        supports_estimation=True,
        model_expression="j-a",
        requires_independent=False,
    )
    if format_version == 1:
        contents = """[STEP1]
FIX = ["J_B"]

[STEP2]
FIT = ["J_B"]
"""
    else:
        contents = """FORMAT_VERSION = 2
[STEP1]
ROLES = [{ FIX = ["J_B"] }]

[STEP2]
ROLES_FROM = "STEP1"
ROLES = [{ FIT = ["J_B"] }]
"""
    method = _write(tmp_path / "method.toml", contents)

    read_method_plan([method]).validate(model)


@pytest.mark.parametrize("format_version", (1, 2))
def test_v1_and_v2_rank_spin_context_before_condition_specificity(
    tmp_path: Path,
    format_version: int,
) -> None:
    model = _parameter_model(
        ParamDefinition(
            "target",
            "TARGET",
            "I3HD1",
            (("h_larmor_frq", 800.0),),
            0.0,
            0.0,
            100.0,
        ),
        ParamDefinition("companion", "VALUE", "I3CD1", (), 1.0, 0.0, 10.0),
        ParamDefinition(
            "global",
            "VALUE",
            "",
            (("h_larmor_frq", 800.0),),
            2.0,
            0.0,
            10.0,
        ),
    )
    roles = (
        'CONSTRAINTS = ["[TARGET] = [VALUE]"]'
        if format_version == 1
        else 'ROLES = [{ CONSTRAIN = ["[TARGET] = [VALUE]"] }]'
    )
    version = "" if format_version == 1 else "FORMAT_VERSION = 2\n"
    method = _write(tmp_path / "method.toml", f"{version}[STEP]\n{roles}\n")

    read_method_plan([method]).validate(model)


def test_v1_and_v2_grid_statistics_normalize_to_the_same_step_semantics(
    tmp_path: Path,
) -> None:
    v1 = _write(
        tmp_path / "v1.toml",
        """
[STEP]
GRID = ["[PB] = (0.1, 0.2)"]
STATISTICS = { MC = 2 }
""",
    )
    v2 = _write(
        tmp_path / "v2.toml",
        """
FORMAT_VERSION = 2
[STEP]

[STEP.SEARCH.GRID]
AXES = ["[PB] = values(0.1, 0.2)"]

[STEP.STATISTICS.MC]
REPLICATES = 2
SEED = 0
""",
    )

    v1_step = read_method_plan([v1]).steps[0]
    v2_step = read_method_plan([v2]).steps[0]

    assert v1_step == v2_step


def test_v2_preserves_broad_to_specific_order_and_step_local_selection(
    tmp_path: Path,
) -> None:
    method = _write(
        tmp_path / "method.toml",
        """
FORMAT_VERSION = 2

[FIRST]
INCLUDE = [24]
ROLES = [
  { FIX = ["DW_AB"] },
  { FIT = ["DW_AB, NUC->G24N-H"] },
  { CONSTRAIN = ["[DW_AB, NUC->G31N-H] = [DW_AB, NUC->G24N-H]"] },
]

[SECOND]
ROLES_FROM = "FIRST"
""",
    )

    plan = read_method_plan([method])
    first, second = plan.steps

    assert tuple(type(action) for action in first.role_actions) == (
        FixAction,
        FitAction,
        ConstrainAction,
    )
    assert first.selection == ProfileSelection(include=("24",), exclude=None)
    assert second.roles_from == "FIRST"
    assert second.selection == ProfileSelection()


def test_v2_selector_is_strict_while_v1_preserves_partial_legacy_parsing(
    tmp_path: Path,
) -> None:
    v1 = _write(
        tmp_path / "legacy.toml",
        '[STEP]\nFIT = ["PB, TYPO->ignored"]\n',
    )
    v2 = _write(
        tmp_path / "strict.toml",
        """
FORMAT_VERSION = 2
[STEP]
ROLES = [{ FIT = ["PB, TYPO->ignored"] }]
""",
    )

    legacy_selector = read_method_plan([v1]).steps[0].role_actions[0].selectors[0]
    assert legacy_selector.render() == "PB"

    with pytest.raises(MethodFormatError) as error:
        read_method_plan([v2])

    diagnostic = str(error.value)
    assert str(v2) in diagnostic
    assert "[STEP]" in diagnostic
    assert "ROLES[0].FIT[0]" in diagnostic
    assert "TYPO->ignored" in diagnostic


def test_constraint_equations_compile_to_typed_ast_and_render_canonically(
    tmp_path: Path,
) -> None:
    method = _write(
        tmp_path / "constraint.toml",
        """
FORMAT_VERSION = 2
[STEP]
ROLES = [
  { CONSTRAIN = [
    "[r2_b, b0->800.03mhz] = -[r1_a] + 2 * ([r2_a] - 1) / 3"
  ] },
]
""",
    )

    action = read_method_plan([method]).steps[0].role_actions[0]
    constraint = action.constraints[0]

    assert (
        constraint.render()
        == "[R2_B, B0->800.0MHz] = -[R1_A] + 2.0 * ([R2_A] - 1.0) / 3.0"
    )


def test_v1_grid_tuple_and_v2_values_normalize_to_the_same_search(
    tmp_path: Path,
) -> None:
    v1 = _write(
        tmp_path / "grid-v1.toml",
        '[STEP]\nGRID = ["[PB] = (0.01, 0.02, 0.05)"]\n',
    )
    v2 = _write(
        tmp_path / "grid-v2.toml",
        """
FORMAT_VERSION = 2
[STEP.SEARCH.GRID]
AXES = ["[PB] = values(0.01, 0.02, 0.05)"]
""",
    )

    assert (
        read_method_plan([v1]).steps[0].search == read_method_plan([v2]).steps[0].search
    )
    assert (
        read_method_plan([v2]).steps[0].search.axes[0].render()
        == "[PB] = values(0.01, 0.02, 0.05)"
    )


def test_v2_de_parses_only_seeded_two_argument_coordinate_ranges(
    tmp_path: Path,
) -> None:
    method = _write(
        tmp_path / "de.toml",
        """
FORMAT_VERSION = 2
[STEP.SEARCH.DE]
SEED = 597
COORDINATES = [
  "[PB] = log(0.01, 0.20)",
  "[KEX_AB] = lin(100, 5000)",
]
""",
    )

    search = read_method_plan([method]).steps[0].search

    assert isinstance(search, DeSearch)
    assert search.seed == 597
    assert tuple(coordinate.render() for coordinate in search.coordinates) == (
        "[PB] = log(0.01, 0.2)",
        "[KEX_AB] = lin(100.0, 5000.0)",
    )

    invalid = _write(
        tmp_path / "invalid-de.toml",
        """
FORMAT_VERSION = 2
[STEP.SEARCH.DE]
SEED = 1
COORDINATES = ["[PB] = values(0.01, 0.20)"]
""",
    )
    with pytest.raises(MethodFormatError, match="DE accepts only lin"):
        read_method_plan([invalid])


def test_v2_compact_and_expanded_statistics_normalize_to_typed_requests(
    tmp_path: Path,
) -> None:
    compact = _write(
        tmp_path / "compact.toml",
        """
FORMAT_VERSION = 2
[STEP]
STATISTICS = { MC = 100, MCMC = 5000 }
""",
    )
    expanded = _write(
        tmp_path / "expanded.toml",
        """
FORMAT_VERSION = 2
[STEP.STATISTICS.MC]
REPLICATES = 100

[STEP.STATISTICS.MCMC]
STEPS = 5000
""",
    )

    assert (
        read_method_plan([compact]).steps[0].statistics
        == read_method_plan([expanded]).steps[0].statistics
    )

    seeded = _write(
        tmp_path / "seeded.toml",
        """
FORMAT_VERSION = 2
[STEP.STATISTICS.BS]
REPLICATES = 20
SEED = 7

[STEP.STATISTICS.BSN]
REPLICATES = 30

[STEP.STATISTICS.MCMC]
STEPS = 5000
BURN = 1000
SEED = 9
""",
    )
    statistics = read_method_plan([seeded]).steps[0].statistics
    assert statistics is not None
    assert statistics.bs is not None
    assert (statistics.bs.replicates, statistics.bs.seed) == (20, 7)
    assert statistics.bsn is not None
    assert (statistics.bsn.replicates, statistics.bsn.seed) == (30, None)
    assert statistics.mcmc is not None
    assert (
        statistics.mcmc.steps,
        statistics.mcmc.burn,
        statistics.mcmc.seed,
    ) == (5000, 1000, 9)

    invalid = _write(
        tmp_path / "invalid-mcmc.toml",
        """
FORMAT_VERSION = 2
[STEP.STATISTICS.MCMC]
STEPS = 5000
THIN = 2
""",
    )
    with pytest.raises(MethodFormatError, match="THIN"):
        read_method_plan([invalid])


def test_multifile_loading_rejects_mixed_formats_and_duplicate_v2_steps(
    tmp_path: Path,
) -> None:
    v1 = _write(tmp_path / "v1.toml", "[ONE]\n")
    v2 = _write(tmp_path / "v2.toml", "FORMAT_VERSION = 2\n[TWO]\n")

    with pytest.raises(MethodFormatError, match="mixed v1/v2"):
        read_method_plan([v1, v2])

    duplicate = _write(
        tmp_path / "duplicate.toml",
        'FORMAT_VERSION = 2\n[TWO]\nROLES = [{ FIX = ["PB"] }]\n',
    )
    with pytest.raises(MethodFormatError, match="Duplicate v2 step.*TWO"):
        read_method_plan([v2, duplicate])

    non_integer_version = _write(
        tmp_path / "float-version.toml", "FORMAT_VERSION = 2.0\n[STEP]\n"
    )
    with pytest.raises(MethodFormatError, match="supports only version 2"):
        read_method_plan([non_integer_version])


def test_v2_rejects_invalid_role_actions_and_future_inheritance(
    tmp_path: Path,
) -> None:
    invalid_action = _write(
        tmp_path / "invalid-action.toml",
        """
FORMAT_VERSION = 2
[STEP]
ROLES = [{ FIX = ["PB"], FIT = ["KEX_AB"] }]
""",
    )
    with pytest.raises(MethodFormatError, match="exactly one.*FIX.*FIT"):
        read_method_plan([invalid_action])

    empty_roles = _write(
        tmp_path / "empty-roles.toml",
        "FORMAT_VERSION = 2\n[STEP]\nROLES = []\n",
    )
    with pytest.raises(MethodFormatError, match="at least one action"):
        read_method_plan([empty_roles])

    future = _write(
        tmp_path / "future.toml",
        """
FORMAT_VERSION = 2
[FIRST]
ROLES_FROM = "SECOND"
[SECOND]
""",
    )
    with pytest.raises(MethodFormatError, match="earlier step"):
        read_method_plan([future])

    fitmethod = _write(
        tmp_path / "fitmethod.toml",
        'FORMAT_VERSION = 2\n[STEP]\nFITMETHOD = "trf"\n',
    )
    with pytest.raises(MethodFormatError, match="V2 has no FITMETHOD.*TRF is implicit"):
        read_method_plan([fitmethod])


def test_v1_synthesizes_role_and_selection_carry_forward_and_keeps_mcmc_controls(
    tmp_path: Path,
) -> None:
    method = _write(
        tmp_path / "legacy-multistep.toml",
        """
[FIRST]
INCLUDE = [24]
FITMETHOD = "least_squares"
FIX = ["PB"]
[FIRST.STATISTICS.MCMC]
STEPS = 5000
BURN = "AUTO"
THIN = 5
WALKERS = 32
SEED = 9
WORKERS = 2

[SECOND]
FIT = ["PB"]
""",
    )

    first, second = read_method_plan([method]).steps

    assert second.roles_from == "FIRST"
    assert second.selection == first.selection
    assert first.statistics is not None
    assert first.statistics.mcmc is not None
    assert first.statistics.mcmc.burn is None
    assert first.statistics.mcmc.thin == 5
    assert first.statistics.mcmc.walkers == 32
    assert first.statistics.mcmc.seed == 9
    assert first.statistics.mcmc.workers == 2
    assert "# V1_ONLY_THIN = 5" in read_method_plan([method]).render()


def test_v1_adapter_resolves_resampling_seed_zero_but_not_mcmc_seed(
    tmp_path: Path,
) -> None:
    method = _write(
        tmp_path / "legacy-statistics.toml",
        """
[STEP]
STATISTICS = { MC = 2, BS = 3, BSN = 4, MCMC = 5 }
""",
    )

    statistics = read_method_plan([method]).steps[0].statistics

    assert statistics is not None
    assert statistics.mc is not None and statistics.mc.seed == 0
    assert statistics.bs is not None and statistics.bs.seed == 0
    assert statistics.bsn is not None and statistics.bsn.seed == 0
    assert statistics.mcmc is not None and statistics.mcmc.seed is None


def test_v1_adapter_preserves_values_coerced_by_the_legacy_schema(
    tmp_path: Path,
) -> None:
    method = _write(
        tmp_path / "legacy-coercions.toml",
        """
[STEP]
INCLUDE = [24]
STATISTICS = { MC = "10", MCMC = { STEPS = "5000", THIN = "5" } }
""",
    )

    step = read_method_plan([method]).steps[0]

    assert step.selection == ProfileSelection(include=("24",), exclude=None)
    assert step.statistics is not None
    assert step.statistics.mc is not None
    assert step.statistics.mc.replicates == 10
    assert step.statistics.mcmc is not None
    assert (step.statistics.mcmc.steps, step.statistics.mcmc.thin) == (5000, 5)


def test_canonical_rendering_is_deterministic_valid_v2_toml(
    tmp_path: Path,
) -> None:
    source = _write(
        tmp_path / "render-source.toml",
        """
FORMAT_VERSION = 2
[STEP]
INCLUDE = [24]
ROLES = [{ FIT = ["pb", "kex_ab"] }]
STATISTICS = { MC = 10 }
[STEP.SEARCH.GRID]
AXES = ["[PB] = LOG(0.01, 0.2, 3)"]
""",
    )
    plan = read_method_plan([source])

    rendered = plan.render()
    normalized = _write(tmp_path / "normalized.toml", rendered)

    assert rendered == plan.render()
    assert read_method_plan([normalized]).steps == plan.steps
    assert "FORMAT_VERSION = 2" in rendered
    assert '{ FIT = ["PB", "KEX_AB"] }' in rendered
    assert '"[PB] = log(0.01, 0.2, 3)"' in rendered
    assert "[STEP.STATISTICS.MC]\nREPLICATES = 10" in rendered


def test_de_validation_requires_one_final_independent_fit_coordinate(
    tmp_path: Path,
) -> None:
    model = _parameter_model(
        ParamDefinition("pb-600", "PB", "", (("h_larmor_frq", 600.0),), 0.1, 0.0, 1.0),
        ParamDefinition("pb-800", "PB", "", (("h_larmor_frq", 800.0),), 0.1, 0.0, 1.0),
    )
    broad = _write(
        tmp_path / "broad-de.toml",
        """
FORMAT_VERSION = 2
[STEP]
ROLES = [{ FIT = ["PB"] }]
[STEP.SEARCH.DE]
SEED = 1
COORDINATES = ["[PB] = lin(0.01, 0.2)"]
""",
    )
    with pytest.raises(MethodFormatError, match="exactly one.*matched 2"):
        read_method_plan([broad]).validate(model)

    exact = _write(
        tmp_path / "exact-de.toml",
        """
FORMAT_VERSION = 2
[STEP]
ROLES = [{ FIT = ["PB"] }]
[STEP.SEARCH.DE]
SEED = 1
COORDINATES = ["[PB, B0->600MHz] = lin(0.01, 0.2)"]
""",
    )
    read_method_plan([exact]).validate(model)


def test_de_validation_rejects_ranges_outside_physical_bounds(tmp_path: Path) -> None:
    model = _parameter_model(
        ParamDefinition("pb", "PB", "", (), 0.1, 0.0, 1.0),
    )
    method = _write(
        tmp_path / "out-of-bounds-de.toml",
        """
FORMAT_VERSION = 2
[STEP]
ROLES = [{ FIT = ["PB"] }]
[STEP.SEARCH.DE]
SEED = 597
COORDINATES = ["[PB] = lin(0.1, 1.1)"]
""",
    )

    with pytest.raises(MethodFormatError, match="outside physical bounds"):
        read_method_plan([method]).validate(model)


def test_representative_multistep_v1_and_v2_workflows_are_structurally_equivalent(
    tmp_path: Path,
) -> None:
    v1 = _write(
        tmp_path / "workflow-v1.toml",
        """
[FIRST]
INCLUDE = [24]
CONSTRAINTS = ["[R2_B] = [R2_A]"]
FIX = ["PB"]
GRID = ["[KEX_AB] = (100, 500)"]

[SECOND]
FIT = ["PB", "KEX_AB"]
STATISTICS = { MC = 10, BS = 20, BSN = 30, MCMC = 5000 }
""",
    )
    v2 = _write(
        tmp_path / "workflow-v2.toml",
        """
FORMAT_VERSION = 2
[FIRST]
INCLUDE = [24]
ROLES = [
  { CONSTRAIN = ["[R2_B] = [R2_A]"] },
  { FIX = ["PB"] },
]
[FIRST.SEARCH.GRID]
AXES = ["[KEX_AB] = values(100, 500)"]

[SECOND]
INCLUDE = [24]
ROLES_FROM = "FIRST"
ROLES = [{ FIT = ["PB", "KEX_AB"] }]

[SECOND.STATISTICS.MC]
REPLICATES = 10
SEED = 0

[SECOND.STATISTICS.BS]
REPLICATES = 20
SEED = 0

[SECOND.STATISTICS.BSN]
REPLICATES = 30
SEED = 0

[SECOND.STATISTICS.MCMC]
STEPS = 5000
""",
    )

    assert read_method_plan([v1]).steps == read_method_plan([v2]).steps


def test_all_shipped_methods_use_canonical_v2_and_parse() -> None:
    methods = sorted((ROOT / "examples").glob("*/*/Methods/*.toml"))

    assert methods
    for method in methods:
        plan = read_method_plan([method])
        assert plan.format_origin is FormatOrigin.V2
        assert plan.render() == method.read_text(encoding="utf-8")

    de_example = next(path for path in methods if path.name == "method_de.toml")
    assert isinstance(read_method_plan([de_example]).steps[0].search, DeSearch)


def test_method_plan_is_deeply_immutable(tmp_path: Path) -> None:
    method = _write(
        tmp_path / "immutable.toml",
        'FORMAT_VERSION = 2\n[STEP]\nROLES = [{ FIT = ["PB"] }]\n',
    )
    plan = read_method_plan([method])

    with pytest.raises(FrozenInstanceError):
        plan.steps[0].name = "CHANGED"  # type: ignore[misc]


def test_constraint_validation_rejects_final_dependency_cycles(tmp_path: Path) -> None:
    model = _parameter_model(
        ParamDefinition("pb", "PB", "", (), 0.1, 0.0, 1.0),
        ParamDefinition("kex", "KEX_AB", "", (), 500.0, 0.0, 5000.0),
    )
    method = _write(
        tmp_path / "cycle.toml",
        """
FORMAT_VERSION = 2
[STEP]
ROLES = [{ CONSTRAIN = [
  "[PB] = [KEX_AB]",
  "[KEX_AB] = [PB]",
] }]
""",
    )

    with pytest.raises(MethodFormatError, match="cycle.*pb.*kex"):
        read_method_plan([method]).validate(model)

    overridden = _write(
        tmp_path / "cycle-overridden.toml",
        method.read_text(encoding="utf-8").replace(
            "] }]",
            '] }, { FIT = ["PB", "KEX_AB"] }]',
        ),
    )
    read_method_plan([overridden]).validate(model)


def test_constraint_resolution_errors_point_to_the_reference_span(
    tmp_path: Path,
) -> None:
    model = _parameter_model(
        ParamDefinition("pb", "PB", "", (), 0.1, 0.0, 1.0),
    )
    equation = "[PB] = [MISSING]"
    method = _write(
        tmp_path / "missing-reference.toml",
        f'FORMAT_VERSION = 2\n[STEP]\nROLES = [{{ CONSTRAIN = ["{equation}"] }}]\n',
    )

    with pytest.raises(MethodFormatError) as error:
        read_method_plan([method]).validate(model)

    assert error.value.source.start == equation.index("MISSING")
    assert error.value.source.end == equation.index("MISSING") + len("MISSING")


def test_fit_cannot_override_a_structurally_unestimable_parameter(
    tmp_path: Path,
) -> None:
    model = _parameter_model(
        ParamDefinition("pb", "PB", "", (), 0.1, 0.0, 1.0),
        supports_estimation=False,
        model_expression="1.0",
    )
    method = _write(
        tmp_path / "unestimable.toml",
        'FORMAT_VERSION = 2\n[STEP]\nROLES = [{ FIT = ["PB"] }]\n',
    )

    with pytest.raises(MethodFormatError, match="do not support estimation"):
        read_method_plan([method]).validate(model)


def test_fit_cannot_override_an_independent_unsupported_parameter(
    tmp_path: Path,
) -> None:
    model = _parameter_model(
        ParamDefinition("pb", "PB", "", (), 0.1, 0.0, 1.0),
        supports_estimation=False,
        requires_independent=True,
    )
    method = _write(
        tmp_path / "unsupported.toml",
        'FORMAT_VERSION = 2\n[STEP]\nROLES = [{ FIT = ["PB"] }]\n',
    )

    with pytest.raises(MethodFormatError, match="do not support estimation"):
        read_method_plan([method]).validate(model)


@pytest.mark.parametrize(
    "selector, message",
    (
        ("PB, B0->600MHz, B0->800MHz", "Duplicate selector qualifier B0"),
        ("PB, [P]->1.0MHz", "Invalid value or unit"),
        ("PB trailing", "Invalid parameter-name token"),
        ("PB, D2O->1.0", "between 0 and 1"),
        ("PB, NUC->G24N-H-C-X", "Invalid NUC spin-system selector"),
        ("PB, NUC->G24N/H", "Invalid NUC spin-system selector"),
    ),
)
def test_v2_selector_diagnostics_reject_the_complete_invalid_input(
    tmp_path: Path,
    selector: str,
    message: str,
) -> None:
    method = _write(
        tmp_path / "invalid-selector.toml",
        f'FORMAT_VERSION = 2\n[STEP]\nROLES = [{{ FIT = ["{selector}"] }}]\n',
    )

    with pytest.raises(MethodFormatError, match=message):
        read_method_plan([method])


@pytest.mark.parametrize(
    "spin_system",
    ("G23C1'", "G23H1'", "G23H5'1", "G23C1'-H1'"),
)
def test_v2_nuc_selector_accepts_and_round_trips_prime_atom_names(
    tmp_path: Path,
    spin_system: str,
) -> None:
    method = _write(
        tmp_path / "prime-nucleus.toml",
        f"""FORMAT_VERSION = 2
[STEP]
ROLES = [{{ FIT = ["DW_AB, NUC->{spin_system}"] }}]
""",
    )

    plan = read_method_plan([method])
    selector = plan.steps[0].role_actions[0].selectors[0]
    normalized = _write(tmp_path / "prime-nucleus-normalized.toml", plan.render())
    round_tripped = read_method_plan([normalized]).steps[0].role_actions[0].selectors[0]

    assert selector.render() == f"DW_AB, NUC->{spin_system}"
    assert round_tripped == selector


def test_selector_rendering_preserves_round_trip_condition_values(
    tmp_path: Path,
) -> None:
    method = _write(
        tmp_path / "selector-rendering.toml",
        """
FORMAT_VERSION = 2
[STEP]
ROLES = [{ FIT = ["PB, [P]->0.00123456789M, [L]->2.3456789e-6M, D2O->0.123456789"] }]
""",
    )

    selector = read_method_plan([method]).steps[0].role_actions[0].selectors[0]

    assert selector.render() == (
        "PB, [P]->0.00123456789M, [L]->2.3456789e-06M, D2O->0.123456789"
    )


def test_oversized_v2_numeric_literals_report_method_errors(tmp_path: Path) -> None:
    count = "9" * 5000
    method = _write(
        tmp_path / "oversized-count.toml",
        f'FORMAT_VERSION = 2\n[STEP.SEARCH.GRID]\nAXES = ["[PB] = lin(0, 1, {count})"]\n',
    )

    with pytest.raises(MethodFormatError, match="Search arguments must be finite"):
        read_method_plan([method])


def test_constraint_and_search_errors_retain_expression_spans(tmp_path: Path) -> None:
    constraint = "[PB] = 2 ** [KEX_AB]"
    constraint_method = _write(
        tmp_path / "constraint-span.toml",
        f'FORMAT_VERSION = 2\n[STEP]\nROLES = [{{ CONSTRAIN = ["{constraint}"] }}]\n',
    )
    with pytest.raises(MethodFormatError) as constraint_error:
        read_method_plan([constraint_method])
    assert constraint_error.value.source.start == constraint.index("2")
    assert constraint_error.value.source.end == len(constraint)

    search = "  [PB] = spline(0, 1, 3)  "
    search_method = _write(
        tmp_path / "search-span.toml",
        f'FORMAT_VERSION = 2\n[STEP.SEARCH.GRID]\nAXES = ["{search}"]\n',
    )
    with pytest.raises(MethodFormatError) as search_error:
        read_method_plan([search_method])
    assert search_error.value.source.start == len(search) - len(search.lstrip())
    assert search_error.value.source.end == len(search.rstrip())


def test_grid_keeps_lin_log_values_and_broad_to_specific_declaration_order(
    tmp_path: Path,
) -> None:
    method = _write(
        tmp_path / "all-grid-forms.toml",
        """
FORMAT_VERSION = 2
[STEP.SEARCH.GRID]
AXES = [
  "[DW_AB] = lin(0, 10, 3)",
  "[KEX_AB] = log(100, 5000, 4)",
  "[DW_AB, NUC->G24N-H] = values(1, 2.5)",
]
""",
    )

    search = read_method_plan([method]).steps[0].search

    assert tuple(axis.render() for axis in search.axes) == (
        "[DW_AB] = lin(0.0, 10.0, 3)",
        "[KEX_AB] = log(100.0, 5000.0, 4)",
        "[DW_AB, NUC->G24N-H] = values(1.0, 2.5)",
    )


def test_v1_duplicate_steps_keep_first_position_and_replace_the_complete_step(
    tmp_path: Path,
) -> None:
    first = _write(
        tmp_path / "first.toml",
        '[A]\nFIT = ["PB"]\n[B]\nFIT = ["KEX_AB"]\n',
    )
    replacement = _write(
        tmp_path / "replacement.toml",
        '[A]\nFIX = ["PB"]\n',
    )

    plan = read_method_plan([first, replacement])

    assert tuple(step.name for step in plan.steps) == ("A", "B")
    assert isinstance(plan.steps[0].role_actions[0], FixAction)
    assert plan.steps[1].roles_from == "A"


def test_grid_validation_applies_broad_to_specific_overrides_before_bounds(
    tmp_path: Path,
) -> None:
    model = _parameter_model(
        ParamDefinition(
            "dw-600", "DW_AB", "", (("h_larmor_frq", 600.0),), 1.0, 0.0, 10.0
        ),
        ParamDefinition(
            "dw-800", "DW_AB", "", (("h_larmor_frq", 800.0),), 1.0, 0.0, 2.0
        ),
    )
    method = _write(
        tmp_path / "grid-overrides.toml",
        """
FORMAT_VERSION = 2
[STEP.SEARCH.GRID]
AXES = [
  "[DW_AB] = lin(0, 10, 3)",
  "[DW_AB, B0->800MHz] = lin(0, 2, 3)",
]
""",
    )

    read_method_plan([method]).validate(model)


def test_grid_resolution_projects_broad_rules_to_active_fit_scope_and_overrides(
    tmp_path: Path,
) -> None:
    model = _parameter_model(
        ParamDefinition("dw-15", "DW_AB", "15N", (), 1.0, 0.0, 10.0),
        ParamDefinition("dw-31", "DW_AB", "31N", (), 1.0, 0.0, 10.0),
        ParamDefinition("dw-99", "DW_AB", "99N", (), 1.0, 0.0, 10.0),
    )
    method = _write(
        tmp_path / "active-grid.toml",
        """
FORMAT_VERSION = 2
[STEP.SEARCH.GRID]
AXES = [
  "[DW_AB] = lin(0, 10, 3)",
  "[DW_AB, NUC->31N] = values(2, 4)",
]
""",
    )
    search = read_method_plan([method]).steps[0].search
    assert search is not None

    resolved = resolve_grid_axes(
        search,  # ty: ignore[invalid-argument-type]
        model,
        active_scope_ids=("dw-15", "dw-31"),
        final_fit_ids=("dw-15", "dw-31"),
    )

    assert tuple(axis.param_id for axis in resolved) == ("dw-15", "dw-31")
    assert resolved[0].values == (0.0, 5.0, 10.0)
    assert resolved[1].values == (2.0, 4.0)
    assert tuple(axis.declaration_ordinal for axis in resolved) == (0, 1)


def test_grid_resolution_rejects_inactive_non_fit_and_out_of_bounds_targets(
    tmp_path: Path,
) -> None:
    model = _parameter_model(
        ParamDefinition("dw-15", "DW_AB", "15N", (), 1.0, 0.0, 10.0),
        ParamDefinition("dw-31", "DW_AB", "31N", (), 1.0, 0.0, 10.0),
    )

    def search(text: str):
        method = _write(
            tmp_path / "runtime-grid.toml",
            f'FORMAT_VERSION = 2\n[STEP.SEARCH.GRID]\nAXES = ["{text}"]\n',
        )
        return read_method_plan([method]).steps[0].search

    with pytest.raises(MethodFormatError, match="no applicable coordinate"):
        resolve_grid_axes(
            search("[DW_AB, NUC->31N] = values(2)"),  # ty: ignore[invalid-argument-type]
            model,
            active_scope_ids=("dw-15",),
            final_fit_ids=("dw-15",),
        )
    with pytest.raises(MethodFormatError, match="final independent FIT"):
        resolve_grid_axes(
            search("[DW_AB, NUC->15N] = values(2)"),  # ty: ignore[invalid-argument-type]
            model,
            active_scope_ids=("dw-15",),
            final_fit_ids=(),
        )
    with pytest.raises(MethodFormatError, match="outside physical bounds"):
        resolve_grid_axes(
            search("[DW_AB, NUC->15N] = values(11)"),  # ty: ignore[invalid-argument-type]
            model,
            active_scope_ids=("dw-15",),
            final_fit_ids=("dw-15",),
        )


def test_broad_grid_rule_projects_away_active_non_fit_matches(tmp_path: Path) -> None:
    model = _parameter_model(
        ParamDefinition("dw-fit", "DW_AB", "15N", (), 1.0, 0.0, 10.0),
        ParamDefinition("dw-fix", "DW_AB", "31N", (), 1.0, 0.0, 10.0),
    )
    method = _write(
        tmp_path / "mixed-active-grid.toml",
        """FORMAT_VERSION = 2
[STEP]
ROLES = [
  { FIX = ["DW_AB"] },
  { FIT = ["DW_AB, NUC->15N"] },
]
[STEP.SEARCH.GRID]
AXES = ["[DW_AB] = values(2, 4)"]
""",
    )
    plan = read_method_plan([method])
    plan.validate(model)
    search = plan.steps[0].search
    assert search is not None

    resolved = resolve_grid_axes(
        search,  # ty: ignore[invalid-argument-type]
        model,
        active_scope_ids=("dw-fit", "dw-fix"),
        final_fit_ids=("dw-fit",),
    )

    assert tuple(axis.param_id for axis in resolved) == ("dw-fit",)


def test_later_specific_grid_rule_replaces_broad_values_before_bounds_check(
    tmp_path: Path,
) -> None:
    model = _parameter_model(
        ParamDefinition("dw-wide", "DW_AB", "15N", (), 1.0, 0.0, 10.0),
        ParamDefinition("dw-narrow", "DW_AB", "31N", (), 1.0, 0.0, 5.0),
    )
    method = _write(
        tmp_path / "override-before-bounds.toml",
        """FORMAT_VERSION = 2
[STEP.SEARCH.GRID]
AXES = [
  "[DW_AB] = values(8)",
  "[DW_AB, NUC->31N] = values(4)",
]
""",
    )
    search = read_method_plan([method]).steps[0].search
    assert search is not None

    resolved = resolve_grid_axes(
        search,  # ty: ignore[invalid-argument-type]
        model,
        active_scope_ids=("dw-wide", "dw-narrow"),
        final_fit_ids=("dw-wide", "dw-narrow"),
    )

    assert tuple((axis.param_id, axis.values) for axis in resolved) == (
        ("dw-wide", (8.0,)),
        ("dw-narrow", (4.0,)),
    )


def test_ordered_complete_roles_validate_fix_fit_constrain_exceptions(
    tmp_path: Path,
) -> None:
    model = _parameter_model(
        ParamDefinition("dw-g24", "DW_AB", "G24N-H", (), 1.0, -10.0, 10.0),
        ParamDefinition("dw-g31", "DW_AB", "G31N-H", (), 1.0, -10.0, 10.0),
        ParamDefinition("dw-g40", "DW_AB", "G40N-H", (), 1.0, -10.0, 10.0),
    )
    method = _write(
        tmp_path / "ordered-roles.toml",
        """
FORMAT_VERSION = 2
[STEP]
ROLES = [
  { FIX = ["DW_AB"] },
  { FIT = ["DW_AB, NUC->G24N-H"] },
  { CONSTRAIN = ["[DW_AB, NUC->G31N-H] = [DW_AB, NUC->G24N-H]"] },
]
[STEP.SEARCH.GRID]
AXES = ["[DW_AB, NUC->G24N-H] = values(1)"]
""",
    )

    read_method_plan([method]).validate(model)

    constrained_search = _write(
        tmp_path / "constrained-search.toml",
        method.read_text(encoding="utf-8").replace(
            'AXES = ["[DW_AB, NUC->G24N-H] = values(1)"]',
            'AXES = ["[DW_AB, NUC->G31N-H] = values(1)"]',
        ),
    )
    with pytest.raises(MethodFormatError, match="not a final independent FIT"):
        read_method_plan([constrained_search]).validate(model)

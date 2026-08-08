"""Behavioral tests for native method-scoped parameterization (#584).

The primary seam is compilation from the sealed parameter model, a current
Analysis Values snapshot, a method declaration, and the selected profile scope.
"""

from __future__ import annotations

import dataclasses
from argparse import Namespace
from pathlib import Path

import pytest

from chemex import chemex as chemex_module
from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Method, Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.experiments.builder import build_experiments
from chemex.models.factory import model_factory
from chemex.parameters.parameterization import (
    AmbiguousParameterReferenceError,
    ConstraintCycleError,
    ConstraintDomainError,
    ConstraintEvaluationError,
    ConstraintProgramMismatchError,
    ConstraintSelfReferenceError,
    IncompatibleParameterizationInputError,
    IncompleteParameterDependenciesError,
    NonFiniteParameterValueError,
    NoParameterMatchError,
    ParameterDeclaration,
    ParameterRole,
    SealedParameterDeclarations,
    SealedParameterModel,
    UnsupportedConstraintExpressionError,
    compile_active_parameterization,
)
from chemex.parameters.sealed import (
    ParamConfig,
    ParamDefinition,
    SealedConfiguration,
    SealedDefinitions,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.userfunctions import user_function_registry
from chemex.parameters.values import AnalysisValuesSnapshot
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parent.parent
DCEST_EXPERIMENT = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Experiments/3hz.toml"
DCEST_PARAMETERS = (
    ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Parameters/parameters.toml"
)
DCEST_METHOD = ROOT / "examples/Experiments/DCEST_15N_HD_EXCH/Methods/method.toml"
MF_EXPERIMENT = ROOT / "examples/Experiments/CEST_15N_TR/Experiments/trosy.toml"
MF_PARAMETERS = ROOT / "examples/Experiments/CEST_15N_TR/Parameters/parameters.toml"
MF_METHOD = ROOT / "examples/Experiments/CEST_15N_TR/Methods/method.toml"
FIT_EXPERIMENT = ROOT / "examples/Experiments/RELAXATION_HZNZ/Experiments/800mhz.toml"
FIT_PARAMETERS = (
    ROOT / "examples/Experiments/RELAXATION_HZNZ/Parameters/parameters.toml"
)
FIT_METHOD = ROOT / "examples/Experiments/RELAXATION_HZNZ/Methods/method.toml"


def _build_dcest_session() -> tuple[AnalysisSession, set[str]]:
    session = AnalysisSession.create()
    session.set_model("2st_hd")
    experiments = build_experiments(
        [DCEST_EXPERIMENT],
        Selection(include=[SpinSystem.from_name("1N")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([DCEST_PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    return session, experiments.param_ids


def _native_fixture(
    declarations: tuple[ParameterDeclaration, ...],
    *,
    values: dict[str, float] | None = None,
    definitions: tuple[ParamDefinition, ...] | None = None,
    model_name: str = "2st",
    occurrence_identity: str = "occurrence-a",
) -> tuple[SealedParameterModel, AnalysisValuesSnapshot]:
    if definitions is None:
        definitions = tuple(
            ParamDefinition(
                item.param_id,
                item.param_id.removeprefix("__"),
                "",
                (),
                1.0,
                -float("inf"),
                float("inf"),
            )
            for item in declarations
        )
    sealed_definitions = SealedDefinitions(
        definitions,
        {item.param_id: index for index, item in enumerate(definitions)},
    )
    current = values or {item.param_id: 1.0 for item in declarations}
    configs = tuple(
        ParamConfig(item.param_id, current[item.param_id], -float("inf"), float("inf"))
        for item in declarations
    )
    configuration = SealedConfiguration(
        configs,
        {item.param_id: index for index, item in enumerate(configs)},
        definitions_identity=sealed_definitions.identity,
    )
    sealed_declarations = SealedParameterDeclarations(declarations)
    parameter_model = SealedParameterModel(
        model_name,
        sealed_definitions,
        configuration,
        sealed_declarations,
    )
    snapshot = AnalysisValuesSnapshot(
        occurrence_identity=occurrence_identity,
        model_identity=f"{model_name}|test",
        definitions_identity=sealed_definitions.identity,
        configuration_identity=configuration.identity,
        revision=0,
        _items=tuple(current.items()),
    )
    return parameter_model, snapshot


def test_shipped_method_compiles_roles_and_resolves_without_mutation() -> None:
    session, required_ids = _build_dcest_session()
    methods = read_methods([DCEST_METHOD])
    method = methods["STEP1"]
    before = session.analysis_values.snapshot()
    legacy_expressions = {
        param_id: parameter.expr
        for param_id, parameter in session.parameters.database._parameters.items()
    }

    parameterization = session.try_compile_parameterization(method, required_ids)
    assert parameterization is not None
    resolved = session.last_resolved_parameter_values
    assert resolved is not None

    definitions = session.parameter_factory.sealed_definitions
    assert definitions is not None
    r1_a = next(
        item.param_id
        for item in definitions
        if item.name == "R1_A" and item.spin_system_name == "1N"
    )
    r1_b = next(
        item.param_id
        for item in definitions
        if item.name == "R1_B" and item.spin_system_name == "1N"
    )
    d2o = next(item.param_id for item in definitions if item.name == "D2O")

    assert parameterization.role(r1_b) is ParameterRole.DERIVED
    assert parameterization.role(d2o) is ParameterRole.FIT
    assert resolved[r1_b] == pytest.approx(0.5 * resolved[r1_a])
    assert session.analysis_values.snapshot() == before
    assert (
        session.parameters.database._parameters[r1_b].expr == legacy_expressions[r1_b]
    )

    step2 = session.try_compile_parameterization(methods["STEP2"], required_ids)
    assert step2 is not None
    assert step2.role(d2o) is ParameterRole.FIX
    assert session.analysis_values.snapshot() == before


def test_constraint_chain_is_deterministic_and_freshly_reresolves() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False),
        ParameterDeclaration("__C", False),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 3.0, "__B": 20.0, "__C": 30.0},
    )
    method = Method(
        constraints=("[B] = [A] * 2", "[C] = [B] + 1"),
        fit=("A",),
    )

    first = compile_active_parameterization(parameter_model, snapshot, method, {"__C"})
    second = compile_active_parameterization(parameter_model, snapshot, method, {"__C"})
    frame = first.frame_from_snapshot(snapshot)

    assert first.program.fingerprint == second.program.fingerprint
    assert first.scope_ids == ("__A", "__B", "__C")
    assert first.program.evaluation_order == ("__B", "__C")
    assert first.resolve(frame)["__C"] == pytest.approx(7.0)
    assert first.resolve(frame.with_updates({"__A": 4.0}))["__C"] == pytest.approx(9.0)
    assert snapshot["__B"] == 20.0
    with pytest.raises(dataclasses.FrozenInstanceError):
        first.source_revision = 1  # ty: ignore[invalid-assignment]


def test_legacy_constraints_fix_fit_precedence_uses_final_fit_role() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 3.0, "__B": 20.0},
    )
    method = Method(
        constraints=("[B] = [A] * 2",),
        fix=("B",),
        fit=("B",),
    )

    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        method,
        {"__B"},
    )

    assert parameterization.role("__B") is ParameterRole.FIT
    assert parameterization.scope_ids == ("__B",)
    assert parameterization.derived_ids == ()
    assert (
        parameterization.resolve(parameterization.frame_from_snapshot(snapshot))["__B"]
        == 20.0
    )


def test_later_constraint_declaration_wins_for_the_same_target() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 3.0, "__B": 20.0},
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(
            constraints=(
                "[B] = [A] * 2",
                "[B] = [A] * 3",
            )
        ),
        {"__B"},
    )

    resolved = parameterization.resolve(parameterization.frame_from_snapshot(snapshot))

    assert resolved["__B"] == 9.0
    with pytest.raises(TypeError):
        resolved["__B"] = 7.0  # ty: ignore[invalid-assignment]


def test_shadowed_constraint_still_receives_public_syntax_validation() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False),
    )
    parameter_model, snapshot = _native_fixture(declarations)

    with pytest.raises(UnsupportedConstraintExpressionError):
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(
                constraints=("[B] = max([A], 1)",),
                fit=("B",),
            ),
            {"__B"},
        )


def test_real_model_free_scientific_expression_matches_legacy_resolution() -> None:
    session = AnalysisSession.create()
    session.set_model("2st.mf")
    experiments = build_experiments(
        [MF_EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G23N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([MF_PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    snapshot = session.analysis_values.snapshot()
    configuration = session.parameter_factory.sealed_configuration
    assert configuration is not None
    assert tuple(snapshot) == tuple(item.param_id for item in configuration)
    assert all(
        snapshot[item.param_id] == item.effective_value
        for item in configuration
        if item.effective_value is not None
    )
    assert any(item.effective_value is None for item in configuration)
    assert not snapshot.occurrence_identity.startswith("bootstrap:")

    method = read_methods([MF_METHOD])["DEFAULT"]
    parameterization = session.compile_parameterization(method, experiments.param_ids)
    resolved = parameterization.resolve(parameterization.frame_from_snapshot(snapshot))
    legacy = session.parameters.build_lmfit_params(experiments.param_ids).valuesdict()
    definitions = session.parameter_factory.sealed_definitions
    assert definitions is not None
    derived_rate_ids = tuple(
        item.param_id
        for item in definitions
        if item.param_id in parameterization.derived_ids and item.name.startswith("R2")
    )

    assert derived_rate_ids
    pb = next(item.param_id for item in definitions if item.name == "PB")
    kex_ab = next(item.param_id for item in definitions if item.name == "KEX_AB")
    assert parameterization.role(pb) is ParameterRole.FIX
    assert parameterization.role(kex_ab) is ParameterRole.FIX
    for param_id in derived_rate_ids:
        assert resolved[param_id] == pytest.approx(
            legacy[param_id], rel=1e-13, abs=1e-13
        )


@pytest.mark.parametrize("model_name", ("2st_binding", "2st_eyring"))
def test_registered_scientific_model_expression_compiles_through_restricted_language(
    model_name: str,
) -> None:
    AnalysisSession.create()
    conditions = Conditions(
        h_larmor_frq=800.0,
        temperature=25.0,
        p_total=1e-3,
        l_total=2e-3,
        d2o=0.1,
    )
    settings = model_factory.create(model_name, conditions)
    name_map = {name: f"__{name.upper()}" for name in settings}
    declarations = tuple(
        ParameterDeclaration(
            name_map[name],
            setting.vary,
            "" if setting.vary else setting.expr.format_map(name_map),
        )
        for name, setting in settings.items()
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        model_name=model_name,
    )

    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        set(name_map.values()),
    )

    assert parameterization.scope_ids == tuple(name_map.values())


def test_supported_arithmetic_has_explicit_precedence_and_unary_semantics() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 6.0, "__B": 0.0},
    )
    method = Method(constraints=("[B] = -([A] - 2) / +2",))
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        method,
        {"__B"},
    )

    resolved = parameterization.resolve(parameterization.frame_from_snapshot(snapshot))

    assert resolved["__B"] == -2.0


def test_no_match_failure_has_stable_type_code_and_context() -> None:
    parameter_model, snapshot = _native_fixture((ParameterDeclaration("__A", False),))

    with pytest.raises(NoParameterMatchError) as raised:
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(fit=("MISSING",)),
            {"__A"},
        )

    assert raised.value.code == "no_match"
    assert raised.value.context["selector"] == "missing"


def test_contextually_incomparable_reference_is_ambiguity() -> None:
    definitions = (
        ParamDefinition("__A_SPIN", "A", "G1N-HN", (), 1.0, -10.0, 10.0),
        ParamDefinition(
            "__A_FIELD", "A", "", (("h_larmor_frq", 800.0),), 1.0, -10.0, 10.0
        ),
        ParamDefinition(
            "__T", "T", "G1N-HN", (("h_larmor_frq", 800.0),), 1.0, -10.0, 10.0
        ),
    )
    declarations = tuple(
        ParameterDeclaration(item.param_id, False) for item in definitions
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        definitions=definitions,
    )

    with pytest.raises(AmbiguousParameterReferenceError) as raised:
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(constraints=("[T] = [A]",)),
            {"__T"},
        )

    assert raised.value.code == "ambiguity"
    assert raised.value.context["candidate_ids"] == ("__A_SPIN", "__A_FIELD")


def test_direct_self_reference_is_rejected_before_evaluation() -> None:
    parameter_model, snapshot = _native_fixture((ParameterDeclaration("__A", False),))

    with pytest.raises(ConstraintSelfReferenceError) as raised:
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(constraints=("[A] = [A]",)),
            {"__A"},
        )

    assert raised.value.code == "self_reference"
    assert raised.value.context["target_id"] == "__A"


def test_indirect_constraint_cycle_reports_exact_ids_and_provenance() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False),
        ParameterDeclaration("__C", False),
    )
    parameter_model, snapshot = _native_fixture(declarations)

    with pytest.raises(ConstraintCycleError) as raised:
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(
                constraints=("[A] = [B]", "[B] = [A]", "[C] = [A]"),
            ),
            {"__C"},
        )

    assert raised.value.code == "cycle"
    assert raised.value.context["param_ids"] == ("__A", "__B")
    assert raised.value.context["constraints"] == (
        ("__A", "method-rule:0", "[b]"),
        ("__B", "method-rule:1", "[a]"),
    )


def test_arithmetic_domain_failure_does_not_return_partial_values() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 0.0, "__B": 10.0},
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(constraints=("[B] = 1 / [A]",)),
        {"__B"},
    )

    with pytest.raises(ConstraintDomainError) as raised:
        parameterization.resolve(parameterization.frame_from_snapshot(snapshot))

    assert raised.value.code == "domain_error"
    assert raised.value.context["target_id"] == "__B"


def test_scientific_function_evaluation_failure_is_typed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def explode(_value: float) -> float:
        msg = "backend exploded"
        raise RuntimeError(msg)

    monkeypatch.setitem(
        user_function_registry.user_function_registry,
        "test_model",
        {"explode": explode},
    )
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False, "explode(__A)"),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        model_name="test_model",
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        {"__B"},
    )

    with pytest.raises(ConstraintEvaluationError) as raised:
        parameterization.resolve(parameterization.frame_from_snapshot(snapshot))

    assert raised.value.code == "evaluation_error"
    assert raised.value.context["function_id"] == "explode"


def test_malformed_scientific_function_component_is_typed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def malformed(_value: float) -> dict[str, object]:
        return {"result": [1.0]}

    monkeypatch.setitem(
        user_function_registry.user_function_registry,
        "test_model",
        {"malformed": malformed},
    )
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration(
            "__B",
            False,
            'malformed(__A)["result"] + 1',
        ),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        model_name="test_model",
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        {"__B"},
    )

    with pytest.raises(ConstraintEvaluationError) as raised:
        parameterization.resolve(parameterization.frame_from_snapshot(snapshot))

    assert raised.value.code == "evaluation_error"
    assert raised.value.context["target_id"] == "__B"
    assert raised.value.context["operator"] == "add"


def test_non_finite_derived_value_has_its_own_failure_category() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 1e308, "__B": 0.0},
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(constraints=("[B] = [A] * [A]",)),
        {"__B"},
    )

    with pytest.raises(NonFiniteParameterValueError) as raised:
        parameterization.resolve(parameterization.frame_from_snapshot(snapshot))

    assert raised.value.code == "non_finite"
    assert raised.value.context["param_id"] == "__B"


def test_unknown_model_dependency_is_incomplete_dependency_failure() -> None:
    declarations = (ParameterDeclaration("__A", False, "__MISSING + 1"),)
    parameter_model, snapshot = _native_fixture(declarations)

    with pytest.raises(IncompleteParameterDependenciesError) as raised:
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(),
            {"__A"},
        )

    assert raised.value.code == "incomplete_dependencies"
    assert raised.value.context["dependency_id"] == "__MISSING"


def test_other_occurrence_snapshot_is_incompatible_even_for_same_model() -> None:
    parameter_model, snapshot = _native_fixture((ParameterDeclaration("__A", False),))
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        {"__A"},
    )
    other_occurrence = dataclasses.replace(
        snapshot,
        occurrence_identity="occurrence-b",
    )

    with pytest.raises(IncompatibleParameterizationInputError) as raised:
        parameterization.frame_from_snapshot(other_occurrence)

    assert raised.value.code == "incompatible_input"


def test_incomplete_and_program_mismatched_frames_are_distinct_failures() -> None:
    parameter_model, snapshot = _native_fixture((ParameterDeclaration("__A", False),))
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        {"__A"},
    )
    frame = parameterization.frame_from_snapshot(snapshot)

    with pytest.raises(IncompleteParameterDependenciesError):
        parameterization.resolve(dataclasses.replace(frame, _items=()))
    with pytest.raises(ConstraintProgramMismatchError) as raised:
        parameterization.resolve(
            dataclasses.replace(frame, program_fingerprint="different-program")
        )

    assert raised.value.code == "program_mismatch"


def test_native_compilation_failure_cannot_veto_real_legacy_fit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = AnalysisSession.create()

    def fail_preview(*_args: object, **_kwargs: object) -> None:
        msg = "native preview failed"
        raise RuntimeError(msg)

    monkeypatch.setattr(session, "compile_parameterization", fail_preview)
    args = Namespace(
        commands="fit",
        model="2st",
        include=[SpinSystem.from_name("G2N-HN")],
        exclude=None,
        experiments=[FIT_EXPERIMENT],
        parameters=[FIT_PARAMETERS],
        method=[FIT_METHOD],
        output=tmp_path,
        plot="nothing",
        workers=1,
        native_threads=1,
    )

    chemex_module.run(args, session=session)

    assert list((tmp_path / "Data").rglob("*.dat"))
    assert list((tmp_path / "Parameters").rglob("*.toml"))
    assert isinstance(session.last_parameterization_failure, RuntimeError)
    assert session.last_parameterization is None
    assert session.last_resolved_parameter_values is None

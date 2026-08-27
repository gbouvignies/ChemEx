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
from chemex.configuration.method_execution import normalize_methods_for_execution
from chemex.configuration.method_plan import (
    ConstrainAction,
    Constraint,
    FitAction,
    FixAction,
    FormatOrigin,
    MethodPlan,
    ParameterSelector,
    SourceRef,
    StepPlan,
)
from chemex.configuration.method_plan import (
    LiteralExpression as MethodLiteralExpression,
)
from chemex.configuration.method_plan import (
    SelectorExpression as MethodSelectorExpression,
)
from chemex.configuration.methods import (
    Method,
    Selection,
    read_method_plan,
)
from chemex.configuration.parameters import read_defaults
from chemex.experiments.builder import build_experiments
from chemex.models.factory import model_factory
from chemex.parameters.database import ParameterIndex
from chemex.parameters.name import ParamName
from chemex.parameters.parameterization import (
    ActiveParameterization,
    AmbiguousParameterReferenceError,
    ConstraintCycleError,
    ConstraintDomainError,
    ConstraintEvaluationError,
    ConstraintProgramMismatchError,
    ConstraintSelfReferenceError,
    IncompatibleParameterizationInputError,
    IncompleteParameterDependenciesError,
    ModelDerivationOverrideError,
    NonFiniteParameterValueError,
    NoParameterMatchError,
    ParameterDeclaration,
    ParameterDeclarationContribution,
    ParameterRole,
    ScientificFunctionBinder,
    SealedParameterDeclarations,
    SealedParameterModel,
    UnsupportedConstraintExpressionError,
    compile_active_parameterization,
    compile_active_parameterization_from_actions,
    seal_parameter_declarations,
)
from chemex.parameters.sealed import (
    ParamConfig,
    ParamDefinition,
    SealedConfiguration,
    SealedDefinitions,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.userfunctions import user_function_registry
from chemex.parameters.values import AnalysisValues, AnalysisValuesSnapshot
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
BINDING_ROOT = ROOT / "examples/Combinations/2stBinding"
BINDING_EXPERIMENTS = tuple(sorted((BINDING_ROOT / "Experiments").glob("*.toml")))
BINDING_PARAMETERS = BINDING_ROOT / "Parameters/params.toml"
BINDING_METHOD = BINDING_ROOT / "Methods/method.toml"
METHYL_ROOT = ROOT / "examples/Combinations/CPMG_CH3_1H_DQ_TQ"
METHYL_EXPERIMENTS = tuple(sorted((METHYL_ROOT / "Experiments").glob("*.toml")))
METHYL_PARAMETERS = METHYL_ROOT / "Parameters/parameters.toml"
METHYL_METHOD = METHYL_ROOT / "Methods/method.toml"
_METHOD_SOURCE = SourceRef(Path("method.toml"), "STEP", "ROLES")

_BINDING_STEP1_R2_B_IDS = {
    "__R2_B_486N_800_0MHZ",
    "__R2_B_488N_800_0MHZ",
    "__R2_B_489N_800_0MHZ",
    "__R2_B_489NE_800_0MHZ",
    "__R2_B_491N_800_0MHZ",
    "__R2_B_492N_800_0MHZ",
    "__R2_B_493N_800_0MHZ",
    "__R2_B_494N_800_0MHZ",
    "__R2_B_495N_800_0MHZ",
    "__R2_B_496N_800_0MHZ",
    "__R2_B_497N_800_0MHZ",
    "__R2_B_498N_800_0MHZ",
    "__R2_B_499N_800_0MHZ",
}
_BINDING_STEP2_R2_B_IDS = {
    f"__R2_B_{spin}_800_0MHZ"
    for spin in (
        "480N",
        "481N",
        "482N",
        "483N",
        "484N",
        "485N",
        "486N",
        "488N",
        "489N",
        "489NE",
        "491N",
        "492N",
        "493N",
        "494N",
        "495N",
        "496N",
        "497N",
        "498N",
        "499N",
        "500N",
        "501N",
        "502N",
        "503N",
        "504N",
        "505N",
        "506N",
        "509N",
        "510N",
        "511N",
        "512N",
        "513N",
        "514N",
        "515N",
        "516N",
        "518N",
        "520N",
        "521N",
        "522N",
        "523N",
        "524N",
        "526N",
        "527N",
        "528N",
        "529N",
        "530N",
        "532N",
        "533N",
        "535N",
        "537N",
        "538N",
        "540N",
        "541N",
        "543N",
        "544N",
        "545N",
        "546N",
        "547N",
        "550N",
        "551N",
        "552N",
        "553N",
        "554N",
        "555N",
        "556N",
        "557N",
        "558N",
        "559N",
        "560N",
        "566N",
        "567N",
        "568N",
        "569N",
        "570N",
        "571N",
        "572N",
        "573N",
        "574N",
        "575N",
    )
}


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


def _build_mf_session() -> tuple[AnalysisSession, set[str]]:
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
    return session, experiments.param_ids


def _build_fit_session() -> tuple[AnalysisSession, set[str]]:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [FIT_EXPERIMENT],
        Selection(include=[SpinSystem.from_name("G2N-HN")], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([FIT_PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    return session, experiments.param_ids


def _native_fixture(
    declarations: tuple[ParameterDeclaration, ...],
    *,
    values: dict[str, float] | None = None,
    definitions: tuple[ParamDefinition, ...] | None = None,
    configuration_bounds: dict[str, tuple[float, float]] | None = None,
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
    bounds = configuration_bounds or {}
    configs = tuple(
        ParamConfig(
            item.param_id,
            current[item.param_id],
            *bounds.get(item.param_id, (-float("inf"), float("inf"))),
        )
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
        f"{model_name}|test",
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


def _compile_context_constraint(
    definitions: tuple[ParamDefinition, ...],
    expression: str,
    target_id: str = "__TARGET",
) -> ActiveParameterization:
    declarations = tuple(
        ParameterDeclaration(definition.param_id, False) for definition in definitions
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        definitions=definitions,
        values={
            definition.param_id: definition.default_value for definition in definitions
        },
    )
    return compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(constraints=(expression,)),
        {target_id},
    )


def test_canonical_ordered_actions_compile_complete_replacement_roles() -> None:
    declarations = (
        ParameterDeclaration("__DW_AB_24N", True),
        ParameterDeclaration("__DW_AB_31N", True),
        ParameterDeclaration("__DW_AB_40N", True),
    )
    definitions = tuple(
        ParamDefinition(
            declaration.param_id,
            "DW_AB",
            spin_system,
            (),
            1.0,
            -float("inf"),
            float("inf"),
        )
        for declaration, spin_system in zip(
            declarations,
            ("24N", "31N", "40N"),
            strict=True,
        )
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        definitions=definitions,
    )
    broad = ParameterSelector("DW_AB")
    subset_a = ParameterSelector("DW_AB", spin_system="24N")
    subset_b = ParameterSelector("DW_AB", spin_system="31N")
    actions = (
        FixAction((broad,), _METHOD_SOURCE),
        FitAction((subset_a,), _METHOD_SOURCE),
        ConstrainAction(
            (
                Constraint(
                    subset_b,
                    MethodSelectorExpression(subset_a),
                    _METHOD_SOURCE,
                ),
            ),
            _METHOD_SOURCE,
        ),
    )

    parameterization = compile_active_parameterization_from_actions(
        parameter_model,
        snapshot,
        actions,
        set(parameter_model.declarations),
    )

    assert parameterization.role("__DW_AB_24N") is ParameterRole.FIT
    assert parameterization.role("__DW_AB_31N") is ParameterRole.DERIVED
    assert parameterization.role("__DW_AB_40N") is ParameterRole.FIX


@pytest.mark.parametrize("replacement", (FitAction, FixAction))
def test_canonical_fit_or_fix_removes_an_earlier_constraint(
    replacement: type[FitAction] | type[FixAction],
) -> None:
    declarations = (ParameterDeclaration("__PB", True),)
    parameter_model, snapshot = _native_fixture(declarations)
    selector = ParameterSelector("PB")
    actions = (
        ConstrainAction(
            (
                Constraint(
                    selector,
                    MethodLiteralExpression(0.25),
                    _METHOD_SOURCE,
                ),
            ),
            _METHOD_SOURCE,
        ),
        replacement((selector,), _METHOD_SOURCE),
    )

    parameterization = compile_active_parameterization_from_actions(
        parameter_model,
        snapshot,
        actions,
        {"__PB"},
    )

    expected = ParameterRole.FIT if replacement is FitAction else ParameterRole.FIX
    assert parameterization.role("__PB") is expected
    assert parameterization.program.constraints == ()


def test_method_plan_resolves_explicit_and_baseline_role_actions_immutably() -> None:
    selector = ParameterSelector("PB")
    fixed = FixAction((selector,), _METHOD_SOURCE)
    fitted = FitAction((selector,), _METHOD_SOURCE)
    plan = MethodPlan(
        FormatOrigin.V2,
        (
            StepPlan("FIRST", role_actions=(fixed,)),
            StepPlan("INHERITED", roles_from="FIRST", role_actions=(fitted,)),
            StepPlan("BASELINE", role_actions=(fitted,)),
        ),
    )

    effective = plan.effective_role_actions()

    assert effective == {
        "FIRST": (fixed,),
        "INHERITED": (fixed, fitted),
        "BASELINE": (fitted,),
    }


def test_inherited_canonical_action_outside_selected_scope_is_inert() -> None:
    declarations = (
        ParameterDeclaration("__DW_AB_24N", True),
        ParameterDeclaration("__DW_AB_31N", True),
    )
    definitions = (
        ParamDefinition("__DW_AB_24N", "DW_AB", "24N", (), 1.0, -10.0, 10.0),
        ParamDefinition("__DW_AB_31N", "DW_AB", "31N", (), 1.0, -10.0, 10.0),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        definitions=definitions,
    )
    inherited = FixAction(
        (ParameterSelector("DW_AB", spin_system="31N"),),
        _METHOD_SOURCE,
    )

    parameterization = compile_active_parameterization_from_actions(
        parameter_model,
        snapshot,
        (inherited,),
        {"__DW_AB_24N"},
    )

    assert parameterization.scope_ids == ("__DW_AB_24N",)
    assert parameterization.role("__DW_AB_24N") is ParameterRole.FIT


def test_shipped_method_compiles_roles_and_resolves_without_mutation() -> None:
    session, required_ids = _build_dcest_session()
    plan = read_method_plan([DCEST_METHOD])
    effective_actions = plan.effective_role_actions()
    before = session.analysis_values.snapshot()
    legacy_expressions = {
        param_id: parameter.expr
        for param_id, parameter in session.parameters.database._parameters.items()
    }

    parameterization = session.compile_parameterization_from_actions(
        effective_actions["STEP1"], required_ids
    )
    resolved = parameterization.resolve(parameterization.frame_from_snapshot(before))

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

    step2 = session.compile_parameterization_from_actions(
        effective_actions["STEP2"], required_ids
    )
    assert step2.role(d2o) is ParameterRole.FIX
    assert session.analysis_values.snapshot() == before


@pytest.mark.parametrize(
    "method",
    (
        Method(fit=("PA",)),
        Method(fix=("PA",)),
        Method(constraints=("[PA] = 0.5",)),
    ),
    ids=("fit", "fix", "constraint"),
)
def test_shipped_2st_hd_population_derivation_cannot_be_overridden(
    method: Method,
) -> None:
    session, required_ids = _build_dcest_session()
    definitions = session.parameter_factory.sealed_definitions
    assert definitions is not None
    pa_id = next(item.param_id for item in definitions if item.name == "PA")
    baseline = session.compile_parameterization(Method(), required_ids)

    assert baseline.role(pa_id) is ParameterRole.DERIVED

    with pytest.raises(ModelDerivationOverrideError) as raised:
        session.compile_parameterization(method, required_ids)

    assert raised.value.code == "model_derivation_override"
    assert raised.value.context["param_ids"] == (pa_id,)


def test_sealing_retains_model_derivation_when_estimation_is_supported() -> None:
    definition = ParamDefinition(
        "__PA",
        "PA",
        "",
        (),
        1.0,
        -float("inf"),
        float("inf"),
    )
    definitions = SealedDefinitions((definition,), {definition.param_id: 0})
    declarations = seal_parameter_declarations(
        definitions,
        {
            definition.param_id: (
                ParameterDeclarationContribution(
                    definition.param_id,
                    False,
                    "1.0 - __PB",
                    "model",
                    True,
                ),
                ParameterDeclarationContribution(
                    definition.param_id,
                    True,
                    "",
                    "experiment",
                ),
            )
        },
    )

    assert declarations[definition.param_id].supports_estimation
    assert declarations[definition.param_id].model_expression == "1.0 - __PB"


def test_non_model_owned_baseline_expression_can_compile_as_fit() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", True, "__A", model_owned=False),
    )
    parameter_model, snapshot = _native_fixture(declarations)

    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(fit=("B",)),
        {"__B"},
    )

    declaration = parameter_model.declarations["__B"]
    assert declaration.model_expression == "__A"
    assert declaration.supports_estimation
    assert not declaration.model_owned
    assert parameterization.role("__B") is ParameterRole.FIT


def test_non_model_owned_baseline_expression_remains_active_constraint() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", True, "__A", model_owned=False),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 3.0, "__B": 20.0},
    )

    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(constraints=("[B] = [A]",)),
        {"__B"},
    )
    resolved = parameterization.resolve(parameterization.frame_from_snapshot(snapshot))

    assert parameterization.role("__B") is ParameterRole.DERIVED
    assert resolved["__B"] == pytest.approx(3.0)


def test_model_owned_expression_cannot_compile_as_fit() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", True, "__A", model_owned=True),
    )
    parameter_model, snapshot = _native_fixture(declarations)

    with pytest.raises(ModelDerivationOverrideError):
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(fit=("B",)),
            {"__B"},
        )


@pytest.mark.parametrize(
    "method",
    (
        Method(fit=("PA",)),
        Method(fix=("PA",)),
        Method(constraints=("[PA] = 0.5",)),
    ),
)
def test_product_role_application_rejects_model_owned_override(method: Method) -> None:
    session, required_ids = _build_dcest_session()

    with pytest.raises(ModelDerivationOverrideError):
        session.compile_parameterization(method, required_ids)


def test_non_estimable_declaration_cannot_compile_as_fit() -> None:
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False, "__A", model_owned=False),
    )
    parameter_model, snapshot = _native_fixture(declarations)

    with pytest.raises(IncompatibleParameterizationInputError) as raised:
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(fit=("B",)),
            {"__B"},
        )

    assert raised.value.context["param_ids"] == ("__B",)


def test_binding_current_roles_compile_all_estimable_r2_b_coordinates() -> None:
    session = AnalysisSession.create()
    session.set_model("2st_binding")
    experiments = build_experiments(
        BINDING_EXPERIMENTS,
        Selection(include=None, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([BINDING_PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    plan = read_method_plan([BINDING_METHOD])
    _plan, operational = normalize_methods_for_execution(plan)
    effective_actions = plan.effective_role_actions()

    experiments.select(operational["STEP1"].selection)
    step1 = session.compile_parameterization_from_actions(
        effective_actions["STEP1"],
        experiments.param_ids,
    )
    step1_fit_ids = {
        param_id
        for param_id in step1.independent_ids
        if step1.role(param_id) is ParameterRole.FIT
    }
    step1_r2_b_ids = {
        param_id for param_id in step1_fit_ids if param_id.startswith("__R2_B_")
    }

    assert len(step1_fit_ids) == 53
    assert step1_r2_b_ids == _BINDING_STEP1_R2_B_IDS
    for param_id in _BINDING_STEP1_R2_B_IDS:
        declaration = session.parameter_factory.sealed_parameter_model.declarations[
            param_id
        ]
        assert declaration.model_expression
        assert declaration.supports_estimation
        assert not declaration.model_owned

    experiments.select(operational["STEP2"].selection)
    step2 = session.compile_parameterization_from_actions(
        effective_actions["STEP2"],
        experiments.param_ids,
    )
    step2_fit_ids = {
        param_id
        for param_id in step2.independent_ids
        if step2.role(param_id) is ParameterRole.FIT
    }
    step2_r2_b_ids = {
        param_id for param_id in step2_fit_ids if param_id.startswith("__R2_B_")
    }

    assert len(step2_fit_ids) == 312
    assert len(_BINDING_STEP2_R2_B_IDS) == 78
    assert step2_r2_b_ids == _BINDING_STEP2_R2_B_IDS


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
    native_session, required_ids = _build_mf_session()
    snapshot = native_session.analysis_values.snapshot()
    configuration = native_session.parameter_factory.sealed_configuration
    assert configuration is not None
    assert tuple(snapshot) == tuple(item.param_id for item in configuration)
    assert all(
        snapshot[item.param_id] == item.effective_value
        for item in configuration
        if item.effective_value is not None
    )
    assert any(item.effective_value is None for item in configuration)
    assert not snapshot.occurrence_identity.startswith("bootstrap:")

    plan = read_method_plan([MF_METHOD])
    parameterization = native_session.compile_parameterization_from_actions(
        plan.effective_role_actions()["DEFAULT"], required_ids
    )
    definitions = native_session.parameter_factory.sealed_definitions
    assert definitions is not None
    pb = next(item.param_id for item in definitions if item.name == "PB")
    kex_ab = next(item.param_id for item in definitions if item.name == "KEX_AB")
    assert parameterization.role(pb) is ParameterRole.FIX
    assert parameterization.role(kex_ab) is ParameterRole.FIX
    assert any(
        item.param_id in parameterization.derived_ids and item.name.startswith("R2")
        for item in definitions
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


def test_finite_derived_value_outside_configured_bounds_commits() -> None:
    declarations = (
        ParameterDeclaration("__R2", False),
        ParameterDeclaration("__R1", False),
        ParameterDeclaration("__R2A", False, "__R2 - __R1", model_owned=True),
    )
    parameter_model, _snapshot = _native_fixture(
        declarations,
        values={"__R2": 1.0, "__R1": 2.0, "__R2A": 12.0},
        configuration_bounds={"__R2A": (0.0, float("inf"))},
    )
    values = AnalysisValues()
    values.initialize(parameter_model.model_identity, parameter_model.configuration)
    initial = values.snapshot()
    parameterization = compile_active_parameterization(
        parameter_model,
        initial,
        Method(),
        {"__R2A"},
    )

    resolved = parameterization.resolve(parameterization.frame_from_snapshot(initial))
    committed = values.commit(
        dict(resolved.items()),
        expected=initial,
        scope=parameterization.scope_ids,
    )

    assert parameterization.role("__R2A") is ParameterRole.DERIVED
    assert resolved["__R2A"] == -1.0
    assert committed["__R2A"] == -1.0
    assert committed.revision == 1


def test_derived_continuity_value_becomes_the_next_independent_value() -> None:
    declarations = (ParameterDeclaration("__A", True),)
    parameter_model, _snapshot = _native_fixture(
        declarations,
        values={"__A": 0.5},
        configuration_bounds={"__A": (0.0, 1.0)},
    )
    values = AnalysisValues()
    values.initialize(parameter_model.model_identity, parameter_model.configuration)
    initial = values.snapshot()
    derived = compile_active_parameterization(
        parameter_model,
        initial,
        Method(constraints=("[A] = 2.0",)),
        {"__A"},
    )
    resolved = derived.resolve(derived.frame_from_snapshot(initial))
    committed = values.commit(
        dict(resolved.items()),
        expected=initial,
        scope=derived.scope_ids,
    )

    independent = compile_active_parameterization(
        parameter_model,
        committed,
        Method(),
        {"__A"},
    )
    next_frame = independent.frame_from_snapshot(committed)
    next_values = dict(next_frame.ordered_items())

    assert derived.role("__A") is ParameterRole.DERIVED
    assert independent.role("__A") is ParameterRole.FIT
    assert next_values["__A"] == 2.0
    assert next_values["__A"] != parameter_model.configuration["__A"].effective_value


@pytest.mark.parametrize(
    ("literal", "expected"),
    (
        ("10", 10.0),
        ("10.25", 10.25),
        (".25", 0.25),
        ("10.", 10.0),
        ("1.25e-2", 0.0125),
        ("2E+3", 2000.0),
    ),
)
def test_public_constraints_accept_legacy_decimal_numeric_syntax(
    literal: str,
    expected: float,
) -> None:
    parameter_model, snapshot = _native_fixture(
        (ParameterDeclaration("__A", False),),
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(constraints=(f"[A] = {literal}",)),
        {"__A"},
    )

    resolved = parameterization.resolve(parameterization.frame_from_snapshot(snapshot))

    assert resolved["__A"] == expected


@pytest.mark.parametrize("literal", ("0x10", "0b10", "0o10", "1_0", "1j"))
def test_public_constraints_reject_python_only_numeric_literals(literal: str) -> None:
    parameter_model, snapshot = _native_fixture(
        (ParameterDeclaration("__A", False),),
    )

    with pytest.raises(UnsupportedConstraintExpressionError) as raised:
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(constraints=(f"[A] = {literal}",)),
            {"__A"},
        )

    assert raised.value.code == "unsupported_expression"


def test_oversized_public_literal_is_a_typed_non_finite_failure() -> None:
    parameter_model, snapshot = _native_fixture(
        (ParameterDeclaration("__A", False),),
    )

    with pytest.raises(NonFiniteParameterValueError) as raised:
        compile_active_parameterization(
            parameter_model,
            snapshot,
            Method(constraints=(f"[A] = {'9' * 400}",)),
            {"__A"},
        )

    assert raised.value.code == "non_finite"


def test_oversized_independent_frame_value_is_a_typed_non_finite_failure() -> None:
    parameter_model, snapshot = _native_fixture(
        (ParameterDeclaration("__A", False),),
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        {"__A"},
    )
    frame = parameterization.frame_from_snapshot(snapshot)

    with pytest.raises(NonFiniteParameterValueError) as raised:
        parameterization.resolve(frame.with_updates({"__A": 10**400}))

    assert raised.value.code == "non_finite"


def test_foreign_model_snapshot_is_rejected_before_program_compilation() -> None:
    parameter_model, snapshot = _native_fixture(
        (ParameterDeclaration("__A", False),),
    )

    with pytest.raises(IncompatibleParameterizationInputError) as raised:
        compile_active_parameterization(
            parameter_model,
            dataclasses.replace(snapshot, model_identity="foreign-model"),
            Method(),
            {"__A"},
        )

    assert raised.value.code == "incompatible_input"
    assert raised.value.context["actual"][0] == "foreign-model"


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


@pytest.mark.parametrize(
    ("method", "required_ids"),
    (
        (Method(fit=("R2",)), {"__R2_A"}),
        (Method(fix=("R2",)), {"__R2_A"}),
        (Method(constraints=("[R2] = 1.0",)), {"__R2_A"}),
        (Method(constraints=("[B] = [R2]",)), {"__B"}),
    ),
    ids=("fit-target", "fix-target", "constraint-target", "constraint-reference"),
)
def test_parameter_selectors_use_legacy_exact_name_matching(
    method: Method,
    required_ids: set[str],
) -> None:
    definitions = (
        ParamDefinition("__R2_A", "R2_A", "", (), 1.0, -10.0, 10.0),
        ParamDefinition("__B", "B", "", (), 1.0, -10.0, 10.0),
    )
    declarations = tuple(
        ParameterDeclaration(item.param_id, False) for item in definitions
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        definitions=definitions,
    )
    legacy_index = ParameterIndex()
    legacy_index.add(ParamName.from_section("R2_A"))

    assert legacy_index.get_matching_ids(ParamName.from_section("R2")) == set()

    with pytest.raises(NoParameterMatchError):
        compile_active_parameterization(
            parameter_model,
            snapshot,
            method,
            required_ids,
        )


def test_exact_spin_reference_outranks_condition_specific_global() -> None:
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

    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(constraints=("[T] = [A]",)),
        {"__T"},
    )

    assert parameterization.program.constraints[0].dependencies == ("__A_SPIN",)


def test_shipped_methyl_constraint_resolves_companion_spin_for_both_states() -> None:
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        METHYL_EXPERIMENTS,
        Selection(
            include=[
                SpinSystem.from_name("I3HD1-CD1"),
                SpinSystem.from_name("V75HG1-CG1"),
                SpinSystem.from_name("V75HG2-CG2"),
            ],
            exclude=None,
        ),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([METHYL_PARAMETERS]))
    assert session.try_build_analysis_values(), repr(
        session.parameter_factory.native_construction_error
    )
    plan = read_method_plan([METHYL_METHOD])
    session.validate_method_plan(plan)
    parameterization = session.compile_parameterization_from_actions(
        plan.effective_role_actions()["STEP1"],
        experiments.param_ids,
    )

    dependencies = {
        constraint.target_id: constraint.dependencies
        for constraint in parameterization.program.constraints
    }
    assert dependencies["__R2ADQ_A_I3HD1_700_2MHZ"] == (
        "__R2DQ_A_I3HD1_700_2MHZ",
        "__R1_A_I3CD1_700_2MHZ",
    )
    assert dependencies["__R2ADQ_B_I3HD1_700_2MHZ"] == (
        "__R2DQ_B_I3HD1_700_2MHZ",
        "__R1_B_I3CD1_700_2MHZ",
    )
    assert dependencies["__R2ADQ_A_V75HG1_700_2MHZ"][-1] == ("__R1_A_V75CG1_700_2MHZ")
    assert dependencies["__R2ADQ_A_V75HG2_700_2MHZ"][-1] == ("__R1_A_V75CG2_700_2MHZ")
    assert dependencies["__R2ADQ_B_V75HG1_700_2MHZ"][-1] == ("__R1_B_V75CG1_700_2MHZ")
    assert dependencies["__R2ADQ_B_V75HG2_700_2MHZ"][-1] == ("__R1_B_V75CG2_700_2MHZ")
    resolved = parameterization.resolve(
        parameterization.frame_from_snapshot(session.analysis_values.snapshot())
    )
    assert resolved["__R2ADQ_A_I3HD1_700_2MHZ"] == pytest.approx(
        resolved["__R2DQ_A_I3HD1_700_2MHZ"] - resolved["__R1_A_I3CD1_700_2MHZ"] / 3.0
    )
    assert resolved["__R2ADQ_B_I3HD1_700_2MHZ"] == pytest.approx(
        resolved["__R2DQ_B_I3HD1_700_2MHZ"] - resolved["__R1_B_I3CD1_700_2MHZ"] / 3.0
    )


def test_exact_spin_reference_outranks_same_group_companion() -> None:
    parameterization = _compile_context_constraint(
        (
            ParamDefinition("__TARGET", "TARGET", "I3HD1", (), 0.0, -10.0, 10.0),
            ParamDefinition("__VALUE_EXACT", "VALUE", "I3HD1", (), 2.0, -10.0, 10.0),
            ParamDefinition(
                "__VALUE_COMPANION", "VALUE", "I3CD1", (), 3.0, -10.0, 10.0
            ),
        ),
        "[TARGET] = [VALUE]",
    )

    assert parameterization.program.constraints[0].dependencies == ("__VALUE_EXACT",)


def test_explicit_spin_selector_remains_authoritative() -> None:
    parameterization = _compile_context_constraint(
        (
            ParamDefinition("__TARGET", "TARGET", "I3HD1", (), 0.0, -10.0, 10.0),
            ParamDefinition("__VALUE_EXACT", "VALUE", "I3HD1", (), 2.0, -10.0, 10.0),
            ParamDefinition(
                "__VALUE_COMPANION", "VALUE", "I3CD1", (), 3.0, -10.0, 10.0
            ),
        ),
        "[TARGET] = [VALUE, NUC->I3CD1]",
    )

    assert parameterization.program.constraints[0].dependencies == (
        "__VALUE_COMPANION",
    )


def test_same_group_companion_reference_outranks_global_fallback() -> None:
    parameterization = _compile_context_constraint(
        (
            ParamDefinition("__TARGET", "TARGET", "I3HD1", (), 0.0, -10.0, 10.0),
            ParamDefinition(
                "__VALUE_COMPANION", "VALUE", "I3CD1", (), 3.0, -10.0, 10.0
            ),
            ParamDefinition("__VALUE_GLOBAL", "VALUE", "", (), 4.0, -10.0, 10.0),
        ),
        "[TARGET] = [VALUE]",
    )

    assert parameterization.program.constraints[0].dependencies == (
        "__VALUE_COMPANION",
    )


def test_same_group_companion_outranks_condition_specific_global() -> None:
    parameterization = _compile_context_constraint(
        (
            ParamDefinition(
                "__TARGET",
                "TARGET",
                "I3HD1",
                (("h_larmor_frq", 800.0),),
                0.0,
                -10.0,
                10.0,
            ),
            ParamDefinition(
                "__VALUE_COMPANION", "VALUE", "I3CD1", (), 3.0, -10.0, 10.0
            ),
            ParamDefinition(
                "__VALUE_GLOBAL",
                "VALUE",
                "",
                (("h_larmor_frq", 800.0),),
                4.0,
                -10.0,
                10.0,
            ),
        ),
        "[TARGET] = [VALUE]",
    )

    assert parameterization.program.constraints[0].dependencies == (
        "__VALUE_COMPANION",
    )


def test_unrelated_group_reference_is_rejected() -> None:
    with pytest.raises(NoParameterMatchError, match="context-compatible"):
        _compile_context_constraint(
            (
                ParamDefinition("__TARGET", "TARGET", "I3HD1", (), 0.0, -10.0, 10.0),
                ParamDefinition(
                    "__VALUE_OTHER", "VALUE", "L7CD2", (), 3.0, -10.0, 10.0
                ),
            ),
            "[TARGET] = [VALUE]",
        )


def test_equally_specific_same_group_companions_remain_ambiguous() -> None:
    with pytest.raises(AmbiguousParameterReferenceError) as raised:
        _compile_context_constraint(
            (
                ParamDefinition("__TARGET", "TARGET", "I3HD1", (), 0.0, -10.0, 10.0),
                ParamDefinition("__VALUE_CD1", "VALUE", "I3CD1", (), 2.0, -10.0, 10.0),
                ParamDefinition("__VALUE_MD1", "VALUE", "I3MD1", (), 3.0, -10.0, 10.0),
            ),
            "[TARGET] = [VALUE]",
        )

    assert raised.value.context["candidate_ids"] == (
        "__VALUE_CD1",
        "__VALUE_MD1",
    )


def test_condition_specificity_is_preserved_for_same_group_companions() -> None:
    parameterization = _compile_context_constraint(
        (
            ParamDefinition(
                "__TARGET",
                "TARGET",
                "I3HD1",
                (("h_larmor_frq", 800.0),),
                0.0,
                -10.0,
                10.0,
            ),
            ParamDefinition(
                "__VALUE_800",
                "VALUE",
                "I3CD1",
                (("h_larmor_frq", 800.0),),
                2.0,
                -10.0,
                10.0,
            ),
            ParamDefinition(
                "__VALUE_GLOBAL_CONDITION",
                "VALUE",
                "I3CD1",
                (),
                3.0,
                -10.0,
                10.0,
            ),
            ParamDefinition(
                "__VALUE_600",
                "VALUE",
                "I3CD1",
                (("h_larmor_frq", 600.0),),
                4.0,
                -10.0,
                10.0,
            ),
        ),
        "[TARGET] = [VALUE]",
    )

    assert parameterization.program.constraints[0].dependencies == ("__VALUE_800",)


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
    def explode(_value: float) -> dict[str, float]:
        msg = "backend exploded"
        raise RuntimeError(msg)

    monkeypatch.setitem(
        user_function_registry.user_function_registry,
        "test_model",
        {"explode": explode},
    )
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False, 'explode(__A)["result"]'),
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


def test_scientific_function_identity_tracks_same_metadata_implementation_changes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def implementation_one(value: float) -> dict[str, float]:
        return {"result": value * 2.0}

    def implementation_two(value: float) -> dict[str, float]:
        return {"result": value * 3.0}

    for implementation in (implementation_one, implementation_two):
        monkeypatch.setattr(implementation, "__name__", "same_metadata")
        monkeypatch.setattr(implementation, "__qualname__", "same_metadata")

    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration(
            "__B",
            False,
            'same_metadata(__A)["result"]',
        ),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 2.0, "__B": 0.0},
        model_name="test_model",
    )
    monkeypatch.setitem(
        user_function_registry.user_function_registry,
        "test_model",
        {"same_metadata": implementation_one},
    )
    first = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        {"__B"},
    )
    first_frame = first.frame_from_snapshot(snapshot)

    monkeypatch.setitem(
        user_function_registry.user_function_registry,
        "test_model",
        {"same_metadata": implementation_two},
    )
    second = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        {"__B"},
    )

    assert first.binder.identity != second.binder.identity
    assert first.program.fingerprint != second.program.fingerprint
    assert first.resolve(first_frame)["__B"] == 4.0
    with pytest.raises(ConstraintProgramMismatchError):
        second.resolve(first_frame)


def test_mutated_scientific_function_is_rejected_at_fresh_binding_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class MutableScientificFunction:
        def __init__(self) -> None:
            self.scale = 2.0

        def __call__(self, value: float) -> dict[str, float]:
            return {"result": value * self.scale}

    function = MutableScientificFunction()
    monkeypatch.setitem(
        user_function_registry.user_function_registry,
        "test_model",
        {"mutable": function},
    )
    declarations = (
        ParameterDeclaration("__A", False),
        ParameterDeclaration("__B", False, 'mutable(__A)["result"]'),
    )
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 2.0, "__B": 0.0},
        model_name="test_model",
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        {"__B"},
    )
    frame = parameterization.frame_from_snapshot(snapshot)
    assert parameterization.resolve(frame)["__B"] == 4.0

    function.scale = 3.0

    with pytest.raises(ConstraintProgramMismatchError) as raised:
        dataclasses.replace(
            parameterization,
            occurrence_identity="fresh-binding-occurrence",
        )

    assert raised.value.code == "program_mismatch"
    assert raised.value.context["function_id"] == "mutable"


def test_valid_scientific_functions_validate_once_per_bound_occurrence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    validation_calls = 0
    original_validate = ScientificFunctionBinder.validate_implementations

    def counted_validate(binder: ScientificFunctionBinder) -> None:
        nonlocal validation_calls
        validation_calls += 1
        original_validate(binder)

    monkeypatch.setattr(
        ScientificFunctionBinder,
        "validate_implementations",
        counted_validate,
    )
    declarations = (ParameterDeclaration("__A", False),)
    parameter_model, snapshot = _native_fixture(
        declarations,
        values={"__A": 2.0},
    )
    parameterization = compile_active_parameterization(
        parameter_model,
        snapshot,
        Method(),
        {"__A"},
    )
    frame = parameterization.frame_from_snapshot(snapshot)
    for _ in range(200):
        parameterization.resolve(frame)

    assert validation_calls == 1


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
    assert raised.value.context["param_id"] == "__B"


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


def test_native_compilation_failure_cannot_fall_back_to_legacy_fit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = AnalysisSession.create()
    preview_calls = 0

    def fail_preview(*_args: object, **_kwargs: object) -> None:
        nonlocal preview_calls
        preview_calls += 1
        msg = "native compilation failed"
        raise RuntimeError(msg)

    monkeypatch.setattr(session, "compile_parameterization_from_actions", fail_preview)
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

    with pytest.raises(RuntimeError, match="native compilation failed"):
        chemex_module.run(args, session=session)

    assert not (tmp_path / "Data").exists()
    assert not (tmp_path / "Parameters").exists()
    assert preview_calls == 1

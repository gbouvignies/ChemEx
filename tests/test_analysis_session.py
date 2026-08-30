from __future__ import annotations

from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemex import chemex as chemex_module
from chemex.configuration.method_input import prepare_method_plan
from chemex.configuration.method_plan import (
    FitAction,
    FixAction,
    FormatOrigin,
    GridAxis,
    GridSearch,
    GridValues,
    McmcRequest,
    MethodPlan,
    ParameterSelector,
    ProfileSelection,
    ResamplingRequest,
    SourceRef,
    StatisticsPlan,
    StepPlan,
)
from chemex.configuration.methods import McmcSettings, Method, Selection, Statistics
from chemex.experiments import builder as builder_module
from chemex.nmr.basis import Basis
from chemex.optimize import fitting as fitting_module
from chemex.optimize import method_plan_execution as method_execution_module
from chemex.optimize import resampling as resampling_module
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting
from chemex.printers import parameters as parameter_printer_module
from chemex.runtime import AnalysisSession, ExecutionSettings
from chemex.runtime import session as session_module


class StubModel:
    def __init__(self) -> None:
        self.reset_calls = 0
        self.model_names: list[str] = []

    def reset(self) -> None:
        self.reset_calls += 1

    def set_model(self, name: str) -> None:
        self.model_names.append(name)


class StubParameters:
    def __init__(self) -> None:
        self.reset_calls = 0
        self.defaults_calls: list[object] = []
        self.fix_all_calls = 0
        self.sort_calls = 0

    def reset(self) -> None:
        self.reset_calls += 1

    def set_defaults(self, defaults: object) -> None:
        self.defaults_calls.append(defaults)

    def fix_all(self) -> None:
        self.fix_all_calls += 1

    def sort(self) -> None:
        self.sort_calls += 1

    def get_parameters(self, _param_ids: object) -> dict[str, ParamSetting]:
        return {}


class StubParameterization:
    independent_ids: tuple[str, ...] = ()

    @staticmethod
    def frame_from_snapshot(snapshot: object) -> object:
        return snapshot

    @staticmethod
    def resolve(_frame: object) -> dict[str, float]:
        return {}


class StubParameterFactory:
    def __init__(self) -> None:
        self.clear_calls = 0
        self.seal_definition_calls = 0
        self.seal_configuration_calls = 0
        self.native_sealing_succeeds = True
        self.sealed_parameter_model = SimpleNamespace(declarations={}, definitions=())

    def clear_cache(self) -> None:
        self.clear_calls += 1

    def reset(self) -> None:
        self.clear_calls += 1

    def try_seal_definitions(self) -> bool:
        self.seal_definition_calls += 1
        return self.native_sealing_succeeds

    def try_seal_configuration(self) -> bool:
        self.seal_configuration_calls += 1
        return self.native_sealing_succeeds


class StubSession:
    def __init__(self) -> None:
        self.parameters = StubParameters()
        self.model = SimpleNamespace(spec=object())
        self.model_names: list[str] = []
        self.parameter_factory = StubParameterFactory()
        self.analysis_values = SimpleNamespace(snapshot=lambda: object())
        self.execution = ExecutionSettings()
        self.build_analysis_values_calls = 0

    def set_model(self, name: str) -> None:
        self.model_names.append(name)

    def try_build_analysis_values(self) -> bool:
        self.build_analysis_values_calls += 1
        return self.parameter_factory.try_seal_configuration()

    def resolve_current_values(self, _required_ids: set[str]) -> dict[str, float]:
        return {}

    def compile_parameterization(
        self,
        _method: Method,
        _required_ids: set[str],
    ) -> StubParameterization:
        return StubParameterization()

    def validate_method_plan(self, _plan: MethodPlan) -> None:
        return None


class FakeExperiments:
    def __init__(self) -> None:
        self.filtered = 0
        self.selections: list[ProfileSelection] = []
        self.param_ids: set[str] = set()

    def filter(self) -> None:
        self.filtered += 1

    def filter_from_values(self, _values: object) -> None:
        self.filtered += 1

    def select_profiles(self, selection: ProfileSelection) -> None:
        self.selections.append(selection)

    def __bool__(self) -> bool:
        return True


class EmptyAfterSelectExperiments(FakeExperiments):
    def __init__(self) -> None:
        super().__init__()
        self.has_profiles = True

    def select_profiles(self, selection: ProfileSelection) -> None:
        super().select_profiles(selection)
        self.has_profiles = False

    def __bool__(self) -> bool:
        return self.has_profiles


class DummyExperiment:
    def __init__(self, filename: Path) -> None:
        self.filename = filename

    def __len__(self) -> int:
        return 1


EXPECTED_EXPERIMENT_COUNT = 2


def make_args(command: str) -> Namespace:
    args = Namespace(
        commands=command,
        model="2st",
        include=None,
        exclude=None,
        experiments=[Path("experiment.toml")],
        parameters=[Path("parameters.toml")],
        method=None,
        output=Path("Output"),
        plot="normal",
    )
    if command == "fit":
        args.workers = "auto"
        args.native_threads = "auto"
    return args


def test_analysis_session_lifecycle_calls_reset_hooks(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = StubModel()
    parameters = StubParameters()
    parameter_factory = StubParameterFactory()
    plugin_calls: list[str] = []

    monkeypatch.setattr(
        session_module,
        "ensure_plugins_registered",
        lambda: plugin_calls.append("plugins"),
    )

    session = AnalysisSession(
        model=model,
        parameters=parameters,
        parameter_factory=parameter_factory,
    )

    session.reset()
    session.set_model("3st")

    np.testing.assert_equal(plugin_calls, ["plugins"])
    np.testing.assert_equal(parameter_factory.clear_calls, 2)
    np.testing.assert_equal(parameters.reset_calls, 1)
    np.testing.assert_equal(model.reset_calls, 1)
    np.testing.assert_equal(model.model_names, ["3st"])


def test_analysis_sessions_own_distinct_runtime_state() -> None:
    session_a = AnalysisSession.create()
    session_b = AnalysisSession.create()

    session_a.set_model("3st")
    session_b.set_model("2st")

    basis_a = Basis(type="ixy", model=session_a.model.spec)
    basis_b = Basis(type="ixy", model=session_b.model.spec)

    param = ParamSetting(ParamName("PB"), value=0.15, vary=False)
    session_a.parameters.add_multiple({param.id_: param})

    np.testing.assert_equal(int(session_a.parameters is session_b.parameters), 0)
    np.testing.assert_equal(
        int(session_a.parameter_factory is session_b.parameter_factory),
        0,
    )
    np.testing.assert_equal(
        int(session_a.analysis_values is session_b.analysis_values), 0
    )
    np.testing.assert_equal(basis_a.model.states, "abc")
    np.testing.assert_equal(basis_b.model.states, "ab")
    np.testing.assert_equal(
        session_b.parameters.get_parameters([param.id_]),
        {},
    )


def test_ensure_plugins_registered_runs_once(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []

    session_module.reset_plugin_registration()
    monkeypatch.setattr(
        session_module,
        "register_kinetic_settings",
        lambda: calls.append("models"),
    )
    monkeypatch.setattr(
        session_module,
        "register_experiments",
        lambda: calls.append("experiments"),
    )

    session_module.ensure_plugins_registered()
    session_module.ensure_plugins_registered()

    np.testing.assert_equal(calls, ["models", "experiments"])


def test_build_experiments_uses_session_construction_state(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    opened: list[Path] = []
    build_calls: list[dict[str, object]] = []

    monkeypatch.setattr(
        builder_module.experiment_types,
        "open",
        lambda filename: (
            opened.append(filename)
            or SimpleNamespace(
                filename=filename,
                experiment_type_name="test",
            )
        ),
    )

    def fake_build(source: SimpleNamespace, **kwargs: object) -> SimpleNamespace:
        build_calls.append(kwargs)
        return SimpleNamespace(
            experiment=DummyExperiment(source.filename),
            notices=(),
        )

    monkeypatch.setattr(
        builder_module.experiment_types,
        "build",
        fake_build,
    )

    experiments = builder_module.build_experiments(
        [Path("a.toml"), Path("b.toml")],
        Selection(include=None, exclude=None),
        session=session,
    )

    np.testing.assert_equal(session.parameters.sort_calls, 1)
    np.testing.assert_equal(len(list(experiments)), EXPECTED_EXPERIMENT_COUNT)
    np.testing.assert_equal(opened, [Path("a.toml"), Path("b.toml")])
    assert all(call["parameters"] is session.parameter_factory for call in build_calls)
    assert all(call["model"] is session.model.spec for call in build_calls)
    np.testing.assert_equal(session.parameter_factory.seal_definition_calls, 1)


def test_build_experiments_seals_empty_definition_set() -> None:
    session = StubSession()

    experiments = builder_module.build_experiments(
        None,
        Selection(include=None, exclude=None),
        session=session,
    )

    assert not experiments
    np.testing.assert_equal(session.parameter_factory.seal_definition_calls, 1)


def test_run_uses_explicit_session_for_fit_flow(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    defaults = object()
    session = StubSession()
    experiments = FakeExperiments()
    recorded: dict[str, object] = {}
    recorded_env: dict[str, str] = {}
    experiment_path = tmp_path / "experiment.toml"
    parameter_path = tmp_path / "parameters.toml"
    experiment_bytes = b'[experiment]\r\nname = "captured-before-parse"\r\n'
    experiment_path.write_bytes(experiment_bytes)
    parameter_path.write_text("[GLOBAL]\nPB = 0.1\n", encoding="utf-8")

    def fake_build_experiments(
        filenames: list[Path] | None,
        selection: Selection,
        *,
        session: StubSession | None = None,
    ) -> FakeExperiments:
        experiment_path.write_text("[changed]\nvalue = true\n", encoding="utf-8")
        recorded["build"] = (filenames, selection, session)
        return experiments

    def fake_execute_method_plan(
        experiments_arg: FakeExperiments,
        methods: MethodPlan,
        path: Path,
        plot_level: str,
        *,
        session: StubSession,
        run_info: object,
    ) -> None:
        recorded["execute_method_plan"] = (
            experiments_arg,
            methods,
            path,
            plot_level,
            session,
            run_info,
        )

    monkeypatch.setattr(chemex_module, "build_experiments", fake_build_experiments)
    monkeypatch.setattr(chemex_module, "read_defaults", lambda _filenames: defaults)
    monkeypatch.setattr(
        chemex_module,
        "execute_method_plan",
        fake_execute_method_plan,
    )
    fake_run_info = SimpleNamespace(
        write_outcome=lambda status, _snapshot, **_kwargs: recorded.setdefault(
            "outcome", status
        )
    )

    def fake_write_run_info(*_args, **kwargs):
        recorded["captured_inputs"] = kwargs["input_files"]
        return recorded.setdefault("run_info", fake_run_info)

    monkeypatch.setattr(chemex_module, "write_run_info", fake_write_run_info)
    monkeypatch.setattr(chemex_module.os, "environ", recorded_env)

    args = make_args("fit")
    args.experiments = [experiment_path]
    args.parameters = [parameter_path]
    args.workers = 3
    args.native_threads = 2
    chemex_module.run(args, session=session)

    np.testing.assert_equal(session.model_names, ["2st"])
    np.testing.assert_equal(
        session.execution, ExecutionSettings(workers=3, native_threads=2)
    )
    np.testing.assert_equal(recorded_env["OMP_NUM_THREADS"], "2")
    np.testing.assert_equal(recorded_env["VECLIB_MAXIMUM_THREADS"], "2")
    np.testing.assert_equal(session.parameters.defaults_calls, [defaults])
    np.testing.assert_equal(session.parameter_factory.seal_configuration_calls, 1)
    np.testing.assert_equal(session.build_analysis_values_calls, 1)
    np.testing.assert_equal(experiments.filtered, 1)
    np.testing.assert_equal(recorded["build"][2], session)
    assert recorded["execute_method_plan"][4] is session
    assert recorded["execute_method_plan"][5] is fake_run_info
    assert recorded["run_info"] is fake_run_info
    captured_inputs = recorded["captured_inputs"]
    assert captured_inputs[0].content == experiment_bytes
    assert recorded["outcome"] == "complete"


def test_run_fails_closed_when_native_configuration_is_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    session.parameter_factory.native_sealing_succeeds = False
    experiments = FakeExperiments()
    fit_ran: list[bool] = []

    monkeypatch.setattr(
        chemex_module,
        "build_experiments",
        lambda *_args, **_kwargs: experiments,
    )
    monkeypatch.setattr(chemex_module, "read_defaults", lambda _filenames: object())
    monkeypatch.setattr(chemex_module, "capture_input_files", lambda _args: ())
    monkeypatch.setattr(
        chemex_module,
        "execute_method_plan",
        lambda *_args, **_kwargs: fit_ran.append(True),
    )

    with pytest.raises(RuntimeError, match="Native parameter initialization failed"):
        chemex_module.run(make_args("fit"), session=session)

    assert fit_ran == []
    assert session.parameter_factory.seal_configuration_calls == 1
    assert session.build_analysis_values_calls == 1


def test_run_fit_prepares_method_plan_before_run_side_effects(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    events: list[str] = []

    def reject_method(*_args: object) -> MethodPlan:
        events.append("prepare")
        raise ValueError("invalid method plan")

    monkeypatch.setattr(chemex_module, "prepare_method_plan", reject_method)
    monkeypatch.setattr(
        chemex_module,
        "write_run_info",
        lambda *_args, **_kwargs: events.append("run_info"),
    )
    monkeypatch.setattr(
        chemex_module,
        "invalidate_planned_outputs",
        lambda *_args: events.append("invalidate"),
    )

    with pytest.raises(ValueError, match="invalid method plan"):
        chemex_module.run_fit(
            make_args("fit"),
            FakeExperiments(),
            session,  # type: ignore[arg-type]
            input_files=(),
            methods={"": Method()},
        )

    assert events == ["prepare"]


def test_run_uses_explicit_session_for_simulation_flow(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    defaults = object()
    session = StubSession()
    experiments = FakeExperiments()
    recorded: dict[str, object] = {}

    def fake_build_experiments(
        filenames: list[Path] | None,
        selection: Selection,
        *,
        session: StubSession | None = None,
    ) -> FakeExperiments:
        recorded["build"] = (filenames, selection, session)
        return experiments

    def fake_execute_simulation(
        experiments_arg: FakeExperiments,
        path: Path,
        *,
        parameter_values: object,
        parameter_model: object,
        parameterization: object,
        plot: bool = False,
    ) -> None:
        recorded["simulation"] = (
            experiments_arg,
            path,
            parameter_values,
            parameter_model,
            parameterization,
            plot,
        )

    monkeypatch.setattr(chemex_module, "build_experiments", fake_build_experiments)
    monkeypatch.setattr(chemex_module, "read_defaults", lambda _filenames: defaults)
    monkeypatch.setattr(chemex_module, "execute_simulation", fake_execute_simulation)

    chemex_module.run(make_args("simulate"), session=session)

    np.testing.assert_equal(session.model_names, ["2st"])
    np.testing.assert_equal(session.parameters.defaults_calls, [defaults])
    np.testing.assert_equal(session.parameter_factory.seal_configuration_calls, 1)
    np.testing.assert_equal(session.build_analysis_values_calls, 1)
    np.testing.assert_equal(session.parameters.fix_all_calls, 0)
    np.testing.assert_equal(recorded["build"][2], session)
    simulation = recorded["simulation"]
    assert simulation[:3] == (experiments, Path("Output"), {})
    assert simulation[3] is session.parameter_factory.sealed_parameter_model
    assert isinstance(simulation[4], StubParameterization)
    assert simulation[5] is True


def test_main_bootstraps_plugins_for_non_run_commands(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[str] = []

    class StubParser:
        def parse_args(self) -> Namespace:
            return Namespace(
                analysis_command=False,
                func=lambda _args: calls.append("func"),
            )

    def build_parser() -> StubParser:
        return StubParser()

    monkeypatch.setattr(chemex_module, "print_logo", lambda: calls.append("logo"))
    monkeypatch.setattr(
        chemex_module,
        "ensure_plugins_registered",
        lambda: calls.append("plugins"),
    )
    monkeypatch.setattr(chemex_module, "build_parser", build_parser)

    chemex_module.main()

    np.testing.assert_equal(calls, ["logo", "plugins", "func"])


def test_main_dispatches_analysis_commands(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []

    class StubParser:
        def parse_args(self) -> Namespace:
            return Namespace(
                analysis_command=True, func=lambda _args: calls.append("func")
            )

    def build_parser() -> StubParser:
        return StubParser()

    monkeypatch.setattr(chemex_module, "print_logo", lambda: calls.append("logo"))
    monkeypatch.setattr(
        chemex_module,
        "ensure_plugins_registered",
        lambda: calls.append("plugins"),
    )
    monkeypatch.setattr(chemex_module, "build_parser", build_parser)
    monkeypatch.setattr(
        chemex_module,
        "run",
        lambda _args, **_kwargs: calls.append("run"),
    )

    chemex_module.main()

    np.testing.assert_equal(calls, ["logo", "plugins", "run"])


def test_run_methods_skips_fit_when_selection_removes_all_profiles(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    experiments = EmptyAfterSelectExperiments()
    calls: list[str] = []

    monkeypatch.setattr(
        method_execution_module, "print_no_data", lambda: calls.append("no_data")
    )

    def fail_if_called(*_args, **_kwargs) -> None:
        pytest.fail("native fitting should not run when no profiles remain selected")

    monkeypatch.setattr(
        method_execution_module,
        "run_native_deterministic",
        fail_if_called,
    )

    fitting_module.run_methods(
        experiments,
        {"": Method(include=["1H"])},
        Path("Output"),
        "normal",
        session=session,
    )

    np.testing.assert_equal(calls, ["no_data"])
    np.testing.assert_equal(len(experiments.selections), 1)


def test_legacy_methods_compatibility_keeps_implicit_selection_inheritance(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(MethodPlan, "validate", lambda *_args: None)
    plan = prepare_method_plan(
        {
            "FIRST": Method(include=["1H"]),
            "SECOND": Method(),
        },
        object(),  # type: ignore[arg-type]
    )

    assert plan.steps[1].selection == plan.steps[0].selection


def test_programmatic_v1_methods_normalize_grid_and_statistics_canonically(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(MethodPlan, "validate", lambda *_args: None)
    plan = prepare_method_plan(
        {
            "STEP": Method(
                fix=["PB"],
                grid=["[PB] = (0.1, 0.2)"],
                statistics=Statistics(
                    mc=2,
                    mcmc=McmcSettings(
                        steps=10,
                        burn=2,
                        thin=2,
                        walkers=4,
                        seed=9,
                        workers=2,
                    ),
                ),
            )
        },
        object(),  # type: ignore[arg-type]
    )

    step = plan.steps[0]
    assert isinstance(step.search, GridSearch)
    assert step.statistics is not None
    assert step.statistics.mc == ResamplingRequest(2, seed=0)
    assert step.statistics.mcmc == McmcRequest(
        steps=10,
        burn=2,
        thin=2,
        walkers=4,
        seed=9,
        workers=2,
    )
    assert prepare_method_plan(plan, object()) is plan  # type: ignore[arg-type]


def test_v2_omitted_selection_restores_the_step_local_all_profiles_baseline(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    experiments = FakeExperiments()
    plan = MethodPlan(
        FormatOrigin.V2,
        (
            StepPlan("FIRST", selection=ProfileSelection(("1H",), None)),
            StepPlan("SECOND"),
        ),
    )
    session.compile_parameterization_from_actions = lambda *_args: object()  # type: ignore[attr-defined]
    monkeypatch.setattr(
        method_execution_module,
        "run_native_deterministic",
        lambda *_args, **_kwargs: None,
    )

    fitting_module.run_methods(
        experiments,
        plan,
        Path("Output"),
        "normal",
        session=session,
    )

    assert experiments.selections == [
        ProfileSelection(("1H",), None),
        ProfileSelection(),
    ]


@pytest.mark.parametrize("origin", tuple(FormatOrigin))
def test_canonical_grid_runs_requested_statistics_independently_of_origin(
    monkeypatch: pytest.MonkeyPatch,
    origin: FormatOrigin,
) -> None:
    session = StubSession()
    experiments = FakeExperiments()
    source = SourceRef(Path("method.toml"), "STEP", "SEARCH.GRID.AXES", 0)
    selector = ParameterSelector("PB", source=source)
    plan = MethodPlan(
        origin,
        (
            StepPlan(
                "STEP",
                search=GridSearch(
                    (GridAxis(selector, GridValues((0.1, 0.2)), source),)
                ),
                statistics=StatisticsPlan(mc=ResamplingRequest(1, seed=7)),
            ),
        ),
    )
    accepted_fit = object()
    statistics_calls: list[tuple[object, object]] = []
    session.compile_parameterization_from_actions = lambda *_args: object()  # type: ignore[attr-defined]
    monkeypatch.setattr(
        method_execution_module,
        "run_native_deterministic",
        lambda *_args, **_kwargs: accepted_fit,
    )
    monkeypatch.setattr(
        method_execution_module,
        "_run_requested_native_statistics",
        lambda _experiments, _path, statistics, fit, **_kwargs: statistics_calls.append(
            (statistics, fit)
        ),
    )

    method_execution_module.execute_method_plan(
        experiments,
        plan,
        Path("Output"),
        "normal",
        session=session,
    )

    assert statistics_calls == [(plan.steps[0].statistics, accepted_fit)]


def test_v2_omitted_resampling_seed_is_resolved_once_for_the_occurrence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    resolved_seed = 0x1234_5678_90AB_CDEF
    calls: list[int] = []
    monkeypatch.setattr(
        method_execution_module.secrets,
        "randbits",
        lambda _bits: resolved_seed,
    )
    monkeypatch.setattr(
        method_execution_module,
        "run_native_resampling",
        lambda *_args, root_seed, **_kwargs: calls.append(root_seed),
    )
    session = StubSession()
    snapshot = object()
    session.analysis_values = SimpleNamespace(snapshot=lambda: snapshot)  # type: ignore[attr-defined]

    method_execution_module._run_native_statistics(
        FakeExperiments(),
        Path("Output"),
        StatisticsPlan(mc=ResamplingRequest(1)),
        object(),  # type: ignore[arg-type]
        session=session,  # type: ignore[arg-type]
    )

    assert calls == [resolved_seed]


def test_executor_delegates_canonical_statistics_in_fixed_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    snapshot = object()
    session.analysis_values = SimpleNamespace(snapshot=lambda: snapshot)  # type: ignore[attr-defined]
    parameterization = object()
    session.compile_parameterization_from_actions = (  # type: ignore[attr-defined]
        lambda *_args: parameterization
    )
    experiments = FakeExperiments()
    fit = object()
    requests = (
        ResamplingRequest(2, seed=11),
        ResamplingRequest(3, seed=12),
        ResamplingRequest(4, seed=13),
    )
    mcmc = McmcRequest(steps=10, burn=2, seed=14, thin=2, walkers=4, workers=2)
    plan = MethodPlan(
        FormatOrigin.V2,
        (
            StepPlan(
                "STEP",
                statistics=StatisticsPlan(*requests, mcmc),
            ),
        ),
    )
    calls: list[tuple[object, ...]] = []

    def deterministic(
        _experiments: object,
        path: Path,
        plot: str,
        *,
        parameterization: object,
        search: object,
        **_kwargs: object,
    ) -> object:
        calls.append(("deterministic", path, plot, parameterization, search))
        return fit

    def resampling(
        _experiments: object,
        path: Path,
        kind: str,
        request: ResamplingRequest,
        actual_fit: object,
        **_kwargs: object,
    ) -> None:
        calls.append((kind, path, request, actual_fit))

    def run_mcmc(
        _experiments: object,
        actual_fit: object,
        request: McmcRequest,
        path: Path,
        **_kwargs: object,
    ) -> None:
        calls.append(("mcmc", path, request, actual_fit))

    monkeypatch.setattr(
        method_execution_module,
        "run_native_deterministic",
        deterministic,
    )
    monkeypatch.setattr(method_execution_module, "run_native_resampling", resampling)
    monkeypatch.setattr(method_execution_module, "run_native_mcmc", run_mcmc)

    method_execution_module.execute_method_plan(
        experiments,
        plan,
        Path("Output"),
        "normal",
        session=session,
    )

    assert calls == [
        ("deterministic", Path("Output"), "normal", parameterization, None),
        ("mc", Path("Output"), requests[0], fit),
        ("bs", Path("Output"), requests[1], fit),
        ("bsn", Path("Output"), requests[2], fit),
        ("mcmc", Path("Output"), mcmc, fit),
    ]


def test_failed_plan_reinvocation_starts_fresh_from_committed_scientific_state(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    experiments = FakeExperiments()
    experiments.param_ids = {"__PB"}
    source = SourceRef(Path("method.toml"), "STEP", "ROLES")
    fixed = FixAction((ParameterSelector("PB"),), source)
    fitted = FitAction((ParameterSelector("PB"),), source)
    plan = MethodPlan(
        FormatOrigin.V2,
        (
            StepPlan("FIRST", role_actions=(fixed,)),
            StepPlan("SECOND", roles_from="FIRST", role_actions=(fitted,)),
            StepPlan("THIRD", role_actions=(fixed,)),
        ),
    )
    compiled: list[tuple[object, ...]] = []
    starting_revisions: list[int] = []
    committed_revision = 0
    execution_count = 0

    def compile_actions(actions: tuple[object, ...], _required: set[str]) -> object:
        compiled.append(actions)
        return actions

    def fail_once(*_args: object, **_kwargs: object) -> None:
        nonlocal committed_revision, execution_count
        starting_revisions.append(committed_revision)
        execution_count += 1
        if execution_count == 2:
            raise RuntimeError("injected second-step failure")
        committed_revision += 1

    session.compile_parameterization_from_actions = compile_actions  # type: ignore[attr-defined]
    monkeypatch.setattr(
        method_execution_module,
        "run_native_deterministic",
        fail_once,
    )

    with pytest.raises(RuntimeError, match="second-step failure"):
        method_execution_module.execute_method_plan(
            experiments, plan, Path("Output"), "normal", session=session
        )
    method_execution_module.execute_method_plan(
        experiments, plan, Path("Output"), "normal", session=session
    )

    assert compiled == [
        (fixed,),
        (fixed, fitted),
        (fixed,),
        (fixed, fitted),
        (fixed,),
    ]
    assert starting_revisions == [0, 1, 1, 2, 3]
    assert committed_revision == 4


@pytest.mark.parametrize("origin", tuple(FormatOrigin))
def test_run_methods_compiles_each_canonical_step_without_origin_or_store_state(
    monkeypatch: pytest.MonkeyPatch,
    origin: FormatOrigin,
) -> None:
    session = StubSession()
    experiments = FakeExperiments()
    experiments.param_ids = {"__PB"}
    source = SourceRef(Path("method.toml"), "STEP", "ROLES")
    selector = ParameterSelector("PB")
    fixed = FixAction((selector,), source)
    fitted = FitAction((selector,), source)
    plan = MethodPlan(
        origin,
        (
            StepPlan("FIRST", role_actions=(fixed,)),
            StepPlan("SECOND", roles_from="FIRST", role_actions=(fitted,)),
            StepPlan("THIRD", role_actions=(fitted,)),
        ),
    )
    compiled: list[tuple[object, ...]] = []
    executed: list[tuple[str, object]] = []

    def compile_actions(
        actions: tuple[object, ...],
        required_ids: set[str],
    ) -> object:
        compiled.append(actions)
        assert required_ids == {"__PB"}
        return actions

    session.compile_parameterization_from_actions = compile_actions  # type: ignore[attr-defined]

    def run_deterministic(
        _experiments: object,
        path: Path,
        _plot: str,
        *,
        session: object,
        parameterization: object,
        search: object,
        run_info: object,
        step_name: str,
    ) -> None:
        assert session is not None
        assert search is None
        assert run_info is None
        assert step_name in {"FIRST", "SECOND", "THIRD"}
        executed.append((path.name, parameterization))

    monkeypatch.setattr(
        method_execution_module,
        "run_native_deterministic",
        run_deterministic,
    )

    method_execution_module.execute_method_plan(
        experiments,
        plan,
        Path("Output"),
        "normal",
        session=session,
    )

    assert compiled == [(fixed,), (fixed, fitted), (fitted,)]
    assert executed == [
        ("FIRST", (fixed,)),
        ("SECOND", (fixed, fitted)),
        ("THIRD", (fitted,)),
    ]
    assert session.parameters.fix_all_calls == 0


def test_resampling_summary_and_correlations_are_written(tmp_path: Path) -> None:
    samples = np.array([[0.1, 200.0], [0.3, 300.0], [0.5, 400.0]])
    parameter_names = ("[PB]", "[KEX_AB]")

    resampling_module._write_resampling_summary(
        tmp_path,
        parameter_names=parameter_names,
        samples=samples,
    )
    resampling_module._write_resampling_correlations(
        tmp_path,
        parameter_names=parameter_names,
        samples=samples,
    )

    summary = (tmp_path / "summary.toml").read_text(encoding="utf-8")
    correlations = (tmp_path / "correlations.tsv").read_text(encoding="utf-8")

    assert '["PB"]' in summary
    assert 'interval = "95% percentile"' in summary
    assert "sample_count = 3" in summary
    assert "median = 3.00000e-01" in summary
    assert "[PB]" in correlations
    assert "[KEX_AB]" in correlations
    assert "1.00000e+00" in correlations


def test_write_file_rejects_unknown_parameter_status(tmp_path: Path) -> None:
    parameter = ParamSetting(ParamName("PB"), value=0.15, vary=False)
    report_parameter = parameter_printer_module.ReportParameter(
        parameter.id_,
        parameter.param_name,
        parameter.value,
    )
    parameters = parameter_printer_module.GlobalLocalParameters(
        {parameter.param_name: report_parameter},
        {},
    )

    with pytest.raises(ValueError, match="Unknown parameter status"):
        parameter_printer_module.write_file(
            parameters,
            "typo",
            tmp_path,
        )


def test_constraint_output_replaces_prefix_colliding_parameter_ids_atomically() -> None:
    rendered = parameter_printer_module._replace_parameter_ids(
        "__A_SPIN + __A",
        {
            "__A": ParamName("A"),
            "__A_SPIN": ParamName("A_SPIN"),
        },
    )

    assert rendered == "[A_SPIN] + [A]"

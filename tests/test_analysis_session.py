from __future__ import annotations

from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemex import chemex as chemex_module
from chemex.configuration.methods import Method, Selection
from chemex.experiments import builder as builder_module
from chemex.nmr.basis import Basis
from chemex.optimize import fitting as fitting_module
from chemex.optimize import helper as helper_module
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
        self.status_calls: list[Method] = []

    def reset(self) -> None:
        self.reset_calls += 1

    def set_defaults(self, defaults: object) -> None:
        self.defaults_calls.append(defaults)

    def fix_all(self) -> None:
        self.fix_all_calls += 1

    def sort(self) -> None:
        self.sort_calls += 1

    def set_parameter_status(self, method: Method) -> None:
        self.status_calls.append(method)


class StubParameterFactory:
    def __init__(self) -> None:
        self.clear_calls = 0
        self.seal_definition_calls = 0
        self.seal_configuration_calls = 0
        self.native_sealing_succeeds = True

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
        self.execution = ExecutionSettings()
        self.build_analysis_values_calls = 0

    def set_model(self, name: str) -> None:
        self.model_names.append(name)

    def try_build_analysis_values(self) -> bool:
        self.build_analysis_values_calls += 1
        return self.parameter_factory.try_seal_configuration()

    def try_compile_parameterization(
        self,
        _method: object,
        _required_ids: set[str],
    ) -> None:
        return None

    def resolve_current_values(self, _required_ids: set[str]) -> dict[str, float]:
        return {}


class WriterParameterStore:
    def __init__(self, parameters: dict[str, ParamSetting]) -> None:
        self.parameters = parameters

    def get_parameters(self, param_ids: object) -> dict[str, ParamSetting]:
        return {
            param_id: parameter
            for param_id, parameter in self.parameters.items()
            if param_id in set(param_ids)
        }


class StatisticsParameterStore(StubParameters):
    def get_parameters(self, param_ids: object) -> dict[str, ParamSetting]:
        parameters = {
            "__PB": ParamSetting(ParamName("PB"), value=0.15, vary=True),
            "__KEX_AB": ParamSetting(ParamName("KEX_AB"), value=250.0, vary=True),
        }
        return {
            param_id: parameter
            for param_id, parameter in parameters.items()
            if param_id in set(param_ids)
        }


class FakeExperiments:
    def __init__(self, parameter_store: object | None = None) -> None:
        self.filtered = 0
        self.selections: list[Selection] = []
        self.parameter_store = parameter_store
        self.param_ids: set[str] = set()

    def filter(self) -> None:
        self.filtered += 1

    def filter_from_values(self, _values: object) -> None:
        self.filtered += 1

    def select(self, selection: Selection) -> None:
        self.selections.append(selection)

    def __bool__(self) -> bool:
        return True


class EmptyAfterSelectExperiments(FakeExperiments):
    def __init__(self, parameter_store: object | None = None) -> None:
        super().__init__(parameter_store=parameter_store)
        self.has_profiles = True

    def select(self, selection: Selection) -> None:
        super().select(selection)
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


def test_build_experiments_uses_session_parameter_store(
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
) -> None:
    defaults = object()
    session = StubSession()
    experiments = FakeExperiments(parameter_store=session.parameters)
    recorded: dict[str, object] = {}
    recorded_env: dict[str, str] = {}

    def fake_build_experiments(
        filenames: list[Path] | None,
        selection: Selection,
        *,
        session: StubSession | None = None,
    ) -> FakeExperiments:
        recorded["build"] = (filenames, selection, session)
        return experiments

    def fake_run_methods(
        experiments_arg: FakeExperiments,
        methods: dict[str, Method],
        path: Path,
        plot_level: str,
        *,
        session: StubSession,
    ) -> None:
        recorded["run_methods"] = (
            experiments_arg,
            methods,
            path,
            plot_level,
            session,
        )

    monkeypatch.setattr(chemex_module, "build_experiments", fake_build_experiments)
    monkeypatch.setattr(chemex_module, "read_defaults", lambda _filenames: defaults)
    monkeypatch.setattr(chemex_module, "run_methods", fake_run_methods)
    monkeypatch.setattr(
        chemex_module,
        "write_run_info",
        lambda *_args, **_kwargs: recorded.setdefault("run_info", True),
    )
    monkeypatch.setattr(
        chemex_module,
        "write_run_outcome",
        lambda _path, status, **_kwargs: recorded.setdefault("outcome", status),
    )
    monkeypatch.setattr(chemex_module.os, "environ", recorded_env)

    args = make_args("fit")
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
    assert recorded["run_methods"][4] is session
    assert recorded["run_info"] is True
    assert recorded["outcome"] == "complete"


def test_run_fails_closed_when_native_configuration_is_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    session.parameter_factory.native_sealing_succeeds = False
    experiments = FakeExperiments(parameter_store=session.parameters)
    fit_ran: list[bool] = []

    monkeypatch.setattr(
        chemex_module,
        "build_experiments",
        lambda *_args, **_kwargs: experiments,
    )
    monkeypatch.setattr(chemex_module, "read_defaults", lambda _filenames: object())
    monkeypatch.setattr(
        chemex_module,
        "run_methods",
        lambda *_args, **_kwargs: fit_ran.append(True),
    )
    monkeypatch.setattr(chemex_module, "write_run_info", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        chemex_module,
        "write_run_outcome",
        lambda *_args, **_kwargs: None,
    )

    with pytest.raises(RuntimeError, match="Native parameter initialization failed"):
        chemex_module.run(make_args("fit"), session=session)

    assert fit_ran == []
    assert session.parameter_factory.seal_configuration_calls == 1
    assert session.build_analysis_values_calls == 1


def test_run_uses_explicit_session_for_simulation_flow(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    defaults = object()
    session = StubSession()
    experiments = FakeExperiments(parameter_store=session.parameters)
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
        plot: bool = False,
    ) -> None:
        recorded["simulation"] = (
            experiments_arg,
            path,
            parameter_values,
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
    np.testing.assert_equal(session.parameters.fix_all_calls, 1)
    np.testing.assert_equal(recorded["build"][2], session)
    assert recorded["simulation"] == (experiments, Path("Output"), {}, True)


def test_run_rejects_experiments_built_under_a_different_session(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    mismatched_store = StubParameters()
    experiments = FakeExperiments(parameter_store=mismatched_store)

    def fake_build_experiments(
        filenames: list[Path] | None,
        selection: Selection,
        *,
        session: StubSession | None = None,
    ) -> FakeExperiments:
        return experiments

    def fail_run_methods(*_args, **_kwargs) -> None:
        pytest.fail("run_methods should not run for a mismatched parameter store")

    def fail_execute_simulation(*_args, **_kwargs) -> None:
        pytest.fail(
            "execute_simulation should not run for a mismatched parameter store"
        )

    monkeypatch.setattr(chemex_module, "build_experiments", fake_build_experiments)
    monkeypatch.setattr(chemex_module, "run_methods", fail_run_methods)
    monkeypatch.setattr(chemex_module, "execute_simulation", fail_execute_simulation)

    with pytest.raises(ValueError, match="does not match the active session"):
        chemex_module.run(make_args("fit"), session=session)


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
    experiments = EmptyAfterSelectExperiments(parameter_store=session.parameters)
    calls: list[str] = []

    monkeypatch.setattr(
        fitting_module, "print_no_data", lambda: calls.append("no_data")
    )

    def fail_if_called(*_args, **_kwargs) -> None:
        pytest.fail("native fitting should not run when no profiles remain selected")

    monkeypatch.setattr(fitting_module, "run_native_deterministic", fail_if_called)

    fitting_module.run_methods(
        experiments,
        {"": Method(include=["1H"])},
        Path("Output"),
        "normal",
        session=session,
    )

    np.testing.assert_equal(calls, ["no_data"])
    np.testing.assert_equal(len(experiments.selections), 1)
    np.testing.assert_equal(session.parameters.status_calls, [])


def test_resampling_summary_and_correlations_are_written(tmp_path: Path) -> None:
    store = WriterParameterStore(
        {
            "__PB": ParamSetting(ParamName("PB"), value=0.15, vary=True),
            "__KEX_AB": ParamSetting(ParamName("KEX_AB"), value=250.0, vary=True),
        },
    )
    samples = np.array([[0.1, 200.0], [0.3, 300.0], [0.5, 400.0]])
    parameter_names = resampling_module._format_parameter_names(
        ["__PB", "__KEX_AB"],
        store,
    )

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


def test_execute_post_fit_writes_parameters_from_experiments_store(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    param = ParamSetting(ParamName("PB"), value=0.15, vary=False)
    experiments = SimpleNamespace(
        param_ids=[param.id_],
        parameter_store=WriterParameterStore({param.id_: param}),
        write=lambda _path: None,
    )
    monkeypatch.setattr(helper_module, "print_writing_results", lambda _path: None)
    monkeypatch.setattr(
        helper_module,
        "_write_statistics",
        lambda *_args, **_kwargs: None,
    )

    helper_module.execute_post_fit(
        experiments,
        tmp_path,
        residuals=np.array([], dtype=float),
        nvarys=0,
    )

    np.testing.assert_equal(int((tmp_path / "Parameters" / "fixed.toml").exists()), 1)


def test_write_file_rejects_unknown_parameter_status(tmp_path: Path) -> None:
    parameter = ParamSetting(ParamName("PB"), value=0.15, vary=False)
    parameters = parameter_printer_module.GlobalLocalParameters(
        {parameter.param_name: parameter},
        {},
    )

    with pytest.raises(ValueError, match="Unknown parameter status"):
        parameter_printer_module.write_file(
            parameters,
            "typo",
            tmp_path,
            WriterParameterStore({}),
        )

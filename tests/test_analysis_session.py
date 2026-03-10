from __future__ import annotations

from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemex import chemex as chemex_module
from chemex.configuration.methods import Method, Selection
from chemex.experiments import builder as builder_module
from chemex.optimize import fitting as fitting_module
from chemex.optimize import helper as helper_module
from chemex.parameters.name import ParamName
from chemex.parameters.setting import ParamSetting
from chemex.printers import parameters as parameter_printer_module
from chemex.runtime import AnalysisSession
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

    def reset_parameters(self) -> None:
        self.reset_calls += 1

    def set_param_defaults(self, defaults: object) -> None:
        self.defaults_calls.append(defaults)

    def fix_all_parameters(self) -> None:
        self.fix_all_calls += 1

    def sort_parameters(self) -> None:
        self.sort_calls += 1

    def set_parameter_status(self, method: Method) -> None:
        self.status_calls.append(method)


class StubSession:
    def __init__(self) -> None:
        self.parameters = StubParameters()
        self.model_names: list[str] = []

    def set_model(self, name: str) -> None:
        self.model_names.append(name)


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
    def build_lmfit_params(self, param_ids: object) -> dict[str, SimpleNamespace]:
        _ = param_ids
        return {"__PB": SimpleNamespace(name="PB", vary=True)}


class FakeExperiments:
    def __init__(self) -> None:
        self.filtered = 0
        self.selections: list[Selection] = []

    def filter(self) -> None:
        self.filtered += 1

    def select(self, selection: Selection) -> None:
        self.selections.append(selection)

    def __bool__(self) -> bool:
        return True


class DummyExperiment:
    def __init__(self, filename: Path) -> None:
        self.filename = filename


EXPECTED_EXPERIMENT_COUNT = 2


def make_args(command: str) -> Namespace:
    return Namespace(
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


def test_analysis_session_lifecycle_calls_reset_hooks(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clear_calls: list[str] = []
    model = StubModel()
    parameters = StubParameters()

    monkeypatch.setattr(
        session_module,
        "clear_settings_cache",
        lambda: clear_calls.append("clear"),
    )
    monkeypatch.setattr(
        session_module,
        "ensure_plugins_registered",
        lambda: clear_calls.append("plugins"),
    )

    session = AnalysisSession(model=model, parameters=parameters)

    session.reset()
    session.set_model("3st")

    np.testing.assert_equal(clear_calls, ["clear", "plugins", "clear"])
    np.testing.assert_equal(parameters.reset_calls, 1)
    np.testing.assert_equal(model.reset_calls, 1)
    np.testing.assert_equal(model.model_names, ["3st"])


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

    monkeypatch.setattr(
        builder_module,
        "build_experiment",
        lambda filename, _selection: DummyExperiment(filename),
    )

    def should_not_be_called() -> None:
        msg = "global parameter store should not be used"
        raise AssertionError(msg)

    monkeypatch.setattr(
        builder_module.database,
        "sort_parameters",
        should_not_be_called,
    )

    experiments = builder_module.build_experiments(
        [Path("a.toml"), Path("b.toml")],
        Selection(include=None, exclude=None),
        session=session,
    )

    np.testing.assert_equal(session.parameters.sort_calls, 1)
    np.testing.assert_equal(len(list(experiments)), EXPECTED_EXPERIMENT_COUNT)


def test_run_uses_explicit_session_for_fit_flow(
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

    def fake_run_methods(
        experiments_arg: FakeExperiments,
        methods: dict[str, Method],
        path: Path,
        plot_level: str,
        *,
        session: StubSession | None = None,
    ) -> None:
        recorded["run_methods"] = (experiments_arg, methods, path, plot_level, session)

    monkeypatch.setattr(chemex_module, "build_experiments", fake_build_experiments)
    monkeypatch.setattr(chemex_module, "read_defaults", lambda _filenames: defaults)
    monkeypatch.setattr(chemex_module, "run_methods", fake_run_methods)

    chemex_module.run(make_args("fit"), session=session)

    np.testing.assert_equal(session.model_names, ["2st"])
    np.testing.assert_equal(session.parameters.defaults_calls, [defaults])
    np.testing.assert_equal(experiments.filtered, 1)
    np.testing.assert_equal(recorded["build"][2], session)
    np.testing.assert_equal(recorded["run_methods"][4], session)


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
        plot: bool = False,
        session: StubSession | None = None,
    ) -> None:
        recorded["simulation"] = (experiments_arg, path, plot, session)

    monkeypatch.setattr(chemex_module, "build_experiments", fake_build_experiments)
    monkeypatch.setattr(chemex_module, "read_defaults", lambda _filenames: defaults)
    monkeypatch.setattr(chemex_module, "execute_simulation", fake_execute_simulation)

    chemex_module.run(make_args("simulate"), session=session)

    np.testing.assert_equal(session.model_names, ["2st"])
    np.testing.assert_equal(session.parameters.defaults_calls, [defaults])
    np.testing.assert_equal(session.parameters.fix_all_calls, 1)
    np.testing.assert_equal(recorded["build"][2], session)
    np.testing.assert_equal(recorded["simulation"][3], session)


def test_main_bootstraps_plugins_for_non_run_commands(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[str] = []

    class StubParser:
        def parse_args(self) -> Namespace:
            return Namespace(func=lambda _args: calls.append("func"))

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


def test_run_methods_passes_session_to_fit_groups(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    experiments = FakeExperiments()
    recorded: dict[str, object] = {}

    def fake_fit_groups(  # noqa: PLR0913
        experiments_arg: FakeExperiments,
        path: Path,
        plot: str,
        fitmethod: str,
        statistics: object,
        *,
        session: StubSession | None = None,
    ) -> None:
        recorded["fit_groups"] = (
            experiments_arg,
            path,
            plot,
            fitmethod,
            statistics,
            session,
        )

    monkeypatch.setattr(fitting_module, "_fit_groups", fake_fit_groups)

    fitting_module.run_methods(
        experiments,
        {"": Method()},
        Path("Output"),
        "normal",
        session=session,
    )

    np.testing.assert_equal(len(experiments.selections), 1)
    np.testing.assert_equal(len(session.parameters.status_calls), 1)
    np.testing.assert_equal(recorded["fit_groups"][5], session)


def test_run_statistics_uses_session_for_header(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    recorded: dict[str, object] = {}
    session = StubSession()
    session.parameters = StatisticsParameterStore()
    experiments = SimpleNamespace(param_ids=["__PB"])

    monkeypatch.setattr(fitting_module, "track", lambda _iterable, **_kwargs: [])

    def fake_print_header(
        grid: object,
        *,
        session: StubSession | None = None,
    ) -> str:
        recorded["grid"] = grid
        recorded["session"] = session
        return "# header\n"

    monkeypatch.setattr(fitting_module, "print_header", fake_print_header)

    fitting_module._run_statistics(  # noqa: SLF001
        experiments,
        tmp_path,
        "leastsq",
        fitting_module.Statistics(mc=1),
        session=session,
    )

    np.testing.assert_equal(recorded["grid"], ["PB"])
    np.testing.assert_equal(recorded["session"], session)


def test_execute_post_fit_writes_parameters_from_session_store(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    param = ParamSetting(ParamName("PB"), value=0.15, vary=False)
    session = StubSession()
    session.parameters = WriterParameterStore({param.id_: param})
    experiments = SimpleNamespace(param_ids=[param.id_], write=lambda _path: None)

    def should_not_be_called(_param_ids: object) -> None:
        msg = "global parameter store should not be used"
        raise AssertionError(msg)

    monkeypatch.setattr(
        parameter_printer_module.database,
        "get_parameters",
        should_not_be_called,
    )
    monkeypatch.setattr(helper_module, "print_writing_results", lambda _path: None)
    monkeypatch.setattr(
        helper_module,
        "_write_statistics",
        lambda *_args, **_kwargs: None,
    )

    helper_module.execute_post_fit(experiments, tmp_path, session=session)

    np.testing.assert_equal(int((tmp_path / "Parameters" / "fixed.toml").exists()), 1)


def test_write_file_rejects_unknown_parameter_status(tmp_path: Path) -> None:
    parameter = ParamSetting(ParamName("PB"), value=0.15, vary=False)
    parameters = parameter_printer_module.GlobalLocalParameters(
        {parameter.param_name: parameter},
        {},
    )

    with pytest.raises(ValueError, match="Unknown parameter status"):
        parameter_printer_module.write_file(parameters, "typo", tmp_path)

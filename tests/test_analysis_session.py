from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

from chemex import chemex as chemex_module
from chemex.configuration.methods import Method, Selection
from chemex.experiments import builder as builder_module
from chemex.optimize import fitting as fitting_module
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

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


class StubParameterFactory:
    def __init__(self) -> None:
        self.clear_calls = 0

    def clear_cache(self) -> None:
        self.clear_calls += 1


class StubSession:
    def __init__(self) -> None:
        self.parameters = StubParameters()
        self.model_names: list[str] = []
        self.parameter_factory = StubParameterFactory()

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
        return {"__PB": SimpleNamespace(name="__PB", vary=True)}

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

    def filter(self) -> None:
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

    monkeypatch.setattr(
        builder_module,
        "build_experiment",
        lambda filename, _selection, *, session=None: (  # noqa: ARG005
            DummyExperiment(filename)
        ),
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
    experiments = FakeExperiments(parameter_store=session.parameters)
    recorded: dict[str, object] = {}

    def fake_fit_groups(
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


def test_run_methods_skips_fit_when_selection_removes_all_profiles(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = StubSession()
    experiments = EmptyAfterSelectExperiments(parameter_store=session.parameters)
    calls: list[str] = []

    monkeypatch.setattr(fitting_module, "print_no_data", lambda: calls.append("no_data"))

    def fail_if_called(*_args, **_kwargs) -> None:
        pytest.fail("_fit_groups should not run when no profiles remain selected")

    monkeypatch.setattr(fitting_module, "_fit_groups", fail_if_called)

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


def test_run_statistics_uses_session_for_header(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    session = StubSession()
    session.parameters = StatisticsParameterStore()
    experiments = SimpleNamespace(
        param_ids=["__PB"],
        parameter_store=session.parameters,
    )

    monkeypatch.setattr(resampling_module, "track", lambda _iterable, **_kwargs: [])

    fitting_module._run_statistics(  # noqa: SLF001
        experiments,
        tmp_path,
        "leastsq",
        fitting_module.Statistics(mc=1),
    )

    assert not (tmp_path / "monte_carlo.out").exists()
    assert (
        tmp_path / "Statistics" / "MonteCarlo" / "samples.tsv"
    ).read_text(encoding="utf-8") == "[PB]\tchisqr\n"
    assert (tmp_path / "Statistics" / "MonteCarlo" / "summary.toml").exists()
    assert (tmp_path / "Statistics" / "MonteCarlo" / "correlations.tsv").exists()
    assert (tmp_path / "Statistics" / "MonteCarlo" / "plots.pdf").stat().st_size > 0
    diagnostics = (
        tmp_path / "Statistics" / "MonteCarlo" / "diagnostics.toml"
    ).read_text(encoding="utf-8")
    assert 'method = "Monte Carlo"' in diagnostics
    assert "requested_samples = 1" in diagnostics
    assert "completed_samples = 0" in diagnostics
    assert 'samples_file = "samples.tsv"' in diagnostics
    assert 'summary_file = "summary.toml"' in diagnostics
    assert 'correlations_file = "correlations.tsv"' in diagnostics
    assert 'plots_file = "plots.pdf"' in diagnostics


def test_run_statistics_dispatches_mcmc_without_resampling(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    recorded: dict[str, object] = {}
    session = StubSession()
    session.parameters = StatisticsParameterStore()
    experiments = SimpleNamespace(
        param_ids=["__PB"],
        parameter_store=session.parameters,
    )

    def fail_generate_exp_for_statistics(*_args, **_kwargs) -> None:
        pytest.fail("MCMC should not generate resampled experiments")

    def fake_run_mcmc(
        experiments_arg: object,
        params_arg: object,
        settings_arg: object,
        path_arg: Path,
    ) -> None:
        recorded["experiments"] = experiments_arg
        recorded["params"] = params_arg
        recorded["settings"] = settings_arg
        recorded["path"] = path_arg

    monkeypatch.setattr(
        resampling_module,
        "generate_exp_for_statistics",
        fail_generate_exp_for_statistics,
    )
    monkeypatch.setattr(fitting_module, "run_mcmc", fake_run_mcmc)

    fitting_module._run_statistics(  # noqa: SLF001
        experiments,
        tmp_path,
        "leastsq",
        fitting_module.Statistics(mcmc=100),
    )

    np.testing.assert_equal(recorded["experiments"], experiments)
    assert "__PB" in recorded["params"]
    assert recorded["settings"].steps == 100
    np.testing.assert_equal(recorded["path"], tmp_path)


def test_resampling_summary_and_correlations_are_written(tmp_path: Path) -> None:
    store = WriterParameterStore(
        {
            "__PB": ParamSetting(ParamName("PB"), value=0.15, vary=True),
            "__KEX_AB": ParamSetting(ParamName("KEX_AB"), value=250.0, vary=True),
        },
    )
    samples = np.array([[0.1, 200.0], [0.3, 300.0], [0.5, 400.0]])
    parameter_names = resampling_module._format_parameter_names(  # noqa: SLF001
        ["__PB", "__KEX_AB"],
        store,
    )

    resampling_module._write_resampling_summary(  # noqa: SLF001
        tmp_path,
        parameter_names=parameter_names,
        samples=samples,
    )
    resampling_module._write_resampling_correlations(  # noqa: SLF001
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


def test_execute_post_fit_writes_parameters_from_session_store(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    param = ParamSetting(ParamName("PB"), value=0.15, vary=False)
    session = StubSession()
    session.parameters = WriterParameterStore({param.id_: param})
    experiments = SimpleNamespace(
        param_ids=[param.id_],
        parameter_store=session.parameters,
        write=lambda _path: None,
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
        parameter_printer_module.write_file(
            parameters,
            "typo",
            tmp_path,
            WriterParameterStore({}),
        )

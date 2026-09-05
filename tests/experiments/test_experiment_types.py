from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

import chemex.experiments.builder as builder_module
from chemex.configuration.methods import Selection
from chemex.containers.dataset import DatasetLoadError
from chemex.containers.experiment import NoDuplicateNoiseNotice
from chemex.experiments import catalog, experiment_types
from chemex.experiments.catalog import cest_15n, wip
from chemex.experiments.experiment_types import (
    ExperimentConfigurationError,
    ExperimentDataError,
    ExperimentFileError,
    ExperimentNameError,
    ExperimentSource,
    ExperimentTomlError,
    InvalidExperimentSourceError,
    UnknownExperimentTypeError,
)
from chemex.experiments.loader import (
    import_module,
    iter_experiment_modules,
    register_experiments,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.plotters.cest import CestPlotter
from chemex.plotters.cpmg import CpmgPlotter
from chemex.plotters.exsy import EXSYPlotter
from chemex.plotters.relaxation import RelaxationPlotter
from chemex.plotters.shift import ShiftPlotter
from chemex.printers.data import (
    CestPrinter,
    CpmgPrinter,
    EXSYPrinter,
    RelaxationPrinter,
    ShiftPrinter,
)
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parents[2]
ALL_PROFILES = Selection(include=None, exclude=None)


def _write_relaxation_input(tmp_path: Path, data_filename: str) -> Path:
    filename = tmp_path / "experiment.toml"
    filename.write_text(
        f"""[experiment]
name = "relaxation_nz"

[conditions]
h_larmor_frq = 600.0

[data]
path = "."

[data.profiles]
G2N = "{data_filename}"
""",
        encoding="utf-8",
    )
    return filename


def _discovered_experiment_type_names() -> set[str]:
    module_names = [
        *iter_experiment_modules(catalog),
        *iter_experiment_modules(wip),
    ]
    return {
        import_module(module_name).EXPERIMENT_TYPE.name for module_name in module_names
    }


def test_all_discovered_experiment_types_register_through_public_surface() -> None:
    expected_names = _discovered_experiment_type_names()

    register_experiments()
    registered = experiment_types.registered_experiment_types()

    assert len(expected_names) == 37
    assert set(registered) == expected_names
    assert all(name == adapter.name for name, adapter in registered.items())


def test_every_registered_experiment_type_can_be_opened(tmp_path: Path) -> None:
    register_experiments()

    for index, name in enumerate(experiment_types.registered_experiment_types()):
        filename = tmp_path / f"experiment-{index}.toml"
        filename.write_text(
            f'[experiment]\nname = "{name}"\n',
            encoding="utf-8",
        )

        source = experiment_types.open(filename)

        assert source.filename == filename
        assert source.experiment_type_name == name


def test_registration_is_idempotent_only_for_the_same_adapter() -> None:
    register_experiments()
    adapter = cest_15n.EXPERIMENT_TYPE

    experiment_types.register_experiment_type(adapter)
    experiment_types.register_experiment_type(adapter)

    with pytest.raises(ValueError, match="different adapter"):
        experiment_types.register_experiment_type(replace(adapter))


def test_opened_source_does_not_alias_parsed_toml(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    register_experiments()
    filename = ROOT / "examples/Experiments/CEST_15N/Experiments/13hz.toml"
    raw_config = experiment_types.load_toml(filename)
    monkeypatch.setattr(experiment_types, "load_toml", lambda _filename: raw_config)
    source = experiment_types.open(filename)
    raw_config["experiment"]["name"] = "unknown"
    raw_config["data"]["profiles"].clear()
    session = AnalysisSession.create()
    session.set_model("2st")

    result = experiment_types.build(
        source,
        selection=Selection(
            include=[SpinSystem.from_name("13N")],
            exclude=None,
        ),
        model=session.model.spec,
        parameters=session.parameter_factory,
    )

    assert source.filename == filename
    assert source.experiment_type_name == "cest_15n"
    assert len(result.experiment.profiles) == 1


def test_experiment_source_cannot_be_constructed_with_arbitrary_state() -> None:
    with pytest.raises(TypeError):
        ExperimentSource(
            Path("forged.toml"),
            "cest_15n",
            {"experiment": {"name": "cest_15n"}},
        )


def test_opened_source_is_immutable(tmp_path: Path) -> None:
    register_experiments()
    filename = tmp_path / "experiment.toml"
    filename.write_text('[experiment]\nname = "cest_15n"\n', encoding="utf-8")
    source = experiment_types.open(filename)

    with pytest.raises(AttributeError, match="immutable"):
        source._filename = tmp_path / "forged.toml"

    assert source.filename == filename
    assert source.experiment_type_name == "cest_15n"


def test_forged_experiment_source_is_a_typed_build_failure() -> None:
    session = AnalysisSession.create()
    session.set_model("2st")
    forged = object.__new__(ExperimentSource)

    with pytest.raises(InvalidExperimentSourceError):
        experiment_types.build(
            forged,
            selection=ALL_PROFILES,
            model=session.model.spec,
            parameters=session.parameter_factory,
        )


def test_stale_experiment_source_is_a_typed_build_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    register_experiments()
    filename = tmp_path / "experiment.toml"
    filename.write_text('[experiment]\nname = "cest_15n"\n', encoding="utf-8")
    source = experiment_types.open(filename)
    session = AnalysisSession.create()
    session.set_model("2st")
    registry_without_source_adapter = dict(
        experiment_types.registered_experiment_types(),
    )
    registry_without_source_adapter.pop(source.experiment_type_name)
    monkeypatch.setattr(
        experiment_types,
        "_REGISTRY",
        registry_without_source_adapter,
    )

    with pytest.raises(InvalidExperimentSourceError) as error_info:
        experiment_types.build(
            source,
            selection=ALL_PROFILES,
            model=session.model.spec,
            parameters=session.parameter_factory,
        )

    assert error_info.value.filename == filename
    assert "Stale ExperimentSource" in error_info.value.explanation


def test_open_reports_expected_input_failures(tmp_path: Path) -> None:
    with pytest.raises(ExperimentFileError):
        experiment_types.open(tmp_path / "missing.toml")


def test_open_reports_invalid_utf8_as_typed_experiment_toml_error(
    tmp_path: Path,
) -> None:
    filename = tmp_path / "experiment.toml"
    filename.write_bytes(b"[experiment]\nname = \xff\n")

    with pytest.raises(ExperimentTomlError) as error_info:
        experiment_types.open(filename)

    assert isinstance(error_info.value.__cause__, UnicodeDecodeError)

    invalid = tmp_path / "invalid.toml"
    invalid.write_text("[experiment\n", encoding="utf-8")
    with pytest.raises(ExperimentTomlError):
        experiment_types.open(invalid)

    missing_name = tmp_path / "missing-name.toml"
    missing_name.write_text("[experiment]\n", encoding="utf-8")
    with pytest.raises(ExperimentNameError):
        experiment_types.open(missing_name)

    unknown = tmp_path / "unknown.toml"
    unknown.write_text('[experiment]\nname = "unknown"\n', encoding="utf-8")
    with pytest.raises(UnknownExperimentTypeError):
        experiment_types.open(unknown)


def test_build_reports_configuration_failure(tmp_path: Path) -> None:
    register_experiments()
    filename = tmp_path / "incomplete.toml"
    filename.write_text('[experiment]\nname = "cest_15n"\n', encoding="utf-8")
    source = experiment_types.open(filename)
    session = AnalysisSession.create()
    session.set_model("2st")

    with pytest.raises(ExperimentConfigurationError):
        experiment_types.build(
            source,
            selection=ALL_PROFILES,
            model=session.model.spec,
            parameters=session.parameter_factory,
        )


def test_build_reports_missing_dataset_file(tmp_path: Path) -> None:
    register_experiments()
    filename = _write_relaxation_input(tmp_path, "missing.dat")
    session = AnalysisSession.create()
    session.set_model("2st")

    with pytest.raises(ExperimentDataError):
        experiment_types.build(
            experiment_types.open(filename),
            selection=ALL_PROFILES,
            model=session.model.spec,
            parameters=session.parameter_factory,
        )


def test_build_reports_malformed_dataset_content(tmp_path: Path) -> None:
    register_experiments()
    data_file = tmp_path / "malformed.dat"
    data_file.write_text("not-a-number 1.0 0.1\n", encoding="utf-8")
    filename = _write_relaxation_input(tmp_path, data_file.name)
    session = AnalysisSession.create()
    session.set_model("2st")

    with pytest.raises(ExperimentDataError) as error_info:
        experiment_types.build(
            experiment_types.open(filename),
            selection=ALL_PROFILES,
            model=session.model.spec,
            parameters=session.parameter_factory,
        )

    assert isinstance(error_info.value.error, DatasetLoadError)
    assert error_info.value.error.filename == data_file
    assert isinstance(error_info.value.__cause__, DatasetLoadError)
    assert isinstance(error_info.value.error.__cause__, ValueError)


def test_experiment_builder_preserves_typed_dataset_failure(
    tmp_path: Path,
) -> None:
    data_file = tmp_path / "malformed.dat"
    data_file.write_text("not-a-number 1.0 0.1\n", encoding="utf-8")
    filename = _write_relaxation_input(tmp_path, data_file.name)
    session = AnalysisSession.create()
    session.set_model("2st")

    with pytest.raises(ExperimentDataError) as error_info:
        builder_module.build_experiments(
            [filename],
            ALL_PROFILES,
            session=session,
        )

    assert isinstance(error_info.value.error, DatasetLoadError)
    assert error_info.value.error.filename == data_file


@pytest.mark.parametrize(
    ("relative_path", "expected_name", "printer_type", "plotter_type"),
    [
        (
            "examples/Experiments/CPMG_15N_IP/Experiments/500mhz.toml",
            "cpmg_15n_ip",
            CpmgPrinter,
            CpmgPlotter,
        ),
        (
            "examples/Experiments/CEST_15N/Experiments/13hz.toml",
            "cest_15n",
            CestPrinter,
            CestPlotter,
        ),
        (
            "examples/Experiments/RELAXATION_NZ/Experiments/800mhz.toml",
            "relaxation_nz",
            RelaxationPrinter,
            RelaxationPlotter,
        ),
        (
            "examples/Combinations/Shifts/Experiments/sqmq_500.toml",
            "shift_15n_sqmq",
            ShiftPrinter,
            ShiftPlotter,
        ),
    ],
)
def test_representative_family_builds_complete_experiment(
    relative_path: str,
    expected_name: str,
    printer_type: type,
    plotter_type: type,
) -> None:
    register_experiments()
    session = AnalysisSession.create()
    session.set_model("2st")
    source = experiment_types.open(ROOT / relative_path)

    result = experiment_types.build(
        source,
        selection=ALL_PROFILES,
        model=session.model.spec,
        parameters=session.parameter_factory,
    )

    assert result.experiment.name == expected_name
    assert result.experiment.profiles
    assert isinstance(result.experiment.printer, printer_type)
    assert isinstance(result.experiment.plotter, plotter_type)
    assert not hasattr(result.experiment, "parameter_store")


def test_exsy_explicit_support_builds_complete_experiment(tmp_path: Path) -> None:
    register_experiments()
    data_file = tmp_path / "profile.dat"
    data_file.write_text(
        "0.0 a a 1.0 0.1\n0.1 a b 0.5 0.1\n",
        encoding="utf-8",
    )
    input_file = tmp_path / "experiment.toml"
    input_file.write_text(
        """[experiment]
name = "noesyfpgpph19"

[conditions]
h_larmor_frq = 600.0

[data]
path = "."
error = "duplicates"

[data.profiles]
G2N = "profile.dat"
""",
        encoding="utf-8",
    )
    session = AnalysisSession.create()
    session.set_model("2st")

    result = experiment_types.build(
        experiment_types.open(input_file),
        selection=ALL_PROFILES,
        model=session.model.spec,
        parameters=session.parameter_factory,
    )

    assert result.experiment.name == "noesyfpgpph19"
    assert len(result.experiment.profiles) == 1
    assert isinstance(result.experiment.printer, EXSYPrinter)
    assert isinstance(result.experiment.plotter, EXSYPlotter)
    assert isinstance(result.notices[0], NoDuplicateNoiseNotice)
    assert result.experiment.estimate_noise("duplicates") is None
    assert isinstance(result.experiment.noise_notices[0], NoDuplicateNoiseNotice)


def test_selection_precedes_profile_and_parameter_construction(
    tmp_path: Path,
) -> None:
    data_file = tmp_path / "profile.dat"
    data_file.write_text("0.0 1.0 0.1\n", encoding="utf-8")
    filename = _write_relaxation_input(tmp_path, data_file.name)
    spin_system = SpinSystem.from_name("G2N")
    included_session = AnalysisSession.create()
    included_session.set_model("2st")
    included = experiment_types.build(
        experiment_types.open(filename),
        selection=Selection(
            include=[spin_system],
            exclude=None,
        ),
        model=included_session.model.spec,
        parameters=included_session.parameter_factory,
    )
    expected_parameter_ids = included.experiment.profiles[0].param_ids
    excluded_session = AnalysisSession.create()
    excluded_session.set_model("2st")

    result = experiment_types.build(
        experiment_types.open(filename),
        selection=Selection(
            include=None,
            exclude=[spin_system],
        ),
        model=excluded_session.model.spec,
        parameters=excluded_session.parameter_factory,
    )

    assert not result.experiment.profiles
    assert excluded_session.parameters.get_parameters(expected_parameter_ids) == {}


def test_cpmg_parameter_ids_and_construction_identity_remain_stable() -> None:
    register_experiments()
    filename = ROOT / "examples/Experiments/CPMG_15N_IP/Experiments/500mhz.toml"
    selection = Selection(
        include=[SpinSystem.from_name("10N")],
        exclude=None,
    )
    session = AnalysisSession.create()
    session.set_model("2st")
    source = experiment_types.open(filename)

    first = experiment_types.build(
        source,
        selection=selection,
        model=session.model.spec,
        parameters=session.parameter_factory,
    )
    expected_ids = {
        "__CS_A_10N",
        "__CS_B_10N",
        "__KAB",
        "__KBA",
        "__PA",
        "__PB",
        "__R1_A_10N_500_0MHZ",
        "__R1_B_10N_500_0MHZ",
        "__R2_A_10N_500_0MHZ",
        "__R2_B_10N_500_0MHZ",
    }
    first_profile = first.experiment.profiles[0]
    assert first_profile.param_ids == expected_ids

    parameter_id = "__R2_A_10N_500_0MHZ"
    parameter = session.parameters.get_parameters([parameter_id])[parameter_id]
    second = experiment_types.build(
        source,
        selection=selection,
        model=session.model.spec,
        parameters=session.parameter_factory,
    )
    rebuilt_parameter = session.parameters.get_parameters([parameter_id])[parameter_id]

    assert second.experiment.profiles[0].param_ids == expected_ids
    assert rebuilt_parameter is parameter
    assert rebuilt_parameter.value == parameter.value


def test_programming_error_from_catalog_adapter_propagates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = AnalysisSession.create()
    session.set_model("2st")

    def fail_profile_calculation(
        _config: object,
        _spin_system: SpinSystem,
    ) -> None:
        message = "broken catalog adapter"
        raise RuntimeError(message)

    adapter = replace(
        cest_15n.EXPERIMENT_TYPE,
        create_profile_calculation=fail_profile_calculation,
    )
    monkeypatch.setattr(experiment_types, "_REGISTRY", {adapter.name: adapter})
    source = experiment_types.open(
        ROOT / "examples/Experiments/CEST_15N/Experiments/13hz.toml",
    )

    with pytest.raises(RuntimeError, match="broken catalog adapter"):
        experiment_types.build(
            source,
            selection=ALL_PROFILES,
            model=session.model.spec,
            parameters=session.parameter_factory,
        )

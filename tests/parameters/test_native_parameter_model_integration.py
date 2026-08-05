"""Integration coverage for the sealed native parameter model.

These tests construct real, representative shipped experiment and parameter
inputs (`examples/Experiments/CPMG_15N_IP`) through the legacy-authoritative
`experiment_types.build`/`ParameterFactory` pipeline, with a
`ParameterDefinitionsCollector` attached, and check that the resulting sealed
`ParameterDefinitions`/`ParameterConfiguration` carry the same stable
identities, values, bounds, and scopes as the legacy catalogs -- for both
ordinary and model-free analysis construction.
"""

from __future__ import annotations

from pathlib import Path

from chemex.configuration.methods import Selection
from chemex.configuration.parameters import read_defaults
from chemex.experiments import experiment_types
from chemex.experiments.loader import register_experiments
from chemex.parameters.configuration import build_parameter_configuration
from chemex.parameters.definitions import ParameterDefinitionsCollector
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parents[2]
EXPERIMENT_FILE = ROOT / "examples/Experiments/CPMG_15N_IP/Experiments/500mhz.toml"
PARAMETER_FILE = ROOT / "examples/Experiments/CPMG_15N_IP/Parameters/parameters.toml"
INCLUDED_RESIDUES = ("1N", "10N")


def _build_representative_session(model_name: str) -> tuple[AnalysisSession, object]:
    register_experiments()
    session = AnalysisSession.create()
    session.set_model(model_name)
    session.parameter_factory.definitions = ParameterDefinitionsCollector()

    selection = Selection(
        include=[SpinSystem.from_name(name) for name in INCLUDED_RESIDUES],
        exclude=None,
    )
    source = experiment_types.open(EXPERIMENT_FILE)
    result = experiment_types.build(
        source,
        selection=selection,
        model=session.model.spec,
        parameters=session.parameter_factory,
    )
    assert result.experiment.profiles

    return session, session.parameter_factory.definitions


def test_native_definitions_match_legacy_catalogs_before_defaults() -> None:
    for model_name in ("2st", "2st.mf"):
        session, collector = _build_representative_session(model_name)

        # Snapshot the legacy-authoritative catalogs before parameter-file
        # defaults are applied, so both native and legacy reflect only the
        # model/spin-system construction.
        ordinary_catalog = session.parameters.get_catalog(model_free=False)
        model_free_catalog = session.parameters.get_catalog(model_free=True)

        ordinary_defs, model_free_defs = collector.seal()

        assert ordinary_defs.ids
        assert model_free_defs.ids
        assert list(ordinary_defs.ids) == sorted(ordinary_defs.ids)
        assert list(model_free_defs.ids) == sorted(model_free_defs.ids)

        legacy_ordinary = ordinary_catalog.get_parameters(ordinary_defs.ids)
        for param_id, definition in ordinary_defs.items():
            legacy = legacy_ordinary[param_id]
            assert definition.param_name == legacy.param_name
            assert definition.value == legacy.value
            assert definition.min == legacy.min
            assert definition.max == legacy.max

        legacy_model_free = model_free_catalog.get_parameters(model_free_defs.ids)
        for param_id, definition in model_free_defs.items():
            legacy = legacy_model_free[param_id]
            assert definition.param_name == legacy.param_name
            assert definition.value == legacy.value
            assert definition.min == legacy.min
            assert definition.max == legacy.max

        # `vary`/`expr` are intentionally excluded from the definition/legacy
        # comparison above: the legacy catalog already reflects each
        # experiment's `to_be_fitted` runtime selection (`_set_to_fit`
        # forces `vary=True`/`expr=""`), while the sealed native definition
        # keeps the model/spin-system construction's own default -- a
        # representative CPMG rate demonstrates the intentional divergence.
        r2_rate_id = next(
            param_id
            for param_id, definition in ordinary_defs.items()
            if definition.param_name.name == "R2_A"
        )
        assert legacy_ordinary[r2_rate_id].vary is True
        assert ordinary_defs[r2_rate_id].vary is False

        # A representative model-free quantity is a canonical definition in
        # `.mf` mode and a auxiliary-initialization-only construction input
        # otherwise; either way it must be part of the sealed model-free
        # catalog.
        assert any(
            definition.param_name.name == "TAUC_A"
            for definition in model_free_defs.values()
        )


def test_native_configuration_matches_legacy_effective_values_ordinary_mode() -> None:
    session, collector = _build_representative_session("2st")
    ordinary_defs, model_free_defs = collector.seal()

    defaults = read_defaults([PARAMETER_FILE])
    session.parameters.set_defaults(defaults)

    ordinary_catalog = session.parameters.get_catalog(model_free=False)
    model_free_catalog = session.parameters.get_catalog(model_free=True)

    ordinary_config = build_parameter_configuration(
        ordinary_defs,
        ordinary_catalog,
        defaults,
    )
    model_free_config = build_parameter_configuration(
        model_free_defs,
        model_free_catalog,
        defaults,
    )

    for param_id in ordinary_defs:
        assert ordinary_config[param_id].value == ordinary_catalog.get_value(param_id)
    for param_id in model_free_defs:
        value = model_free_config[param_id].value
        assert value == model_free_catalog.get_value(param_id)

    # `[GLOBAL] PB = 0.07` and per-residue `[CS_A] 1N/10N` overrides from the
    # shipped parameter file must be reflected, with provenance recorded.
    pb_id = next(
        param_id
        for param_id, definition in ordinary_defs.items()
        if definition.param_name.name == "PB" and not definition.param_name.spin_system
    )
    assert ordinary_config[pb_id].value == 0.07
    assert ordinary_config[pb_id].overridden is True

    for residue in INCLUDED_RESIDUES:
        cs_id = next(
            param_id
            for param_id, definition in ordinary_defs.items()
            if definition.param_name.name == "CS_A"
            and str(definition.param_name.spin_system) == residue
        )
        assert ordinary_config[cs_id].overridden is True


def test_native_configuration_matches_legacy_effective_values_model_free_mode() -> None:
    session, collector = _build_representative_session("2st.mf")
    ordinary_defs, model_free_defs = collector.seal()

    defaults = read_defaults([PARAMETER_FILE])
    session.parameters.set_defaults(defaults)

    model_free_catalog = session.parameters.get_catalog(model_free=True)
    model_free_config = build_parameter_configuration(
        model_free_defs,
        model_free_catalog,
        defaults,
    )

    for param_id in model_free_defs:
        assert model_free_config[param_id].value == model_free_catalog.get_value(
            param_id,
        )

    tauc_id = next(
        param_id
        for param_id, definition in model_free_defs.items()
        if definition.param_name.name == "TAUC_A"
    )
    assert model_free_config[tauc_id].value == 3.0
    assert model_free_config[tauc_id].overridden is True

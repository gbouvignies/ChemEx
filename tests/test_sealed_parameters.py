"""Tests for sealed canonical parameter definitions and configuration (#582).

Primary seam: real construction path through AnalysisSession, ParameterFactory,
and shipped experiment inputs. Verified against the legacy-authoritative catalog.
"""

from __future__ import annotations

from dataclasses import FrozenInstanceError, fields
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.configuration.methods import Selection
from chemex.configuration.parameters import DefaultSetting, read_defaults
from chemex.experiments.builder import build_experiments
from chemex.nmr.basis import Basis
from chemex.parameters.name import ParamName
from chemex.parameters.sealed import (
    ConfigurationMismatchError,
    DefinitionConflictError,
    ParamConfig,
    ParamDefinition,
    SealedDefinitions,
    build_sealed_configuration,
    extract_condition_entries,
)
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

ROOT = Path(__file__).parents[1]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_cpmg_config(
    h_larmor_frq: float = 800.0,
    profiles: dict[str, str] | None = None,
) -> SimpleNamespace:
    """Minimal config object mimicking ExperimentConfiguration for CPMG 15N IP."""
    if profiles is None:
        profiles = {}
    conditions = Conditions(h_larmor_frq=h_larmor_frq, label=("2h",))
    return SimpleNamespace(
        conditions=conditions,
        to_be_fitted=SimpleNamespace(rates=["r2"], model_free=[]),
    )


def _build_ordinary_session_with_two_residues() -> tuple[
    AnalysisSession, dict[str, str], dict[str, str]
]:
    """Build a 2st ordinary session with two spin systems through the real path."""
    session = AnalysisSession.create()
    session.set_model("2st")

    config = _make_cpmg_config()
    basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
    spin_a = SpinSystem.from_name("G23N-HN")
    spin_b = SpinSystem.from_name("A45N-HN")

    ids_a = session.parameter_factory.create_parameters(
        config,
        basis=basis,
        spin_system=spin_a,
    )
    ids_b = session.parameter_factory.create_parameters(
        config,
        basis=basis,
        spin_system=spin_b,
    )

    return session, ids_a, ids_b


# ---------------------------------------------------------------------------
# Test 1: Ordinary analysis — definitions match legacy catalog
# ---------------------------------------------------------------------------


class TestOrdinaryDefinitionsMatchLegacy:
    """Verify that ordinary-mode sealed definitions match the legacy catalog."""

    def test_sealed_definitions_are_produced(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        assert session.parameter_factory.sealed_definitions is not None

    def test_definition_param_ids_match_legacy_catalog(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()

        definitions = session.parameter_factory.sealed_definitions
        legacy_params = session.parameters.database._parameters

        definition_ids = {d.param_id for d in definitions}
        legacy_ids = set(legacy_params.keys())

        assert definition_ids == legacy_ids

    def test_definition_fields_match_legacy_param_settings(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()

        definitions = session.parameter_factory.sealed_definitions
        legacy_params = session.parameters.database._parameters

        for defn in definitions:
            legacy = legacy_params[defn.param_id]
            assert defn.name == legacy.param_name.name, (
                f"{defn.param_id}: name mismatch"
            )
            assert defn.spin_system_name == str(legacy.param_name.spin_system), (
                f"{defn.param_id}: spin_system mismatch"
            )
            assert defn.default_value == legacy.value or (
                defn.default_value is None and legacy.value is None
            ), (
                f"{defn.param_id}: default_value mismatch {defn.default_value} vs {legacy.value}"
            )

    def test_definition_bounds_match_legacy(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()

        definitions = session.parameter_factory.sealed_definitions
        legacy_params = session.parameters.database._parameters

        for defn in definitions:
            legacy = legacy_params[defn.param_id]
            np.testing.assert_equal(defn.lower_bound, legacy.min)
            np.testing.assert_equal(defn.upper_bound, legacy.max)

    def test_native_type_fields_are_scoped_to_definition_and_configuration(
        self,
    ) -> None:
        assert [item.name for item in fields(ParamDefinition)] == [
            "param_id",
            "name",
            "spin_system_name",
            "condition_entries",
            "default_value",
            "lower_bound",
            "upper_bound",
        ]
        assert [item.name for item in fields(ParamConfig)] == [
            "param_id",
            "effective_value",
            "lower_bound",
            "upper_bound",
        ]


# ---------------------------------------------------------------------------
# Test 2: Model-free analysis — complete definition set
# ---------------------------------------------------------------------------


class TestModelFreeDefinitionsComplete:
    """Model-free definitions include MF quantities and rate identities."""

    def test_model_free_definitions_match_active_catalog(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st.mf")

        config = SimpleNamespace(
            conditions=Conditions(h_larmor_frq=800.0, label=("2h",)),
            to_be_fitted=SimpleNamespace(rates=[], model_free=["tauc", "s2"]),
        )
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
        spin = SpinSystem.from_name("G23N-HN")

        session.parameter_factory.create_parameters(
            config,
            basis=basis,
            spin_system=spin,
        )
        session.parameter_factory.seal_definitions()

        definitions = session.parameter_factory.sealed_definitions
        # In model-free mode, database property dispatches to _database_mf
        legacy_mf_params = session.parameters.database._parameters

        definition_ids = {d.param_id for d in definitions}
        legacy_ids = set(legacy_mf_params.keys())

        assert definition_ids == legacy_ids

    def test_model_free_definitions_include_tauc_and_s2(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st.mf")

        config = SimpleNamespace(
            conditions=Conditions(h_larmor_frq=800.0, label=("2h",)),
            to_be_fitted=SimpleNamespace(rates=[], model_free=["tauc", "s2"]),
        )
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
        spin = SpinSystem.from_name("G23N-HN")

        session.parameter_factory.create_parameters(
            config,
            basis=basis,
            spin_system=spin,
        )
        session.parameter_factory.seal_definitions()

        definitions = session.parameter_factory.sealed_definitions
        names = {d.name for d in definitions}

        assert "TAUC_A" in names
        assert "S2_A" in names

    def test_model_free_definitions_include_rate_identities(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st.mf")

        config = SimpleNamespace(
            conditions=Conditions(h_larmor_frq=800.0, label=("2h",)),
            to_be_fitted=SimpleNamespace(rates=[], model_free=["tauc", "s2"]),
        )
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
        spin = SpinSystem.from_name("G23N-HN")

        session.parameter_factory.create_parameters(
            config,
            basis=basis,
            spin_system=spin,
        )
        session.parameter_factory.seal_definitions()

        definitions = session.parameter_factory.sealed_definitions
        names = {d.name for d in definitions}

        # Rate identities should be present even in model-free mode
        r2_names = {n for n in names if n.startswith("R2")}
        assert len(r2_names) > 0, "Rate identities missing from model-free definitions"

    def test_ordinary_mode_excludes_model_free_only_quantities(self) -> None:
        """Model-free auxiliary quantities (TAUC, S2, KHH) that exist only in
        parameters_mf must not appear in ordinary canonical definitions."""
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()

        definitions = session.parameter_factory.sealed_definitions
        definition_ids = {d.param_id for d in definitions}

        # These should match the ordinary catalog, not the model-free catalog
        ordinary_ids = set(session.parameters._database._parameters.keys())
        mf_only_ids = (
            set(session.parameters._database_mf._parameters.keys()) - ordinary_ids
        )

        assert not (definition_ids & mf_only_ids), (
            f"Ordinary definitions contain model-free-only IDs: "
            f"{definition_ids & mf_only_ids}"
        )


# ---------------------------------------------------------------------------
# Test 3: Equivalent duplicate contributions are deduplicated
# ---------------------------------------------------------------------------


class TestEquivalentDuplicateDeduplication:
    """Multiple spin systems contributing the same global parameter."""

    def test_global_params_deduplicated(self) -> None:
        session, ids_a, ids_b = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()

        definitions = session.parameter_factory.sealed_definitions

        # kex_ab is global — contributed by both spin systems
        kex_ids = [d.param_id for d in definitions if d.name == "KEX_AB"]
        assert len(kex_ids) == 1, f"Expected 1 KEX_AB definition, got {len(kex_ids)}"

    def test_residue_specific_params_not_deduplicated(self) -> None:
        session, ids_a, ids_b = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()

        definitions = session.parameter_factory.sealed_definitions

        # Chemical shifts are residue-specific
        cs_defs = [d for d in definitions if d.name == "CS_A"]
        assert len(cs_defs) >= 2, f"Expected ≥2 CS_A definitions, got {len(cs_defs)}"


# ---------------------------------------------------------------------------
# Test 4: Conflicting contributions raise DefinitionConflictError
# ---------------------------------------------------------------------------


class TestConflictingContributions:
    """Conflicting definition fields for the same param_id must be rejected."""

    def test_conflicting_default_value_raises(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st")

        # Create a definition contribution manually
        defn_a = ParamDefinition(
            param_id="__KEX_AB_T_25_0C",
            name="KEX_AB",
            spin_system_name="",
            condition_entries=(),
            default_value=200.0,
            lower_bound=0.0,
            upper_bound=float("inf"),
        )
        defn_b = ParamDefinition(
            param_id="__KEX_AB_T_25_0C",
            name="KEX_AB",
            spin_system_name="",
            condition_entries=(),
            default_value=500.0,  # CONFLICT
            lower_bound=0.0,
            upper_bound=float("inf"),
        )

        session.parameter_factory._definition_contributions["__KEX_AB_T_25_0C"] = [
            defn_a,
            defn_b,
        ]

        with pytest.raises(DefinitionConflictError) as exc_info:
            session.parameter_factory.seal_definitions()

        assert exc_info.value.param_id == "__KEX_AB_T_25_0C"
        assert exc_info.value.field == "default_value"

    def test_conflicting_bounds_raises(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st")

        defn_a = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(),
            default_value=0.05,
            lower_bound=0.0,
            upper_bound=1.0,
        )
        defn_b = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(),
            default_value=0.05,
            lower_bound=-1.0,  # CONFLICT
            upper_bound=1.0,
        )

        session.parameter_factory._definition_contributions["__PB"] = [defn_a, defn_b]

        with pytest.raises(DefinitionConflictError) as exc_info:
            session.parameter_factory.seal_definitions()

        assert exc_info.value.param_id == "__PB"
        assert exc_info.value.field == "lower_bound"

    def test_reports_every_conflicting_field_and_all_contributors(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st")
        defn_a = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(),
            default_value=0.05,
            lower_bound=0.0,
            upper_bound=1.0,
        )
        defn_b = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(),
            default_value=0.10,
            lower_bound=-1.0,
            upper_bound=1.0,
        )
        session.parameter_factory._definition_contributions["__PB"] = [
            (defn_a, "experiment-a"),  # type: ignore[list-item]
            (defn_b, "experiment-b"),  # type: ignore[list-item]
        ]

        with pytest.raises(DefinitionConflictError) as exc_info:
            session.parameter_factory.seal_definitions()

        assert exc_info.value.conflicting_fields == (
            "default_value",
            "lower_bound",
        )
        assert exc_info.value.contributors == ("experiment-a", "experiment-b")
        assert "experiment-a" in str(exc_info.value)
        assert "experiment-b" in str(exc_info.value)


# ---------------------------------------------------------------------------
# Test 5: Deep public immutability
# ---------------------------------------------------------------------------


class TestDeepImmutability:
    """All sealed objects must be structurally immutable through public API."""

    def test_param_definition_fields_frozen(self) -> None:
        defn = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(("temperature", 25.0),),
            default_value=0.05,
            lower_bound=0.0,
            upper_bound=1.0,
        )
        with pytest.raises(FrozenInstanceError):
            defn.default_value = 0.1  # type: ignore[misc]

    def test_param_definition_condition_entries_immutable(self) -> None:
        defn = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(("temperature", 25.0),),
            default_value=0.05,
            lower_bound=0.0,
            upper_bound=1.0,
        )
        with pytest.raises(TypeError):
            defn.condition_entries[0] = ("h_larmor_frq", 600.0)  # type: ignore[index]

    def test_mutable_constructor_inputs_are_copied(self) -> None:
        entries = [["temperature", 25.0]]
        defn = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=entries,  # type: ignore[arg-type]
            default_value=0.05,
            lower_bound=0.0,
            upper_bound=1.0,
        )
        index = {"__PB": 0}
        definitions = SealedDefinitions(_definitions=(defn,), _index=index)

        entries[0][1] = 30.0
        index["__FAKE"] = 0

        assert defn.condition_entries == (("temperature", 25.0),)
        assert "__FAKE" not in definitions

    def test_param_config_fields_frozen(self) -> None:
        cfg = ParamConfig(
            param_id="__PB",
            effective_value=0.05,
            lower_bound=0.0,
            upper_bound=1.0,
        )
        with pytest.raises(FrozenInstanceError):
            cfg.effective_value = 0.1  # type: ignore[misc]

    def test_sealed_definitions_cannot_be_mutated(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        definitions = session.parameter_factory.sealed_definitions

        with pytest.raises((FrozenInstanceError, TypeError, AttributeError)):
            definitions._definitions = ()  # type: ignore[misc]

    def test_sealed_configuration_cannot_be_mutated(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        session.parameters.set_defaults([])
        session.parameter_factory.seal_configuration()
        config = session.parameter_factory.sealed_configuration

        with pytest.raises((FrozenInstanceError, TypeError, AttributeError)):
            config._configs = ()  # type: ignore[misc]

    def test_lookup_indexes_cannot_be_mutated(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        definitions = session.parameter_factory.sealed_definitions

        with pytest.raises(TypeError):
            definitions._index["__FAKE"] = 0  # type: ignore[index]

        session.parameters.set_defaults([])
        session.parameter_factory.seal_configuration()
        configuration = session.parameter_factory.sealed_configuration

        with pytest.raises(TypeError):
            configuration._index["__FAKE"] = 0  # type: ignore[index]


# ---------------------------------------------------------------------------
# Test 6: Canonical ordering is deterministic
# ---------------------------------------------------------------------------


class TestCanonicalOrdering:
    """Ordering must be deterministic regardless of contribution order."""

    def test_ordering_independent_of_spin_system_order(self) -> None:
        # Build with G23N first, A45N second
        session_ab, _, _ = _build_ordinary_session_with_two_residues()
        session_ab.parameter_factory.seal_definitions()
        ids_ab = [d.param_id for d in session_ab.parameter_factory.sealed_definitions]

        # Build with A45N first, G23N second (reversed order)
        session_ba = AnalysisSession.create()
        session_ba.set_model("2st")
        config = _make_cpmg_config()
        basis = Basis(type="ixy", spin_system="nh", model=session_ba.model.spec)
        spin_a = SpinSystem.from_name("A45N-HN")
        spin_b = SpinSystem.from_name("G23N-HN")

        session_ba.parameter_factory.create_parameters(
            config,
            basis=basis,
            spin_system=spin_a,
        )
        session_ba.parameter_factory.create_parameters(
            config,
            basis=basis,
            spin_system=spin_b,
        )
        session_ba.parameter_factory.seal_definitions()
        ids_ba = [d.param_id for d in session_ba.parameter_factory.sealed_definitions]

        assert ids_ab == ids_ba

    def test_residues_and_conditions_use_semantic_numeric_order(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st")
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
        for field, residue in (
            (1000.0, "G10N-HN"),
            (800.0, "A2N-HN"),
        ):
            session.parameter_factory.create_parameters(
                _make_cpmg_config(h_larmor_frq=field),
                basis=basis,
                spin_system=SpinSystem.from_name(residue),
            )
        session.parameter_factory.seal_definitions()
        definitions = session.parameter_factory.sealed_definitions

        cs_residues = [
            defn.spin_system_name for defn in definitions if defn.name == "CS_A"
        ]
        r2_fields = [
            dict(defn.condition_entries)["h_larmor_frq"]
            for defn in definitions
            if defn.name == "R2_A"
        ]

        assert cs_residues == ["A2N", "G10N"]
        assert r2_fields == [800.0, 1000.0]


class TestSealingLifecycle:
    def test_model_cannot_change_after_definition_construction(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        definitions_identity = session.parameter_factory.sealed_definitions.identity

        with pytest.raises(RuntimeError, match="reset the analysis session"):
            session.set_model("3st")

        assert session.model.name == "2st"
        assert (
            session.parameter_factory.sealed_definitions.identity
            == definitions_identity
        )

    def test_contributions_are_rejected_after_definition_sealing(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        config = _make_cpmg_config()
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)

        with pytest.raises(RuntimeError, match="definitions.*sealed"):
            session.parameter_factory.create_parameters(
                config,
                basis=basis,
                spin_system=SpinSystem.from_name("L67N-HN"),
            )

    def test_definition_sealing_cannot_be_repeated(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()

        with pytest.raises(RuntimeError, match="definitions.*sealed"):
            session.parameter_factory.seal_definitions()

    def test_configuration_requires_defaults_initialization(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()

        with pytest.raises(RuntimeError, match="defaults"):
            session.parameter_factory.seal_configuration()

    def test_configuration_sealing_cannot_be_repeated(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        session.parameters.set_defaults([])
        session.parameter_factory.seal_configuration()

        with pytest.raises(RuntimeError, match="configuration.*sealed"):
            session.parameter_factory.seal_configuration()

    def test_defaults_cannot_be_applied_repeatedly(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        session.parameters.set_defaults([])

        with pytest.raises(RuntimeError, match="defaults.*already.*applied"):
            session.parameters.set_defaults([])

    def test_defaults_cannot_change_after_configuration_sealing(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        session.parameters.set_defaults([])
        session.parameter_factory.seal_configuration()
        identity = session.parameter_factory.sealed_configuration.identity

        with pytest.raises(RuntimeError, match="configuration.*sealed"):
            session.parameters.set_defaults(
                [(ParamName.from_section("PB"), DefaultSetting(0.20))]
            )

        assert session.parameter_factory.sealed_configuration.identity == identity


# ---------------------------------------------------------------------------
# Test 7: Configuration matches post-defaults legacy state
# ---------------------------------------------------------------------------


class TestConfigurationMatchesLegacy:
    """SealedConfiguration must match the legacy catalog after set_defaults()."""

    def test_configuration_matches_legacy_catalog(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()

        # Apply defaults (empty in this case)
        session.parameters.set_defaults([])
        session.parameter_factory.seal_configuration()

        config = session.parameter_factory.sealed_configuration
        legacy_params = session.parameters.database._parameters

        for cfg in config:
            legacy = legacy_params[cfg.param_id]
            assert cfg.effective_value == legacy.value or (
                cfg.effective_value is None and legacy.value is None
            ), f"{cfg.param_id}: value mismatch"
            np.testing.assert_equal(cfg.lower_bound, legacy.min)
            np.testing.assert_equal(cfg.upper_bound, legacy.max)

    def test_configuration_param_ids_match_definitions(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        session.parameters.set_defaults([])
        session.parameter_factory.seal_configuration()

        def_ids = {d.param_id for d in session.parameter_factory.sealed_definitions}
        cfg_ids = {c.param_id for c in session.parameter_factory.sealed_configuration}

        assert def_ids == cfg_ids

    def test_missing_definition_configuration_is_rejected(self) -> None:
        definition = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(),
            default_value=0.05,
            lower_bound=0.0,
            upper_bound=1.0,
        )
        definitions = SealedDefinitions(
            _definitions=(definition,),
            _index={"__PB": 0},
        )

        with pytest.raises(ConfigurationMismatchError) as exc_info:
            build_sealed_configuration(definitions, {})

        assert exc_info.value.missing == ("__PB",)

    def test_configuration_captures_effective_overrides(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        defaults = [
            (
                ParamName.from_section("PB"),
                DefaultSetting(0.08, min=0.01, max=0.20),
            ),
        ]
        session.parameters.set_defaults(defaults)
        session.parameter_factory.seal_configuration()

        config = session.parameter_factory.sealed_configuration
        pb_id = next(defn.param_id for defn in config if defn.param_id == "__PB")
        pb = config[pb_id]

        assert pb.effective_value == 0.08
        assert pb.lower_bound == 0.01
        assert pb.upper_bound == 0.20

    def test_definitions_and_configuration_have_deterministic_identities(self) -> None:
        session_a, _, _ = _build_ordinary_session_with_two_residues()
        session_a.parameter_factory.seal_definitions()
        session_a.parameters.set_defaults([])
        session_a.parameter_factory.seal_configuration()

        session_b = AnalysisSession.create()
        session_b.set_model("2st")
        config = _make_cpmg_config()
        basis = Basis(type="ixy", spin_system="nh", model=session_b.model.spec)
        for residue in ("A45N-HN", "G23N-HN"):
            session_b.parameter_factory.create_parameters(
                config,
                basis=basis,
                spin_system=SpinSystem.from_name(residue),
            )
        session_b.parameter_factory.seal_definitions()
        session_b.parameters.set_defaults([])
        session_b.parameter_factory.seal_configuration()

        assert session_a.parameter_factory.sealed_definitions.identity == (
            session_b.parameter_factory.sealed_definitions.identity
        )
        assert session_a.parameter_factory.sealed_configuration.identity == (
            session_b.parameter_factory.sealed_configuration.identity
        )
        assert (
            session_a.parameter_factory.sealed_configuration.definitions_identity
            == session_a.parameter_factory.sealed_definitions.identity
        )


class TestShippedInputParity:
    @pytest.mark.parametrize(
        ("model", "experiment", "parameters"),
        [
            (
                "2st",
                "examples/Experiments/CPMG_15N_IP/Experiments/500mhz.toml",
                "examples/Experiments/CPMG_15N_IP/Parameters/parameters.toml",
            ),
            (
                "2st.mf",
                "examples/Experiments/CEST_15N_TR/Experiments/trosy.toml",
                "examples/Experiments/CEST_15N_TR/Parameters/parameters.toml",
            ),
        ],
    )
    def test_real_inputs_match_active_legacy_catalog(
        self,
        model: str,
        experiment: str,
        parameters: str,
    ) -> None:
        session = AnalysisSession.create()
        session.set_model(model)
        build_experiments(
            [ROOT / experiment],
            Selection(include=None, exclude=None),
            session=session,
        )
        session.parameters.set_defaults(read_defaults([ROOT / parameters]))
        session.parameter_factory.seal_configuration()

        definitions = session.parameter_factory.sealed_definitions
        configuration = session.parameter_factory.sealed_configuration
        legacy = session.parameters.database._parameters

        assert [defn.param_id for defn in definitions] == list(legacy)
        assert [item.param_id for item in configuration] == list(legacy)
        assert all(
            str(ROOT / experiment) in contribution.contributor
            for contributions in session.parameter_factory._definition_contributions.values()
            for contribution in contributions
        )
        for defn, item in zip(definitions, configuration, strict=True):
            setting = legacy[defn.param_id]
            assert item.effective_value == setting.value
            assert item.lower_bound == setting.min
            assert item.upper_bound == setting.max


# ---------------------------------------------------------------------------
# Test: extract_condition_entries helper
# ---------------------------------------------------------------------------


class TestExtractConditionEntries:
    """Verify the condition extraction helper produces sorted immutable tuples."""

    def test_extracts_set_fields(self) -> None:
        conditions = Conditions(h_larmor_frq=800.0, temperature=25.0)
        entries = extract_condition_entries(conditions)

        assert entries == (("h_larmor_frq", 800.0), ("temperature", 25.0))

    def test_omits_none_fields(self) -> None:
        conditions = Conditions(temperature=25.0)
        entries = extract_condition_entries(conditions)

        assert entries == (("temperature", 25.0),)

    def test_empty_conditions(self) -> None:
        conditions = Conditions()
        entries = extract_condition_entries(conditions)

        assert entries == ()

    def test_result_is_deeply_immutable(self) -> None:
        conditions = Conditions(h_larmor_frq=800.0, temperature=25.0)
        entries = extract_condition_entries(conditions)

        assert isinstance(entries, tuple)
        assert all(isinstance(e, tuple) for e in entries)

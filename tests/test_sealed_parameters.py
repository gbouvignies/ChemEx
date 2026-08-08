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
    InvalidConfigurationError,
    InvalidDefinitionError,
    ParamConfig,
    ParamDefinition,
    SealedDefinitions,
    build_sealed_configuration,
    extract_condition_entries,
)
from chemex.parameters.setting import ParamSetting
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
# Definition defaults and physical bounds are valid before sealing
# ---------------------------------------------------------------------------


class TestDefinitionValidation:
    @pytest.mark.parametrize(
        ("field", "default_value", "lower_bound", "upper_bound"),
        [
            ("default_value", np.nan, 0.0, 1.0),
            ("lower_bound", 0.5, np.nan, 1.0),
            ("upper_bound", 0.5, 0.0, np.nan),
            ("lower_bound", 0.5, 0.8, 0.2),
            ("default_value", -0.1, 0.0, 1.0),
            ("default_value", 1.1, 0.0, 1.0),
        ],
    )
    def test_invalid_definition_cannot_be_sealed(
        self,
        field: str,
        default_value: float,
        lower_bound: float,
        upper_bound: float,
    ) -> None:
        session = AnalysisSession.create()
        definition = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(),
            default_value=default_value,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
        )
        session.parameter_factory._definition_contributions["__PB"] = [definition]

        with pytest.raises(InvalidDefinitionError) as exc_info:
            session.parameter_factory.seal_definitions()

        assert exc_info.value.field == field
        assert session.parameter_factory.sealed_definitions is None

    @pytest.mark.parametrize(
        ("lower_bound", "upper_bound"),
        [
            (-np.inf, np.inf),
            (0.0, np.inf),
            (-np.inf, 0.0),
        ],
    )
    def test_definition_accepts_supported_infinite_bounds(
        self,
        lower_bound: float,
        upper_bound: float,
    ) -> None:
        session = AnalysisSession.create()
        definition = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(),
            default_value=0.0,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
        )
        session.parameter_factory._definition_contributions["__PB"] = [definition]

        session.parameter_factory.seal_definitions()

        assert session.parameter_factory.sealed_definitions["__PB"] == definition


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

    def test_real_construction_detects_conflict_before_legacy_collapse(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st_binding")
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
        spin_system = SpinSystem.from_name("A2N-HN")

        for contributor, p_total in (
            ("experiment-a", 1.0001e-3),
            ("experiment-b", 1.0002e-3),
        ):
            config = SimpleNamespace(
                conditions=Conditions(p_total=p_total, l_total=2.0e-3),
                to_be_fitted=SimpleNamespace(rates=[], model_free=[]),
            )
            session.parameter_factory.create_parameters(
                config,
                basis=basis,
                spin_system=spin_system,
                contributor=contributor,
            )

        assert session.parameter_factory.try_seal_definitions() is False
        error = session.parameter_factory.native_construction_error
        assert isinstance(error, DefinitionConflictError)
        assert error.conflicting_fields == ("condition_entries",)
        assert {source.split(";", maxsplit=1)[0] for source in error.contributors} == {
            "experiment-a",
            "experiment-b",
        }
        assert error.param_id in session.parameters.database._parameters

    def test_definition_collection_failure_preserves_legacy_parameters(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        session = AnalysisSession.create()
        session.set_model("2st")
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)

        def fail_native_collection(*_args: object, **_kwargs: object) -> None:
            msg = "native definition collection failed"
            raise RuntimeError(msg)

        monkeypatch.setattr(
            session.parameter_factory,
            "_collect_definitions",
            fail_native_collection,
        )

        name_map = session.parameter_factory.create_parameters(
            _make_cpmg_config(),
            basis=basis,
            spin_system=SpinSystem.from_name("A2N-HN"),
        )

        assert name_map
        legacy_parameters = session.parameters.get_parameters(name_map.values())
        assert set(name_map.values()) <= set(legacy_parameters)
        assert session.parameter_factory.try_seal_definitions() is False
        assert isinstance(
            session.parameter_factory.native_construction_error,
            RuntimeError,
        )
        with pytest.raises(RuntimeError, match="reset the analysis session"):
            session.set_model("3st")
        assert session.model.name == "2st"


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

    def test_residue_order_matches_legacy_for_multi_digit_numbers(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st")
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
        for residue in ("A10N-HN", "A2N-HN"):
            session.parameter_factory.create_parameters(
                _make_cpmg_config(),
                basis=basis,
                spin_system=SpinSystem.from_name(residue),
            )
        session.parameters.sort()
        session.parameter_factory.seal_definitions()
        definitions = session.parameter_factory.sealed_definitions

        cs_residues = [
            defn.spin_system_name for defn in definitions if defn.name == "CS_A"
        ]

        assert cs_residues == ["A2N", "A10N"]
        assert [defn.param_id for defn in definitions if defn.name == "CS_A"] == [
            param_id
            for param_id, setting in session.parameters.database._parameters.items()
            if setting.param_name.name == "CS_A"
        ]

    def test_condition_order_matches_legacy_for_different_digit_widths(self) -> None:
        session = AnalysisSession.create()
        session.set_model("2st")
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
        spin_system = SpinSystem.from_name("A2N-HN")
        for field in (800.0, 1000.0):
            session.parameter_factory.create_parameters(
                _make_cpmg_config(h_larmor_frq=field),
                basis=basis,
                spin_system=spin_system,
            )
        session.parameters.sort()
        session.parameter_factory.seal_definitions()
        definitions = session.parameter_factory.sealed_definitions

        r2_ids = [defn.param_id for defn in definitions if defn.name == "R2_A"]

        assert r2_ids == [
            "__R2_A_A2N_1000_0MHZ",
            "__R2_A_A2N_800_0MHZ",
        ]
        assert r2_ids == [
            param_id
            for param_id, setting in session.parameters.database._parameters.items()
            if setting.param_name.name == "R2_A"
        ]


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

    def test_reset_clears_and_reconstructs_sealed_parameter_state(self) -> None:
        session, _, _ = _build_ordinary_session_with_two_residues()
        session.parameter_factory.seal_definitions()
        session.parameters.set_defaults([])
        session.parameter_factory.seal_configuration()
        definitions_identity = session.parameter_factory.sealed_definitions.identity
        configuration_identity = session.parameter_factory.sealed_configuration.identity

        session.reset()

        assert session.parameter_factory.sealed_definitions is None
        assert session.parameter_factory.sealed_configuration is None
        assert session.parameter_factory.native_construction_error is None

        session.set_model("2st")
        config = _make_cpmg_config()
        basis = Basis(type="ixy", spin_system="nh", model=session.model.spec)
        for residue in ("G23N-HN", "A45N-HN"):
            session.parameter_factory.create_parameters(
                config,
                basis=basis,
                spin_system=SpinSystem.from_name(residue),
            )
        session.parameter_factory.seal_definitions()
        session.parameters.set_defaults([])
        session.parameter_factory.seal_configuration()

        assert session.parameter_factory.sealed_definitions.identity == (
            definitions_identity
        )
        assert session.parameter_factory.sealed_configuration.identity == (
            configuration_identity
        )


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

    def test_configuration_rejects_reversed_bounds(self) -> None:
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
        parameter = ParamSetting(
            ParamName("PB"),
            value=0.5,
            min=0.8,
            max=0.2,
        )

        with pytest.raises(InvalidConfigurationError, match="lower_bound"):
            build_sealed_configuration(definitions, {"__PB": parameter})

    @pytest.mark.parametrize(
        ("field", "value", "lower_bound", "upper_bound"),
        [
            ("effective_value", np.nan, 0.0, 1.0),
            ("lower_bound", 0.5, np.nan, 1.0),
            ("upper_bound", 0.5, 0.0, np.nan),
        ],
    )
    def test_configuration_rejects_nan_values_and_bounds(
        self,
        field: str,
        value: float,
        lower_bound: float,
        upper_bound: float,
    ) -> None:
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
        parameter = ParamSetting(
            ParamName("PB"),
            value=value,
            min=lower_bound,
            max=upper_bound,
        )

        with pytest.raises(InvalidConfigurationError) as exc_info:
            build_sealed_configuration(definitions, {"__PB": parameter})

        assert exc_info.value.field == field

    @pytest.mark.parametrize("value", [-0.1, 1.1])
    def test_configuration_rejects_initial_value_outside_bounds(
        self,
        value: float,
    ) -> None:
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
        parameter = ParamSetting(
            ParamName("PB"),
            value=value,
            min=0.0,
            max=1.0,
        )

        with pytest.raises(
            InvalidConfigurationError,
            match="effective_value",
        ):
            build_sealed_configuration(definitions, {"__PB": parameter})

    @pytest.mark.parametrize(
        ("lower_bound", "upper_bound"),
        [
            (-np.inf, np.inf),
            (0.0, np.inf),
            (-np.inf, 0.0),
        ],
    )
    def test_configuration_accepts_supported_infinite_bounds(
        self,
        lower_bound: float,
        upper_bound: float,
    ) -> None:
        definition = ParamDefinition(
            param_id="__PB",
            name="PB",
            spin_system_name="",
            condition_entries=(),
            default_value=0.0,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
        )
        definitions = SealedDefinitions(
            _definitions=(definition,),
            _index={"__PB": 0},
        )
        parameter = ParamSetting(
            ParamName("PB"),
            value=0.0,
            min=lower_bound,
            max=upper_bound,
        )

        configuration = build_sealed_configuration(
            definitions,
            {"__PB": parameter},
        )

        assert configuration["__PB"] == ParamConfig(
            param_id="__PB",
            effective_value=0.0,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
        )

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

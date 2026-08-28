from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from chemex.configuration.base import ExperimentConfiguration
from chemex.configuration.conditions import Conditions
from chemex.models.factory import model_factory
from chemex.nmr.basis import Basis
from chemex.parameters.database import ParameterStore
from chemex.parameters.parameterization import (
    ParameterDeclarationContribution,
    SealedParameterModel,
    seal_parameter_declarations,
)
from chemex.parameters.relaxation import RelaxationPsdBlock, SealedRelaxationDomains
from chemex.parameters.sealed import (
    DefinitionContribution,
    ParamDefinition,
    SealedConfiguration,
    SealedDefinitions,
    build_sealed_configuration,
    canonicalize_definitions,
    extract_condition_entries,
)
from chemex.parameters.setting import LocalSettings, Parameters, ParamSetting
from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spins import build_spin_param_settings

ConfigConditionsType = ExperimentConfiguration[Any, Conditions, Any]


def _set_to_fit(
    parameters: Parameters,
    name_map: dict[str, str],
    fitted: list[str],
) -> None:
    pool = set(name_map)
    for fitted_name in fitted:
        selection = {name for name in pool if name.startswith(fitted_name)}
        for name in selection:
            param_id = name_map[name]
            parameters[param_id].vary = True
            parameters[param_id].expr = ""
        pool -= selection


def _build_parameters(
    settings: LocalSettings,
    spin_system: SpinSystem,
    conditions: Conditions,
) -> tuple[dict[str, str], Parameters]:
    param_names = {
        local_name: setting.name_setting.get_param_name(spin_system, conditions)
        for local_name, setting in settings.items()
    }

    name_map = {
        local_name: param_name.id_ for local_name, param_name in param_names.items()
    }

    parameters: Parameters = {}
    for local_name, setting in settings.items():
        param_id = name_map[local_name]
        param_name = param_names[local_name]
        expression = setting.expr.format_map(name_map)
        parameters[param_id] = ParamSetting(
            param_name=param_name,
            value=setting.value,
            min=setting.min,
            max=setting.max,
            vary=setting.vary,
            expr=expression,
        )

    return name_map, parameters


@dataclass
class ParameterFactory:
    parameter_store: ParameterStore
    _settings_cache: dict[
        tuple[str, bool, bool, bool, Basis, Conditions],
        tuple[LocalSettings, LocalSettings, frozenset[str]],
    ] = field(default_factory=dict)
    _definition_contributions: dict[str, list[DefinitionContribution]] = field(
        default_factory=dict,
    )
    _declaration_contributions: dict[
        str,
        list[ParameterDeclarationContribution],
    ] = field(default_factory=dict)
    _model_free_definition_contributions: dict[
        str,
        list[DefinitionContribution],
    ] = field(default_factory=dict)
    _model_free_declaration_contributions: dict[
        str,
        list[ParameterDeclarationContribution],
    ] = field(default_factory=dict)
    _relaxation_blocks: set[RelaxationPsdBlock] = field(default_factory=set)
    _model_free_relaxation_blocks: set[RelaxationPsdBlock] = field(
        default_factory=set,
    )
    _sealed_definitions: SealedDefinitions | None = field(default=None, repr=False)
    _sealed_configuration: SealedConfiguration | None = field(default=None, repr=False)
    _sealed_parameter_model: SealedParameterModel | None = field(
        default=None,
        repr=False,
    )
    _native_construction_error: Exception | None = field(
        default=None,
        init=False,
        repr=False,
    )

    @property
    def sealed_definitions(self) -> SealedDefinitions | None:
        """Return the immutable definitions once construction has sealed them."""
        return self._sealed_definitions

    @property
    def sealed_configuration(self) -> SealedConfiguration | None:
        """Return immutable per-analysis configuration once it has been sealed."""
        return self._sealed_configuration

    @property
    def sealed_parameter_model(self) -> SealedParameterModel | None:
        """Return the complete immutable native parameter model once sealed."""
        return self._sealed_parameter_model

    @property
    def native_construction_error(self) -> Exception | None:
        """Return the failure that prevented native parameter construction."""
        return self._native_construction_error

    def disable_native_candidate(self, error: Exception) -> None:
        """Record the first native construction failure for fail-closed reporting."""
        if self._native_construction_error is None:
            self._native_construction_error = error

    def _build_settings(
        self,
        basis: Basis,
        conditions: Conditions,
    ) -> tuple[LocalSettings, LocalSettings, frozenset[str]]:
        key = (
            basis.model.name,
            basis.model.model_free,
            basis.model.temp_coef,
            basis.model.residue_specific,
            basis,
            conditions,
        )
        if key not in self._settings_cache:
            settings_spins, settings_spins_mf = build_spin_param_settings(
                basis,
                conditions,
            )
            settings_kinetics = model_factory.create_for_model(basis.model, conditions)
            settings = settings_kinetics | settings_spins
            settings_mf = settings_kinetics | settings_spins_mf
            self._settings_cache[key] = (
                settings,
                settings_mf,
                frozenset(settings_kinetics),
            )
        return self._settings_cache[key]

    def _collect_definitions(
        self,
        parameters: Parameters,
        *,
        contributor: str,
        contributions: dict[str, list[DefinitionContribution]] | None = None,
    ) -> None:
        """Extract ParamDefinition contributions from parameters before legacy merge."""
        target = (
            self._definition_contributions if contributions is None else contributions
        )
        for param_id, param_setting in parameters.items():
            contribution = DefinitionContribution(
                definition=ParamDefinition(
                    param_id=param_id,
                    name=param_setting.param_name.name,
                    spin_system_name=str(param_setting.param_name.spin_system),
                    condition_entries=extract_condition_entries(
                        param_setting.param_name.conditions,
                    ),
                    default_value=param_setting.value,
                    lower_bound=param_setting.min,
                    upper_bound=param_setting.max,
                ),
                contributor=contributor,
            )
            target.setdefault(param_id, []).append(contribution)

    @staticmethod
    def _relaxation_blocks_for_basis(
        basis: Basis,
        name_map: dict[str, str],
    ) -> tuple[RelaxationPsdBlock, ...]:
        return tuple(
            RelaxationPsdBlock(
                f"{basis.spin_system}:{specification.family}:{specification.state}",
                specification.state,
                tuple(name_map[item] for item in specification.diagonal_names),
                tuple(
                    (row, column, name_map[item])
                    for row, column, item in specification.off_diagonal_names
                ),
            )
            for specification in basis.relaxation_psd_specs
        )

    def _collect_declarations(
        self,
        parameters: Parameters,
        *,
        contributor: str,
        model_owned_ids: frozenset[str],
        supports_estimation_ids: set[str],
        default_fit_ids: set[str],
        contributions: (
            dict[str, list[ParameterDeclarationContribution]] | None
        ) = None,
    ) -> None:
        target = (
            self._declaration_contributions if contributions is None else contributions
        )
        for param_id, parameter in parameters.items():
            contribution = ParameterDeclarationContribution(
                param_id=param_id,
                supports_estimation=param_id in supports_estimation_ids,
                model_expression=parameter.expr,
                contributor=contributor,
                model_owned=param_id in model_owned_ids and bool(parameter.expr),
                requires_independent=not bool(parameter.expr),
                fits_by_default=param_id in default_fit_ids,
            )
            target.setdefault(param_id, []).append(contribution)

    def create_parameters(
        self,
        config: ConfigConditionsType,
        *,
        basis: Basis,
        spin_system: SpinSystem,
        contributor: str | None = None,
    ) -> dict[str, str]:
        if self._sealed_definitions is not None:
            msg = "Parameter definitions are already sealed; contributions are closed"
            raise RuntimeError(msg)

        settings, settings_mf, kinetic_names = self._build_settings(
            basis,
            config.conditions,
        )

        name_map, parameters = _build_parameters(
            settings,
            spin_system,
            config.conditions,
        )

        name_map_mf, parameters_mf = _build_parameters(
            settings_mf,
            spin_system,
            config.conditions,
        )

        conditions = extract_condition_entries(config.conditions)
        construction_context = (
            f"model={self.parameter_store.model.name!r}, basis={basis.type!r}, "
            f"profile={str(spin_system)!r}, conditions={conditions!r}"
        )
        if contributor is not None:
            construction_context = f"{contributor}; {construction_context}"
        _set_to_fit(parameters, name_map, config.to_be_fitted.rates)
        _set_to_fit(parameters_mf, name_map_mf, config.to_be_fitted.model_free)
        default_fit_ids = {
            param_id for param_id, parameter in parameters.items() if parameter.vary
        }
        model_free_default_fit_ids = {
            param_id for param_id, parameter in parameters_mf.items() if parameter.vary
        }
        supports_estimation_ids = default_fit_ids | {
            name_map[name]
            for name, setting in settings.items()
            if setting.supports_estimation
        }
        model_free_supports_estimation_ids = model_free_default_fit_ids | {
            name_map_mf[name]
            for name, setting in settings_mf.items()
            if setting.supports_estimation
        }

        if self._native_construction_error is None:
            native_parameters = (
                parameters_mf if self.parameter_store.model.model_free else parameters
            )
            native_name_map = (
                name_map_mf if self.parameter_store.model.model_free else name_map
            )
            model_owned_ids = frozenset(
                native_name_map[name]
                for name in kinetic_names
                if name in native_name_map
            )
            try:
                self._collect_definitions(
                    native_parameters,
                    contributor=construction_context,
                )
                self._collect_declarations(
                    native_parameters,
                    contributor=construction_context,
                    model_owned_ids=model_owned_ids,
                    supports_estimation_ids=(
                        model_free_supports_estimation_ids
                        if self.parameter_store.model.model_free
                        else supports_estimation_ids
                    ),
                    default_fit_ids=(
                        model_free_default_fit_ids
                        if self.parameter_store.model.model_free
                        else default_fit_ids
                    ),
                )
                self._relaxation_blocks.update(
                    self._relaxation_blocks_for_basis(basis, native_name_map)
                )
                if not self.parameter_store.model.model_free:
                    model_free_owned_ids = frozenset(
                        name_map_mf[name]
                        for name in kinetic_names
                        if name in name_map_mf
                    )
                    self._collect_definitions(
                        parameters_mf,
                        contributor=construction_context,
                        contributions=self._model_free_definition_contributions,
                    )
                    self._collect_declarations(
                        parameters_mf,
                        contributor=construction_context,
                        model_owned_ids=model_free_owned_ids,
                        supports_estimation_ids=model_free_supports_estimation_ids,
                        default_fit_ids=model_free_default_fit_ids,
                        contributions=self._model_free_declaration_contributions,
                    )
                    self._model_free_relaxation_blocks.update(
                        self._relaxation_blocks_for_basis(basis, name_map_mf)
                    )
            except Exception as error:  # noqa: BLE001 - checkpoint-1 isolation boundary
                self._native_construction_error = error

        self.parameter_store.add_multiple(parameters)
        self.parameter_store.add_multiple_mf(parameters_mf)

        selection = set(name_map) & basis.required_names

        return {local_name: name_map[local_name] for local_name in selection}

    def seal_definitions(self) -> None:
        """Canonicalize accumulated definitions into an immutable sealed collection.

        Must be called after all experiments have contributed their parameters.
        Equivalent duplicate contributions are collapsed; conflicting fields
        raise DefinitionConflictError.

        """
        if self._sealed_definitions is not None:
            msg = "Parameter definitions are already sealed"
            raise RuntimeError(msg)
        self._sealed_definitions = canonicalize_definitions(
            self._definition_contributions,
        )

    def seal_configuration(self) -> None:
        """Build sealed configuration from post-defaults legacy catalog state.

        Must be called after set_defaults() has completed. Requires
        seal_definitions() to have been called first.

        """
        if self._sealed_definitions is None:
            msg = "seal_definitions() must be called before seal_configuration()"
            raise RuntimeError(msg)
        if self._sealed_configuration is not None:
            msg = "Parameter configuration is already sealed"
            raise RuntimeError(msg)
        if not self.parameter_store.defaults_applied:
            msg = "Parameter defaults must be applied before configuration sealing"
            raise RuntimeError(msg)

        configuration = build_sealed_configuration(
            self._sealed_definitions,
            self.parameter_store.database._parameters,
        )
        declarations = seal_parameter_declarations(
            self._sealed_definitions,
            self._declaration_contributions,
        )
        self.parameter_store.lock_configuration()
        self._sealed_configuration = configuration
        self._sealed_parameter_model = SealedParameterModel(
            model_name=self.parameter_store.model.name,
            model_identity=self.parameter_store.model.identity,
            definitions=self._sealed_definitions,
            configuration=configuration,
            declarations=declarations,
            relaxation_domains=SealedRelaxationDomains(tuple(self._relaxation_blocks)),
        )

    def build_model_free_parameter_model(self) -> SealedParameterModel | None:
        """Build the native model-free bootstrap used by ordinary models."""
        if self.parameter_store.model.model_free:
            return None
        definitions = canonicalize_definitions(
            self._model_free_definition_contributions,
        )
        configuration = build_sealed_configuration(
            definitions,
            self.parameter_store._database_mf._parameters,
        )
        declarations = seal_parameter_declarations(
            definitions,
            self._model_free_declaration_contributions,
        )
        return SealedParameterModel(
            model_name=self.parameter_store.model.name,
            model_identity=self.parameter_store.model.identity,
            definitions=definitions,
            configuration=configuration,
            declarations=declarations,
            relaxation_domains=SealedRelaxationDomains(
                tuple(self._model_free_relaxation_blocks)
            ),
        )

    def try_seal_definitions(self) -> bool:
        """Attempt to seal definitions and retain the first construction failure."""
        if self._native_construction_error is not None:
            return False
        try:
            self.seal_definitions()
        except Exception as error:  # noqa: BLE001 - checkpoint-1 isolation boundary
            self._native_construction_error = error
            return False
        return True

    def try_seal_configuration(self) -> bool:
        """Attempt to seal configuration and retain the first construction failure."""
        if self._native_construction_error is not None:
            return False
        try:
            self.seal_configuration()
        except Exception as error:  # noqa: BLE001 - checkpoint-1 isolation boundary
            self._native_construction_error = error
            return False
        return True

    def clear_cache(self) -> None:
        if (
            self._definition_contributions
            or self._sealed_definitions is not None
            or self._native_construction_error is not None
        ):
            msg = (
                "Parameter construction has started; reset the analysis session "
                "before changing its model"
            )
            raise RuntimeError(msg)
        self._clear_state()

    def reset(self) -> None:
        """Clear all construction state as part of a full analysis reset."""
        self._clear_state()

    def _clear_state(self) -> None:
        self._settings_cache.clear()
        self._definition_contributions.clear()
        self._declaration_contributions.clear()
        self._model_free_definition_contributions.clear()
        self._model_free_declaration_contributions.clear()
        self._relaxation_blocks.clear()
        self._model_free_relaxation_blocks.clear()
        self._sealed_definitions = None
        self._sealed_configuration = None
        self._sealed_parameter_model = None
        self._native_construction_error = None


def create_parameters(
    config: ConfigConditionsType,
    *,
    basis: Basis,
    spin_system: SpinSystem,
    parameter_factory: ParameterFactory,
    contributor: str | None = None,
) -> dict[str, str]:
    return parameter_factory.create_parameters(
        config,
        basis=basis,
        spin_system=spin_system,
        contributor=contributor,
    )

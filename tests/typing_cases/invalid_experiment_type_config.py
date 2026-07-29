from typing import Self

from chemex.configuration.base import ToBeFitted
from chemex.configuration.conditions import Conditions
from chemex.experiments.catalog.cpmg_15n_ip import Cpmg15NIpSettings
from chemex.experiments.experiment_types import ProfileCalculation, cpmg_type
from chemex.models.model import ModelSpec
from chemex.parameters.spin_system import SpinSystem


class IncompatibleConfig:
    """Satisfies the other build fields but omits generic data settings."""

    model: ModelSpec
    conditions: Conditions
    experiment: Cpmg15NIpSettings

    @classmethod
    def model_validate(
        cls,
        obj: object,
        *,
        context: object | None = None,
    ) -> Self:
        raise NotImplementedError(obj, context)

    @property
    def to_be_fitted(self) -> ToBeFitted:
        return ToBeFitted()


def create_profile_calculation(
    config: IncompatibleConfig,
    spin_system: SpinSystem,
) -> ProfileCalculation:
    raise NotImplementedError(config, spin_system)


INVALID_TYPE = cpmg_type(
    name="typing.invalid_config",
    config_type=IncompatibleConfig,
)(create_profile_calculation)

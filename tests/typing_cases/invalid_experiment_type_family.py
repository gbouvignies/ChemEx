from chemex.experiments.catalog.cpmg_15n_ip import Cpmg15NIpConfig
from chemex.experiments.catalog.shift_15n_sq import (
    create_profile_calculation,
)
from chemex.experiments.experiment_types import cpmg_type

INVALID_TYPE = cpmg_type(
    name="typing.invalid",
    config_type=Cpmg15NIpConfig,
)(create_profile_calculation)

from chemex.experiments.catalog.cest_15n import (
    Cest15NConfig,
)
from chemex.experiments.catalog.cest_15n import (
    create_profile_calculation as create_cest_profile,
)
from chemex.experiments.catalog.cpmg_15n_ip import (
    Cpmg15NIpConfig,
)
from chemex.experiments.catalog.cpmg_15n_ip import (
    create_profile_calculation as create_cpmg_profile,
)
from chemex.experiments.catalog.relaxation_nz import (
    RelaxationNzConfig,
)
from chemex.experiments.catalog.relaxation_nz import (
    create_profile_calculation as create_relaxation_profile,
)
from chemex.experiments.catalog.shift_15n_sq import (
    Shift15NSqConfig,
)
from chemex.experiments.catalog.shift_15n_sq import (
    create_profile_calculation as create_shift_profile,
)
from chemex.experiments.experiment_types import (
    cest_type,
    cpmg_type,
    relaxation_type,
    shift_type,
)

CPMG_TYPE = cpmg_type(
    name="typing.cpmg",
    config_type=Cpmg15NIpConfig,
)(create_cpmg_profile)

CEST_TYPE = cest_type(
    name="typing.cest",
    config_type=Cest15NConfig,
)(create_cest_profile)

RELAXATION_TYPE = relaxation_type(
    name="typing.relaxation",
    config_type=RelaxationNzConfig,
)(create_relaxation_profile)

SHIFT_TYPE = shift_type(
    name="typing.shift",
    config_type=Shift15NSqConfig,
)(create_shift_profile)

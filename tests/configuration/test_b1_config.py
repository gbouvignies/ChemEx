"""Test B1InhomogeneityMixin configuration parsing."""

from typing import ClassVar

import numpy as np
import pytest
from pydantic import ValidationError

from chemex.configuration.experiment import B1InhomogeneityMixin
from chemex.experiments.catalog.dcest_13c import DCest13CSettings
from chemex.experiments.catalog.dcest_15n import DCest15NSettings
from chemex.nmr.distributions.dephasing import DephasingConfig
from chemex.nmr.distributions.gaussian import GaussianDistributionConfig


class DummyB1Settings(B1InhomogeneityMixin):
    pass


class DummyDefaultGaussianB1Settings(B1InhomogeneityMixin):
    legacy_b1_inh_scale_default: ClassVar[float | None] = 0.2
    legacy_b1_inh_res_default: ClassVar[int | None] = 7


class DummyDefaultDephasingB1Settings(B1InhomogeneityMixin):
    legacy_b1_inh_scale_default: ClassVar[float | None] = np.inf
    legacy_b1_inh_res_default: ClassVar[int | None] = 11


def test_b1_inhomogeneity_mixin_parses_distribution_to_typed_config():
    settings = DummyB1Settings.model_validate(
        {
            "b1_frq": 15.0,
            "b1_distribution": {
                "type": "gaussian",
                "scale": 0.1,
                "res": 5,
            },
        }
    )

    assert settings.b1_distribution is not None
    assert settings.b1_distribution.type == "gaussian"

    distribution = settings.get_b1_distribution()
    assert len(distribution.values) == 5
    assert np.isclose(distribution.weights.sum(), 1.0)


def test_b1_inhomogeneity_mixin_defaults_distribution_type_to_gaussian():
    settings = DummyB1Settings.model_validate(
        {
            "b1_frq": 15.0,
            "b1_distribution": {
                "scale": 0.1,
                "res": 5,
            },
        }
    )

    assert settings.b1_distribution is not None
    assert settings.b1_distribution.type == "gaussian"


def test_b1_inhomogeneity_mixin_without_distribution_uses_single_point():
    settings = DummyB1Settings.model_validate({"b1_frq": 15.0})

    distribution = settings.get_b1_distribution()

    np.testing.assert_allclose(distribution.values, [15.0])
    np.testing.assert_allclose(distribution.weights, [1.0])


def test_b1_inhomogeneity_mixin_rejects_unknown_distribution_type():
    with pytest.raises(ValidationError, match="Unknown B1 distribution type"):
        DummyB1Settings.model_validate(
            {
                "b1_frq": 15.0,
                "b1_distribution": {"type": "nonexistent"},
            }
        )


def test_b1_inhomogeneity_mixin_normalizes_legacy_gaussian_fields():
    settings = DummyB1Settings.model_validate(
        {
            "b1_frq": 15.0,
            "b1_inh_scale": 0.2,
            "b1_inh_res": 7,
        }
    )

    assert isinstance(settings.b1_distribution, GaussianDistributionConfig)
    assert settings.b1_distribution.scale == 0.2
    assert settings.b1_distribution.res == 7
    assert settings.b1_inh_scale == 0.2
    assert settings.b1_inh_res == 7


def test_b1_inhomogeneity_mixin_normalizes_legacy_dephasing_scale():
    settings = DummyB1Settings.model_validate(
        {
            "b1_frq": 15.0,
            "b1_inh_scale": np.inf,
        }
    )

    assert isinstance(settings.b1_distribution, DephasingConfig)


def test_b1_inhomogeneity_mixin_normalizes_defaulted_legacy_gaussian_fields():
    settings = DummyDefaultGaussianB1Settings.model_validate({"b1_frq": 15.0})

    assert isinstance(settings.b1_distribution, GaussianDistributionConfig)
    assert settings.b1_distribution.scale == 0.2
    assert settings.b1_distribution.res == 7
    assert settings.b1_inh_scale == 0.2
    assert settings.b1_inh_res == 7


def test_b1_inhomogeneity_mixin_normalizes_defaulted_legacy_fields_with_init():
    settings = DummyDefaultGaussianB1Settings(b1_frq=15.0)

    assert isinstance(settings.b1_distribution, GaussianDistributionConfig)
    assert settings.b1_distribution.scale == 0.2
    assert settings.b1_distribution.res == 7


def test_b1_inhomogeneity_mixin_normalizes_defaulted_legacy_dephasing_fields():
    settings = DummyDefaultDephasingB1Settings.model_validate({"b1_frq": 15.0})

    assert isinstance(settings.b1_distribution, DephasingConfig)
    assert settings.b1_inh_scale == np.inf
    assert settings.b1_inh_res == 11


def test_b1_inhomogeneity_mixin_rejects_legacy_and_typed_distribution_together():
    with pytest.raises(
        ValidationError,
        match="Use either 'b1_distribution' or the legacy",
    ):
        DummyB1Settings.model_validate(
            {
                "b1_frq": 15.0,
                "b1_distribution": {"type": "gaussian", "scale": 0.1, "res": 5},
                "b1_inh_scale": 0.2,
            }
        )


def test_b1_inhomogeneity_mixin_rejects_unknown_fields():
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        DummyB1Settings.model_validate({"b1_frq": 15.0, "b1_frqq": 12.0})


@pytest.mark.parametrize(
    ("settings_cls", "experiment_name"),
    [
        (DCest13CSettings, "dcest_13c"),
        (DCest15NSettings, "dcest_15n"),
    ],
)
def test_dcest_legacy_b1_frq_alias_maps_only_to_b1_eff(settings_cls, experiment_name):
    settings = settings_cls.model_validate(
        {
            "name": experiment_name,
            "time_t1": 0.1,
            "sw": 100.0,
            "carrier": 120.0,
            "pw90": 25.0e-6,
            "b1_frq": 2000.0,
        }
    )

    assert settings.b1_eff == 2000.0
    assert settings.b1_frq is None
    settings.model_name = "2st"
    assert settings.b1_eff == 2000.0
    assert settings.b1_frq is None


@pytest.mark.parametrize(
    ("settings_cls", "experiment_name"),
    [
        (DCest13CSettings, "dcest_13c"),
        (DCest15NSettings, "dcest_15n"),
    ],
)
def test_dcest_rejects_conflicting_b1_eff_and_legacy_b1_frq(
    settings_cls, experiment_name
):
    with pytest.raises(
        ValidationError,
        match="Use either 'b1_eff' or its legacy alias 'b1_frq'",
    ):
        settings_cls.model_validate(
            {
                "name": experiment_name,
                "time_t1": 0.1,
                "sw": 100.0,
                "carrier": 120.0,
                "pw90": 25.0e-6,
                "b1_eff": 1500.0,
                "b1_frq": 2000.0,
            }
        )

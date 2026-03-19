from __future__ import annotations

import numpy as np

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.constants import Distribution
from chemex.nmr.distributions.custom import CustomDistributionConfig
from chemex.nmr.distributions.gaussian import GaussianDistributionConfig
from chemex.nmr.distributions.gaussian import generate as generate_gaussian
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.parameters.spin_system import SpinSystem


def make_liouvillian() -> LiouvillianIS:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    return LiouvillianIS(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )


def test_b1_i_updates_preserve_active_distribution_model() -> None:
    liouvillian = make_liouvillian()
    distribution = CustomDistributionConfig(
        scales=[0.8, 1.0, 1.2],
        weights=[0.2, 0.6, 0.2],
    )

    liouvillian.set_b1_i_inhomogeneity(10.0, distribution)
    liouvillian.b1_i = 20.0

    np.testing.assert_allclose(liouvillian.b1_i_dist.values, [16.0, 20.0, 24.0])
    np.testing.assert_allclose(liouvillian.b1_i_dist.weights, [0.2, 0.6, 0.2])


def test_set_b1_i_distribution_rescales_sampled_profile() -> None:
    liouvillian = make_liouvillian()
    distribution = Distribution(
        np.array([8.0, 10.0, 12.0]),
        np.array([0.2, 0.6, 0.2]),
    )

    liouvillian.set_b1_i_distribution(distribution, nominal=10.0)
    liouvillian.b1_i = 20.0

    np.testing.assert_allclose(liouvillian.b1_i_dist.values, [16.0, 20.0, 24.0])
    np.testing.assert_allclose(liouvillian.b1_i_dist.weights, [0.2, 0.6, 0.2])


def test_set_b1_i_inhomogeneity_switches_back_to_gaussian_model() -> None:
    liouvillian = make_liouvillian()
    liouvillian.set_b1_i_inhomogeneity(
        10.0,
        CustomDistributionConfig(
            scales=[0.85, 1.0, 1.15],
            weights=[0.2, 0.5, 0.3],
        ),
    )

    liouvillian.set_b1_i_inhomogeneity(
        10.0,
        GaussianDistributionConfig(scale=0.1, res=5),
    )

    expected = generate_gaussian(10.0, scale=0.1, res=5)
    np.testing.assert_allclose(liouvillian.b1_i_dist.values, expected.values)
    np.testing.assert_allclose(liouvillian.b1_i_dist.weights, expected.weights)


def test_set_b1_i_inhomogeneity_replaces_custom_distribution_with_gaussian() -> None:
    liouvillian = make_liouvillian()
    liouvillian.set_b1_i_inhomogeneity(
        10.0,
        CustomDistributionConfig(
            scales=[0.8, 1.0, 1.2],
            weights=[0.2, 0.6, 0.2],
        ),
    )
    liouvillian.set_b1_i_inhomogeneity(
        10.0,
        GaussianDistributionConfig(scale=0.15, res=7),
    )

    expected = generate_gaussian(10.0, scale=0.15, res=7)
    np.testing.assert_allclose(liouvillian.b1_i_dist.values, expected.values)
    np.testing.assert_allclose(liouvillian.b1_i_dist.weights, expected.weights)

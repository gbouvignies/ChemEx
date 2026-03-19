from __future__ import annotations

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.constants import Distribution
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.nmr.spectrometer import Spectrometer
from chemex.parameters.spin_system import SpinSystem


def make_spectrometer() -> Spectrometer:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    liouvillian = LiouvillianIS(
        SpinSystem(name="G23N-HN"),
        basis,
        Conditions(h_larmor_frq=600.0),
    )
    liouvillian.update({"r2_i_a": 4.0, "cs_i_a": 1.0})
    return Spectrometer(liouvillian)


def test_calculate_shifts_returns_liouvillian_eigenvalue_imaginary_parts() -> None:
    spectrometer = make_spectrometer()
    expected = np.linalg.eigvals(
        np.squeeze(spectrometer.liouvillian.l_free).astype(np.complex128),
    ).imag

    shifts = spectrometer.calculate_shifts()

    np.testing.assert_allclose(np.sort(shifts), np.sort(expected))


def test_calculate_shifts_rejects_multi_point_distributions() -> None:
    spectrometer = make_spectrometer()
    spectrometer.jeff_i = Distribution(
        np.array([-15.0, 0.0, 15.0]),
        np.array([0.2, 0.6, 0.2]),
    )

    with pytest.raises(
        ValueError,
        match="requires single-point B1 and Jeff distributions",
    ):
        spectrometer.calculate_shifts()

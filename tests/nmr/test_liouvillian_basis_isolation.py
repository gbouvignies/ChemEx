from __future__ import annotations

import numpy as np
import pytest

from chemex.configuration.conditions import Conditions
from chemex.models.model import ModelSpec
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS
from chemex.parameters.spin_system import SpinSystem


def test_basis_reference_matrices_are_read_only() -> None:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())

    with pytest.raises(ValueError, match="read-only"):
        basis.matrices["cs_i_a"][0, 1] = 1.0


def test_shared_basis_does_not_leak_scaling_between_liouvillians() -> None:
    basis = Basis(type="ixyz", spin_system="nh", model=ModelSpec())
    spin_system = SpinSystem(name="G23N-HN")
    conditions_600 = Conditions(h_larmor_frq=600.0)
    conditions_800 = Conditions(h_larmor_frq=800.0)
    par_values = {"cs_i_a": 1.0}

    liouvillian_600 = LiouvillianIS(spin_system, basis, conditions_600)
    liouvillian_600.update(par_values)
    reference_liouvillian = liouvillian_600.l_free.copy()
    reference_matrix = basis.matrices["cs_i_a"].copy()

    liouvillian_800 = LiouvillianIS(spin_system, basis, conditions_800)
    liouvillian_800.update(par_values)
    liouvillian_600.update(par_values)

    np.testing.assert_allclose(basis.matrices["cs_i_a"], reference_matrix)
    np.testing.assert_allclose(liouvillian_600.l_free, reference_liouvillian)
    assert not np.allclose(liouvillian_600.l_free, liouvillian_800.l_free)

from __future__ import annotations

import numpy as np
import pytest

from chemex.experiments.catalog import cest_13c, dcest_15n
from chemex.nmr.constants import J_EFF, get_multiplet
from chemex.parameters.spin_system import SpinSystem


def test_get_multiplet_normalizes_weights_for_singlet_doublet_and_triplet() -> None:
    singlet = get_multiplet("m", "ce")
    np.testing.assert_allclose(singlet.values, [0.0])
    np.testing.assert_allclose(singlet.weights, [1.0])

    doublet = get_multiplet("a", "cb")
    np.testing.assert_allclose(doublet.values, [-17.5, 17.5])
    np.testing.assert_allclose(doublet.weights, [0.5, 0.5])

    triplet = get_multiplet("f", "cd1")
    order = np.argsort(triplet.values)
    np.testing.assert_allclose(triplet.values[order], [-57.0, 0.0, 57.0])
    np.testing.assert_allclose(triplet.weights[order], [0.25, 0.5, 0.25])


def test_all_j_eff_multiplets_have_unit_total_weight_and_zero_mean() -> None:
    for symbol, nuclei in J_EFF.items():
        for nucleus in nuclei:
            multiplet = get_multiplet(symbol, nucleus)
            assert multiplet.weights.sum() == pytest.approx(1.0)
            assert float(np.sum(multiplet.values * multiplet.weights)) == pytest.approx(
                0.0,
                abs=1e-12,
            )


@pytest.mark.parametrize(
    ("symbol", "expected_nuclei"),
    (
        (
            "f",
            {
                "cg",
                "cd1",
                "cd2",
                "cd",
                "cqd",
                "ce1",
                "ce2",
                "ce",
                "cqe",
                "cz",
            },
        ),
        (
            "y",
            {
                "cg",
                "cd1",
                "cd2",
                "cd",
                "cqd",
                "ce1",
                "ce2",
                "ce",
                "cqe",
                "cz",
            },
        ),
        (
            "w",
            {
                "cg",
                "cd1",
                "cd2",
                "ce2",
                "ce3",
                "cz2",
                "cz3",
                "ch2",
                "ne1",
            },
        ),
        ("u", {"n", "c", "ca", "cb", "c1'", "c6", "n3"}),
    ),
)
def test_j_eff_contains_aromatic_and_selenocysteine_entries(
    symbol: str,
    expected_nuclei: set[str],
) -> None:
    assert expected_nuclei <= set(J_EFF[symbol])


def test_aromatic_j_eff_couplings_satisfy_bond_reciprocity() -> None:
    # CB-CG aliphatic-aromatic bond (52.5 Hz)
    for symbol in ("f", "y", "w"):
        assert 52.5 in J_EFF[symbol]["cb"]
        assert 52.5 in J_EFF[symbol]["cg"]

    # Tyr phenolic C-OH increases CE-CZ coupling (66.0 Hz)
    assert J_EFF["y"]["ce1"] == (57.0, 66.0)
    assert J_EFF["y"]["ce2"] == (57.0, 66.0)
    assert J_EFF["y"]["cz"] == (66.0, 66.0)

    # Trp indole ring couplings match NE1 (-14.4, -14.4)
    assert J_EFF["w"]["ne1"] == (-14.4, -14.4)
    assert -14.4 in J_EFF["w"]["cd1"]
    assert -14.4 in J_EFF["w"]["ce2"]
    assert 66.0 in J_EFF["w"]["cg"]
    assert 66.0 in J_EFF["w"]["cd1"]
    assert 53.0 in J_EFF["w"]["cg"]
    assert 53.0 in J_EFF["w"]["cd2"]
    assert 60.0 in J_EFF["w"]["ce2"]
    assert 60.0 in J_EFF["w"]["cz2"]


@pytest.mark.parametrize("spin_name", ("F10CD1-HD1", "Y20CE1-HE1", "W30CH2-HH2"))
def test_cest_13c_spectrometer_builds_for_uniformly_13c_labeled_aromatics(
    spin_name: str,
) -> None:
    config = cest_13c.Cest13CConfig.model_validate(
        {
            "experiment": {
                "name": "cest_13c",
                "time_t1": 0.3,
                "carrier": 128.0,
                "b1_frq": 25.0,
            },
            "conditions": {
                "h_larmor_frq": 600.0,
                "temperature": 25.0,
                "label": ["13c"],
            },
            "data": {},
        }
    )

    spectrometer = cest_13c.build_spectrometer(config, SpinSystem(name=spin_name))

    assert spectrometer.jeff_i.weights.sum() == pytest.approx(1.0)
    assert spectrometer.jeff_i.values.size >= 3


def test_dcest_15n_spectrometer_builds_for_uniformly_13c_labeled_sec_amide() -> None:
    config = dcest_15n.DCest15NConfig.model_validate(
        {
            "experiment": {
                "name": "dcest_15n",
                "time_t1": 0.4,
                "carrier": 120.0,
                "b1_eff": 15.0,
                "pw90": 40.0e-6,
                "sw": 400.0,
            },
            "conditions": {
                "h_larmor_frq": 600.0,
                "temperature": 25.0,
                "label": ["13c"],
            },
            "data": {},
        }
    )

    spectrometer = dcest_15n.build_spectrometer(config, SpinSystem(name="U5N-HN"))

    assert spectrometer.jeff_i.weights.sum() == pytest.approx(1.0)
    assert spectrometer.jeff_i.values.size > 1

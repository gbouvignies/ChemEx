from __future__ import annotations

import pickle

from chemex.parameters.spin_system import SpinSystem
from chemex.parameters.spin_system.atom import Atom
from chemex.parameters.spin_system.constants import STANDARD_ATOM_NAMES
from chemex.parameters.spin_system.group import Group


def test_spin_system_name_parser_does_not_mutate_input_dict() -> None:
    raw = {"name": "G23N-HN"}

    spin_system = SpinSystem.model_validate(raw)

    assert spin_system.name == "G23N-H"
    assert raw == {"name": "G23N-HN"}


def test_group_pickle_rebuilds_search_keys() -> None:
    group = Group("G23")

    restored = pickle.loads(pickle.dumps(group))  # noqa: S301

    assert restored == group
    assert restored in restored.search_keys


def test_atom_pickle_rebuilds_search_keys() -> None:
    atom = Atom("HN")

    restored = pickle.loads(pickle.dumps(atom))  # noqa: S301

    assert restored == atom
    assert restored in restored.search_keys
    assert restored.nucleus in restored.search_keys


def test_spin_system_pickle_rebuilds_nested_search_keys() -> None:
    spin_system = SpinSystem.from_name("G23N-HN")

    _ = spin_system.search_keys
    _ = spin_system.groups
    restored = pickle.loads(pickle.dumps(spin_system))  # noqa: S301

    assert restored == spin_system
    assert restored.search_keys == spin_system.search_keys


def test_strict_spin_system_parser_accepts_every_registered_atom_name() -> None:
    for atom_name in STANDARD_ATOM_NAMES:
        spin_system = SpinSystem.from_name_strict(f"G23{atom_name}")

        assert spin_system.atoms["i"].name == atom_name


def test_spin_system_group_site_identifies_cross_nucleus_companions() -> None:
    target = SpinSystem.from_name("V75HG1")

    assert target.shares_group_site(SpinSystem.from_name("V75CG1"))
    assert not target.shares_group_site(SpinSystem.from_name("V75CG2"))
    assert not target.shares_group_site(SpinSystem.from_name("V71CG1"))

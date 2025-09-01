from __future__ import annotations

from functools import total_ordering
from re import search

from chemex.parameters.spin_system.atom import Atom
from chemex.parameters.spin_system.constants import STANDARD_ATOM_NAMES
from chemex.parameters.spin_system.group import Group


@total_ordering
class Spin:
    def __init__(self, name: str, group_for_completion: Group | None = None) -> None:
        self.group, self.atom = self.split_group_atom(name.strip().upper())
        if not self.group and group_for_completion:
            self.group = group_for_completion
        self.search_keys = self.group.search_keys | self.atom.search_keys

    @staticmethod
    def split_group_atom(name: str) -> tuple[Group, Atom]:
        if name == "?":
            return Group(""), Atom("")
        found_digit = search("[0-9]", name)
        first_digit = found_digit.start() if found_digit else 0
        found_atom = search("[HCNQM]", name[first_digit:])
        if not found_atom:
            if name in STANDARD_ATOM_NAMES:
                return Group(""), Atom(name)
            return Group(name), Atom("")
        atom_index = first_digit + found_atom.start()
        return Group(name[:atom_index]), Atom(name[atom_index:])

    @property
    def name(self) -> str:
        return f"{self.group}{self.atom}"

    def match(self, other: Spin) -> bool:
        return self.group.match(other.group) and self.atom.match(other.atom)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        self_tuple = (self.atom.nucleus.name, self.group, self.atom)
        other_tuple = (other.atom.nucleus.name, other.group, other.atom)
        return self_tuple < other_tuple

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def __bool__(self) -> bool:
        return bool(self.name)

    def __repr__(self) -> str:
        return self.name

    def __str__(self) -> str:
        return self.name

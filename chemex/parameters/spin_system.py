from __future__ import annotations

from collections.abc import Hashable
from collections.abc import Iterable
from collections.abc import Sequence
from dataclasses import dataclass
from dataclasses import field
from enum import auto
from enum import Enum
from functools import cache
from functools import cached_property
from functools import total_ordering
from itertools import chain
from itertools import combinations
from re import search
from typing import Literal
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from chemex.nmr.liouvillian import Basis


SPIN_ALIASES = "isx"


class Nucleus(Enum):
    """Define all the different types of atoms"""

    H1 = auto()
    N15 = auto()
    C13 = auto()
    X = auto()


# Conversion dictionary from atom letter to corresponding Nucleus
STR_TO_NUCLEUS: dict[str, Nucleus] = {
    "H": Nucleus.H1,
    "Q": Nucleus.H1,
    "M": Nucleus.H1,
    "N": Nucleus.N15,
    "C": Nucleus.C13,
    "X": Nucleus.X,
}

# Conversion dictionary from 3-letter to 1-letter amino-acid convention
AAA_TO_A = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


# Conversion dictionary to correct different ways to spell nuclei
CORRECT_ATOM_NAME = {"HN": "H", "C'": "C", "CO": "C"}

# fmt: off
STANDARD_ATOM_NAMES = {
    # Protein
    'C', 'CA', 'CB', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CG', 'CG1', 'CG2',
    'CH2', 'CQD', 'CQE', 'CQG', 'CZ', 'CZ2', 'CZ3', 'H', 'H2', 'H3', 'HA', 'HA2',
    'HA3', 'HB', 'HB1', 'HB2', 'HB3', 'HD', 'HD1', 'HD11', 'HD12', 'HD13', 'HD2',
    'HD21', 'HD22', 'HD23', 'HD3', 'HE', 'HE1', 'HE2', 'HE21', 'HE22', 'HE3', 'HG',
    'HG1', 'HG11', 'HG12', 'HG13', 'HG2', 'HG21', 'HG22', 'HG23', 'HG3', 'HH', 'HH1',
    'HH11', 'HH12', 'HH2', 'HH21', 'HH22', 'HZ', 'HZ1', 'HZ2', 'HZ3', 'MB', 'MD', 'MD1',
    'MD2', 'ME', 'MG', 'MG1', 'MG2', 'MZ', 'N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH',
    'NH1', 'NH2', 'NQH', 'NZ', 'QA', 'QB', 'QD', 'QD2', 'QE', 'QE2', 'QG', 'QG1', 'QH1',
    'QH2', 'QMD', 'QMG', 'QQH', 'QR', 'QZ',
    # DNA/RNA
    "C1'", 'C2', "C2'", "C3'", 'C4', "C4'", 'C5', "C5'", 'C6', 'C7', 'C8',
    'H1', "H1'", 'H2', "H2'", "H2'1", "H2'2", 'H21', 'H22', 'H3', "H3'", "H4'", 'H41', 'H42', 'H5', "H5'1", "H5'2", 'H6', 'H61', 'H62', 'H71', 'H72', 'H73', 'H8',
    'M7', 'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9'}
# fmt: on


@dataclass(order=True)
class Atom:
    name: str
    nucleus: Nucleus = field(init=False, compare=False)
    search_keys: set[Hashable] = field(init=False, default_factory=set, compare=False)

    def __post_init__(self):
        name = self.name.strip().upper()
        self.name = CORRECT_ATOM_NAME.get(name, name)
        self.nucleus = STR_TO_NUCLEUS.get(self.name[:1], Nucleus.X)
        self.search_keys.update({self, self.nucleus})

    def match(self, other: Atom) -> bool:
        return other.name.startswith(self.name)

    def __hash__(self) -> int:
        return hash(self.name)

    def __bool__(self) -> bool:
        return bool(self.name)

    def __str__(self):
        return self.name


@total_ordering
class Group:
    NO_NUMBER: int = -100000000

    def __init__(self, name: str) -> None:
        self.symbol, self.number, self.suffix = self.parse_group(name.strip().upper())
        self.symbol = AAA_TO_A.get(self.symbol, self.symbol)
        self.name = self.get_name()
        self.search_keys: set = {self} if self else set()

    def parse_group(self, name: str) -> tuple[str, int, str]:
        if found := search("[0-9]+", name.strip().upper()):
            return name[: found.start()], int(found.group()), name[found.end() :]
        return name, self.NO_NUMBER, ""

    def get_name(self) -> str:
        number = "" if self.number == self.NO_NUMBER else self.number
        return f"{self.symbol}{number}{self.suffix}"

    def match(self, other: Group) -> bool:
        symbol = other.symbol == self.symbol or not self.symbol
        number = other.number == self.number or self.number == self.NO_NUMBER
        suffix = other.suffix == self.suffix or not self.suffix
        return number and symbol and suffix

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, Group):
            return NotImplemented
        return self.number < other.number

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Group):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def __bool__(self) -> bool:
        return bool(self.name)

    def __str__(self) -> str:
        return self.name


@total_ordering
class Spin:
    def __init__(self, name: str, group_for_completion: Group | None = None) -> None:
        self.group, self.atom = self.split_group_atom(name.strip().upper())
        if not self.group and group_for_completion:
            self.group = group_for_completion
        self.name = self.get_name()
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

    def get_name(self) -> str:
        return f"{self.group}{self.atom}"

    def match(self, other: Spin):
        return self.group.match(other.group) and self.atom.match(other.atom)

    def __lt__(self, other: Spin) -> bool:
        self_tuple = (self.atom.nucleus.name, self.group, self.atom)
        other_tuple = (other.atom.nucleus.name, other.group, other.atom)
        return self_tuple < other_tuple

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Spin):
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


def powerset(iterable):
    "powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))


@total_ordering
class SpinSystem:
    @staticmethod
    @cache
    def _parse_spin_system(name: str) -> dict[str, Spin]:
        if not name:
            return {}
        split = name.split("-")
        spins = {}
        last_group = None
        for short_name, name_spin in zip(SPIN_ALIASES, split):
            spin = Spin(name_spin, last_group)
            spins[short_name] = spin
            last_group = spin.group
        return spins

    @staticmethod
    def _spins2name(spins: Iterable[Spin]) -> str:
        spin_names = []
        last_group: Group = Group("")
        for spin in spins:
            spin_name = str(spin.atom) if spin.group == last_group else str(spin)
            spin_names.append(spin_name)
            last_group = spin.group
        return "-".join(spin_names)

    def __init__(self, name: str | int | None = None):
        if name is None:
            name = ""
        if isinstance(name, int):
            name = str(name)
        self.spins = self._parse_spin_system(name.strip().upper())
        self.name = self._spins2name(self.spins.values())
        self.search_keys: set[Hashable] = set()
        self.search_keys = self.search_keys.union(
            *(spin.search_keys for spin in self.spins.values())
        )

    def __deepcopy__(self, memo) -> SpinSystem:
        return SpinSystem(self.name)

    @cached_property
    def names(self) -> dict[str, str]:
        result = {}
        for alias_set in powerset(SPIN_ALIASES):
            if set(alias_set).issubset(self.spins):
                key = "".join(alias_set)
                name = self._spins2name(self.spins[alias] for alias in alias_set)
                result[key] = name
        return result

    @cached_property
    def groups(self) -> dict[str, Group]:
        return {alias: spin.group for alias, spin in self.spins.items()}

    @cached_property
    def symbols(self) -> dict[str, str]:
        return {alias: group.symbol for alias, group in self.groups.items()}

    @cached_property
    def numbers(self) -> dict[str, int]:
        return {alias: group.number for alias, group in self.groups.items()}

    @cached_property
    def atoms(self) -> dict[str, Atom]:
        return {alias: spin.atom for alias, spin in self.spins.items()}

    @cached_property
    def nuclei(self) -> dict[str, Nucleus]:
        return {alias: atom.nucleus for alias, atom in self.atoms.items()}

    def match(self, other: SpinSystem) -> bool:
        return all(
            spin.match(other_spin)
            for spin, other_spin in zip(self.spins.values(), other.spins.values())
        )

    def part_of(self, selection: Sequence[SpinSystem] | str) -> bool:
        if isinstance(selection, str):
            return selection.lower() in ("all", "*")
        return any(item.match(self) for item in selection)

    def correct(self, basis: Basis) -> SpinSystem:
        spins = []
        last_spin = Spin("")
        for letter, atom in basis.atoms.items():
            spin = self.spins.get(letter, last_spin)
            if not spin.atom.name.startswith(atom.upper()):
                spin.atom = Atom(f"{atom}{spin.atom.name[1:]}")
            last_spin = spin
            spins.append(spin)
        return SpinSystem(self._spins2name(spins))

    def build_sub_spin_system(
        self, spin_system_part: Literal["i", "s", "is", "g", ""]
    ) -> SpinSystem:
        if spin_system_part == "g":
            return SpinSystem(str(self.groups.get("i", "")))
        if spin_system_part == "i":
            return SpinSystem(str(self.spins.get("i", "")))
        if spin_system_part == "s":
            return SpinSystem(str(self.spins.get("s", "")))
        if spin_system_part == "is":
            return self
        return SpinSystem()

    @classmethod
    def validate(cls, v) -> SpinSystem:
        if not isinstance(v, str) and not isinstance(v, int) and v is not None:
            raise TypeError("string required")
        return cls(v)

    @classmethod
    def __get_validators__(cls):
        """Method that allows `pydantic` to validate and convert list
        of strings to list of SpinSystem objects."""
        yield cls.validate

    def __and__(self, other) -> SpinSystem:
        if not isinstance(other, SpinSystem):
            return NotImplemented
        if self == other:
            return SpinSystem(self.name)
        if spins := set(self.spins.values()) & set(other.spins.values()):
            return SpinSystem("-".join(spin.name for spin in spins))
        groups = set(self.groups.values()) & set(other.groups.values())
        if len(groups) == 1:
            return SpinSystem("-".join(group.name for group in groups))
        return SpinSystem("")

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, SpinSystem):
            return NotImplemented
        return tuple(self.spins.values()) < tuple(other.spins.values())

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SpinSystem):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def __bool__(self) -> bool:
        return bool(self.name)

    def __str__(self) -> str:
        return self.name


if __name__ == "__main__":
    print(SpinSystem("G23N-G23HN"), SpinSystem("GLY023N-HN"), SpinSystem("g23N-hN"))
    print(SpinSystem("G23N-G23HN") == SpinSystem("GLY023N-HN"))
    print(SpinSystem("G23N-G23HN").match(SpinSystem("GLY023N-HN")))
    print(SpinSystem("G23N-G23HN") & SpinSystem("G23C"))
    print(SpinSystem("") == SpinSystem())
    group = Group("L99")
    spin = Spin("HD1", group)
    print(f"spin = {spin}, spin.group = {spin.group}, spin.atom = {spin.atom}")
    print(SpinSystem("GLY023N-HN").search_keys)
    print(SpinSystem("GLY023N-HN").names)
    print(SpinSystem("G23C1'").names)

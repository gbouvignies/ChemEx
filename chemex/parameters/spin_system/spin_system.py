from __future__ import annotations

from functools import cache, cached_property, total_ordering
from typing import TYPE_CHECKING, Annotated, Any, Literal, TypeVar

from pydantic import PlainValidator

from .atom import Atom
from .group import Group
from .spin import Spin
from .utilities import powerset

if TYPE_CHECKING:
    from collections.abc import Hashable, Iterable, Sequence

    from chemex.nmr.basis import Basis

    from .nucleus import Nucleus


Self = TypeVar("Self", bound="SpinSystem")

SPIN_ALIASES = "isx"


@total_ordering
class SpinSystem:
    @staticmethod
    @cache
    def _parse_spin_system(name: str) -> dict[str, Spin]:
        if not name:
            return {}
        split = name.split("-")
        spins: dict[str, Spin] = {}
        last_group = None
        for short_name, name_spin in zip(SPIN_ALIASES, split, strict=False):
            spin = Spin(name_spin, last_group)
            spins[short_name] = spin
            last_group = spin.group
        return spins

    @staticmethod
    def _spins2name(spins: Iterable[Spin]) -> str:
        spin_names: list[str] = []
        last_group: Group = Group("")
        for spin in spins:
            spin_name = str(spin.atom) if spin.group == last_group else str(spin)
            spin_names.append(spin_name)
            last_group = spin.group
        return "-".join(spin_names)

    def __init__(self, name: str | int | None = None) -> None:
        if name is None:
            name = ""
        if isinstance(name, int):
            name = str(name)
        self.spins = self._parse_spin_system(name.strip().upper())
        self.name = self._spins2name(self.spins.values())
        self.search_keys: set[Hashable] = set()
        self.search_keys = self.search_keys.union(
            *(spin.search_keys for spin in self.spins.values()),
        )

    def __deepcopy__(self: Self, memo: dict[Any, Any]) -> Self:
        return type(self)(self.name)

    @cached_property
    def names(self) -> dict[str, str]:
        result: dict[str, str] = {}
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

    def match(self: Self, other: Self) -> bool:
        spin_pairs = zip(self.spins.values(), other.spins.values(), strict=False)
        return all(spin.match(other_spin) for spin, other_spin in spin_pairs)

    def part_of(self: Self, selection: Sequence[Self] | str) -> bool:
        if isinstance(selection, str):
            return selection.lower() in ("all", "*")
        return any(item.match(self) for item in selection)

    def correct(self: Self, basis: Basis) -> Self:
        spins: list[Spin] = []
        last_spin = Spin("")
        for letter, atom in basis.atoms.items():
            spin = Spin(self.spins.get(letter, last_spin).name)
            if not spin.atom.name.startswith(atom.upper()):
                spin.atom = Atom(f"{atom}{spin.atom.name[1:]}")
            last_spin = spin
            spins.append(spin)
        return type(self)(self._spins2name(spins))

    def build_sub_spin_system(
        self: Self,
        spin_system_part: Literal["i", "s", "is", "g", ""],
    ) -> Self:
        if spin_system_part == "g":
            return type(self)(str(self.groups.get("i", "")))
        if spin_system_part == "i":
            return type(self)(str(self.spins.get("i", "")))
        if spin_system_part == "s":
            return type(self)(str(self.spins.get("s", "")))
        if spin_system_part == "is":
            return self
        return type(self)()

    def __and__(self: Self, other: object) -> Self:
        if not isinstance(other, type(self)):
            return NotImplemented
        if self == other:
            return type(self)(self.name)
        if spins := set(self.spins.values()) & set(other.spins.values()):
            return type(self)("-".join(spin.name for spin in spins))
        groups = set(self.groups.values()) & set(other.groups.values())
        if len(groups) == 1:
            return type(self)("-".join(group.name for group in groups))
        return type(self)("")

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return tuple(self.spins.values()) < tuple(other.spins.values())

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def __bool__(self) -> bool:
        return bool(self.name)

    def __str__(self) -> str:
        return self.name


# We now create an `Annotated` wrapper that we'll use as the annotation for
# fields on `BaseModel`
PydanticSpinSystem = Annotated[SpinSystem, PlainValidator(lambda x: SpinSystem(x))]


# if __name__ == "__main__":
#     print(SpinSystem("G23N-G23HN"), SpinSystem("GLY023N-HN"), SpinSystem("g23N-hN"))
#     print(SpinSystem("G23N-G23HN") == SpinSystem("GLY023N-HN"))
#     print(SpinSystem("G23N-G23HN").match(SpinSystem("GLY023N-HN")))
#     print(SpinSystem("G23N-G23HN") & SpinSystem("G23C"))
#     print(SpinSystem("") == SpinSystem())
#     group = Group("L99")
#     spin = Spin("HD1", group)
#     print(f"spin = {spin}, spin.group = {spin.group}, spin.atom = {spin.atom}")
#     print(SpinSystem("GLY023N-HN").search_keys)
#     print(SpinSystem("GLY023N-HN").names)
#     print(SpinSystem("G23C1'").names)

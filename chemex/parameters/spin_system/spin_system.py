from __future__ import annotations

from collections.abc import Hashable, Iterable, Sequence
from functools import cache, cached_property, total_ordering
from typing import TYPE_CHECKING, Any, Literal, Self

from pydantic import BaseModel, Field, InstanceOf, computed_field, model_validator

from chemex.parameters.spin_system.atom import Atom
from chemex.parameters.spin_system.group import Group
from chemex.parameters.spin_system.nucleus import Nucleus
from chemex.parameters.spin_system.spin import Spin
from chemex.parameters.spin_system.utilities import powerset

if TYPE_CHECKING:
    from chemex.nmr.basis import Basis

SPIN_ALIASES = "isx"


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


def _spins2name(spins: Iterable[Spin]) -> str:
    spin_names: list[str] = []
    last_group: Group = Group("")
    for spin in spins:
        spin_name = str(spin.atom) if spin.group == last_group else str(spin)
        spin_names.append(spin_name)
        last_group = spin.group
    return "-".join(spin_names)


@total_ordering
class SpinSystem(BaseModel):
    name: str = ""
    spins: dict[str, InstanceOf[Spin]] = Field(default_factory=dict)

    @model_validator(mode="before")
    @classmethod
    def parse_name(cls, model: Any) -> Any:
        if isinstance(model, dict) and "name" in model:
            name = str(model["name"]).strip().upper()
            model["spins"] = _parse_spin_system(name)
            model["name"] = _spins2name(model["spins"].values())
        return model

    @classmethod
    def from_name(cls, name: int | str) -> Self:
        if isinstance(name, int):
            name = str(name)
        return cls(name=name)

    @computed_field
    @cached_property
    def search_keys(self) -> set[Hashable]:
        search_keys: set[Hashable] = set()
        return search_keys.union(*(spin.search_keys for spin in self.spins.values()))

    @cached_property
    def names(self) -> dict[str, str]:
        result: dict[str, str] = {}
        for alias_set in powerset(SPIN_ALIASES):
            if set(alias_set).issubset(self.spins):
                key = "".join(alias_set)
                name = _spins2name(self.spins[alias] for alias in alias_set)
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

    def match(self, other: Self) -> bool:
        spin_pairs = zip(self.spins.values(), other.spins.values(), strict=False)
        return all(spin.match(other_spin) for spin, other_spin in spin_pairs)

    def part_of(self, selection: Sequence[Self] | str) -> bool:
        if isinstance(selection, str):
            return selection.lower() in ("all", "*")
        return any(item.match(self) for item in selection)

    def correct(self, basis: Basis) -> Self:
        spins: list[Spin] = []
        last_spin = Spin("")
        for letter, atom in basis.atoms.items():
            spin = Spin(self.spins.get(letter, last_spin).name)
            if not spin.atom.name.startswith(atom.upper()):
                spin.atom = Atom(f"{atom}{spin.atom.name[1:]}")
            last_spin = spin
            spins.append(spin)
        return type(self)(name=_spins2name(spins))

    def build_sub_spin_system(
        self,
        spin_system_part: Literal["i", "s", "is", "g", ""],
    ) -> Self:
        if spin_system_part == "g":
            return type(self)(name=str(self.groups.get("i", "")))
        if spin_system_part == "i":
            return type(self)(name=str(self.spins.get("i", "")))
        if spin_system_part == "s":
            return type(self)(name=str(self.spins.get("s", "")))
        if spin_system_part == "is":
            return self
        return type(self)()

    def __and__(self, other: object) -> Self:
        if not isinstance(other, type(self)):
            return NotImplemented
        if self == other:
            return type(self)(name=self.name)
        if spins := set(self.spins.values()) & set(other.spins.values()):
            return type(self)(name="-".join(spin.name for spin in spins))
        groups = set(self.groups.values()) & set(other.groups.values())
        if len(groups) == 1:
            return type(self)(name="-".join(group.name for group in groups))
        return type(self)()

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


# if __name__ == "__main__":
#     print(
#         SpinSystem(name="G23N-G23HN"),
#         SpinSystem(name="GLY023N-HN"),
#         SpinSystem(name="g23N-hN"),
#     )
#     print(SpinSystem(name="G23N-G23HN") == SpinSystem(name="GLY023N-HN"))
#     print(SpinSystem(name="G23N-G23HN").match(SpinSystem(name="GLY023N-HN")))
#     print(SpinSystem(name="G23N-G23HN") & SpinSystem(name="G23C"))
#     print(SpinSystem(name="") == SpinSystem())
#     group = Group("L99")
#     spin = Spin("HD1", group)
#     print(f"spin = {spin}, spin.group = {spin.group}, spin.atom = {spin.atom}")
#     print(SpinSystem(name="GLY023N-HN").search_keys)
#     print(SpinSystem(name="GLY023N-HN").names)
#     print(SpinSystem(name="G23C1'").names)

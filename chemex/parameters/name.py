from __future__ import annotations

import re
from collections.abc import Hashable
from collections.abc import Iterable
from dataclasses import dataclass
from dataclasses import field
from functools import cached_property
from re import Pattern

from rapidfuzz.process import extractOne

from chemex.configuration.conditions import Conditions
from chemex.parameters.spin_system import SpinSystem

_EXPAND = {"-": "_", "+": "_", ".": "_"}

_FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"

_RE_SECTION = re.compile(
    rf"""
        (^\s*(?P<name>\w+)) |
        (NUC\s*->\s*(?P<spin_system>(\w|-)+)) |
        (T\s*->s*(?P<temperature>{_FLOAT})) |
        (B0\s*->\s*(?P<h_larmor_frq>{_FLOAT})) |
        (\[P]\s*->\s*(?P<p_total>{_FLOAT})) |
        (\[L]\s*->\s*(?P<l_total>{_FLOAT})) |
        (D2O\s*->\s*(?P<d2o>{_FLOAT}))
    """,
    re.IGNORECASE | re.VERBOSE,
)


def _parse(re_to_match: Pattern, text: str) -> dict[str, str]:
    return {
        key: value
        for match in re_to_match.finditer(text)
        for key, value in match.groupdict().items()
        if value is not None
    }


def _multireplace(string: str, replacements: dict[str, str]) -> str:
    """
    Given a string and a replacement map, it returns the replaced string.

    :param str string: string to execute replacements on
    :param dict replacements: replacement dictionary {value to find: value to replace}
    :rtype: str

    """
    # Place longer ones first to keep shorter substrings from matching where the longer
    # ones should take place. For instance given the replacements
    # {'ab': 'AB', 'abc': 'ABC'} against the string 'hey abc', it should produce
    # 'hey ABC' and not 'hey ABc'
    substrings = sorted(replacements, key=len, reverse=True)

    # Create a big OR regex that matches any of the substrings to replace
    regexp = re.compile("|".join(re.escape(substring) for substring in substrings))

    # For each match, look up the new string in the replacements
    return regexp.sub(lambda match: replacements[match.group(0)], string)


def _expand(string: str) -> str:
    return _multireplace(string, _EXPAND)


@dataclass(order=True, unsafe_hash=True)
class ParamName:
    name: str = ""
    spin_system: SpinSystem = field(default_factory=SpinSystem)
    conditions: Conditions = field(default_factory=Conditions.construct)
    search_keys: set[Hashable] = field(init=False, compare=False)

    def __post_init__(self):
        self.name = self.name.strip().upper()
        self.conditions = self.conditions.rounded()
        self.search_keys = (
            self.spin_system.search_keys | self.conditions.search_keys | {self.name}
        )

    @classmethod
    def from_section(cls, section: str = "") -> ParamName:
        parsed = _parse(_RE_SECTION, section.strip(' []"'))
        name = parsed.get("name")
        if name is None:
            name = ""
        spin_system = SpinSystem(parsed.get("spin_system"))
        conditions = Conditions.parse_obj(parsed)
        return ParamName(name, spin_system, conditions)

    @cached_property
    def section(self) -> str:
        parts = (self.name, self.conditions.section)
        return ", ".join(parts).strip(", ").upper()

    @cached_property
    def section_res(self) -> str:
        parts = [self.name]
        if self.spin_system:
            parts.append(f"NUC->{self.spin_system}")
        parts.append(self.conditions.section)
        return ", ".join(parts).strip(", ").upper()

    @cached_property
    def folder(self) -> str:
        parts = (self.name, str(self.spin_system), self.conditions.folder)
        return "_".join(parts).strip("_").upper()

    @cached_property
    def id(self) -> str:
        return f"__{_expand(self.folder)}"

    def match(self, other: ParamName) -> bool:
        if self.name and self.name not in other.name:
            return False
        if not self.spin_system.match(other.spin_system):
            return False
        if not self.conditions.match(other.conditions):
            return False
        return True

    def get_closest_id(self, ids: Iterable[str]) -> str:
        best_match, _, _ = extractOne(self.id, ids)
        return best_match

    def __and__(self, other: ParamName) -> ParamName:
        name = self.name if self.name == other.name else ""
        spin_system = self.spin_system & other.spin_system
        conditions = self.conditions & other.conditions
        return ParamName(name, spin_system, conditions)

    def __bool__(self) -> bool:
        return bool(self.name) or bool(self.spin_system) or bool(self.conditions)

    def __str__(self) -> str:
        return f"[{self.section_res}]"

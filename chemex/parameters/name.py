from __future__ import annotations

import dataclasses as dc
import re
from functools import cached_property
from functools import total_ordering
from typing import Any

from chemex.containers.conditions import Conditions
from chemex.nmr.spin_system import SpinSystem


_DECORATORS = {
    "name": "__n_{}_n__",
    "spin_system": "__r_{}_r__",
    "temperature": "__t_{}_t__",
    "h_larmor_frq": "__b_{}_b__",
    "p_total": "__p_{}_p__",
    "l_total": "__l_{}_l__",
    "d2o": "__d_{}_d__",
}
_EXPAND = {"-": "__minus__", "+": "__plus__", ".": "__point__"}
_RE_FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
_RE_NAME = re.compile(r"^(?P<name>.+?)(_(?P<spin>i|s|is))?(_(?P<state>[a-h]{1,2}))?$")
_RE_SECTION = re.compile(
    r"""
                (^\s*(?P<name>\w+)) |
                (NUC\s*->\s*(?P<spin_system>(\w|-)+)) |
                (T\s*->s*(?P<temperature>{0})) |
                (B0\s*->\s*(?P<h_larmor_frq>{0})) |
                (\[P]\s*->\s*(?P<p_total>{0})) |
                (\[L]\s*->\s*(?P<l_total>{0})) |
                (D2O\s*->\s*(?P<d2o>{0})) |
            """.format(
        _RE_FLOAT
    ),
    re.IGNORECASE | re.VERBOSE,
)


@total_ordering
class ParamName:
    name: str
    spin_system: SpinSystem
    conditions: Conditions

    def __init__(
        self,
        name: str | None = None,
        spin_system: SpinSystem | None = None,
        conditions: Conditions | None = None,
    ):
        if name is None:
            name = ""

        if spin_system is None:
            spin_system = SpinSystem("")

        if conditions is None:
            conditions = Conditions()

        self.name = name.strip().upper()
        self.spin_system = spin_system
        self.conditions = conditions.rounded()

    @classmethod
    def from_dict(cls, dict_: dict[str, Any]) -> ParamName:
        name = dict_.get("name")
        spin_system = SpinSystem(str(dict_.get("spin_system", "")))
        conditions = Conditions.from_dict(dict_)
        return cls(name, spin_system, conditions)

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "spin_system": self.spin_system,
            **dc.asdict(self.conditions),
        }

    @classmethod
    def from_section(cls, section: str = "") -> ParamName:
        parsed = _re_to_dict(_RE_SECTION, section)
        return cls.from_dict(parsed)

    @cached_property
    def full(self) -> str:
        full_name = "".join(
            _DECORATORS[name].format(value)
            for name, value in self.to_dict().items()
            if value
        )
        return _expand(full_name)

    def _to_section_name(self, show_spin_system: bool = False) -> str:
        formatters = {
            "name": "{}",
            "spin_system": "NUC->{}",
            "temperature": "T->{:.1f}C",
            "h_larmor_frq": "B0->{:.1f}MHz",
            "p_total": "[P]->{:e}M",
            "l_total": "[L]->{:e}M",
            "d2o": "D2O->{:.4f}",
        }
        components = [
            formatters[name].format(value)
            for name, value in self.to_dict().items()
            if name in formatters
            and value not in [None, ""]
            and (name != "spin_system" or (show_spin_system and value))
        ]
        return ",".join(components).upper()

    @cached_property
    def section(self) -> str:
        return self._to_section_name()

    @cached_property
    def section_res(self) -> str:
        return self._to_section_name(show_spin_system=True)

    @cached_property
    def folder(self) -> str:
        formatters = {
            "name": "{}",
            "spin_system": "{}",
            "temperature": "{:.1f}C",
            "h_larmor_frq": "{:.1f}MHz",
            "p_total": "P{:e}M",
            "l_total": "L{:e}M",
            "d2o": "D{:.4f}",
        }
        components = [
            formatters[name].format(value)
            for name, value in self.to_dict().items()
            if value and name in formatters
        ]
        return "_".join(components).upper()

    def match(self, other: ParamName) -> bool:
        name_match = self.name == other.name or not self.name
        spin_system_match = self.spin_system.match(other.spin_system)
        conditions_match = self.conditions.match(other.conditions)
        return name_match and spin_system_match and conditions_match

    def __repr__(self) -> str:
        return f"[{self.section_res}]"

    def __hash__(self) -> int:
        return hash(self._members())

    def __eq__(self, other: ParamName) -> bool:
        return self._members() == other._members()

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, ParamName):
            return NotImplemented
        if not self.spin_system and other.spin_system:
            return True
        elif self.spin_system and not other.spin_system:
            return False
        return self._members() < other._members()

    def __and__(self, other: ParamName) -> ParamName:
        both = {
            "name": self.name if self.name == other.name else None,
            "spin_system": self.spin_system & other.spin_system,
            "conditions": self.conditions & other.conditions,
        }
        return ParamName(**both)

    def __bool__(self) -> bool:
        return bool(self.name) or bool(self.spin_system) or bool(self.conditions)

    def _members(self):
        temperature = h_larmor_frq = p_total = l_total = -1e16
        conditions = self.conditions
        if conditions.temperature is not None:
            temperature = conditions.temperature
        if conditions.h_larmor_frq is not None:
            h_larmor_frq = conditions.h_larmor_frq
        if conditions.p_total is not None:
            p_total = conditions.p_total
        if self.conditions.l_total is not None:
            l_total = conditions.l_total
        return (
            self.name,
            temperature,
            h_larmor_frq,
            p_total,
            l_total,
            self.spin_system,
        )


def _expand(string):
    return _multireplace(str(string), _EXPAND)


def _multireplace(string, replacements):
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


def _re_to_dict(re_to_match, text):
    return {
        key: value
        for match in re_to_match.finditer(text)
        for key, value in match.groupdict().items()
        if value is not None
    }


def get_pnames(settings, conditions=None, spin_system=None):
    pnames = {}
    cdict = dc.asdict(conditions)
    for name, setting in settings.items():
        attributes = setting["attributes"]
        cdict_subset = {key: cdict[key] for key in cdict if key in attributes}
        conditions_ = Conditions.from_dict(cdict_subset)
        name_ = _squeeze_name(name, setting.get("ext", ""))
        spin_system_ = None
        if "spin_system" in attributes:
            spin_system_ = _name_to_spin_system(name, spin_system)
        pnames[name] = ParamName(name_, spin_system_, conditions_)
    return pnames


def _squeeze_name(name, ext=""):
    parsed = _re_to_dict(_RE_NAME, name)
    if "spin" in parsed:
        return f"{parsed['name']}{ext}_{parsed['state']}"
    return name


def _name_to_spin_system(name: str, spin_system: SpinSystem) -> SpinSystem:
    """TODO: Refactor this"""
    if not spin_system:
        return spin_system
    parsed = _re_to_dict(_RE_NAME, name)
    spin = parsed.get("spin")
    if spin is None:
        return SpinSystem(spin_system.names.get("i", ""))
    return SpinSystem(spin_system.names.get(spin, ""))

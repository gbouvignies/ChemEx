import dataclasses as dc
import functools as ft
import re

import chemex.nmr.spin_system as cns


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
_COMPRESS = {val: key for key, val in _EXPAND.items()}
_RE_FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
_RE_NAME = re.compile(r"^(?P<name>.+?)(_(?P<spin>i|s|is))?(_(?P<state>[a-h]{1,2}))?$")


@ft.total_ordering
class ParamName:
    def __init__(
        self,
        name=None,
        spin_system=None,
        temperature=None,
        h_larmor_frq=None,
        p_total=None,
        l_total=None,
        d2o=None,
    ):
        if name is None:
            name = ""
        self.name = name.lower()
        self.spin_system = spin_system
        self.temperature = (
            round(float(temperature), 1) if temperature is not None else None
        )
        self.h_larmor_frq = (
            round(float(h_larmor_frq), 1) if h_larmor_frq is not None else None
        )
        self.p_total = float(p_total) if p_total is not None else None
        self.l_total = float(l_total) if l_total is not None else None
        self.d2o = float(d2o) if d2o is not None else None

    @property
    def spin_system(self):
        return self._spin_system.name

    @spin_system.setter
    def spin_system(self, value):
        if isinstance(value, cns.SpinSystem):
            self._spin_system = value
        else:
            self._spin_system = cns.SpinSystem(value)

    def to_dict(self):
        return {
            "name": self.name,
            "spin_system": self.spin_system,
            "temperature": self.temperature,
            "h_larmor_frq": self.h_larmor_frq,
            "p_total": self.p_total,
            "l_total": self.l_total,
            "d2o": self.d2o,
        }

    @classmethod
    def from_full_name(cls, full_name=""):
        parser = re.compile(
            r"""
                (__n_(?P<name>.*?)_n__)?
                (__r_(?P<spin_system>.*?)_r__)?
                (__t_(?P<temperature>.*?)_t__)?
                (__b_(?P<h_larmor_frq>.*?)_b__)?
                (__p_(?P<p_total>.*?)_p__)?
                (__l_(?P<l_total>.*?)_l__)?
                (__d_(?P<d2o>.*?)_d__)?
            """,
            re.IGNORECASE | re.VERBOSE,
        )
        return cls(**_re_to_dict(parser, _compress(full_name)))

    def to_full_name(self):
        full_name = "".join(
            _DECORATORS[name].format(value)
            for name, value in self.to_dict().items()
            if value
        )
        return _expand(full_name)

    @classmethod
    def from_section(cls, section=""):
        parser = re.compile(
            r"""
                (^\s*(?P<name>\w+)) |
                (NUC\s*->\s*(?P<spin_system>(\w|-)+)) |
                (T\s*->s*(?P<temperature>{0})) |
                (B0\s*->\s*(?P<h_larmor_frq>{0})) |
                (\[P\]\s*->\s*(?P<p_total>{0})) |
                (\[L\]\s*->\s*(?P<l_total>{0})) |
                (\[D2O/H2O\]\s*->\s*(?P<d2o>{0})) |
            """.format(
                _RE_FLOAT
            ),
            re.IGNORECASE | re.VERBOSE,
        )
        parsed = _re_to_dict(parser, section)
        return cls(**parsed)

    def to_section_name(self, show_spin_system=False):
        formatters = {
            "name": "{}",
            "spin_system": "NUC->{}",
            "temperature": "T->{:.1f}C",
            "h_larmor_frq": "B0->{:.1f}MHz",
            "p_total": "[P]->{:e}M",
            "l_total": "[L]->{:e}M",
            "d2o": "D2O/H2O->{:.4f}",
        }
        components = []
        for name, value in self.to_dict().items():
            to_be_showed = (
                (name != "spin_system" or show_spin_system)
                and value is not None
                and value != ""
            )
            if to_be_showed:
                components.append(formatters[name].format(value))
        return ", ".join(components).upper()

    def to_folder_name(self):
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
            if value
        ]
        return "-".join(components).upper()

    def match(self, other):
        re_name = self._to_re()
        return re_name.match(other)

    @ft.lru_cache()
    def _to_re(self):
        if self._spin_system is not None:
            spin_system = self._spin_system.to_re().pattern
        else:
            spin_system = ""
        pattern = _get_re_component(self.name, "name")
        pattern += _get_re_component(spin_system, "spin_system")
        pattern += _get_re_component(self.temperature, "temperature")
        pattern += _get_re_component(self.h_larmor_frq, "h_larmor_frq")
        pattern += _get_re_component(self.p_total, "p_total")
        pattern += _get_re_component(self.l_total, "l_total")
        pattern += _get_re_component(self.d2o, "d2o")
        pattern = _clean_re(pattern)
        return re.compile(pattern, re.IGNORECASE)

    def __repr__(self):
        return self.to_section_name(show_spin_system=True)

    def __hash__(self):
        return hash(self._members())

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ParamName):
            return NotImplemented
        return self._members() == other._members()

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, ParamName):
            return NotImplemented
        return self._members() < other._members()

    def __and__(self, other: object) -> "ParamName":
        if not isinstance(other, ParamName):
            return NotImplemented
        both = {
            name: getattr(self, name)
            for name in (
                "name",
                "temperature",
                "h_larmor_frq",
                "p_total",
                "l_total",
                "d2o",
            )
            if getattr(self, name) == getattr(other, name)
        }
        both["spin_system"] = self._spin_system & other._spin_system
        return ParamName(**both)

    def __bool__(self):
        return bool(str(self))

    def _members(self):
        return (
            self.name,
            self.temperature,
            self.h_larmor_frq,
            self.p_total,
            self.l_total,
            self._spin_system,
        )


def _get_re_component(value, kind):
    if not value:
        return "([a-z_][a-z0-9_]*)?"
    if kind != "spin_system":
        value = _expand(value)
    return _DECORATORS[kind].format(value)


def _clean_re(value):
    pattern1 = "([a-z_][a-z0-9_]*)?"
    pattern2 = f"({re.escape(pattern1)})+$"
    pattern3 = f"({re.escape(pattern1)})+"
    value_ = re.sub(pattern2, "", value)
    value_ = re.sub(pattern3, pattern1, value_)
    return value_


def _expand(string):
    return _multireplace(str(string), _EXPAND)


def _compress(string):
    return _multireplace(str(string), _COMPRESS)


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
    regexp = re.compile("|".join([re.escape(substring) for substring in substrings]))

    # For each match, look up the new string in the replacements
    return regexp.sub(lambda match: replacements[match.group(0)], string)


def _re_to_dict(re_to_match, text):
    return {
        key: value
        for match in re_to_match.finditer(text)
        for key, value in match.groupdict().items()
        if value is not None
    }


def remove_state(name):
    parsed = _re_to_dict(_RE_NAME, name)
    name_ = f"{parsed['name']}"
    if "spin" in parsed:
        name_ += f"_{parsed['spin']}"
    return name_


def get_pnames(settings, conditions=None, spin_system=None):
    attributes = dc.asdict(conditions)
    if spin_system is not None:
        attributes["spin_system"] = spin_system
    pnames = {}
    for name, setting in settings.items():
        attributes_ = {key: attributes.get(key) for key in setting["attributes"]}
        name_ = _squeeze_name(name, setting.get("ext", ""))
        spin_system_ = _name_to_spin_system(name, attributes_.get("spin_system"))
        attributes_["spin_system"] = spin_system_
        pnames[name] = ParamName(name_, **attributes_).to_full_name()
    return pnames


def _squeeze_name(name, ext=""):
    parsed = _re_to_dict(_RE_NAME, name)
    if "spin" in parsed:
        return f"{parsed['name']}{ext}_{parsed['state']}"
    return name


def _name_to_spin_system(name, spin_system):
    if spin_system is None:
        return None
    parsed = _re_to_dict(_RE_NAME, name)
    spin = parsed.get("spin")
    if spin is None:
        return spin_system.names["i"]
    return spin_system.names.get(spin)

"""The peaks module contains the code for handling peak assignments and
resonances."""
import functools as ft
import itertools as it
import re
import string


_ALIASES = "isx"
RE_NAME = re.compile(
    r"""
        (^\s*|-)
        (
            (?P<symbol>\D?)              # one letter amino acid (optional)
            0*(?P<number>[0-9]+|[*])     # residue number
            (?P<suffix>[abd-gi-mopr-z]*) # suffix (optional)
        )?
        (?P<nucleus>                     # nucleus name (e.g., CA, HG, ...)
            (?P<atom>[hncq])             # nucleus type
            [a-z0-9]*                    # nucleus name - nucleus type
        )?
    """,
    re.IGNORECASE | re.VERBOSE,
)


@ft.total_ordering
class SpinSystem:
    def __init__(self, name=None):
        if name is None:
            name = ""
        self.name = str(name)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        spin_list = _name_to_spins(str(value).upper())
        self._spins = dict(zip(_ALIASES, spin_list))
        self._name = _spins_to_name(self._spins.values())

    @property
    def names(self):
        result = {}
        for aliases in _powerset(_ALIASES):
            if set(aliases).issubset(self._spins):
                key = "".join(aliases)
                name = _spins_to_name(self._spins[alias] for alias in aliases)
                result[key] = name
        return result

    @property
    def groups(self):
        return {key: spin["group"] for key, spin in self._spins.items()}

    @property
    def symbols(self):
        return {key: spin["symbol"] for key, spin in self._spins.items()}

    @property
    def atoms(self):
        return {key: spin["atom"] for key, spin in self._spins.items()}

    @property
    def nuclei(self):
        return {key: spin["nucleus"] for key, spin in self._spins.items()}

    @property
    def numbers(self):
        return {
            key: int(spin["number"])
            for key, spin in self._spins.items()
            if spin["number"]
        }

    def to_re(self):
        return _spins_to_re(self._spins)

    def part_of(self, selection):
        selection_ = (SpinSystem(item) for item in selection)
        return any(self & name == name for name in selection_)

    def __and__(self, other):
        if isinstance(other, str):
            other = SpinSystem(other)
        if self.name == other.name:
            return SpinSystem(self.name)
        names = set(self.names.values()) & set(other.names.values())
        if names:
            return SpinSystem("-".join(names))
        groups = set(self.groups.values()) & set(other.groups.values())
        if len(groups) == 1:
            return SpinSystem(groups.pop())
        numbers = set(self.numbers.values()) & set(other.numbers.values())
        if len(numbers) == 1:
            return SpinSystem(numbers.pop())
        return SpinSystem()

    def __repr__(self):
        return str(self.name).upper()

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        if isinstance(other, str):
            other = SpinSystem(other)
        self_tuple = tuple(
            zip(self.atoms.values(), self.numbers.values(), self.nuclei.values())
        )
        other_tuple = tuple(
            zip(other.atoms.values(), other.numbers.values(), other.nuclei.values())
        )
        return self_tuple < other_tuple

    def __len__(self):
        return len(self._spins)


def _name_to_spins(name):
    """Get spins from an assignment."""
    spins = []
    last_spin = {}
    for match in re.finditer(RE_NAME, name):
        spin = match.groupdict()
        for key, value in spin.items():
            if value is None:
                spin[key] = last_spin.get(key, "")
        spin["group"] = "{symbol}{number}{suffix}".format_map(spin)
        spin["name"] = "{group}{nucleus}".format_map(spin)
        spins.append(spin)
        last_spin = spin
    return spins


def _spins_to_name(spins):
    """Get assignment from resonances."""
    parts = []
    last_spin = {}
    for spin in spins:
        if spin["group"] != last_spin.get("group", ""):
            parts.append("{group}{nucleus}".format_map(spin))
        else:
            parts.append("{nucleus}".format_map(spin))
        last_spin = spin
    return "-".join(parts)


def _spins_to_re(spins):
    if _is_empty(spins):
        return re.compile("")
    last_spin = {}
    re_expr = ""
    for index, spin in enumerate(spins.values()):
        if index > 0:
            re_expr += "(-|__minus__)"
        if spin["group"] != last_spin.get("group"):
            re_expr += f"(?P<group{index}>"
            re_expr += _to_re(spin["symbol"], r"(\D?)")
            re_expr += _to_re(spin["number"], r"([0-9]+)")
            re_expr += _to_re(spin["suffix"], r"([abd-gi-mopr-z]*)")
            re_expr += ")"
        else:
            re_expr += f"(?P=group{index-1})?"
        if spin["nucleus"] != spin["atom"]:
            re_expr += _to_re(spin["nucleus"], r"[hncq][a-z0-9]*")
        else:
            re_expr += _to_re(spin["atom"], r"[hncq]")
            re_expr += r"[a-z0-9]*"
        last_spin = spin
    return re.compile(re_expr, re.IGNORECASE)


def _is_empty(spins):
    values = {value for spin in spins.values() for value in spin.values()}
    return len(values) == 1 and values.pop() == ""


def _to_re(value, default):
    if value in ("", "*"):
        return default
    return value


def _powerset(iterable):
    """powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return it.chain.from_iterable(it.combinations(s, r + 1) for r in range(len(s)))


def get_state_names(state_nb):
    return string.ascii_lowercase[:state_nb]

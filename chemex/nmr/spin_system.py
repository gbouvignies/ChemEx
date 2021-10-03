"""The peaks module contains the code for handling peak assignments and
resonances."""
import functools as ft
import itertools as it
import re


_ALIASES = "isx"
RE_NAME = re.compile(
    r"""
        (^\s*|-)
        (
            (?P<symbol>(\D?|\D{3}?))              # one letter amino acid (optional)
            0*(?P<number>[0-9]+|[*])     # residue number
            (?P<suffix>[abd-gi-mopr-wyz]*) # suffix (optional)
        )?
        (?P<nucleus>                     # nucleus name (e.g., CA, HG, ...)
            (?P<atom>[hncqx])             # nucleus type
            [a-z0-9]*                    # nucleus name - nucleus type
        )
    """,
    re.IGNORECASE | re.VERBOSE,
)

RE_NAME_GROUP = re.compile(
    r"""
        (^\s*|-)
        (
            (?P<symbol>(\D?|\D{3}?))              # one letter amino acid (optional)
            0*(?P<number>[0-9]+|[*])     # residue number
            (?P<suffix>[abd-gi-mopr-wyz]*) # suffix (optional)
        )
    """,
    re.IGNORECASE | re.VERBOSE,
)

_AA_CODE = {
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


@ft.total_ordering
class SpinSystem:
    def __init__(self, name=None):
        if name is None:
            name = ""
        self.name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._spins = _name_to_spins(str(value).upper())
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
    def xnames(self):
        result = {}
        for aliases in _powerset(_ALIASES):
            if set(aliases).issubset(self._spins):
                key = "".join(aliases)
                name = _spins_to_xname(self._spins[alias] for alias in aliases)
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

    def correct(self, basis):

        if not self:
            return self

        atoms = {
            letter: atom for letter, atom in basis.atoms.items() if letter in basis.type
        }
        spins = {}
        for letter, atom in atoms.items():
            for spin in self._spins.values():
                if spin["atom"].upper() == atom.upper() and letter in basis.type:
                    spins[letter] = spin
                    break

        for letter, atom in atoms.items():
            if letter in basis.type and letter not in spins:
                spins[letter] = spins["i"].copy()
                spins[letter]["nucleus"] = f"{atom}{spins[letter]['nucleus'][1:]}"

        return SpinSystem(_spins_to_name(spins.values()))

    def __and__(self, other):
        if isinstance(other, str):
            other = SpinSystem(other)
        if self.name == other.name:
            return SpinSystem(self.name)
        names = set(self.names.values()) & set(other.names.values())
        if names:
            return SpinSystem("-".join(names))
        xnames = set(self.xnames.values()) & set(other.xnames.values())
        if xnames:
            return SpinSystem("-".join(xnames))
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
        if isinstance(other, str):
            other = SpinSystem(other)
        if isinstance(other, SpinSystem):
            return self.name == other.name
        return NotImplemented

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

    def __bool__(self):
        return bool(self._spins)


def _name_to_spins(name):
    """Get spins from an assignment."""
    spins = []
    last_spin = {}
    re_name = RE_NAME if re.match(RE_NAME, name) else RE_NAME_GROUP
    for match in re.finditer(re_name, name):
        spin = {k: "" for k in ("symbol", "number", "suffix", "nucleus")}
        spin.update(match.groupdict())
        if not any(spin.values()):
            continue
        if spin["symbol"] is not None and spin["symbol"].upper() in _AA_CODE:
            spin["symbol"] = _AA_CODE[spin["symbol"]]
        if spin["nucleus"] is not None and spin["nucleus"].upper() == "HN":
            spin["nucleus"] = "H"
        for key, value in spin.items():
            if value is None:
                spin[key] = last_spin.get(key, "")
        spin["group"] = "{symbol}{number}{suffix}".format_map(spin)
        spin["name"] = "{group}{nucleus}".format_map(spin)
        spins.append(spin)
        last_spin = spin
    return dict(zip(_ALIASES, spins))


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


def _spins_to_xname(spins):
    """Get assignment from resonances."""
    parts = []
    last_spin = {}
    for spin in spins:
        nucleus = spin["nucleus"]
        xnucleus = f"X{nucleus[1:]}" if nucleus[1:] else nucleus
        if spin["group"] != last_spin.get("group", ""):
            parts.append(f"{spin['group']}{xnucleus}")
        else:
            parts.append(f"{xnucleus}")
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
        if spin["group"] == last_spin.get("group"):
            re_expr += f"(?P=group{index-1})?"
        else:
            re_expr += f"(?P<group{index}>"
            re_expr += _to_re(spin["symbol"], r"(\D?)")
            re_expr += _to_re(spin["number"], r"([0-9]+)")
            re_expr += _to_re(spin["suffix"], r"([abd-gi-mopr-z]*)")
            re_expr += ")"
        if spin["nucleus"] == spin["atom"]:
            re_expr += _to_re(spin["atom"], r"[hncq]")
            re_expr += r"[a-z0-9]*"
        else:
            re_expr += _to_re(spin["nucleus"], r"[hncq][a-z0-9]*")
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

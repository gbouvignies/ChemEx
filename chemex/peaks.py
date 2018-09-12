"""The peaks module contains the code for handling peak assignments and
resonances."""

import collections
import functools
import re

re_peak_name = re.compile(
    """
        (^\s*|\-)
        (                                # group name
            (?P<symbol>\D?)              # one letter amino acid (optional)
            0*(?P<number>[0-9]+)         # residue number
            (?P<suffix>[abd-gi-mopr-z]*) # suffix (optional)
        )?
        (?P<nucleus>                     # nucleus name (e.g., CA, HG, ...)
            (?P<atom>[hncq])             # nucleus type
            [a-z0-9]*                    # nucleus name - nucleus type
        )?
    """,
    re.IGNORECASE | re.VERBOSE,
)

spin_names = ("i", "s", "x")


@functools.total_ordering
class Peak(object):
    """Peak class."""

    def __init__(self, assignment=None):

        if assignment is None:
            assignment = ""

        self.assignment = assignment

    @property
    def assignment(self):
        return self._assignment

    @assignment.setter
    def assignment(self, value):
        resonance_list = get_resonances(value.upper())
        self._resonances = collections.OrderedDict(
            (spin, res) for spin, res in zip("isx", resonance_list)
        )
        self._assignment = get_assignment(self._resonances.values())

    @property
    def names(self):
        return collections.OrderedDict(
            (
                ("i", self._resonances.get("i", {}).get("name", "")),
                ("s", self._resonances.get("s", {}).get("name", "")),
                (
                    "is",
                    get_assignment(
                        self._resonances[_] for _ in "is" if _ in self._resonances
                    ),
                ),
            )
        )

    @property
    def symbols(self):
        return collections.OrderedDict(
            (
                ("i", self._resonances.get("i", {}).get("symbol", "")),
                ("s", self._resonances.get("s", {}).get("symbol", "")),
            )
        )

    @property
    def atoms(self):
        return {key: res.get("atom", "") for key, res in self._resonances.items()}

    @property
    def nuclei(self):
        return {key: res.get("nucleus", "") for key, res in self._resonances.items()}

    def __repr__(self):
        return self.assignment

    def __eq__(self, other):
        return self.assignment == other.assignment

    def __hash__(self):
        return hash(self.assignment)

    def __lt__(self, other):
        if isinstance(other, str):
            other = Peak(other)

        self_tuple = tuple(
            (resonance["nucleus"], int(resonance["number"]))
            for resonance in self._resonances.values()
            if resonance["number"]
        )

        other_tuple = tuple(
            (resonance["nucleus"], int(resonance["number"]))
            for resonance in other._resonances.values()
            if resonance["number"]
        )

        return self_tuple < other_tuple

    def __len__(self):
        return len(
            [resonance for resonance in self._resonances.values() if resonance["name"]]
        )

    def intersection(self, other):
        """TODO: method docstring."""

        assignments = {self.assignment, other.assignment}

        names = {resonance["name"] for resonance in self._resonances.values()}
        names.update({resonance["name"] for resonance in other._resonances.values()})

        groups = {resonance["group"] for resonance in self._resonances.values()}
        groups.update({resonance["group"] for resonance in other._resonances.values()})

        if len(assignments) == 1:
            assignment_new = assignments.pop()
        elif len(names) == 1:
            assignment_new = names.pop()
        elif len(groups) == 1:
            assignment_new = groups.pop()
        else:
            assignment_new = ""

        return Peak(assignment_new)


def get_resonances(assignment):
    """Get resonances from an assignment."""

    resonances = []
    last_resonance = None

    for match in re.finditer(re_peak_name, assignment.lower()):

        resonance = match.groupdict()

        for key, value in resonance.items():
            if value is None:
                if last_resonance is not None:
                    resonance[key] = last_resonance[key]
                else:
                    resonance[key] = ""

        resonance["group"] = "".join(
            resonance.get(_, "") for _ in ("symbol", "number", "suffix")
        )
        resonance["name"] = "".join(resonance.get(_, "") for _ in ("group", "nucleus"))

        resonances.append(resonance)

        last_resonance = resonance

    return resonances


def get_assignment(resonances):
    """Get assignment from resonances."""
    parts = []
    last_group = None

    for resonance in resonances:
        if resonance["name"]:
            part = []
            if resonance["group"] is not None and resonance["group"] != last_group:
                part.append(resonance["group"])
            part.append(resonance["nucleus"])
            parts.append("".join(part))
            last_group = resonance["group"]

    return "-".join(parts)

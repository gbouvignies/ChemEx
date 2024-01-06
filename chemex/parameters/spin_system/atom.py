"""The atom module defines the Atom class, used to identify atoms."""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass, field

from .constants import CORRECT_ATOM_NAME
from .nucleus import Nucleus, str2nucleus


@dataclass(order=True)
class Atom:
    """Represents an atom in a molecular structure.

    This class encapsulates details about an atom, including its name,
    associated nucleus type, and search keys for matching and identification.

    Attributes:
        name (str): The name of the atom.
        nucleus (Nucleus): The type of nucleus associated with the atom, determined
                           automatically.
        search_keys (set[Hashable]): Set of keys used for searching or matching the
                                     atom.
    """

    name: str
    nucleus: Nucleus = field(init=False, compare=False)
    search_keys: set[Hashable] = field(init=False, default_factory=set, compare=False)

    def __post_init__(self) -> None:
        """Initializes the Atom instance, correcting the atom name and determining the nucleus type."""
        name = self.name.strip().upper()
        self.name = CORRECT_ATOM_NAME.get(name, name)
        self.nucleus = str2nucleus(self.name[:1])
        self.search_keys.update({self, self.nucleus})

    def match(self, other: Atom) -> bool:
        """Checks if this atom matches another atom based on the name.

        Args:
            other (Atom): The other atom to compare with.

        Returns:
            bool: True if the other atom's name starts with the name of this atom, False otherwise.
        """
        return other.name.startswith(self.name)

    def __hash__(self) -> int:
        """Generates a hash value for an Atom instance.

        Returns:
            int: The hash value of the atom.
        """
        return hash(self.name)

    def __bool__(self) -> bool:
        """Determines the boolean value of the Atom instance.

        Returns:
            bool: True if the atom has a valid name, False otherwise.
        """
        return bool(self.name)

    def __str__(self) -> str:
        """Returns the string representation of the Atom instance.

        Returns:
            str: The name of the atom.
        """
        return self.name

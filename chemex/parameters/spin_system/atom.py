"""The atom module defines the Atom class, used to identify atoms."""

from __future__ import annotations

from collections.abc import Hashable
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Self

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

    def __deepcopy__(self, memo: dict[int, Self]) -> Self:
        """Creates a deep copy of the Atom instance.

        Args:
            memo (dict[int, Self]): A dictionary of memoized objects.

        Returns:
            Self: A deep copy of the Atom instance.
        """
        if id(self) in memo:
            return memo[id(self)]

        # Create a new instance of Atom without calling __init__
        cls = self.__class__
        new_atom = cls.__new__(cls)

        # Copy all attributes to the new instance
        new_atom.name = deepcopy(self.name, memo)
        new_atom.nucleus = deepcopy(self.nucleus, memo)

        # Copy search_keys excluding self
        new_search_keys = deepcopy(
            {key for key in self.search_keys if key is not self}, memo
        )
        new_atom.search_keys = new_search_keys

        # Add the new_atom to its own search_keys set
        new_atom.search_keys.add(new_atom)

        # Memoize the new instance
        memo[id(self)] = new_atom

        return new_atom

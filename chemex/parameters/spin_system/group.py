from __future__ import annotations

from collections.abc import Hashable
from copy import deepcopy
from functools import cache, total_ordering
from re import search
from typing import Self

from .constants import AAA_TO_A

NO_NUMBER = -100000000


@cache
def parse_group(name: str) -> tuple[str, int, str]:
    """Parses the group name into symbol, number, and suffix.

    Args:
        name (str): The name of the group to be parsed.

    Returns:
        tuple[str, int, str]: A tuple containing the symbol, number, and suffix.
    """
    if found := search("[0-9]+", name.strip().upper()):
        symbol = name[: found.start()]
        corrected_symbol = AAA_TO_A.get(symbol, symbol)
        number = int(found.group())
        suffix = name[found.end() :]
        return corrected_symbol, number, suffix
    return name, NO_NUMBER, ""


@total_ordering
class Group:
    """Represents a molecular group in a spin system.

    This class encapsulates the information about an amino acid residue or nucleic
    acid base, including its symbol, sequence number, and suffix for identifying
    different states. It provides methods for parsing group names and comparing group
    objects.

    Attributes:
        symbol (str): The symbol of the group.
        number (int): The number associated with the group, if any.
        suffix (str): Any suffix associated with the group name.
        search_keys (set[Hashable]): A set of keys used for searching or matching the
                                     group.
    """

    def __init__(self, name: str) -> None:
        """Initializes a Group object with a given name.

        The constructor parses the input name into symbol, number, and suffix
        components.

        Args:
            name (str): The name of the group to be parsed.
        """
        self.symbol, self.number, self.suffix = parse_group(name.strip().upper())
        self.search_keys: set[Hashable] = {self} if self else set()

    @property
    def name(self) -> str:
        """Represents the full name of the group.

        Returns:
            str: The full name of the group, combining symbol, number, and suffix.
        """
        number = "" if self.number == NO_NUMBER else self.number
        return f"{self.symbol}{number}{self.suffix}"

    def match(self, other: Group) -> bool:
        """Determines if this group matches another group.

        Args:
            other (Group): Another group to compare with.

        Returns:
            bool: True if the groups match, False otherwise.
        """
        symbol = other.symbol == self.symbol or not self.symbol
        number = self.number in (other.number, NO_NUMBER)
        suffix = other.suffix == self.suffix or not self.suffix
        return number and symbol and suffix

    def __lt__(self, other: object) -> bool:
        """Less than comparison for ordering of Group objects.

        Args:
            other (object): Another object to compare with.

        Returns:
            bool: True if this group is considered less than the other, False otherwise.
        """
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.number < other.number

    def __eq__(self, other: object) -> bool:
        """Equality comparison for Group objects.

        Args:
            other (object): Another object to compare with.

        Returns:
            bool: True if the groups are equal, False otherwise.
        """
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        """Generates a hash value for a Group object.

        Returns:
            int: The hash value of the group.
        """
        return hash((self.symbol, self.number, self.suffix))

    def __bool__(self) -> bool:
        """Boolean representation of the Group object.

        Returns:
            bool: True if the group has a valid name, False otherwise.
        """
        return bool(self.name)

    def __str__(self) -> str:
        """String representation of the Group object.

        Returns:
            str: The full name of the group.
        """
        return self.name

    def __deepcopy__(self, memo: dict[int, Self]) -> Self:
        """Deep copy of the Group object.

        Args:
            memo (dict[int, object]): A dictionary for tracking copied objects.

        Returns:
            Group: A deep copy of the group.
        """
        if id(self) in memo:
            return memo[id(self)]

        # Create a new instance of Group without calling __init__
        cls = self.__class__
        new_group = cls.__new__(cls)

        # Deep copy all attributes to the new instance
        new_group.symbol = deepcopy(self.symbol, memo)
        new_group.number = deepcopy(self.number, memo)
        new_group.suffix = deepcopy(self.suffix, memo)

        # Copy search_keys excluding self
        new_search_keys = deepcopy(
            {key for key in self.search_keys if key is not self}, memo
        )
        new_group.search_keys = new_search_keys

        # Add the new_group to its own search_keys set
        new_group.search_keys.add(new_group)

        # Memoize the new instance
        memo[id(self)] = new_group

        return new_group

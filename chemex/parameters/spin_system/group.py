from __future__ import annotations

from functools import total_ordering
from re import search
from typing import TYPE_CHECKING

from .constants import AAA_TO_A

if TYPE_CHECKING:
    from collections.abc import Hashable


NO_NUMBER = -100000000


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
        self.symbol, self.number, self.suffix = self.parse_group(name.strip().upper())
        self.symbol = AAA_TO_A.get(self.symbol, self.symbol)
        self.search_keys: set[Hashable] = {self} if self else set()

    def parse_group(self, name: str) -> tuple[str, int, str]:
        """Parses the group name into symbol, number, and suffix.

        Args:
            name (str): The name of the group to be parsed.

        Returns:
            tuple[str, int, str]: A tuple containing the symbol, number, and suffix.
        """
        if found := search("[0-9]+", name.strip().upper()):
            return name[: found.start()], int(found.group()), name[found.end() :]
        return name, NO_NUMBER, ""

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
        return hash(self.name)

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

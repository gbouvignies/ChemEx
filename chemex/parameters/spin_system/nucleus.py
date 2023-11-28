from enum import Enum, auto


class Nucleus(Enum):
    """Enumeration of different types of atomic nuclei.

    This Enum class represents different types of atomic nuclei commonly
    used in NMR spectroscopy.

    Attributes:
        H1: Represents Hydrogen-1 nucleus.
        N15: Represents Nitrogen-15 nucleus.
        C13: Represents Carbon-13 nucleus.
        X: Placeholder for unknown or unspecified nucleus type.
    """

    H1 = auto()
    N15 = auto()
    C13 = auto()
    X = auto()


# Dictionary mapping atom letter abbreviations to corresponding Nucleus enum values.
STR_TO_NUCLEUS: dict[str, Nucleus] = {
    "H": Nucleus.H1,
    "Q": Nucleus.H1,  # Q, M are aliases for Hydrogen-1
    "M": Nucleus.H1,
    "N": Nucleus.N15,
    "C": Nucleus.C13,
    "X": Nucleus.X,
}


def str2nucleus(atom_letter: str) -> Nucleus:
    """Converts a single-letter atom name to the corresponding Nucleus enum.

    This function maps a single-letter abbreviation for an atom (e.g., 'H', 'C')
    to the corresponding Nucleus enumeration member. If the atom letter is not
    recognized, it defaults to 'X', representing an unknown or unspecified nucleus.

    Args:
        atom_letter (str): A single-letter abbreviation of an atom.

    Returns:
        Nucleus: The corresponding Nucleus enumeration member.

    Examples:
        >>> str2nucleus('H')
        <Nucleus.H1: 1>
        >>> str2nucleus('C')
        <Nucleus.C13: 3>
        >>> str2nucleus('Z')  # Not recognized
        <Nucleus.X: 4>
    """
    return STR_TO_NUCLEUS.get(atom_letter, Nucleus.X)

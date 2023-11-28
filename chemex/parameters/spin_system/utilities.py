from collections.abc import Iterable
from itertools import chain, combinations
from typing import TypeVar

T = TypeVar("T")


def powerset(iterable: Iterable[T]) -> Iterable[tuple[T, ...]]:
    """Generates all possible non-empty subsets of the given iterable.

    This function returns the power set of the input iterable, excluding the empty set.
    The subsets are returned as tuples. The order of subsets follows the order of
    elements in the input iterable.

    Args:
        iterable (Iterable[Any]): An iterable for which to generate the power set.

    Returns:
        Iterable[Any]: An iterable of tuples, each representing a non-empty subset of
        the input iterable.

    Examples:
        >>> powerset([1, 2, 3])
        (1,) (2,) (3,) (1, 2) (1, 3) (2, 3) (1, 2, 3)

    Note:
        This implementation uses `itertools.chain` and `itertools.combinations`, and is
        efficient for small to medium-sized iterables.

    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))

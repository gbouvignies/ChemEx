"""Utility module for list and string manipulation in Pydantic models.

This module provides utility functions for ensuring a variable is in list format and
for converting strings to lowercase. It also includes a Pydantic BaseModel subclass
which applies these utilities to convert all keys in a model to lowercase.

Typical usage example:

  variable = "Sample"
  variable_list = ensure_list(variable)
  lower_variable = to_lower("HELLO")
  model_instance = BaseModelLowerCase.parse_obj({"Name": "Alice"})
"""

from typing import Any, TypeVar, overload

T = TypeVar("T")


def ensure_list(variable: T | list[T] | None) -> list[T]:
    """Ensures that the input variable is returned as a list.

    If the input variable is already a list, it is returned as-is.
    If the variable is None, an empty list is returned.
    Otherwise, the variable is wrapped in a list and returned.

    Parameters
    ----------
    variable (T | list[T] | None): The variable to be ensured as a list.

    Returns
    -------
    list[T]: The input variable as a list.

    """
    if isinstance(variable, list):
        return variable
    if variable is None:
        return []
    return [variable]


@overload
def to_lower(string: str) -> str: ...


@overload
def to_lower(string: T) -> T: ...


def to_lower(string: object) -> object:
    """Converts a string to lowercase if it is of type str.

    If the input is not a string, it is returned unchanged.

    Parameters
    ----------
    string (T): The string to be converted to lowercase.

    Returns
    -------
    T: The lowercase string if input is a string; otherwise, the original input.

    """
    if isinstance(string, str):
        return string.lower()
    return string


def key_to_lower(model: dict[str, Any]) -> dict[str, Any]:
    """Converts keys of the model dictionary to lowercase."""
    return {k.lower(): v for k, v in model.items()}

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
from __future__ import annotations

from typing import TypeVar, cast

from pydantic import BaseModel, model_validator

T = TypeVar("T")


def ensure_list(variable: T | list[T] | None) -> T | list[T]:
    """Ensures that the input variable is returned as a list.

    If the input variable is already a list, it is returned as-is.
    If the variable is None, an empty list is returned.
    Otherwise, the variable is wrapped in a list and returned.

    Parameters:
    variable (T | list[T] | None): The variable to be ensured as a list.

    Returns:
    T | list[T]: The input variable as a list.
    """
    if isinstance(variable, list):
        return cast(list[T], variable)
    if variable is None:
        return []
    return [variable]


def to_lower(string: T) -> T:
    """Converts a string to lowercase if it is of type str.

    If the input is not a string, it is returned unchanged.

    Parameters:
    string (T): The string to be converted to lowercase.

    Returns:
    T: The lowercase string if input is a string; otherwise, the original input.
    """
    if isinstance(string, str):
        return string.lower()
    return string


class BaseModelLowerCase(BaseModel):
    """A Pydantic BaseModel class that converts all keys in the model to lowercase.

    This is achieved using a model validator that operates before model initialization.
    The validator applies the `to_lower` function to each key of the model.
    """

    @model_validator(mode="before")
    def key_to_lower(cls, model: dict[str, T]) -> dict[str, T]:
        """Model validator to convert all dictionary keys to lowercase.

        Parameters:
        model (dict[str, T]): The dictionary model with string keys.

        Returns:
        dict[str, T]: The modified dictionary with all keys in lowercase.
        """
        return {to_lower(k): v for k, v in model.items()}

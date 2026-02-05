"""Common Pydantic type aliases with constraints for ChemEx configuration.

This module defines reusable constrained types that ensure consistency across
the codebase and reduce boilerplate in field definitions.
"""

from __future__ import annotations

from typing import Annotated

from pydantic import AfterValidator, Field

# Constants
MAX_TEMPERATURE = 500.0  # Maximum physically reasonable temperature in Kelvin


# Validators for common patterns
def validate_positive(v: float) -> float:
    """Validate that a value is positive (> 0)."""
    if v <= 0:
        msg = f"Must be positive, got {v}"
        raise ValueError(msg)
    return v


def validate_non_negative(v: float) -> float:
    """Validate that a value is non-negative (>= 0)."""
    if v < 0:
        msg = f"Must be non-negative, got {v}"
        raise ValueError(msg)
    return v


def validate_probability(v: float) -> float:
    """Validate that a value is a valid probability [0, 1]."""
    if not 0.0 <= v <= 1.0:
        msg = f"Must be between 0 and 1, got {v}"
        raise ValueError(msg)
    return v


def validate_temperature(v: float) -> float:
    """Validate that a temperature is physically reasonable."""
    if not 0.0 < v < MAX_TEMPERATURE:
        msg = f"Temperature must be between 0 and {MAX_TEMPERATURE} K, got {v}"
        raise ValueError(msg)
    return v


# Basic constrained types
PositiveFloat = Annotated[
    float,
    Field(gt=0.0),
    AfterValidator(validate_positive),
]
"""A positive float value (> 0)."""

NonNegativeFloat = Annotated[
    float,
    Field(ge=0.0),
    AfterValidator(validate_non_negative),
]
"""A non-negative float value (>= 0)."""

PositiveInt = Annotated[
    int,
    Field(gt=0),
]
"""A positive integer value (> 0)."""

NonNegativeInt = Annotated[
    int,
    Field(ge=0),
]
"""A non-negative integer value (>= 0)."""

Probability = Annotated[
    float,
    Field(ge=0.0, le=1.0),
    AfterValidator(validate_probability),
]
"""A probability value between 0 and 1 (inclusive)."""


# Physics and NMR specific types
PulseWidth = Annotated[
    float,
    Field(
        gt=0.0,
        description="Pulse width in seconds (must be positive)",
    ),
]
"""Pulse width in seconds (> 0)."""

Frequency = Annotated[
    float,
    Field(
        description="Frequency in Hz",
    ),
]
"""Frequency in Hz (can be negative for offset frequencies)."""

PositiveFrequency = Annotated[
    float,
    Field(
        gt=0.0,
        description="Frequency in Hz (must be positive)",
    ),
]
"""Positive frequency in Hz (> 0)."""

Temperature = Annotated[
    float,
    Field(
        gt=0.0,
        lt=MAX_TEMPERATURE,
        description=f"Temperature in Kelvin (0-{MAX_TEMPERATURE} K)",
    ),
    AfterValidator(validate_temperature),
]
f"""Temperature in Kelvin, physically reasonable range (0-{MAX_TEMPERATURE} K)."""

RelaxationTime = Annotated[
    float,
    Field(
        gt=0.0,
        description="Relaxation time in seconds (must be positive)",
    ),
]
"""Relaxation time in seconds (> 0)."""

Delay = Annotated[
    float,
    Field(
        ge=0.0,
        description="Time delay in seconds (must be non-negative)",
    ),
]
"""Time delay in seconds (>= 0)."""

B1Field = Annotated[
    float,
    Field(
        ge=0.0,
        description="B1 field strength in Hz (must be non-negative)",
    ),
]
"""B1 field strength in Hz (>= 0). Zero indicates no field applied."""

ChemicalShift = Annotated[
    float,
    Field(
        description="Chemical shift in ppm",
    ),
]
"""Chemical shift in ppm (can be any float value)."""

ExchangeRate = Annotated[
    float,
    Field(
        ge=0.0,
        description="Exchange rate constant in s⁻¹ (must be non-negative)",
    ),
]
"""Exchange rate constant in s⁻¹ (>= 0)."""

Population = Annotated[
    float,
    Field(
        ge=0.0,
        le=1.0,
        description="Population fraction (0-1)",
    ),
    AfterValidator(validate_probability),
]
"""Population fraction between 0 and 1."""


# Optional versions of common types
OptionalPulseWidth = PulseWidth | None
OptionalFrequency = Frequency | None
OptionalPositiveFrequency = PositiveFrequency | None
OptionalTemperature = Temperature | None
OptionalDelay = Delay | None
OptionalB1Field = B1Field | None

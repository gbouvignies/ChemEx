"""Canonical temperature-polynomial semantics for ``.tc`` chemical shifts."""

from __future__ import annotations

import math
import sys

from chemex.parameters.setting import NameSetting, ParamLocalSetting

TREF_DEFAULT = 25.0
TREF_MINIMUM = math.nextafter(-273.15, math.inf)


def _reference_temperature_setting() -> ParamLocalSetting:
    return ParamLocalSetting(
        name_setting=NameSetting(
            "tref",
            "",
            allow_residue_specific=False,
        ),
        value=TREF_DEFAULT,
        min=TREF_MINIMUM,
        max=sys.float_info.max,
        model_owned=True,
    )


def add_reference_shift_polynomial(
    settings: dict[str, ParamLocalSetting],
    *,
    temperature: float,
) -> None:
    """Replace state-A shifts with a centered first-order polynomial."""
    settings.setdefault("tref", _reference_temperature_setting())
    for nucleus in ("i", "s"):
        settings[f"cs0_{nucleus}_a"] = ParamLocalSetting(
            name_setting=NameSetting("cs0_a", nucleus),
            value=0.0,
            min=-100.0,
            max=300.0,
            supports_estimation=True,
        )
        settings[f"cs1_{nucleus}_a"] = ParamLocalSetting(
            name_setting=NameSetting("cs1_a", nucleus),
            value=0.0,
            min=-1.0,
            max=1.0,
            supports_estimation=True,
        )
        shift = settings[f"cs_{nucleus}_a"]
        shift.vary = False
        shift.expr = (
            f"{{cs0_{nucleus}_a}} + ({temperature} - {{tref}}) * {{cs1_{nucleus}_a}}"
        )
        shift.supports_estimation = False
        shift.model_owned = True


def add_shift_difference_polynomial(
    settings: dict[str, ParamLocalSetting],
    *,
    state: str,
    temperature: float,
) -> None:
    """Replace each A-to-X shift difference with the same centered basis."""
    settings.setdefault("tref", _reference_temperature_setting())
    for nucleus in ("i", "s"):
        settings[f"dw0_{nucleus}_a{state}"] = ParamLocalSetting(
            name_setting=NameSetting(f"dw0_a{state}", nucleus),
            value=0.0,
            min=-100.0,
            max=100.0,
            vary=True,
            supports_estimation=True,
        )
        settings[f"dw1_{nucleus}_a{state}"] = ParamLocalSetting(
            name_setting=NameSetting(f"dw1_a{state}", nucleus),
            value=0.0,
            min=-1.0,
            max=1.0,
            vary=True,
            supports_estimation=True,
        )
        difference = settings[f"dw_{nucleus}_a{state}"]
        difference.vary = False
        difference.expr = (
            f"{{dw0_{nucleus}_a{state}}} + ({temperature} - {{tref}})"
            f" * {{dw1_{nucleus}_a{state}}}"
        )
        difference.supports_estimation = False
        difference.model_owned = True


def canonical_control_guidance(name: str) -> str | None:
    """Describe the independent coordinates controlling a derived shift."""
    upper = name.upper()
    if upper == "TREF":
        return "set TREF only as a scalar in [GLOBAL]"
    if upper == "CS_A":
        return "use CS0_A and CS1_A"
    if upper.startswith("DW_A"):
        state_pair = upper.split("_", 1)[1]
        return f"use DW0_{state_pair} and DW1_{state_pair}"
    if upper.startswith("CS_") and len(upper) == 4 and upper[-1] != "A":
        state = upper[-1]
        return (
            f"use CS0_A, CS1_A, DW0_A{state}, and DW1_A{state}; "
            f"CS_{state} is composed from the A-state reference"
        )
    return None

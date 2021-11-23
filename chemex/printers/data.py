from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from chemex.containers.data import Data


class Printer(Protocol):
    header: str = ""

    def print(self, name: str, data: Data) -> str:
        ...


@dataclass
class ShiftPrinter:
    header = "\n".join(
        [
            "#",
            "# Unit: PPB",
            "#",
            f"# {'NAME':>16s} {'SHIFT (EXP)':>17s} {'ERROR (EXP)':>17s} {'SHIFT (CALC)':>17s}\n",
        ]
    )

    def print(self, name: str, data: Data) -> str:
        return f"{name:>18s} {data.exp[0]: 17.3f} {data.err[0]: 17.3f} {data.calc[0]: 17.3f}\n"


@dataclass
class ProfilePrinter:
    first_column_name: str = ""
    first_column_fmt: str = "12.2f"
    header = ""

    def print(self, name: str, data: Data) -> str:
        output = [
            f"[{name}]",
            f"# {self.first_column_name:>12s} {'INTENSITY (EXP)':>17s} {'ERROR (EXP)':>17s} {'INTENSITY (CALC)':>17s}",
        ]

        for metadata, exp, err, calc, mask in zip(
            data.metadata, data.exp, data.err, data.calc, data.mask
        ):
            start = " " if mask else "#"
            end = "" if mask else " # NOT USED IN THE FIT"
            output.append(
                f"{start} {metadata: {self.first_column_fmt}} {exp: 17.8e} {err: 17.8e} {calc: 17.8e}{end}"
            )
        return "\n".join(output) + "\n\n"


@dataclass
class CestPrinter(ProfilePrinter):
    first_column_name: str = "OFFSET (HZ)"
    first_column_fmt: str = "12.2f"


@dataclass
class CpmgPrinter(ProfilePrinter):
    first_column_name: str = "NCYC"
    first_column_fmt: str = "12.0f"


@dataclass
class RelaxationPrinter(ProfilePrinter):
    first_column_name: str = "TIME (S)"
    first_column_fmt: str = "12.2f"

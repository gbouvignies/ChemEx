from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from chemex.containers.data import Data


class Printer(Protocol):
    header: str
    simulation: bool

    def print(self, name: str, data: Data) -> str:
        ...


@dataclass
class ShiftPrinter:
    header = "\n".join(
        [
            "#",
            "# Unit: PPB",
            "#",
            (
                f"# {'NAME':>16s} {'SHIFT (EXP)':>17s}"
                f" {'ERROR (EXP)':>17s} {'SHIFT (CALC)':>17s}\n"
            ),
        ],
    )

    def print(self, name: str, data: Data) -> str:
        return (
            f"{name:>18s} {data.exp[0]: 17.3f} {data.err[0]: 17.3f}"
            f" {data.calc[0]: 17.3f}\n"
        )


@dataclass
class ProfilePrinter:
    first_column_name: str = ""
    first_column_fmt: str = "12.2f"
    header: str = ""
    simulation: bool = False

    def print(self, name: str, data: Data) -> str:
        output: list[list[str]] = []

        output.append([f"[{name}]"])

        header: list[str] = [
            "#",
            f"{self.first_column_name:>12s}",
            f"{'INTENSITY (EXP)':>17s}",
            f"{'ERROR (EXP)':>17s}",
            f"{'INTENSITY (CALC)':>17s}",
        ]
        if self.simulation:
            del header[2:4]
        output.append(header)

        for metadata, exp, err, calc, mask in zip(
            data.metadata,
            data.exp,
            data.err,
            data.calc,
            data.mask,
            strict=True,
        ):
            start = " " if mask else "#"
            end = "" if mask else " # NOT USED IN THE FIT"
            line: list[str] = []
            line.extend(
                [
                    f"{start}",
                    f"{metadata:{self.first_column_fmt}}",
                    f"{exp: 17.8e}",
                    f"{err: 17.8e}",
                    f"{calc: 17.8e}",
                    f"{end}",
                ]
            )
            if self.simulation:
                del line[2:4]
            output.append(line)

        return "\n".join(" ".join(line) for line in output) + "\n\n"


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


@dataclass
class EXSYPrinter(ProfilePrinter):
    first_column_name: str = "TIME (S)"
    first_column_fmt: str = "12.6f"
    header = ""
    simulation: bool = False

    def print(self, name: str, data: Data) -> str:
        output: list[list[str]] = []

        output.append([f"[{name}]"])

        header: list[str] = [
            "#",
            f"{self.first_column_name:>12s}",
            f"{'STATE1':>12s}",
            f"{'STATE2':>12s}",
            f"{'INTENSITY (EXP)':>17s}",
            f"{'ERROR (EXP)':>17s}",
            f"{'INTENSITY (CALC)':>17s}",
        ]
        if self.simulation:
            del header[4:6]
        output.append(header)

        for metadata, exp, err, calc, mask in zip(
            data.metadata,
            data.exp,
            data.err,
            data.calc,
            data.mask,
            strict=True,
        ):
            start = " " if mask else "#"
            end = "" if mask else " # NOT USED IN THE FIT"
            time = metadata["times"]
            state1 = metadata["states1"]
            state2 = metadata["states1"]
            line: list[str] = []
            line.extend(
                [
                    f"{start}",
                    f"{time: {self.first_column_fmt}}",
                    f"{state1:>12s}",
                    f"{state2:>12s}",
                    f"{exp: 17.8e}",
                    f"{err: 17.8e}",
                    f"{calc: 17.8e}",
                    f"{end}",
                ]
            )
            if self.simulation:
                del line[4:6]
            output.append(line)

        return "\n".join(" ".join(line) for line in output) + "\n\n"

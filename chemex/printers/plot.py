from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from chemex.containers.data import Data


class PlotPrinter(Protocol):
    def print_exp(self, name: str, data: Data) -> str:
        ...

    def print_calc(self, name: str, data: Data) -> str:
        ...


@dataclass
class DataPlotPrinter:
    first_column_name: str = ""
    first_column_fmt: str = "12.2f"

    def print_exp(self, name: str, data: Data) -> str:
        output = [
            f"[{name}]",
            f"# {self.first_column_name:>12s} {'INTENSITY (EXP)':>17s} {'ERROR (EXP)':>17s}",
        ]

        for metadata, exp, err, mask in zip(
            data.metadata, data.exp, data.err, data.mask
        ):
            start = " " if mask else "#"
            end = "" if mask else " # NOT USED IN THE FIT"
            output.append(
                f"{start} {metadata: {self.first_column_fmt}}  {exp: 17.8e} {err[0]: 17.8e}{end}"
            )
        return "\n".join(output) + "\n\n"

    def print_calc(self, name: str, data: Data) -> str:
        output = [
            f"[{name}]",
            f"# {self.first_column_name:>12s} {'INTENSITY (CALC)':>17s}",
        ]

        for metadata, calc in zip(data.metadata, data.calc):
            output.append(f"  {metadata: {self.first_column_fmt}}  {calc: 17.8e}")

        return "\n".join(output) + "\n\n"


@dataclass
class CpmgDataPlotPrinter(DataPlotPrinter):
    first_column_name: str = "NU_CPMG (Hz)"
    first_column_fmt: str = "12.2f"

    def print_exp(self, name: str, data: Data) -> str:
        output = [
            f"[{name}]",
            f"# {self.first_column_name:>12s} {'R2 (EXP)':>17s} {'ERROR DOWN (EXP)':>17s} {'ERROR UP (EXP)':>17s}",
        ]

        for metadata, exp, err, mask in zip(
            data.metadata, data.exp, data.err, data.mask
        ):
            start = " " if mask else "#"
            end = "" if mask else " # NOT USED IN THE FIT"
            output.append(
                f"{start} {metadata: {self.first_column_fmt}}  {exp: 17.8e} {err[0]: 17.8e} {err[1]: 17.8e}{end}"
            )
        return "\n".join(output) + "\n\n"

    def print_calc(self, name: str, data: Data) -> str:
        output = [f"[{name}]", f"# {self.first_column_name:>12s} {'R2 (CALC)':>17s}"]

        for metadata, calc in zip(data.metadata, data.calc):
            output.append(f"  {metadata: {self.first_column_fmt}}  {calc: 17.8e}")

        return "\n".join(output) + "\n\n"


data_plot_printers: dict[str, PlotPrinter] = {
    "cest": DataPlotPrinter("OFFSET (PPM)", "12.2f"),
    "cpmg": CpmgDataPlotPrinter(),
    "relaxation": DataPlotPrinter("TIME (S)", "12.2f"),
}

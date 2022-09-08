from __future__ import annotations

from collections import Counter
from inspect import cleandoc
from pathlib import Path

from pydantic import ValidationError
from rich import box
from rich.console import Console
from rich.padding import Padding
from rich.panel import Panel
from rich.rule import Rule
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

from chemex import __version__

console = Console()


LOGO = "\n".join(
    [
        r"   ________                   ______",
        r"  / ____/ /_  ___  ____ ___  / ____/  __",
        r" / /   / __ \/ _ \/ __ `__ \/ __/ | |/_/",
        r"/ /___/ / / /  __/ / / / / / /____>  <",
        r"\____/_/ /_/\___/_/ /_/ /_/_____/_/|_|",
        "",
        "",
    ]
)

EXPERIMENT_NAME = """
[experiment]
name = "experiment_name"
...
"""


def print_logo() -> None:
    logo = Text(LOGO, style="blue")
    description = "Analysis of NMR chemical exchange data\n\n"
    version = "Version: "
    version_number = Text(f"{__version__}", style="red")
    all_text = Text.assemble(logo, description, version, version_number)
    panel = Panel.fit(all_text)
    console.print(panel)


def print_loading_experiments() -> None:
    console.print("\nLoading datasets...", style="bold yellow")


def get_reading_exp_text(filename: Path, type: str = "", total_nb: int = 0) -> Text:

    parts = (
        "  • Reading ",
        (f"{filename}", "green"),
        " (",
        (f"{type}", "blue"),
        ")",
    )

    if total_nb == 0:
        return Text.assemble(*parts, " ...")

    nb_text = Text(f"{total_nb}", style="bold")

    return Text.assemble(*parts, " -> ", nb_text, " profiles")


def print_reading_defaults() -> None:
    console.print("\nReading default parameters...", style="bold yellow")


def print_reading_methods() -> None:
    console.print("\nReading methods...", style="bold yellow")


def print_start_fit() -> None:
    console.print("\nStarting the fits...", style="bold yellow")


def print_running_simulations() -> None:
    console.print("\nRunning simulations...", style="bold yellow")


def print_step_name(name: str, index: int, total: int) -> None:
    text = Text.assemble(
        "Running ",
        (f"{name}", "magenta"),
        f" ({index}/{total})\n",
    )
    console.print(Padding(Rule(text, end="\n\n"), (1, 0, 0, 1)), width=51)


def print_selecting_profiles(selected_nb: int) -> None:
    console.print(
        Text.from_markup(f"  • Selecting profiles -> [blue]{selected_nb}[/] profiles")
    )


def print_no_data() -> None:
    console.print("  • No data to fit")


def print_fitmethod(fit_method: str) -> None:
    console.print(Text.assemble("  • Fit method -> ", (f"{fit_method}", "blue")))


def print_status_changes(
    vary: Counter[str], fixed: Counter[str], constrained: Counter[str]
) -> None:

    if not any([vary, fixed, constrained]):
        return

    console.print(Text("  • Setting parameter status..."))

    table = Table()
    table.add_column("Name")
    table.add_column("Status")
    table.add_column("Updated", justify="right")

    for name, nb in vary.items():
        table.add_row(name, "fitted", str(nb))

    for name, nb in fixed.items():
        table.add_row(name, "fixed", str(nb))

    for name, nb in constrained.items():
        table.add_row(name, "constrained", str(nb))

    table.box = box.SIMPLE_HEAD

    console.print(Padding.indent(table, 3))


def print_minimizing() -> None:
    console.print(Text("  • Running the minimization..."))


def print_running_grid() -> None:
    console.print(Text("  • Running the grid search..."))


def print_running_statistics(name: str) -> None:
    console.print(Text(f"  • Running {name} simulations..."))


def print_section(name: str) -> None:
    console.print(Text(f"  • Section {name}"))


def create_chi2_table() -> Table:
    table = Table()
    table.add_column("Iteration", justify="right")
    table.add_column("χ²", justify="right")
    table.add_column("Reduced χ²", justify="right")
    table.box = box.SIMPLE
    return table


def print_chi2(chisqr: float, redchi: float) -> None:
    console.print()
    console.print(
        Padding.indent(Text.from_markup(f"        χ²: [bold]{chisqr:15.1f}[/]"), 3)
    )
    console.print(
        Padding.indent(Text.from_markup(f"Reduced χ²: [bold]{redchi:15.3f}[/]"), 3)
    )
    console.print()


def print_writing_results(path: Path) -> None:
    console.print(f"  • Writing results in [green]{path}")


def print_making_plots() -> None:
    console.print(Text("  • Making plots..."))


def print_plot_filename(filename: Path, extra: bool = True) -> None:
    text = f"    ‣ [green]{filename}[/]"

    if extra:
        text += " [[green].fit[/], [green].exp[/]]"

    console.print(Text.from_markup(text))


def print_group_name(text: str) -> None:
    console.print(Padding(Rule(text, characters="⋅"), (1, 0, 0, 3)), width=49)


def print_file_not_found(filename: Path) -> None:
    console.print()
    console.print(f"[red]The file '{filename}' is empty or does not exist!")


def print_error(error: FileNotFoundError) -> None:
    console.print()
    console.print(f"[red]Error: {error}")


def print_toml_error(filename: Path, error_message: Exception) -> None:
    console.print()
    console.print(Text.from_markup(f"[red]Error in the TOML file '{filename}'"))
    console.print(
        Text.from_markup(f"[red]-> {error_message}\n"),
    )
    console.print(
        "Please check the website [link]https://toml.io[/] for a complete description of the TOML format."
    )


def print_experiment_name_error(filename: Path) -> None:
    console.print()
    console.print(f"[red]The experiment name is missing from the file '{filename}'\n")
    console.print(
        f"Please make sure that the 'name' field is provided in '{filename}':\n"
    )
    console.print(
        Syntax(EXPERIMENT_NAME, "toml"),
        "\nRun 'chemex info' to obtain the full list of the available experiments.",
    )


def print_pydantic_parsing_error(filename: Path, error: ValidationError) -> None:
    console.print()
    console.print(Text.from_markup(f"[red]Error(s) while parsing '{filename}'"))
    for err in error.errors():
        location = " -> ".join(str(loc) for loc in err["loc"])
        console.print("  • ", location, ":", err["msg"], style="red")


EXPERIMENT_NOT_FOUND = cleandoc(
    """[red bold]EXPERIMENT LOADING FAILED [/]

    '{name}' is not part of our experiment collection!
    Run 'chemex info' to obtain the full list of the available experiments.
"""
)

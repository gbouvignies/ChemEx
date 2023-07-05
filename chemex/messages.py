from __future__ import annotations

from typing import TYPE_CHECKING

from rich import box
from rich.console import Console
from rich.padding import Padding
from rich.panel import Panel
from rich.rule import Rule
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

from chemex import __version__

if TYPE_CHECKING:
    from collections import Counter
    from pathlib import Path

    from pydantic import ValidationError

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


def print_chi2_table_header() -> None:
    header = Text(f"{'Iteration':>9s}  {'χ²':>12s}  {'Reduced χ²':>12s}")
    console.print()
    console.print(Padding.indent(header, 5), style="bold")
    console.print(Padding.indent("─" * 39, 4))


def print_chi2_table_line(iteration: int, chisqr: float, redchi: float) -> None:
    line = Text(f"{iteration:>9d}  {chisqr:>12.1f}  {redchi:>12.3f}")
    console.print(Padding.indent(line, 5))


def print_chi2_table_footer(iteration: int, chisqr: float, redchi: float) -> None:
    footer = Text(f"{iteration:>9d}  {chisqr:>12.1f}  {redchi:>12.3f}")
    console.print(Padding.indent("─" * 39, 4))
    console.print(Padding.indent(footer, 5), style="bold")
    console.print()


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


def print_file_not_found_error(error: FileNotFoundError) -> None:
    console.print()
    console.print(f"[red]Error: {error}")


def print_toml_error(filename: Path, error_message: Exception) -> None:
    console.print()
    console.print(Text.from_markup(f"[red]Error in the TOML file '{filename}'"))
    console.print(
        Text.from_markup(f"[red]-> {error_message}\n"),
    )
    console.print(
        "Please check the website [link]https://toml.io[/]"
        " for a complete description of the TOML format."
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


def print_calculation_stopped_error() -> None:
    console.print()
    console.print("[red] -- Keyboard Interrupt: Calculation stopped --")
    console.print()


def print_plotting_canceled() -> None:
    console.print()
    console.print("[red] -- Keyboard Interrupt: Plotting cancelled --")
    console.print()


def print_value_error() -> None:
    console.print()
    console.print("[red] -- Got a ValueError: Calculation stopped --")
    console.print()


def print_warning_positive_jnh() -> None:
    console.print()
    console.print(
        "[yellow] -- WARNING: Some 1J(NH) couplings are set with positive values --"
    )
    console.print(
        "This can cause the TROSY and anti-TROSY components to be switched in some"
        " experiments."
    )
    console.print()


def print_warning_negative_jch() -> None:
    console.print()
    console.print(
        "[yellow] -- WARNING: Some 1J(CH) couplings are set with negative values --"
    )
    console.print(
        "This can cause the TROSY and anti-TROSY components to be switched in some"
        " experiments."
    )
    console.print()


def print_error_grid_settings(entry: str) -> None:
    console.print()
    console.print("[red] -- ERROR: Error reading grid settings:")
    console.print(Text(f'    "{entry}"', style="red"))
    console.print()
    console.print(
        "Please make sure that the grid settings are provided in the correct format"
    )
    console.print(Text("    [PB] = lin(<min>, <max>, <nb of points>)"))
    console.print(Text("    [PB] = log(<min>, <max>, <nb of points>)"))
    console.print(Text("    [PB] = (<value1>, <value2>, ..., <valuen>)"))
    console.print()


def print_error_constraints(expression: str) -> None:
    console.print()
    console.print(f'[red] -- ERROR: Error reading constraints -> "{expression}" --')
    console.print()


def print_grid_statistic_warning() -> None:
    console.print()
    console.print(
        "[yellow] -- WARNING: The 'GRID' and 'STATISTICS' options are mutually "
        "exclusive. Only the 'GRID' calculation will be run."
    )
    console.print()


def print_model_error(name: str) -> None:
    from chemex.models.factory import model_factory

    console.print()
    console.print(f"[red] -- ERROR: The model '{name}' is not available.")
    console.print()
    console.print("The available models are:")
    for model_name in sorted(model_factory.set):
        console.print(f"    - '{model_name}'")
    console.print()


def print_not_implemented_noise_method_warning(
    filename: Path, kind: str, implemented: tuple[str, ...]
) -> None:
    warning_message = (
        f"[yellow] -- WARNING: Experiment {filename.name}[/yellow]: "
        f"The '{kind}' method is not implemented. "
        f"Please choose one of the following methods: {', '.join(implemented)}. "
        f"Defaulting to 'file' method."
    )
    console.print()
    console.print(warning_message)
    console.print()


def print_no_duplicate_warning(filename: Path) -> None:
    warning_message = (
        f"[yellow] -- WARNING: Experiment {filename.name}[/yellow]: "
        f"Some profiles do not have duplicate points. "
        f"Uncertainties cannot be estimated and will be directly taken from the files. "
        f"Defaulting to 'file' method."
    )
    console.print()
    console.print(warning_message)
    console.print()


FITMETHOD_ERROR_MESSAGE = """\
  - "FITMETHOD" must be one of:
        "leastsq": Levenberg-Marquardt (default)
        "least_squares": least-squares minimization (Trust Region Reflective method)
        "differential_evolution": differential evolution
        "brute": brute force method
        "basinhopping": basinhopping
        "ampgo": Adaptive Memory Programming for Global Optimization
        "nelder": Nelder-Mead
        "lbfgsb": L-BFGS-B
        "powell": Powell
        "cg": Conjugate-Gradient
        "newton": Newton-CG
        "cobyla": Cobyla
        "bfgs": BFGS
        "tnc": Truncated Newton
        "trust-ncg": Newton-CG trust-region
        "trust-exact": nearly exact trust-region
        "trust-krylov": Newton GLTR trust-region
        "trust-constr": trust-region for constrained optimization
        "dogleg": Dog-leg trust-region
        "slsqp": Sequential Linear Squares Programming
        "emcee": Maximum likelihood via Monte-Carlo Markov Chain
        "shgo": Simplicial Homology Global Optimization
        "dual_annealing": Dual Annealing optimization"""

INCLUDE_ERROR_MESSAGE = """\
  - "INCLUDE" must be a list of nuclei to be included or '*' to select all the profiles.
    Example: ['G3H', 'S4H', 'W5H'] or [3, 4, 5]"""

EXCLUDE_ERROR_MESSAGE = """\
  - "EXCLUDE" must be a list of nuclei to be excluded.

    Example: ['G3H', 'S4H', 'W5H'] or [3, 4, 5]"""

FIT_ERROR_MESSAGE = """\
  - "FIT" must be a list of parameters to be fitted.

    Example: ["PB", "KEX_AB"]"""

FIX_ERROR_MESSAGE = """\
  - "FIX" must be a list of parameters to be fixed.

    Example: ["PB", "KEX_AB"]"""

CONSTRAINTS_ERROR_MESSAGE = """\
  - "CONSTRAINTS" must be a list of expressions.

    Example: ["[R2_B] = [R2_A] / 2", "[R1_B] = [R1_A]"]"""

GRID_ERROR_MESSAGE = """\
  - "GRID" must be a list of parameters to be gridded with the associated grid
    definition.

    Example: GRID    = [
        "[KEX_AB] = log(100.0, 600.0, 10)",
        "[PB] = log(0.03, 0.15, 10)",
        "[DW_AB] = lin(0.0, 10.0, 5)",
    ]"""

STATISTICS_ERROR_MESSAGE = """\
  - "STATISTICS" must be a dictionary with keys 'MC', 'BS', 'BSN'

    Example: { "MC"=10 } or { "MC"=10, "BS"=10 }"""

METHOD_ERROR_MESSAGES = {
    "fitmethod": FITMETHOD_ERROR_MESSAGE,
    "include": INCLUDE_ERROR_MESSAGE,
    "exclude": EXCLUDE_ERROR_MESSAGE,
    "fit": FIT_ERROR_MESSAGE,
    "fix": FIX_ERROR_MESSAGE,
    "constraints": CONSTRAINTS_ERROR_MESSAGE,
    "grid": GRID_ERROR_MESSAGE,
    "statistics": STATISTICS_ERROR_MESSAGE,
}


def print_method_error(filename: Path, section: str, options: set[int | str]) -> None:
    console.print()
    console.print(
        f"[red] -- ERROR: Error reading section '[{section}]' of '{filename}' --"
    )
    for option in options:
        console.print(METHOD_ERROR_MESSAGES.get(str(option)))
    console.print()

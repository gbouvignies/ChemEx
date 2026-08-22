"""Module for command-line interface functionalities in ChemEx software package.

This module is an integral part of the ChemEx software, dedicated to enhancing user
interaction with a command-line interface for analyzing NMR chemical exchange data.
It leverages the 'rich' Python library to provide rich text and table formatting,
making the data presentation more readable and engaging. The module includes various
functions for displaying progress, handling errors, and showing results of data
analysis steps such as dataset loading, fitting processes, simulations, and plotting.

Typical usage example:

  from your_module import print_logo, print_loading_experiments
  print_logo()
  print_loading_experiments()
"""

from collections import Counter
from collections.abc import Callable, Mapping
from contextlib import suppress
from dataclasses import replace
from pathlib import Path
from time import monotonic

from pydantic import ValidationError
from rich import box
from rich.console import Console
from rich.live import Live
from rich.padding import Padding
from rich.panel import Panel
from rich.rule import Rule
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

from chemex import __version__
from chemex.optimize.progress import (
    FitProgressContext,
    ProgressEvent,
    ProgressPhase,
    ProgressRateLimiter,
    ProgressUpdate,
)

console = Console()


class MinimizationProgressReporter:
    """Own Rich rendering and UX policy for one visible deterministic fit."""

    def __init__(
        self,
        output_console: Console,
        *,
        interactive: bool,
        retained_observation_count: int,
        controlled_parameter_count: int,
        step_name: str = "",
        group_labels: Mapping[frozenset[str], str] | None = None,
        grid: bool = False,
        enabled: bool = True,
        clock: Callable[[], float] = monotonic,
    ) -> None:
        self._console = output_console
        self._interactive = interactive
        self._retained_observation_count = retained_observation_count
        self._controlled_parameter_count = controlled_parameter_count
        self._step_name = step_name
        self._group_labels = {} if group_labels is None else dict(group_labels)
        self._grid = grid
        self._enabled = enabled
        self._clock = clock
        self._started_at: float | None = None
        self._live: Live | None = None
        self._limiters: dict[FitProgressContext, ProgressRateLimiter] = {}
        self._noninteractive_limiter = ProgressRateLimiter(interactive=False)
        self._noninteractive_started = False
        self._noninteractive_reported_best: dict[FitProgressContext, float] = {}
        self._grid_limiter = ProgressRateLimiter(interactive=interactive)
        self._grid_started = False
        self._completed_by_context: dict[FitProgressContext, int] = {}
        self._finished = False
        self._active = False

    @property
    def is_active(self) -> bool:
        """Whether the reporter currently owns an active rendering scope."""
        return self._active

    def __enter__(self) -> "MinimizationProgressReporter":
        self._active = True
        if not self._enabled:
            return self
        self._started_at = self._clock()
        if self._interactive:
            try:
                self._live = Live(
                    "",
                    console=self._console,
                    auto_refresh=False,
                    transient=True,
                )
                self._live.start(refresh=True)
            except Exception:  # noqa: BLE001 - reporting is non-scientific
                self._disable()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        _exc_value: BaseException | None,
        _traceback: object,
    ) -> bool:
        if self._enabled and exc_type is not None and not self._finished:
            status = (
                "interrupted" if issubclass(exc_type, KeyboardInterrupt) else "failed"
            )
            self._stop_live()
            with suppress(Exception):
                self._render_final(None, status)
        else:
            self._stop_live()
        self._active = False
        return False

    def observe(self, context: FitProgressContext, event: ProgressEvent) -> None:
        """Observe one contextual scalar event without affecting the fit."""
        if not self._enabled:
            return
        try:
            event = replace(event, elapsed_seconds=self._elapsed())
            self._completed_by_context[context] = max(
                event.objective_evaluations_completed,
                self._completed_by_context.get(context, 0),
            )
            update = (
                self._observe_grid(event)
                if self._grid
                else (
                    self._limiter(context).observe(event)
                    if self._interactive
                    else self._observe_noninteractive(context, event)
                )
            )
            if event.phase is ProgressPhase.TERMINATED and context.group_total == 1:
                update = None
            if update is not None:
                self._render(
                    _format_progress(
                        context,
                        update,
                        self._context_label(context),
                        self._step_name,
                    )
                )
        except KeyboardInterrupt:
            self._stop_live()
            raise
        except Exception:  # noqa: BLE001 - reporting failure is non-scientific
            self._disable()

    def finish(
        self,
        *,
        final_chi_square: float | None,
        terminal_status: str,
    ) -> None:
        """Close live rendering and emit one persistent aggregate terminal line."""
        if not self._enabled or self._finished:
            return
        self._finished = True
        self._stop_live()
        try:
            self._render_final(final_chi_square, terminal_status)
        except KeyboardInterrupt:
            raise
        except Exception:  # noqa: BLE001 - reporting is non-scientific
            self._disable()

    def _limiter(self, context: FitProgressContext) -> ProgressRateLimiter:
        return self._limiters.setdefault(
            context,
            ProgressRateLimiter(interactive=self._interactive),
        )

    def _observe_grid(self, event: ProgressEvent) -> ProgressUpdate | None:
        if event.phase is ProgressPhase.STARTED:
            if self._grid_started:
                return None
            self._grid_started = True
            return self._grid_limiter.observe(event)
        if event.phase is ProgressPhase.TERMINATED:
            return None
        probe = replace(
            event,
            current_chi_square=None,
            best_chi_square=None,
        )
        selected = self._grid_limiter.observe(probe)
        return None if selected is None else ProgressUpdate(event, None)

    def _observe_noninteractive(
        self,
        context: FitProgressContext,
        event: ProgressEvent,
    ) -> ProgressUpdate | None:
        phase = (
            ProgressPhase.STARTED
            if not self._noninteractive_started
            else ProgressPhase.EVALUATED
        )
        self._noninteractive_started = True
        elapsed = event.elapsed_seconds
        probe = replace(
            event,
            phase=phase,
            current_chi_square=None,
            best_chi_square=None,
            elapsed_seconds=elapsed,
            terminal_status=None,
        )
        if self._noninteractive_limiter.observe(probe) is None:
            return None
        previous = self._noninteractive_reported_best.get(context)
        best = event.best_chi_square
        relative_change = (
            None
            if best is None or previous is None
            else 0.0
            if previous == 0.0
            else (best - previous) / abs(previous)
        )
        if best is not None:
            self._noninteractive_reported_best[context] = best
        return ProgressUpdate(
            replace(event, phase=phase, elapsed_seconds=elapsed, terminal_status=None),
            relative_change,
        )

    def _context_label(self, context: FitProgressContext) -> str:
        return self._group_labels.get(frozenset(context.controlled_ids), "")

    def _elapsed(self) -> float:
        return (
            0.0
            if self._started_at is None
            else max(0.0, self._clock() - self._started_at)
        )

    def _render(self, line: str) -> None:
        text = Text(line)
        if self._live is not None:
            self._live.update(text, refresh=True)
        else:
            self._console.print(text)

    def _render_final(self, chi_square: float | None, status: str) -> None:
        reduced = (
            None
            if chi_square is None
            else chi_square
            / max(
                1,
                self._retained_observation_count - self._controlled_parameter_count,
            )
        )
        parts = [f"eval {sum(self._completed_by_context.values())} total"]
        if chi_square is not None:
            parts.append(f"final χ² {_format_scalar(chi_square)}")
        if reduced is not None:
            parts.append(f"red. χ² {_format_scalar(reduced)}")
        parts.extend((f"{self._elapsed():.1f} s", status))
        prefix = f"[{self._step_name}] " if self._step_name else ""
        self._console.print(Text(prefix + " · ".join(parts)))

    def _stop_live(self) -> None:
        if self._live is None:
            return
        live, self._live = self._live, None
        with suppress(Exception):
            live.stop()

    def _disable(self) -> None:
        self._enabled = False
        self._stop_live()


class UncertaintyProgressReporter:
    """Render the single automatic post-fit uncertainty phase."""

    def __init__(
        self,
        output_console: Console,
        *,
        clock: Callable[[], float] = monotonic,
    ) -> None:
        self._console = output_console
        self._clock = clock
        self._started_at: float | None = None
        self._finished = False

    def start(self) -> None:
        """Start timing immediately before automatic uncertainty derivation."""
        if self._started_at is not None:
            return
        self._started_at = self._clock()
        with suppress(Exception):
            self._console.print(Text("  • Estimating parameter uncertainties..."))

    def finish(self, status: str) -> None:
        """Emit one concise terminal status without affecting scientific work."""
        if self._finished:
            return
        self._finished = True
        elapsed = (
            0.0
            if self._started_at is None
            else max(0.0, self._clock() - self._started_at)
        )
        with suppress(Exception):
            self._console.print(Text(f"    {status} · {elapsed:.1f} s"))


def _format_scalar(value: float) -> str:
    return f"{value:.6g}"


def _format_progress(
    context: FitProgressContext,
    update: ProgressUpdate,
    group_label: str,
    step_name: str,
) -> str:
    event = update.event
    labels: list[str] = [step_name] if step_name else []
    if context.grid_seed_ordinal is not None and context.grid_seed_total is not None:
        labels.append(
            f"GRID seed {context.grid_seed_ordinal}/{context.grid_seed_total}"
        )
    if context.group_total > 1:
        labels.append(f"group {context.group_ordinal}/{context.group_total}")
        if group_label:
            labels.append(group_label)
    prefix = f"[{' · '.join(labels)}] " if labels else ""
    parts = [
        (
            f"eval {event.objective_evaluations_completed}/"
            f"{event.objective_request_budget}"
        )
    ]
    if event.best_chi_square is None:
        parts.append("starting")
    else:
        parts.append(f"best χ² {_format_scalar(event.best_chi_square)}")
        reduced = event.reduced_chi_square
        if reduced is not None:
            parts.append(f"red. χ² {_format_scalar(reduced)}")
    if update.relative_best_change is not None:
        parts.append(f"Δχ² {100.0 * update.relative_best_change:+.2f}%")
    parts.append(f"{event.elapsed_seconds:.1f} s")
    if event.terminal_status:
        parts.append(event.terminal_status)
    return prefix + " · ".join(parts)


LOGO = r"""
   ________                   ______
  / ____/ /_  ___  ____ ___  / ____/  __
 / /   / __ \/ _ \/ __ `__ \/ __/ | |/_/
/ /___/ / / /  __/ / / / / / /____>  <
\____/_/ /_/\___/_/ /_/ /_/_____/_/|_|


"""
"""ASCII art logo of the ChemEx software."""

EXPERIMENT_NAME = """
[experiment]
name = "experiment_name"
...
"""
"""Template for the experiment name in configuration files."""


def print_logo() -> None:
    """Print the ChemEx software logo with version information."""
    logo = Text(LOGO, style="blue")
    description = "Analysis of NMR chemical exchange data\n\n"
    version = "Version: "
    version_number = Text(f"{__version__}", style="red")
    all_text = Text.assemble(logo, description, version, version_number)
    panel = Panel.fit(all_text)
    console.print(panel)


def print_loading_experiments() -> None:
    """Display a loading message for datasets."""
    console.print("\nLoading datasets...", style="bold yellow")


def get_reading_exp_text(filename: Path, name: str = "", total_nb: int = 0) -> Text:
    """Generate a formatted message for reading an experiment file.

    Args:
        filename (Path): The path of the file being read.
        name (str, optional): The experiment name. Defaults to an empty string.
        total_nb (int, optional): The total number of profiles. Defaults to 0.

    Returns:
        Text: A formatted message indicating the file name, type, and profile count.

    """
    parts = (
        "  • Reading ",
        (f"{filename}", "green"),
        " (",
        (f"{name}", "blue"),
        ")",
    )

    if total_nb == 0:
        return Text.assemble(*parts, " ...")

    nb_text = Text(f"{total_nb}", style="bold")

    return Text.assemble(*parts, " -> ", nb_text, " profiles")


def print_reading_defaults() -> None:
    """Print a message indicating the reading of default parameters."""
    console.print("\nReading default parameters...", style="bold yellow")


def print_reading_methods() -> None:
    """Display a message indicating that methods are being read."""
    console.print("\nReading methods...", style="bold yellow")


def print_start_fit() -> None:
    """Display a message indicating the start of the fitting process."""
    console.print("\nStarting the fits...", style="bold yellow")


def print_running_simulations() -> None:
    """Display a message indicating that simulations are running."""
    console.print("\nRunning simulations...", style="bold yellow")


def print_step_name(name: str, index: int, total: int) -> None:
    """Print the name of the current step in the analysis process.

    Args:
        name (str): Name of the step.
        index (int): Current step index.
        total (int): Total number of steps.

    """
    text = Text.assemble(
        "Running ",
        (f"{name}", "magenta"),
        f" ({index}/{total})\n",
    )
    console.print(Padding(Rule(text, end="\n\n"), (1, 0, 0, 1)), width=51)


def print_selecting_profiles(selected_nb: int) -> None:
    """Display a message indicating the number of profiles selected.

    Args:
        selected_nb (int): The number of selected profiles.

    """
    console.print(
        Text.from_markup(f"  • Selecting profiles -> [blue]{selected_nb}[/] profiles"),
    )


def print_no_data() -> None:
    """Inform the user that no data is available for fitting."""
    console.print("  • No data to fit")


def print_fitmethod(fit_method: str) -> None:
    """Display the chosen method for fitting.

    Args:
        fit_method (str): The method being used for fitting.

    """
    console.print(Text.assemble("  • Fit method -> ", (f"{fit_method}", "blue")))


def print_status_changes(
    vary: Counter[str],
    fixed: Counter[str],
    constrained: Counter[str],
) -> None:
    """Print changes in parameter status like varied, fixed, or constrained.

    Args:
        vary (Counter[str]): Parameters varied during fitting.
        fixed (Counter[str]): Parameters fixed during fitting.
        constrained (Counter[str]): Parameters constrained during fitting.

    """
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
    """Display a message indicating that the minimization process is running."""
    console.print(Text("  • Running the minimization..."))


def print_running_grid() -> None:
    """Inform the user that the grid search is in progress."""
    console.print(Text("  • Running the grid search..."))


def print_running_statistics(name: str) -> None:
    """Display a message indicating that statistical simulations are running.

    Args:
        name (str): The name of the statistical simulation.

    """
    console.print(Text(f"  • Running {name} simulations..."))


def print_section(name: str) -> None:
    """Print the name of the current section in the analysis process.

    Args:
        name (str): Name of the section.

    """
    console.print(Text(f"  • Section {name}"))


def print_chi2_table_header() -> None:
    """Print the header for the chi-squared table."""
    header = Text(f"{'Evaluation':>10s}  {'χ²':>12s}  {'Reduced χ²':>12s}")
    console.print()
    console.print(Padding.indent(header, 5), style="bold")
    console.print(Padding.indent("─" * 39, 4))


def print_chi2_table_line(iteration: int, chisqr: float, redchi: float) -> None:
    """Print a line in the chi-squared table with evaluation and chi-squared values.

    Args:
        iteration (int): Objective evaluation number (legacy keyword name).
        chisqr (float): Chi-squared value.
        redchi (float): Reduced chi-squared value.

    """
    line = Text(f"{iteration:>9d}  {chisqr:>12.1f}  {redchi:>12.3f}")
    console.print(Padding.indent(line, 5))


def print_chi2_table_footer(iteration: int, chisqr: float, redchi: float) -> None:
    """Print the footer for the chi-squared table.

    Args:
        iteration (int): Objective evaluation number (legacy keyword name).
        chisqr (float): Chi-squared value.
        redchi (float): Reduced chi-squared value.

    """
    footer = Text(f"{iteration:>9d}  {chisqr:>12.1f}  {redchi:>12.3f}")
    console.print(Padding.indent("─" * 39, 4))
    console.print(Padding.indent(footer, 5), style="bold")
    console.print()


def print_chi2(chisqr: float, redchi: float) -> None:
    """Display the chi-squared and reduced chi-squared values.

    Args:
        chisqr (float): Chi-squared value.
        redchi (float): Reduced chi-squared value.

    """
    console.print()
    console.print(
        Padding.indent(Text.from_markup(f"        χ²: [bold]{chisqr:15.1f}[/]"), 3),
    )
    console.print(
        Padding.indent(Text.from_markup(f"Reduced χ²: [bold]{redchi:15.3f}[/]"), 3),
    )
    console.print()


def print_writing_results(path: Path) -> None:
    """Inform the user about the location where results are being written.

    Args:
        path (Path): Path to the file or directory where results are saved.

    """
    console.print(f"  • Writing results in [green]{path}")


def print_making_plots() -> None:
    """Display a message indicating that plots are being generated."""
    console.print(Text("  • Making plots..."))


def print_plot_filename(filename: Path, *, extra: bool = True) -> None:
    """Display the filename of the plot being generated.

    Args:
        filename (Path): The path of the plot file.
        extra (bool, optional): Indicates if extra information is included.
                                Defaults to True.

    """
    text = f"    ‣ [green]{filename}[/]"

    if extra:
        text += " [[green].fit[/], [green].exp[/]]"

    console.print(Text.from_markup(text))


def print_group_name(text: str) -> None:
    """Print the name of the group in the analysis process.

    Args:
        text (str): Name of the group.

    """
    console.print(Padding(Rule(text, characters="⋅"), (1, 0, 0, 3)), width=49)


def print_file_not_found(filename: Path) -> None:
    """Inform the user that a specified file was not found.

    Args:
        filename (Path): The path of the file that was not found.

    """
    console.print()
    console.print(f"[red]The file '{filename}' is empty or does not exist!")


def print_file_not_found_error(error: FileNotFoundError) -> None:
    """Display the file not found error message.

    Args:
        error (FileNotFoundError): The FileNotFoundError exception instance.

    """
    console.print()
    console.print(f"[red]Error: {error}")


def print_dataset_error(filename: Path, explanation: str) -> None:
    """Display a concise error for invalid user-supplied dataset content."""
    console.print()
    console.print(f"[red]Error: Invalid data in '{filename}': {explanation}")


def print_toml_error(filename: Path, error_message: Exception) -> None:
    """Display an error related to TOML file parsing.

    Args:
        filename (Path): Path of the TOML file.
        error_message (Exception): The exception instance with the error message.

    """
    console.print()
    console.print(Text.from_markup(f"[red]Error in the TOML file '{filename}'"))
    console.print(
        Text.from_markup(f"[red]-> {error_message}\n"),
    )
    console.print(
        "Please check the website [link]https://toml.io[/]"
        " for a complete description of the TOML format.",
    )


def print_experiment_name_error(filename: Path) -> None:
    """Display an error message when the experiment name is missing from a file.

    Args:
        filename (Path): Path of the file with the missing experiment name.

    """
    console.print()
    console.print(f"[red]The experiment name is missing from the file '{filename}'\n")
    console.print(
        f"Please make sure that the 'name' field is provided in '{filename}':\n",
    )
    console.print(
        Syntax(EXPERIMENT_NAME, "toml"),
        "\nRun 'chemex info' to obtain the full list of the available experiments.",
    )


def print_experiment_type_error(filename: Path, name: str) -> None:
    """Display an error when an input names an unknown Experiment Type."""
    console.print()
    console.print(
        f"[red]Unknown Experiment Type '{name}' in the file '{filename}'.[/]",
    )
    console.print(
        "Run 'chemex info' to obtain the full list of available experiments.",
    )


def print_experiment_source_error(filename: Path | None, explanation: str) -> None:
    """Display an error for an invalid or stale opened Experiment source."""
    console.print()
    location = f" for '{filename}'" if filename is not None else ""
    console.print(f"[red]Invalid Experiment source{location}: {explanation}")


def print_pydantic_parsing_error(filename: Path, error: ValidationError) -> None:
    """Display errors encountered while parsing a file using Pydantic.

    Args:
        filename (Path): Path of the file being parsed.
        error (ValidationError): Instance detailing parsing errors.

    """
    console.print()
    console.print(Text.from_markup(f"[red]Error(s) while parsing '{filename}'"))
    for err in error.errors():
        location = " -> ".join(str(loc) for loc in err["loc"])
        console.print("  • ", location, ":", err["msg"], style="red")


def print_calculation_stopped_error() -> None:
    """Inform the user that the calculation was manually stopped."""
    console.print()
    console.print("[red] -- Keyboard Interrupt: Calculation stopped --")
    console.print()


def print_plotting_canceled() -> None:
    """Display a message indicating that the plotting process was cancelled."""
    console.print()
    console.print("[red] -- Keyboard Interrupt: Plotting cancelled --")
    console.print()


def print_value_error() -> None:
    """Inform the user that a ValueError occurred, stopping the calculation."""
    console.print()
    console.print("[red] -- Got a ValueError: Calculation stopped --")
    console.print()


def print_warning_positive_jnh() -> None:
    """Warn about positive 1J(NH) coupling values."""
    console.print()
    console.print(
        "[yellow] -- WARNING: Some 1J(NH) couplings are set with positive values --",
    )
    console.print(
        "This can cause the TROSY and anti-TROSY components to be switched in some"
        " experiments.",
    )
    console.print()


def print_warning_negative_jch() -> None:
    """Warn about negative 1J(CH) coupling values."""
    console.print()
    console.print(
        "[yellow] -- WARNING: Some 1J(CH) couplings are set with negative values --",
    )
    console.print(
        "This can cause the TROSY and anti-TROSY components to be switched in some"
        " experiments.",
    )
    console.print()


def print_error_grid_settings(entry: str, detail: str | None = None) -> None:
    """Display an error related to grid settings.

    Args:
        entry (str): The problematic entry in the grid settings.
        detail (str | None): Optional explanation of the validation failure.

    """
    console.print()
    console.print("[red] -- ERROR: Error reading grid settings:")
    console.print(Text(f'    "{entry}"', style="red"))
    if detail is not None:
        console.print(Text(f"    {detail}", style="red"))
    console.print()
    console.print(
        "Please make sure that the grid settings are provided in the correct format",
    )
    console.print(Text("    [PB] = lin(<min>, <max>, <nb of points>)"))
    console.print(Text("    [PB] = log(<min>, <max>, <nb of points>)"))
    console.print(Text("    [PB] = (<value1>, <value2>, ..., <valuen>)"))
    console.print()


def print_error_constraints(expression: str, detail: str | None = None) -> None:
    """Display an error message for issues with constraint expressions.

    Args:
        expression (str): The problematic constraint expression.
        detail (str | None): Optional explanation of the validation failure.

    """
    console.print()
    console.print(f'[red] -- ERROR: Error reading constraints -> "{expression}" --')
    if detail is not None:
        console.print(Text(f"    {detail}", style="red"))
    console.print()


def print_grid_statistic_warning() -> None:
    """Warn that 'GRID' and 'STATISTICS' options are mutually exclusive."""
    console.print()
    console.print(
        "[yellow] -- WARNING: The 'GRID' and 'STATISTICS' options are mutually "
        "exclusive. Only the 'GRID' calculation will be run.",
    )
    console.print()


def print_mcmc_no_vary_warning() -> None:
    """Warn that MCMC was requested without fitted parameters."""
    console.print()
    console.print(
        "[yellow] -- WARNING: MCMC was requested, but no parameters are marked "
        "as fitted. Skipping MCMC.",
    )
    console.print()


def print_mcmc_unbounded_warning(parameters: list[str]) -> None:
    """Warn that MCMC was requested with unbounded parameters."""
    console.print()
    console.print(
        "[yellow] -- WARNING: Some fitted parameters do not have finite MCMC bounds.",
    )
    console.print(
        "Uniform priors are inferred from parameter bounds; finite lower and "
        "upper bounds are recommended for MCMC.",
    )
    console.print("Affected parameters:")
    for parameter in parameters:
        console.print(f"    - {parameter}")
    console.print()


def print_model_error(name: str) -> None:
    """Display an error message for unavailable models.

    Args:
        name (str): The name of the model that is not available.

    """
    from chemex.models.factory import model_factory

    console.print()
    console.print(f"[red] -- ERROR: The model '{name}' is not available.")
    console.print()
    console.print("The available models are:")
    for model_name in sorted(model_factory.set):
        console.print(f"    - '{model_name}'")
    console.print()


def print_model_suffix_error(
    name: str,
    unknown_suffixes: set[str],
    supported_suffixes: tuple[str, ...],
) -> None:
    """Display an error message for unsupported model suffixes."""
    suffixes = ", ".join(f"'.{suffix}'" for suffix in sorted(unknown_suffixes))

    console.print()
    console.print(
        f"[red] -- ERROR: The model '{name}' uses unsupported suffixes: {suffixes}.",
    )
    console.print()
    console.print("The supported model suffixes are:")
    for suffix in supported_suffixes:
        console.print(f"    - '.{suffix}'")
    console.print()


def print_not_implemented_noise_method_warning(
    filename: Path,
    kind: str,
    implemented: tuple[str, ...],
) -> None:
    """Warn about unimplemented noise methods for an experiment.

    Args:
        filename (Path): Path of the experiment file.
        kind (str): The kind of noise method that is not implemented.
        implemented (tuple[str, ...]): Tuple of implemented methods.

    """
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
    """Warn about the absence of duplicate points in some profiles.

    Args:
        filename (Path): Path of the experiment file.

    """
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
  - "FITMETHOD" supports only "trf" (the default).
    "least_squares" is accepted temporarily as an alias for "trf".
    Historical optimizer names, including "leastsq" and
    "differential_evolution", are not supported."""

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
  - "STATISTICS" must be a dictionary with keys 'MC', 'BS', 'BSN', 'MCMC'

    Example: { "MC"=10 }, { "MCMC"=5000 }, or { "MC"=10, "BS"=10 }"""

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


def print_wrong_option(option: str) -> str:
    """Display a message for an invalid option in configuration files.

    Args:
        option (str): The option that is invalid or incorrectly formatted.

    Returns:
        str: A message indicating the error with the specific option.

    """
    return (
        f"\n  - '{option.upper()}' is not a valid option.\n"
        "\nPlease check the website "
        "[link]"
        "https://gbouvignies.github.io/ChemEx/docs/user_guide/fitting/method_files"
        "[/] "
        "for a complete list of valid options."
    )


def print_method_error(filename: Path, section: str, options: set[int | str]) -> None:
    """Display errors for invalid methods or options in a specific section of a file.

    Args:
        filename (Path): The file containing the error.
        section (str): The section of the file with the erroneous method or option.
        options (set[int | str]): Set of options or methods that are incorrect.

    """
    console.print()
    console.print(
        f"[red] -- ERROR: Error reading section '[{section}]' of '{filename}' --",
    )
    for option in options:
        console.print(
            METHOD_ERROR_MESSAGES.get(str(option), print_wrong_option(str(option))),
        )
    console.print()

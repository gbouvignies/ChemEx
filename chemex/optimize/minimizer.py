from __future__ import annotations

import sys
from contextlib import nullcontext
from copy import deepcopy
from dataclasses import dataclass
from dataclasses import field

import numpy as np
from lmfit import Minimizer
from lmfit import Parameters
from rich.live import Live
from rich.padding import Padding
from rich.table import Table

from chemex.containers.experiments import Experiments
from chemex.messages import console
from chemex.messages import create_chi2_table


@dataclass
class Reporter:
    last_chisqr: float = 1e32
    threshold: float = -1.0e-3
    table: Table = field(default_factory=create_chi2_table)

    def iter_cb(
        self, params: Parameters, iteration: int, residuals: np.ndarray
    ) -> None:

        chisqr = (residuals**2).sum()
        change = (chisqr - self.last_chisqr) / self.last_chisqr

        if change > self.threshold or iteration < 0:
            return

        self.last_chisqr = chisqr

        ndata = len(residuals)
        nvarys = len(
            [param for param in params.values() if param.vary and not param.expr]
        )
        redchi = chisqr / max(1, ndata - nvarys)

        self.table.add_row(f"{iteration:d}", f"{chisqr:.1f}", f"{redchi:.3f}")


def minimize(
    experiments: Experiments,
    params: Parameters,
    fitmethod: str | None = None,
    verbose: bool = True,
) -> Parameters:

    if fitmethod is None:
        fitmethod = "leastsq"

    kws = {
        "brute": {"keep": "all"},
        "basinhopping": {"niter_success": 10},
    }

    iter_cb = None
    reporter = None
    live = nullcontext()

    if verbose:
        reporter = Reporter()
        iter_cb = reporter.iter_cb
        indented_table = Padding.indent(reporter.table, 3)
        live = Live(
            indented_table,
            console=console,
            refresh_per_second=20,
        )

    minimizer = Minimizer(experiments.residuals, params, iter_cb=iter_cb)

    with live:
        try:
            result = minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
            if verbose and reporter is not None:
                reporter.table.columns[0].footer = f"{result.nfev:d}"
                reporter.table.columns[1].footer = f"{result.chisqr:.1f}"
                reporter.table.columns[2].footer = f"{result.redchi:.3f}"
                reporter.table.show_footer = True
        except KeyboardInterrupt:
            sys.stderr.write("\n -- Keyboard Interrupt: minimization stopped --\n")
            result = minimizer.result
        except ValueError:
            result = minimizer.result
            sys.stderr.write("\n -- Got a ValueError: minimization stopped --\n")

    return deepcopy(result.params)  # type: ignore

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass

import numpy as np
from lmfit import Minimizer
from lmfit import Parameters

from chemex.containers.experiments import Experiments
from chemex.messages import print_calculation_stopped_error
from chemex.messages import print_chi2_table_footer
from chemex.messages import print_chi2_table_header
from chemex.messages import print_chi2_table_line
from chemex.messages import print_value_error


@dataclass
class Reporter:
    last_chisqr: float = 1e32
    threshold: float = -1.0e-3

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

        self.print_line(iteration, chisqr, redchi)

    def print_line(self, iteration: int, chisqr: float, redchi: float) -> None:
        print_chi2_table_line(iteration, chisqr, redchi)

    def print_header(self) -> None:
        print_chi2_table_header()

    def print_footer(self, iteration: int, chisqr: float, redchi: float) -> None:
        print_chi2_table_footer(iteration, chisqr, redchi)


def minimize(
    experiments: Experiments, params: Parameters, fitmethod: str | None = None
) -> Parameters:

    if fitmethod is None:
        fitmethod = "leastsq"

    kws = {
        "brute": {"keep": "all"},
        "basinhopping": {"niter_success": 10},
    }

    minimizer = Minimizer(experiments.residuals, params)
    minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
    return deepcopy(minimizer.result.params)  # type: ignore


def minimize_with_report(
    experiments: Experiments, params: Parameters, fitmethod: str | None = None
) -> Parameters:

    if fitmethod is None:
        fitmethod = "leastsq"

    kws = {
        "brute": {"keep": "all"},
        "basinhopping": {"niter_success": 10},
    }

    reporter = Reporter()

    minimizer = Minimizer(experiments.residuals, params, iter_cb=reporter.iter_cb)

    reporter.print_header()

    try:
        minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
        reporter.print_footer(
            minimizer.result.nfev, minimizer.result.chisqr, minimizer.result.redchi
        )
    except KeyboardInterrupt:
        print_calculation_stopped_error()
    except ValueError:
        print_value_error()

    return deepcopy(minimizer.result.params)  # type: ignore

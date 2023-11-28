from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import TYPE_CHECKING

import lmfit
import numpy as np
from lmfit import Minimizer, Parameters

from chemex.messages import (
    print_calculation_stopped_error,
    print_chi2_table_footer,
    print_chi2_table_header,
    print_chi2_table_line,
    print_value_error,
)
from chemex.optimize.grouping import ParamTree, create_group_tree

if TYPE_CHECKING:
    from chemex.containers.experiments import Experiments
    from chemex.typing import ArrayFloat


@dataclass
class Reporter:
    last_chisqr: float = 1e32
    threshold: float = -1.0e-3

    def iter_cb(
        self,
        params: Parameters,
        iteration: int,
        residuals: ArrayFloat,
        *_args,
        **_kwargs,
    ) -> None:
        chisqr = (residuals**2).sum()
        change = (chisqr - self.last_chisqr) / self.last_chisqr

        if change > self.threshold or iteration < 0:
            return

        self.last_chisqr = chisqr

        ndata = len(residuals)
        nvarys = len(
            [param for param in params.values() if param.vary and not param.expr],
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
    experiments: Experiments,
    params: Parameters,
    fitmethod: str,
) -> Parameters:
    kws = {
        "brute": {"keep": "all"},
        "basinhopping": {"niter_success": 10},
    }

    minimizer = Minimizer(experiments.residuals, params)
    minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
    return deepcopy(minimizer.result.params)


def minimize_with_report(
    experiments: Experiments,
    params: Parameters,
    fitmethod: str,
) -> Parameters:
    kws = {
        "brute": {"keep": "all"},
        "basinhopping": {"niter_success": 10},
        "differential_evolution": {"polish": True, "disp": True, "init": "sobol"},
        "ampgo": {"disp": True},
    }

    reporter = Reporter()

    minimizer = Minimizer(experiments.residuals, params, iter_cb=reporter.iter_cb)

    reporter.print_header()

    try:
        result = minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
        reporter.print_footer(result.nfev, result.chisqr, result.redchi)
    except KeyboardInterrupt:
        print_calculation_stopped_error()
    except ValueError:
        print_value_error()

    return deepcopy(minimizer.result.params)


def residuals_hierarchical(params: Parameters, param_tree: ParamTree) -> ArrayFloat:
    # print(param_tree.ids_to_fit)
    if not param_tree.branches:
        return param_tree.experiments.residuals(params)

    for param_id in param_tree.ids_to_fit:
        params[param_id].vary = False

    residuals_list: list[ArrayFloat] = []

    for branch in param_tree.branches:
        branch_results = lmfit.minimize(
            residuals_hierarchical,
            params,
            args=(branch,),
            method="leastsq",
        )
        residuals_list.append(branch_results.residual)

    residuals = np.concatenate(residuals_list, axis=0)

    for param_id in param_tree.ids_to_fit:
        params[param_id].vary = True

    return residuals


def minimize_hierarchical(experiments: Experiments, params: Parameters, fitmethod: str):
    kws = {
        "brute": {"keep": "all"},
        "basinhopping": {"niter_success": 10},
        "differential_evolution": {"polish": True, "disp": True, "init": "sobol"},
        "ampgo": {"disp": True},
    }

    reporter = Reporter()

    param_tree = create_group_tree(experiments)

    new_params = params.copy()
    for param in new_params.values():
        if param.name not in param_tree.ids_to_fit:
            param.vary = False

    minimizer = Minimizer(
        residuals_hierarchical,
        params,
        fcn_args=(param_tree,),
        iter_cb=reporter.iter_cb,
    )

    reporter.print_header()

    try:
        result = minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
        reporter.print_footer(result.nfev, result.chisqr, result.redchi)
    except KeyboardInterrupt:
        print_calculation_stopped_error()
    except ValueError:
        print_value_error()

    return deepcopy(minimizer.result.params)

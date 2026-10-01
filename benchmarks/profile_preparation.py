"""Reproduce the production request-budget qualification without outer search.

Use the single-native-thread environment in benchmarks/README.md. Historical
SEARCH.DE comparisons remain in basin_discovery_results.json; no DE backend is
retained. The difficult DCEST state is one already evaluated research point.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch

from benchmarks.basin_discovery import build_context, measured, refine

import chemex.optimize.profile_preparation as preparation

POLICIES = ((100, None), (200, None), (500, None), (500, 1000), (600, 3000))


def compare_budgets(case: str) -> list[dict]:
    context = build_context(case)
    root = context.problem
    hold = tuple(key for key, *_ in context.coordinates if not key.startswith("__DW_"))
    mirror = tuple(key for key in root.controlled_ids if key.startswith("__DW_"))
    rows = []
    endpoints = {}
    for multiplier, ceiling in POLICIES:
        with patch.object(
            preparation,
            "preparation_request_budget",
            lambda n, multiplier=multiplier, ceiling=ceiling: min(
                multiplier * (n + 1), ceiling or 100000
            ),
        ):
            result, trace = measured(
                context,
                lambda: preparation.prepare_profile_start(
                    root, hold, mirror, context.parameterization, context.engine
                ),
            )
        reused = result.vector in endpoints
        if not reused:
            endpoints[result.vector] = refine(context, result.vector)
        row = {
            "case": case,
            "multiplier": multiplier,
            "ceiling": ceiling,
            "preparation": trace,
            "fallback_factors": result.fallback_factors,
            "root_fallback": result.root_fallback,
            "final": endpoints[result.vector],
            "reused_final_measurement": reused,
        }
        rows.append(row)
        print(
            case,
            multiplier,
            ceiling,
            trace["requests"],
            row["final"]["chi2"],
            flush=True,
        )
    return rows


def dcest_stress() -> list[dict]:
    retained = json.loads(
        Path(__file__).with_name("basin_discovery_results.json").read_text()
    )
    run = next(item for item in retained["runs"] if item["case"] == "dcest")
    context = build_context("dcest")
    hold = tuple(item[0] for item in run["outer_coordinates"])
    updates = dict(zip(hold, run["points"][1], strict=True))
    root = context.problem.restart_from(
        tuple(
            updates.get(key, value)
            for key, value in zip(
                context.problem.controlled_ids, context.problem.start, strict=True
            )
        )
    )
    mirror = tuple(key for key in root.controlled_ids if key.startswith("__DW_"))
    rows = []
    for cap in (1000, 3000):
        with patch.object(preparation, "preparation_request_budget", return_value=cap):
            result, trace = measured(
                context,
                lambda: preparation.prepare_profile_start(
                    root, hold, mirror, context.parameterization, context.engine
                ),
            )
        rows.append(
            {
                "cap": cap,
                "held": updates,
                "preparation": trace,
                "fallback": result.fallback_factors,
                "root_fallback": result.root_fallback,
                "original_start_retained": result.vector == root.start,
            }
        )
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case", choices=("cpmg", "cest", "dcest", "all"), default="all"
    )
    parser.add_argument("--stress", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = ("cpmg", "cest", "dcest") if args.case == "all" else (args.case,)
    output = {
        "budget_comparison": [row for case in cases for row in compare_budgets(case)]
    }
    if args.stress:
        output["dcest_stress"] = dcest_stress()
    args.output.write_text(json.dumps(output, indent=2) + "\n")

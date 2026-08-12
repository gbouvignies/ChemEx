"""Diagnostic plot rendering for native step-root publications."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from chemex.evaluation.native import EvaluationPlan, EvaluationResult


def write_native_plots(
    path: Path,
    plan: EvaluationPlan,
    result: EvaluationResult,
) -> None:
    """Render one aggregate PDF containing every evaluated profile."""
    path.mkdir()
    with PdfPages(path / "profiles.pdf") as pdf:
        for ordinal, descriptor in enumerate(plan.profiles):
            start = descriptor.observation_offset
            stop = start + descriptor.observation_count
            x = np.arange(descriptor.observation_count)
            retained = np.asarray(descriptor.mask, dtype=np.bool_)
            figure, axis = plt.subplots(figsize=(7.5, 5.0))
            axis.errorbar(
                x[retained],
                np.asarray(descriptor.experimental_values)[retained],
                yerr=np.asarray(descriptor.uncertainties)[retained],
                fmt="o",
                label="experimental (retained)",
            )
            if np.any(~retained):
                axis.scatter(
                    x[~retained],
                    np.asarray(descriptor.experimental_values)[~retained],
                    marker="x",
                    label="experimental (masked)",
                )
            axis.plot(
                x,
                result.normalized_calculations[start:stop],
                "-",
                label="calculated",
            )
            axis.set_title(f"Profile {ordinal:04d}")
            axis.set_xlabel("Observation ordinal")
            axis.set_ylabel("Intensity")
            axis.legend()
            figure.tight_layout()
            pdf.savefig(figure)
            plt.close(figure)

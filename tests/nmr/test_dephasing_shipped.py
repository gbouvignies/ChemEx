"""Qualify complete B1 dephasing on shipped CEST pulse matrices."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chemex.configuration.methods import Selection
from chemex.configuration.parameters import read_defaults
from chemex.experiments.builder import build_experiments
from chemex.parameters.spin_system import SpinSystem
from chemex.runtime import AnalysisSession

EXAMPLES = Path(__file__).parents[2] / "examples/Experiments"


@pytest.mark.parametrize(
    ("family", "experiment_file", "spin", "delay", "consumer"),
    [
        ("CEST_1HN_IP_AP", "30hz.toml", "74H-N", 0.4, "pulse_i"),
        ("CEST_13C_LABEL_CN", "23hz.toml", "L18CB", 0.25, "pulse_i"),
        ("CEST_15N_CW", "dec.toml", "52N-H", 0.5, "pulse_is"),
    ],
)
def test_shipped_dephasing_pulses_match_legacy_underdamped_outputs(
    family: str, experiment_file: str, spin: str, delay: float, consumer: str
) -> None:
    example = EXAMPLES / family
    session = AnalysisSession.create()
    session.set_model("2st")
    experiments = build_experiments(
        [example / "Experiments" / experiment_file],
        Selection(include=[SpinSystem.from_name(spin)], exclude=None),
        session=session,
    )
    session.parameters.set_defaults(
        read_defaults([example / "Parameters/parameters.toml"])
    )
    assert session.try_build_analysis_values()
    profile = next(iter(experiments)).profiles[0]
    profile.update_spectrometer_from_values(
        dict(session.resolve_current_values(experiments.param_ids))
    )
    spectrometer = profile.spectrometer
    spectrometer.offset_i = 0.0
    engine = spectrometer._engine
    assert engine.b1_i_dist.dephasing

    liouv = engine.l_free + engine.l_b1x_i
    if consumer == "pulse_is":
        liouv = liouv + engine.l_b1x_s
        actual = spectrometer.pulse_is(delay, 0.0, 0.0)
    else:
        actual = spectrometer.pulse_i(delay, 0.0)

    eigenvalues, eigenvectors = np.linalg.eig(liouv)
    oscillatory = np.abs(eigenvalues.imag) >= 1.0e-6
    assert oscillatory.any()
    assert np.max(eigenvalues.real[oscillatory]) < -1.0
    # The old multiplier happens to underflow for these damped shipped modes.
    legacy_eigenvalues = np.where(oscillatory, eigenvalues * 1.0e9, eigenvalues)
    legacy_weights = np.exp(delay * legacy_eigenvalues)
    legacy = (eigenvectors * legacy_weights[..., None, :]) @ np.linalg.inv(eigenvectors)

    np.testing.assert_allclose(actual, legacy.real, rtol=0.0, atol=2.0e-13)

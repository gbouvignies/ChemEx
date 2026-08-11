from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chemex.configuration.methods import Selection
from chemex.configuration.parameters import read_defaults
from chemex.experiments.builder import build_experiments
from chemex.experiments.catalog.cpmg_15n_rc import (
    EXPERIMENT_NAME,
    Cpmg15NRcConfig,
    Cpmg15NRcSequence,
)
from chemex.models.model import ModelSpec
from chemex.runtime import AnalysisSession

TIME_T2 = 0.040
NCYC_MAX = 40


def _build_config(model: str = "2st", **overrides: object) -> Cpmg15NRcConfig:
    experiment = {
        "name": EXPERIMENT_NAME,
        "time_t2": TIME_T2,
        "carrier": 118.0,
        "pw90": 40.0e-6,
        "ncyc_max": NCYC_MAX,
    } | overrides
    return Cpmg15NRcConfig.model_validate(
        {
            "experiment": experiment,
            "conditions": {"h_larmor_frq": 600.0, "temperature": 25.0},
            "data": {"path": ".", "profiles": {}},
            "model": ModelSpec(name=model, states="ab"[: 1 if model == "1st" else 2]),
        },
    )


NCYCS = (0, 2, 4, 8, 12, 18, 24, 32, 40)


def _write_input(tmp_path: Path) -> Path:
    (tmp_path / "93N-HN.out").write_text(
        "\n".join(f"{ncyc:6d} 1.0e6 1.0e4" for ncyc in NCYCS) + "\n",
    )
    filename = tmp_path / "experiment.toml"
    filename.write_text(
        f"""[experiment]
name = "{EXPERIMENT_NAME}"
time_t2 = {TIME_T2}
carrier = 118.0
pw90 = 40.0e-6
ncyc_max = {NCYC_MAX}

[conditions]
h_larmor_frq = 600.0
temperature = 25.0

[data]
path = "."

[data.profiles]
93N = "93N-HN.out"
""",
    )
    return filename


def _calculated_profile(tmp_path: Path, model: str) -> np.ndarray:
    filename = _write_input(tmp_path)
    parameters = tmp_path / "parameters.toml"
    parameters.write_text(
        """[GLOBAL]
R2_I_A = 20.0
TAUC_A = 9.0
PB = 0.03
KEX_AB = 800.0
DW_AB = 3.0

[CS_A]
93N = 118.0
""",
    )

    session = AnalysisSession()
    session.set_model(model)
    experiments = build_experiments(
        [filename],
        Selection(include=None, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(read_defaults([parameters]))
    params = session.parameters.build_lmfit_params(experiments.param_ids)
    (experiment,) = list(experiments)
    (profile,) = list(experiment)
    return profile.calculate(params)


def test_settings_expose_the_expected_coherences() -> None:
    settings = _build_config().experiment

    assert settings.start_terms == ["2izsz"]
    assert settings.detection == "[iz_a]"
    assert settings.even_ncycs is True
    assert settings.pw180 == pytest.approx(2.0 * settings.pw90)


def test_reference_plane_is_flagged() -> None:
    metadata = np.array([0.0, 2.0, 40.0])

    assert Cpmg15NRcSequence.is_reference(metadata).tolist() == [True, False, False]


@pytest.mark.parametrize(
    ("ncyc", "expected"),
    [(0.0, 0), (1.0, 0), (2.0, 2), (3.0, 2), (40.0, 40)],
)
def test_odd_ncyc_is_rounded_down_like_the_pulse_sequence(
    ncyc: float,
    expected: int,
) -> None:
    assert Cpmg15NRcSequence.even_ncyc(ncyc) == expected


def test_inter_pulse_delay_matches_the_rounded_cpmg_frequency() -> None:
    settings = _build_config().experiment
    sequence = Cpmg15NRcSequence(settings)

    tau_cps, _, _ = sequence._get_delays(np.array([0.0, 2.0, 40.0]))

    nu_cpmg = 2.0 / TIME_T2
    expected = 1.0 / (4.0 * nu_cpmg) - 0.375 * settings.pw180
    assert tau_cps[2.0] == pytest.approx(expected)
    assert 0.0 not in tau_cps


def test_compensation_delay_vanishes_at_the_largest_ncyc() -> None:
    settings = _build_config().experiment
    sequence = Cpmg15NRcSequence(settings)

    _, deltas, _ = sequence._get_delays(np.array([0.0, NCYC_MAX]))

    assert deltas[float(NCYC_MAX)] == pytest.approx(0.0)
    assert deltas[0.0] == pytest.approx(0.25 * NCYC_MAX * settings.pw180)


def test_baseline_is_flat_without_exchange(tmp_path: Path) -> None:
    """The equal anti-phase / in-phase periods and the R1 compensation delay
    should leave no residual dependence on the CPMG frequency."""
    intensities = _calculated_profile(tmp_path, "1st")
    rates = -np.log(intensities[1:] / intensities[0]) / TIME_T2

    assert float(np.ptp(rates)) < 0.1


def test_dispersion_is_monotonic_with_exchange(tmp_path: Path) -> None:
    intensities = _calculated_profile(tmp_path, "2st")
    rates = -np.log(intensities[1:] / intensities[0]) / TIME_T2

    assert np.all(np.diff(rates) <= 1.0e-6)

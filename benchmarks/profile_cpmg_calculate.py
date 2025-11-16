"""Profile pulse_sequence.calculate to find the actual bottleneck."""

import sys
import time
from functools import reduce
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chemex.models.loader import register_kinetic_settings
from chemex.experiments.loader import register_experiments

register_kinetic_settings()
register_experiments()


def profile_cpmg_calculate():
    """Profile the components of CPMG pulse_sequence.calculate."""
    from chemex.experiments.builder import build_experiments
    from chemex.configuration.methods import Selection
    from chemex.parameters import database

    print("=" * 70)
    print("PROFILING pulse_sequence.calculate")
    print("=" * 70)

    example_dir = Path(__file__).parent.parent / "examples" / "Experiments" / "CPMG_15N_IP_0013"
    selection = Selection(include="*", exclude=None)

    experiments = build_experiments([
        example_dir / "Experiments" / "600mhz.toml",
    ], selection)

    profile = list(experiments)[0].profiles[0]
    params = database.build_lmfit_params(experiments.param_ids)
    profile.update_spectrometer(params)

    spectrometer = profile.spectrometer
    data = profile.data
    settings = profile.pulse_sequence.settings

    print(f"Profile: {profile.spin_system}")
    print(f"ncyc values: {sorted(set(data.metadata))}")
    print(f"Matrix size: {spectrometer.liouvillian.size}x{spectrometer.liouvillian.size}")

    # --- Profile individual operations ---

    ncycs = data.metadata
    tau_cps, deltas, all_delays = profile.pulse_sequence._get_delays(ncycs)

    print(f"\nNumber of unique delays: {len(all_delays)}")
    print(f"Number of unique ncyc: {len(set(ncycs))}")

    # 1. Profile delays calculation (expm)
    iterations = 100
    start = time.perf_counter()
    for _ in range(iterations):
        _ = spectrometer.delays(all_delays)
    elapsed = time.perf_counter() - start
    print(f"\nspectrometer.delays (expm): {elapsed/iterations*1000:.4f} ms")

    # 2. Profile p90_i and p180_i access
    start = time.perf_counter()
    for _ in range(iterations):
        _ = spectrometer.p90_i
        _ = spectrometer.p180_i
    elapsed = time.perf_counter() - start
    print(f"p90_i + p180_i access: {elapsed/iterations*1000:.4f} ms")

    # 3. Profile the loop over ncyc values
    delays = dict(zip(all_delays, spectrometer.delays(all_delays), strict=True))
    d_neg = delays[settings.t_neg]
    d_pos = delays[settings.t_pos]
    d_delta = {ncyc: delays[delay] for ncyc, delay in deltas.items()}
    d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}
    p90 = spectrometer.p90_i
    p180 = spectrometer.p180_i
    start_mag = spectrometer.get_start_magnetization(settings.start_terms)

    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in set(ncycs) - {0.0}:
            phases = profile.pulse_sequence._get_phases(ncyc)
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg = reduce(np.matmul, echo[phases.T])
            _ = spectrometer.detect(
                d_delta[ncyc] @ p90[3] @ d_neg @ cpmg @ d_neg @ p90[1] @ start_mag,
            )
    elapsed = time.perf_counter() - start
    print(f"ncyc loop (matmul + detect): {elapsed/iterations*1000:.4f} ms")

    # 4. Break down the ncyc loop further
    ncyc = list(set(ncycs) - {0.0})[0]  # Take one ncyc value
    phases = profile.pulse_sequence._get_phases(ncyc)

    start = time.perf_counter()
    for _ in range(1000):
        echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
    elapsed = time.perf_counter() - start
    print(f"\nDetailed breakdown (1000 iter):")
    print(f"  echo = d_cp @ p180 @ d_cp: {elapsed/1000*1000:.4f} ms")

    start = time.perf_counter()
    echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
    for _ in range(1000):
        cpmg = reduce(np.matmul, echo[phases.T])
    elapsed = time.perf_counter() - start
    print(f"  reduce(matmul, echo[phases.T]): {elapsed/1000*1000:.4f} ms")
    print(f"    phases.shape: {phases.shape}")
    print(f"    echo.shape: {echo.shape}")
    print(f"    Number of matmuls: {phases.shape[1]}")

    start = time.perf_counter()
    cpmg = reduce(np.matmul, echo[phases.T])
    for _ in range(1000):
        result = d_delta[ncyc] @ p90[3] @ d_neg @ cpmg @ d_neg @ p90[1] @ start_mag
    elapsed = time.perf_counter() - start
    print(f"  Final matmul chain: {elapsed/1000*1000:.4f} ms")

    start = time.perf_counter()
    for _ in range(1000):
        _ = spectrometer.detect(result)
    elapsed = time.perf_counter() - start
    print(f"  spectrometer.detect: {elapsed/1000*1000:.4f} ms")

    # 5. Alternative: Use matrix power instead of reduce
    print("\n--- Optimization Test: matrix_power vs reduce ---")
    echo_single = d_cp[ncyc] @ p180 @ d_cp[ncyc]
    echo_single = echo_single[0]  # Take first phase

    start = time.perf_counter()
    for _ in range(1000):
        _ = np.linalg.matrix_power(echo_single, 10)
    elapsed = time.perf_counter() - start
    print(f"matrix_power (n=10): {elapsed/1000*1000:.4f} ms")

    start = time.perf_counter()
    for _ in range(1000):
        result = echo_single
        for _ in range(9):
            result = result @ echo_single
    elapsed = time.perf_counter() - start
    print(f"explicit loop matmul (n=10): {elapsed/1000*1000:.4f} ms")


if __name__ == "__main__":
    profile_cpmg_calculate()

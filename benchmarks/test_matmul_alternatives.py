"""Test alternatives to reduce(np.matmul) for matrix chain multiplication."""

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


def test_matmul_alternatives():
    """Test different approaches for matrix chain multiplication."""
    from chemex.experiments.builder import build_experiments
    from chemex.configuration.methods import Selection
    from chemex.parameters import database

    print("=" * 70)
    print("MATRIX CHAIN MULTIPLICATION ALTERNATIVES")
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

    ncycs = data.metadata
    tau_cps, deltas, all_delays = profile.pulse_sequence._get_delays(ncycs)

    # Pre-compute delays
    delays = dict(zip(all_delays, spectrometer.delays(all_delays), strict=True))
    d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}
    p180 = spectrometer.p180_i

    unique_ncycs = sorted(set(ncycs) - {0.0})

    # Pre-compute echos and phases
    all_phases = {ncyc: profile.pulse_sequence._get_phases(ncyc) for ncyc in unique_ncycs}
    all_echos = {ncyc: d_cp[ncyc] @ p180 @ d_cp[ncyc] for ncyc in unique_ncycs}

    print(f"Matrix size: {all_echos[unique_ncycs[0]].shape[-2]}x{all_echos[unique_ncycs[0]].shape[-1]}")

    # Test different approaches for a few ncyc values
    test_ncycs = [2, 10, 20, 40]

    for ncyc in test_ncycs:
        if ncyc not in unique_ncycs:
            continue

        echo = all_echos[ncyc]
        phases = all_phases[ncyc]
        matrices = list(echo[phases.T])

        print(f"\nncyc={ncyc} ({len(matrices)} matrices, shape: {matrices[0].shape}):")

        # Method 1: reduce(np.matmul)
        iterations = 1000
        start = time.perf_counter()
        for _ in range(iterations):
            result1 = reduce(np.matmul, echo[phases.T])
        t1 = (time.perf_counter() - start) / iterations * 1000
        print(f"  reduce(np.matmul):     {t1:.4f} ms")

        # Method 2: np.linalg.multi_dot only works for 2D matrices
        # Since our matrices are >2D, we can't use it directly
        print(f"  np.linalg.multi_dot:   N/A (matrices are {matrices[0].ndim}D)")

        # Method 3: Explicit loop (baseline)
        start = time.perf_counter()
        for _ in range(iterations):
            result3 = matrices[0]
            for mat in matrices[1:]:
                result3 = result3 @ mat
        t3 = (time.perf_counter() - start) / iterations * 1000
        print(f"  Explicit loop:         {t3:.4f} ms (speedup: {t1/t3:.2f}x)")
        print(f"    Results match: {np.allclose(result1, result3)}")

        # Method 4: Pre-allocate intermediate results
        start = time.perf_counter()
        for _ in range(iterations):
            result4 = echo[phases.T[0]]
            for phase_idx in phases.T[1:]:
                result4 = result4 @ echo[phase_idx]
        t4 = (time.perf_counter() - start) / iterations * 1000
        print(f"  Direct indexing loop:  {t4:.4f} ms (speedup: {t1/t4:.2f}x)")
        print(f"    Results match: {np.allclose(result1, result4)}")

    # Total across all ncyc values
    print("\n" + "=" * 70)
    print("TOTAL TIME ACROSS ALL NCYC VALUES")
    print("=" * 70)

    iterations = 1000

    # Method 1: reduce
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            _ = reduce(np.matmul, all_echos[ncyc][all_phases[ncyc].T])
    t_reduce = (time.perf_counter() - start) / iterations * 1000
    print(f"reduce(np.matmul):     {t_reduce:.4f} ms")

    # Method 2: direct indexing loop
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            phases_t = all_phases[ncyc].T
            result = all_echos[ncyc][phases_t[0]]
            for phase_idx in phases_t[1:]:
                result = result @ all_echos[ncyc][phase_idx]
    t_direct = (time.perf_counter() - start) / iterations * 1000
    print(f"Direct indexing loop:  {t_direct:.4f} ms (speedup: {t_reduce/t_direct:.2f}x)")

    # Method 3: explicit loop
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            matrices = all_echos[ncyc][all_phases[ncyc].T]
            result = matrices[0]
            for mat in matrices[1:]:
                result = result @ mat
    t_loop = (time.perf_counter() - start) / iterations * 1000
    print(f"Explicit loop:         {t_loop:.4f} ms (speedup: {t_reduce/t_loop:.2f}x)")

    # Calculate impact
    print(f"\nFor 48 profiles × 2000 iterations:")
    time_saved = (t_reduce - min(t_direct, t_loop)) * 48 * 2000 / 1000
    print(f"  Best alternative saves: {time_saved:.1f} seconds")


if __name__ == "__main__":
    test_matmul_alternatives()

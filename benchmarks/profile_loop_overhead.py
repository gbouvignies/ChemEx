"""Profile Python loop overhead vs actual computation in CPMG calculate."""

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


def profile_loop_overhead():
    """Profile the overhead in the CPMG ncyc loop."""
    from chemex.experiments.builder import build_experiments
    from chemex.configuration.methods import Selection
    from chemex.parameters import database

    print("=" * 70)
    print("PROFILING LOOP OVERHEAD IN pulse_sequence.calculate")
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

    print(f"Number of unique ncyc values (excluding 0): {len(set(ncycs) - {0.0})}")

    # Pre-compute delays
    delays = dict(zip(all_delays, spectrometer.delays(all_delays), strict=True))
    d_neg = delays[settings.t_neg]
    d_delta = {ncyc: delays[delay] for ncyc, delay in deltas.items()}
    d_cp = {ncyc: delays[delay] for ncyc, delay in tau_cps.items()}
    p90 = spectrometer.p90_i
    p180 = spectrometer.p180_i
    start_mag = spectrometer.get_start_magnetization(settings.start_terms)

    unique_ncycs = sorted(set(ncycs) - {0.0})

    # 1. Profile entire loop as-is
    iterations = 1000
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            phases = profile.pulse_sequence._get_phases(ncyc)
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg = reduce(np.matmul, echo[phases.T])
            _ = spectrometer.detect(
                d_delta[ncyc] @ p90[3] @ d_neg @ cpmg @ d_neg @ p90[1] @ start_mag,
            )
    t_total = (time.perf_counter() - start) / iterations * 1000
    print(f"\nTotal ncyc loop time: {t_total:.4f} ms")

    # 2. Profile just the Python for-loop overhead
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            pass
    t_loop = (time.perf_counter() - start) / iterations * 1000
    print(f"Python for-loop overhead: {t_loop:.4f} ms ({t_loop/t_total*100:.1f}%)")

    # 3. Profile _get_phases calls
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            _ = profile.pulse_sequence._get_phases(ncyc)
    t_phases = (time.perf_counter() - start) / iterations * 1000
    print(f"_get_phases calls: {t_phases:.4f} ms ({t_phases/t_total*100:.1f}%)")

    # 4. Profile dictionary lookups
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            _ = d_cp[ncyc]
            _ = d_delta[ncyc]
    t_dict = (time.perf_counter() - start) / iterations * 1000
    print(f"Dictionary lookups: {t_dict:.4f} ms ({t_dict/t_total*100:.1f}%)")

    # 5. Profile echo computation
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            _ = d_cp[ncyc] @ p180 @ d_cp[ncyc]
    t_echo = (time.perf_counter() - start) / iterations * 1000
    print(f"Echo computation: {t_echo:.4f} ms ({t_echo/t_total*100:.1f}%)")

    # 6. Profile reduce operation
    echos = {ncyc: d_cp[ncyc] @ p180 @ d_cp[ncyc] for ncyc in unique_ncycs}
    phases_dict = {ncyc: profile.pulse_sequence._get_phases(ncyc) for ncyc in unique_ncycs}

    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            _ = reduce(np.matmul, echos[ncyc][phases_dict[ncyc].T])
    t_reduce = (time.perf_counter() - start) / iterations * 1000
    print(f"reduce(matmul) operation: {t_reduce:.4f} ms ({t_reduce/t_total*100:.1f}%)")

    # 7. Profile final matmul chain
    cpmgs = {ncyc: reduce(np.matmul, echos[ncyc][phases_dict[ncyc].T]) for ncyc in unique_ncycs}
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            _ = d_delta[ncyc] @ p90[3] @ d_neg @ cpmgs[ncyc] @ d_neg @ p90[1] @ start_mag
    t_chain = (time.perf_counter() - start) / iterations * 1000
    print(f"Final matmul chain: {t_chain:.4f} ms ({t_chain/t_total*100:.1f}%)")

    # 8. Profile detect operation
    mags = {ncyc: d_delta[ncyc] @ p90[3] @ d_neg @ cpmgs[ncyc] @ d_neg @ p90[1] @ start_mag for ncyc in unique_ncycs}
    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            _ = spectrometer.detect(mags[ncyc])
    t_detect = (time.perf_counter() - start) / iterations * 1000
    print(f"spectrometer.detect: {t_detect:.4f} ms ({t_detect/t_total*100:.1f}%)")

    # Calculate overhead vs computation
    computation = t_echo + t_reduce + t_chain + t_detect
    overhead = t_loop + t_phases + t_dict
    print(f"\nSummary:")
    print(f"  Pure computation: {computation:.4f} ms ({computation/t_total*100:.1f}%)")
    print(f"  Python overhead: {overhead:.4f} ms ({overhead/t_total*100:.1f}%)")

    # 9. Test vectorized approaches
    print("\n" + "=" * 70)
    print("OPTIMIZATION TEST: Vectorize ncyc loop")
    print("=" * 70)

    # Approach 1: Pre-compute everything
    start = time.perf_counter()
    for _ in range(iterations):
        # Pre-compute all phases
        all_phases = {ncyc: profile.pulse_sequence._get_phases(ncyc) for ncyc in unique_ncycs}
        # Pre-compute all echos
        all_echos = {ncyc: d_cp[ncyc] @ p180 @ d_cp[ncyc] for ncyc in unique_ncycs}
        # Compute cpmg matrices
        all_cpmgs = {ncyc: reduce(np.matmul, all_echos[ncyc][all_phases[ncyc].T]) for ncyc in unique_ncycs}
        # Final computation
        for ncyc in unique_ncycs:
            _ = spectrometer.detect(
                d_delta[ncyc] @ p90[3] @ d_neg @ all_cpmgs[ncyc] @ d_neg @ p90[1] @ start_mag,
            )
    t_precompute = (time.perf_counter() - start) / iterations * 1000
    print(f"Pre-compute approach: {t_precompute:.4f} ms (speedup: {t_total/t_precompute:.2f}x)")

    # Approach 2: Batch the final detection
    start = time.perf_counter()
    for _ in range(iterations):
        # Compute all magnetization vectors
        mags = []
        for ncyc in unique_ncycs:
            phases = profile.pulse_sequence._get_phases(ncyc)
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg = reduce(np.matmul, echo[phases.T])
            mag = d_delta[ncyc] @ p90[3] @ d_neg @ cpmg @ d_neg @ p90[1] @ start_mag
            mags.append(mag)
        # Batch detect (if possible)
        # For now, still loop (detect expects single vector)
        intensities = [spectrometer.detect(mag) for mag in mags]
    t_batch = (time.perf_counter() - start) / iterations * 1000
    print(f"Batch computation: {t_batch:.4f} ms (speedup: {t_total/t_batch:.2f}x)")

    # Approach 3: Cache _get_phases computation
    # The phases for a given ncyc are always the same
    print("\nCaching _get_phases (since it's deterministic):")
    cached_phases = {ncyc: profile.pulse_sequence._get_phases(ncyc) for ncyc in unique_ncycs}

    start = time.perf_counter()
    for _ in range(iterations):
        for ncyc in unique_ncycs:
            phases = cached_phases[ncyc]  # Direct lookup instead of computation
            echo = d_cp[ncyc] @ p180 @ d_cp[ncyc]
            cpmg = reduce(np.matmul, echo[phases.T])
            _ = spectrometer.detect(
                d_delta[ncyc] @ p90[3] @ d_neg @ cpmg @ d_neg @ p90[1] @ start_mag,
            )
    t_cached = (time.perf_counter() - start) / iterations * 1000
    print(f"With cached phases: {t_cached:.4f} ms (speedup: {t_total/t_cached:.2f}x)")

    # Check the impact
    time_saved = (t_total - t_cached) * 48 * 2000 / 1000  # seconds
    print(f"\nFor 48 profiles × 2000 iterations:")
    print(f"  Potential time saved: {time_saved:.1f} seconds")


if __name__ == "__main__":
    profile_loop_overhead()

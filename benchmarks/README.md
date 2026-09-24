# ChemEx benchmarks

These exploratory benchmarks measure three current ChemEx workloads. They are not scientific acceptance tests or CI performance gates. Run them from the repository root after `uv sync --locked`.

Set native thread counts **before** starting Python so NumPy and SciPy see the same settings in every run:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/chemex-benchmark-mpl"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/chemex-benchmark-cache"
export UV_CACHE_DIR="${TMPDIR:-/tmp}/chemex-benchmark-uv-cache"
```

| Command | Input and measurement |
| --- | --- |
| `uv run --no-sync python benchmarks/run_benchmarks.py cpmg-residual --runs 5` | `CPMG_15N_IP_0013`, 600 MHz, L3N, 2st: native residual evaluation. Reports a cold evaluation with a fresh evaluator and a second cached evaluation; excludes experiment construction. |
| `uv run --no-sync python benchmarks/run_benchmarks.py dcest-residual --runs 5` | `DCEST_15N_HD_EXCH`, 3 Hz, 1N, 2st_hd: the same native residual measurements. |
| `uv run --no-sync python benchmarks/run_benchmarks.py cpmg-fit --runs 3 --timeout 180` | `CPMG_15N_IP_0013`, 600 MHz, L3N, 2st: complete CLI fit in a child process, with plotting disabled, one worker, and one native thread. Wall time includes Python startup, fitting, uncertainty, and output. Each run uses a fresh temporary output directory and must publish a complete outcome. |

The runner checks that both its own process and the fit subprocess import `src/chemex` from this checkout. It prints input paths, model, spin system, timing medians, runtime versions, source path, Git revision, SHA-256 fingerprints of the source Python tree and runner, platform, CPU count, and thread environment. The fingerprints identify the measured code even with uncommitted edits or stale installed distribution metadata; the runner fails if the source or runner changes during a measurement. Residual workloads also print the result length and a SHA-256 fingerprint to help detect an unintended workload change. Compare timings only with the same inputs, thread settings, cache state, Python and library versions, and machine. Use more repetitions when investigating small differences; the short commands above are primarily for baseline and smoke measurements.

## Initial local measurements

These are exploratory measurements from **2026-09-24** at revision `16f1bc697d98ab362301d5035c759ca2ce04d43b` on macOS arm64 (18 logical CPUs), Python 3.13.13, NumPy 2.5.1, and SciPy 1.18.0. All five native thread variables in the commands above were set to 1; Matplotlib and XDG caches were in writable temporary directories. The imported `chemex` source was this checkout, version 2026.09.2 in `pyproject.toml`; the existing environment's installed distribution metadata reported 2026.9.1. This mismatch is recorded, not treated as a numerical difference.

| Workload | Repeats | Median cold | Median cached | Median full fit |
| --- | ---: | ---: | ---: | ---: |
| CPMG residual, 1 profile / 29 residuals | 5 | 0.926 ms | 0.0118 ms | — |
| D-CEST residual, 1 profile / 42 residuals | 5 | 3.567 ms | 0.0125 ms | — |
| CPMG fit, 1 profile | 3 | — | — | 1.094 s |

These numbers establish that the retained commands measure working scientific paths. They are not performance thresholds or evidence of a speedup.

The retired `benchmark_large_native_uncertainty.py` attempted a 96-profile CPMG fit with phase timings and a numerical acceptance record. Its instrumentation imports removed optimization modules, so it cannot validate the current implementation. The retained single-profile fit does not replace a full-size covariance or uncertainty performance measurement. Add a focused workload only if that path needs a future performance investigation; use the existing uncertainty tests for current correctness checks.

## Historical results

[`baseline_results.txt`](baseline_results.txt) and [`optimized_results.txt`](optimized_results.txt) are preserved snapshots from **2025-11-14**. They came from the retired experimental suite, include synthetic matrix comparisons and failed imports, and lack enough environment information for a reliable comparison with current ChemEx. Their names describe the historical experiment; neither file is a baseline or an optimization result for the current implementation. The old scripts were retired because they depended on removed modules or duplicated exploratory comparisons that did not measure the current native fitting path. No speedup claimed in those files is assumed here.

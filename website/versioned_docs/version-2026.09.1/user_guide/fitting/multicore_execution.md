---
sidebar_position: 8
---

# Multicore Execution

ChemEx can use multiple CPU cores during the statistics phase of `chemex fit`.
This is most useful for MCMC sampling and for Monte Carlo or bootstrap
uncertainty estimates.

The default command already enables the recommended modern behavior:

```shell
chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -o Output
```

By default, `--workers auto` chooses a conservative number of worker slots, up
to 8 CPUs. This gives typical workstations useful parallelism without starting
an unexpectedly large worker pool. Use `--workers 1` for serial execution, or
`--workers 0` to explicitly use all CPUs visible to ChemEx.

## What Runs in Parallel

The `--workers` option applies to fit statistics:

-   MCMC sampling requested with `STATISTICS = {"MCMC" = ...}`.
-   Monte Carlo, bootstrap, and nucleus-specific bootstrap refits requested with
    `MC`, `BS`, or `BSN`.

The normal deterministic fit still runs as one optimization task. Grid searches,
simulation runs, plotting, and output writing are not controlled by
`--workers`.

Workers are active only while the parallel statistics task is running. Native
MCMC uses worker processes with one isolated native evaluator per process, so
CPU-bound likelihood calculations can execute on separate CPU cores. Native MC,
BS, and BSN use worker threads with isolated native evaluators. The earlier
deterministic fit and later output-writing phase remain serial.

## Command-Line Controls

### `--workers N|auto`

Controls the number of ChemEx worker slots used by fit statistics.

| Value | Meaning |
| ---- | ------- |
| `auto` | Conservative default, capped at 8 workers. |
| `1` | Serial execution. Useful for debugging and reproducibility checks. |
| `N` | Use `N` worker slots for native statistics. |
| `0` | Use all CPUs visible to the current process. |

For long MCMC runs on a dedicated machine, it can be reasonable to set an
explicit value:

```shell
chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method_stat.toml \
           -o OutputStat \
           --workers 10
```

### `--native-threads N|auto`

Controls native numerical library threads, such as BLAS or OpenMP threads.

The default, `--native-threads auto`, leaves native thread settings untouched for
serial runs. When ChemEx starts multiple workers, it sets native numerical
threads to 1 before starting the worker pool; MCMC worker processes inherit that
setting. This avoids oversubscription, where each ChemEx worker also starts many
BLAS or OpenMP threads.

Most users should keep the default. Use an explicit value only when you are
benchmarking a specific machine:

```shell
chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method_stat.toml \
           -o OutputStat \
           --workers 10 \
           --native-threads 1
```

## Method files and workers

Canonical version 2 method files do not contain execution settings. MCMC and
resampling use the command-line `--workers` value; walker topology and native
thread coordination remain ChemEx policy. Method-local MCMC `WORKERS` is
deprecated v1-only syntax and follows the v1 removal window described in the
[Method Files](method_files.md#version-1-compatibility) guide.

## Practical Guidance

Start with the defaults. They are designed to give useful multicore performance
on modern machines without requiring manual tuning.

Use explicit settings when you have a reason:

-   `--workers 1` for serial debugging.
-   `--workers 0` for a dedicated machine where ChemEx may use all CPUs.
-   `--workers N --native-threads 1` for manual benchmarking of long MCMC runs.

Do not expect every fit to scale linearly. Parallel execution helps most when
each likelihood evaluation or refit is expensive enough to dominate process-pool
overhead. Very short MCMC runs or small bootstrap jobs may show little speedup.

Interactive MCMC runs show completed ensemble transitions against the requested
step count, together with elapsed time and an estimated time remaining. One
completed transition advances the display by one step, independent of the
number of walkers. Redirected or otherwise non-interactive output reports a
concise start and terminal line instead of printing one line per step.

## Diagnostics

Statistics diagnostics record the effective worker count. MCMC diagnostics also
record the direct emcee sampler engine, timing information, acceptance
fractions, and autocorrelation diagnostics.

`sampling_seconds` covers sampler execution, worker-pool startup and shutdown,
and structural capture qualification. Posterior selection, summary, and output
timings are recorded separately. Capture qualification checks topology, content
identities, bounds, and evidence lineage; it does not re-evaluate every stored
likelihood after the authoritative native sampler has already evaluated it.

Look at `Statistics/MCMC/diagnostics.toml` after a run to confirm the effective
settings:

```toml
sampler = "emcee via ChemEx direct EnsembleSampler"
engine = "native MCMC"
workers = 8
root_seed = 1234
sampling_seconds = 120.42
```

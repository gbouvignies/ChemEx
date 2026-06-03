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

By default, `--workers auto` chooses a conservative number of worker processes,
up to 8 CPUs. This gives typical workstations useful parallelism without
starting an unexpectedly large number of Python processes. Use `--workers 1` for
serial execution, or `--workers 0` to explicitly use all CPUs visible to ChemEx.

## What Runs in Parallel

The `--workers` option applies to fit statistics:

-   MCMC sampling requested with `STATISTICS = {"MCMC" = ...}`.
-   Monte Carlo, bootstrap, and nucleus-specific bootstrap refits requested with
    `MC`, `BS`, or `BSN`.

The normal deterministic fit still runs as one optimization task. Grid searches,
simulation runs, plotting, and output writing are not controlled by
`--workers`.

Worker processes are created only while the parallel statistics task is running.
For example, Activity Monitor or `top` should show multiple Python processes
during MCMC sampling, but not necessarily during the earlier deterministic fit
or the later output-writing phase.

## Command-Line Controls

### `--workers N|auto`

Controls the number of ChemEx worker processes used by fit statistics.

| Value | Meaning |
| ---- | ------- |
| `auto` | Conservative default, capped at 8 workers. |
| `1` | Serial execution. Useful for debugging and reproducibility checks. |
| `N` | Use `N` worker processes. |
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
serial runs. When ChemEx starts multiple worker processes, it sets native
threads to 1 inside the worker-pool context. This avoids oversubscription, where
each Python worker also starts many native threads.

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

## MCMC Worker Overrides

MCMC settings in a method file may include an optional `WORKERS` value:

```toml
[STEP1.STATISTICS.MCMC]
STEPS = 5000
WORKERS = 4
```

When `WORKERS` is omitted, MCMC inherits the command-line `--workers` setting.
Use the method-file value only when a specific MCMC step needs a different worker
count from the rest of the fit statistics.

Method-file `WORKERS` must be a positive integer. The special `0` and `auto`
values are command-line controls only.

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

## Diagnostics

Statistics diagnostics record the effective worker count. MCMC diagnostics also
record the direct emcee sampler engine, timing information, acceptance
fractions, and autocorrelation diagnostics.

Look at `Statistics/MCMC/diagnostics.toml` after a run to confirm the effective
settings:

```toml
sampler = "emcee via ChemEx direct EnsembleSampler"
workers = 8
sampling_seconds = 120.42
```

---
sidebar_position: 6
---

# Method Files

Method files describe which profiles and parameter roles a fit uses, whether a
deterministic search precedes the local fit, and which statistical analyses run
afterward. Pass one or more method files with `chemex fit -m`.

Version 2 is the canonical format for new method files. Every v2 file starts
with:

```toml title="method.toml"
FORMAT_VERSION = 2

[DEFAULT]
```

The empty step performs one bounded trust-region reflective (TRF) fit using the
parameter model's baseline roles and all profiles. V2 has no `FITMETHOD`; when
neither `SEARCH.GRID` nor `SEARCH.DE` is present, full-coordinate TRF is
implicit.

## The v2 mental model

Each step has six straightforward rules:

1. It starts from baseline parameter roles unless `ROLES_FROM` explicitly names
   an earlier step.
2. `ROLES` actions are applied from top to bottom.
3. A later matching action replaces the parameter's complete earlier role.
4. Numerical starting values automatically come from the latest successfully
   committed fit, independently of `ROLES_FROM`.
5. `SEARCH.GRID` or `SEARCH.DE` chooses an optional deterministic search;
   otherwise ChemEx runs the implicit full TRF directly.
6. `STATISTICS` analyzes the committed deterministic result and never changes
   its central parameter values.

Steps execute in declaration order. Results are written under the step name for
multi-step methods.

## Ordered parameter roles

`ROLES` is an ordered array. Every entry contains exactly one action:

- `FIT` makes matching user-controlled parameters independent fit coordinates.
- `FIX` holds matching user-controlled parameters at their current values.
- `CONSTRAIN` derives each target from its equation.

```toml
FORMAT_VERSION = 2

[STEP]
ROLES = [
  { FIX = ["PB", "KEX_AB"] },
  { FIT = ["DW_AB"] },
]
```

Later actions are intentional overrides, not conflicts. `FIT` or `FIX` after a
constraint removes that constraint; `CONSTRAIN` after `FIT` makes the target
derived.

### Broad rule followed by specific exceptions

The following example fixes every `DW_AB`, fits the G24 value, and constrains
the G31 value to the fitted G24 value. Every other `DW_AB` remains fixed:

```toml
FORMAT_VERSION = 2

[STEP]
ROLES = [
  { FIX = ["DW_AB"] },
  { FIT = ["DW_AB, NUC->G24N-H"] },
  { CONSTRAIN = [
    "[DW_AB, NUC->G31N-H] = [DW_AB, NUC->G24N-H]",
  ] },
]
```

There is no separate precedence rule: the visible list order is the complete
precedence rule.

### Selectors and constraints

Selectors start with a parameter name and may add qualifiers. V2 consumes the
complete selector and rejects unknown, duplicated, or trailing fragments.

```toml
ROLES = [
  { FIT = [
    "R2_A, NUC->G23N, B0->800.13MHz, T->23C",
    "R1_A, B0->800.13MHz",
    "DW_AB, NUC->N",
  ] },
  { CONSTRAIN = [
    "[R1_B] = 0.5 * [R1_A]",
    "[R2_B, NUC->N] = [R2_A, NUC->N]",
  ] },
]
```

Supported qualifiers are `NUC`, `T`, `B0`, `[P]`, `[L]`, and `D2O`. Constraint
right-hand sides use finite numbers, bracketed parameter references,
parentheses, and `+`, `-`, `*`, and `/`.

Model-owned derived parameters remain protected: method actions cannot turn a
structurally derived parameter into an independent fit coordinate.

## Multi-step methods

Role inheritance and numerical value continuity are separate.

```toml
FORMAT_VERSION = 2

[GLOBAL_FIT]
INCLUDE = [15, 31, 33, 34, 37]
ROLES = [{ FIT = ["PB", "KEX_AB"] }]

[LOCAL_FITS]
INCLUDE = "ALL"
ROLES_FROM = "GLOBAL_FIT"
ROLES = [
  { FIX = ["PB", "KEX_AB"] },
  { FIT = ["DW_AB"] },
]
```

`ROLES_FROM = "GLOBAL_FIT"` copies only that earlier step's effective
role/constraint setup. It does not copy values, profile selection, search,
statistics, seeds, or execution settings. It must name one unique earlier step.

Whether or not `ROLES_FROM` is present, `LOCAL_FITS` starts numerically from the
latest values committed by `GLOBAL_FIT`. Omitting `ROLES_FROM` starts from
baseline roles, **not** baseline numerical values.

`INCLUDE` and `EXCLUDE` are step-local. Omitting both selects all profiles for
that step; selection never inherits. They accept spin-system or group names,
residue numbers, and `"ALL"`/`"*"`:

```toml
INCLUDE = ["G2", "A4", "C5", "H6"]
EXCLUDE = [12, 18]
```

## Choosing the deterministic search

| Method structure | Use it when | What ChemEx does |
| --- | --- | --- |
| No `SEARCH` | The committed start is suitable | Runs one full-coordinate TRF |
| `SEARCH.GRID` | A profiled chi-square landscape is scientifically useful | Holds the GRID coordinates at each point, optimizes the other fitted coordinates, and commits one coherent joint grid solution |
| `SEARCH.DE` | A few basin-defining coordinates need broad stochastic coverage | Searches only those coordinates, then seeds exactly one normal full-coordinate TRF |

GRID and DE never change parameter roles. Their targets must resolve to final
independent `FIT` coordinates.

### GRID

```toml
FORMAT_VERSION = 2

[STEP]
ROLES = [{ FIT = ["PB", "KEX_AB", "DW_AB"] }]

[STEP.SEARCH.GRID]
AXES = [
  "[PB] = values(0.01, 0.02, 0.05, 0.10, 0.20)",
  "[KEX_AB] = log(100.0, 5000.0, 12)",
  "[DW_AB] = lin(0.0, 10.0, 5)",
]
```

GRID accepts only `lin(low, high, count)`, `log(low, high, count)`, and
`values(v1, v2, ...)`. Broad axes may expand to several final fitted
coordinates. Axis declarations apply top-to-bottom, so a later specific axis
replaces an earlier broad axis for the same concrete coordinate.

GRID calculates a profiled chi-square surface. At each grid point, the declared
GRID coordinates are held exactly at their declared values while every other
final independent `FIT` coordinate is optimized by the normal native TRF
solver. Derived and fixed parameters keep their normal constraint and role
semantics. When no non-GRID fitted coordinate remains, ChemEx evaluates the
point directly without starting a zero-variable optimizer.

Broad selectors are resolved only against the current step's active final
`FIT` scope. Exact objective factorization avoids Cartesian products of
unrelated local axes. For example, factors that share global `KEX_AB` and `PB`
but each own a residue-local `DW_AB` are profiled separately and then reduced
and summed exactly. ChemEx selects one coherent joint assignment: one shared
global grid point, its corresponding best local point in every factor, and the
nuisance-parameter solutions attached to those points. The complete root state
is freshly evaluated and accepted before one atomic commit. Factor results and
independent one-dimensional marginal minima never commit.

If the feasible domain contains an affine restriction, ChemEx conservatively
uses one exact unfactorized profiled grid. This can be slower, but prevents a
cross-coordinate domain restriction from being mistaken for independent
factor objectives.

GRID does not release its own axes through an automatic final TRF. A later
method step can do so explicitly with its normal roles, starting from the
committed GRID state. Because a discrete GRID coordinate was not continuously
optimized in that step, ChemEx withholds deterministic covariance errors for a
GRID step. Requested MC, BS, BSN, or MCMC analyses still start from the accepted
GRID state; a following ordinary Direct step is the normal route to local
covariance uncertainty.

### Selected-coordinate DE

```toml
FORMAT_VERSION = 2

[STEP]
ROLES = [{ FIT = ["PB", "KEX_AB", "DW_AB"] }]

[STEP.SEARCH.DE]
SEED = 597
COORDINATES = [
  "[PB] = log(0.001, 0.20)",
  "[KEX_AB] = log(100.0, 5000.0)",
]
```

Each DE coordinate must resolve to exactly one final fitted coordinate and use
`lin(low, high)` or `log(low, high)`. Ranges must be finite, ordered, within the
physical parameter bounds, and positive for logarithmic searches. `SEED` is a
required unsigned 64-bit integer.

During DE, other fitted coordinates such as `DW_AB` remain at the committed
step-start values. The best eligible DE candidate initializes one fresh TRF
that releases the complete final `FIT` set. Only that TRF can be accepted,
committed, or used for uncertainty and statistics. DE failure does not fall back
to direct TRF.

Use GRID to inspect exact profiled chi-square landscapes. Use DE when a small
set of coordinates needs stochastic basin exploration before one complete
continuous fit. See the shipped
three-state DCEST example at
`examples/Experiments/DCEST_15N_3States/Methods/method_de.toml`.

## Statistical analyses

Statistics run once after successful deterministic acceptance and commit. They
do not mutate the committed central fit.

The compact form requests MC, bootstrap, nucleus-specific bootstrap, and MCMC:

```toml
FORMAT_VERSION = 2

[STEP]
STATISTICS = { MC = 100, BS = 100, BSN = 100, MCMC = 5000 }
```

Expanded MC, BS, and BSN requests add an optional reproducibility seed:

```toml
[STEP.STATISTICS.MC]
REPLICATES = 100
SEED = 680

[STEP.STATISTICS.BS]
REPLICATES = 100
SEED = 681

[STEP.STATISTICS.BSN]
REPLICATES = 100
SEED = 682
```

When a seed is omitted, ChemEx resolves one unsigned 64-bit root seed and records
it before the analysis starts.

Expanded MCMC exposes only `STEPS`, optional integer `BURN`, and optional
unsigned 64-bit `SEED`:

```toml
[STEP.STATISTICS.MCMC]
STEPS = 20000
BURN = 4000
SEED = 680
```

`STEPS` counts raw post-initialization ensemble iterations. Omitting `BURN`
invokes ChemEx's automatic burn/convergence policy. If that policy cannot
establish a defensible retained window, ChemEx preserves raw or partial chain
evidence and diagnostics but publishes no authoritative posterior summary; it
does not silently fall back to zero burn. A numeric `BURN` fixes the discard
window but does not suppress convergence, autocorrelation, ESS, or MCSE
diagnostics.

Walker topology, initialization, proposal policy, and parallel execution are
ChemEx policy. Configure parallelism globally with `chemex fit --workers`; v2
has no method-local `WORKERS`, `WALKERS`, `THIN`, `POLICY`, `MODE`, backend
options, or parameter-update behavior.

MC, BS, and BSN write complete summaries only if every requested replicate
succeeds. MCMC likewise withholds authoritative summary products from an
incomplete chain. Successful evidence and diagnostics remain available, and no
statistical analysis replaces deterministic fitted values or covariance errors.

## Version 1 compatibility

Version 1 is supported only for the frozen deprecation window; do not use it for
new method files. If v2 first ships in non-patch release **R**, v1 remains
supported but deprecated in **R** and the following non-patch release **R+1**,
then is removed in **R+2**. Patch releases never remove v1. Calendar versions
are assigned during release planning.

The same boundary applies to all v1-only behavior: `FITMETHOD`, the
`least_squares` alias, bare GRID tuples, permissive legacy selectors, implicit
cross-step role/profile carry-over, `BURN = "AUTO"`, `THIN`, explicit
`WALKERS`, method-local MCMC `WORKERS`, and other legacy MCMC controls.

During that window, omitting `FORMAT_VERSION` selects the frozen v1 adapter. V1
continues to apply role buckets in `CONSTRAINTS` → `FIX` → `FIT` order and to
carry roles and profile selection implicitly between steps. Mixed v1/v2 input
in one invocation is rejected. Duplicate v2 step names are errors and are never
partially merged.

### V1-to-v2 role equivalents

This v1 step:

```toml
[STEP]
CONSTRAINTS = ["[R2_B] = [R2_A]"]
FIX = ["PB"]
FIT = ["KEX_AB"]
```

has the following v2 equivalent because the ordered actions preserve v1's
bucket precedence:

```toml
FORMAT_VERSION = 2

[STEP]
ROLES = [
  { CONSTRAIN = ["[R2_B] = [R2_A]"] },
  { FIX = ["PB"] },
  { FIT = ["KEX_AB"] },
]
```

For later steps, add `ROLES_FROM = "PREVIOUS_STEP"` only when the old file
relied on role/constraint carry-over. Repeat `INCLUDE` or `EXCLUDE` when the old
file relied on selection carry-over. Do not add either merely to continue fitted
values; committed numerical values continue automatically.

Move a v1 `GRID = [...]` list under `[STEP.SEARCH.GRID]` as `AXES = [...]`, and
replace bare tuples such as `"[PB] = (0.01, 0.05, 0.10)"` with
`"[PB] = values(0.01, 0.05, 0.10)"`. Remove `FITMETHOD`; implicit TRF is the
only v2 local fit.

## Outputs

GRID results are written under `Grid/`. Selected-coordinate DE is an initializer
and creates no separate result tree. MC, BS, BSN, and MCMC outputs are written
under their corresponding `Statistics/` subdirectories. See
[Outputs](outputs.mdx) for the complete layout and diagnostics.

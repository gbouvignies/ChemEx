# Issue #712: basin discovery on current ChemEx

This is the retained historical research report. The production follow-up and
request-budget qualification are in [profile_preparation.md](profile_preparation.md).
The original DE implementation and runnable comparison are preserved at research
commit `88d8a5bd`; current production retires DE.

Research date: 2026-10-01. This is exploratory evidence and a prototype, **not a
production search feature, numerical oracle, or performance gate**.

## Recommendation

Do not make profiled Sobol sampling the routine next step. These measurements
support a smaller intervention: optional, exact factor-local nuisance profiling
at the existing coupling start, with targeted DW branch trials, followed by the
unchanged complete native TRF. Retain ordinary TRF as the default. Use explicit
outer sampling only when that simpler procedure demonstrably misses useful
basins. The extra Sobol points did not improve the best final endpoint in any
of the three measured problems.

Profiling helps CEST, branch trials help CPMG, and ordinary TRF wins the DCEST
example. The combined profiled sampling/branch strategy is **not generally
superior** to simpler approaches. In particular, a two-start full TRF is already
sufficient for this CPMG scope. Independent factor trials remain preferable to
a root Cartesian product when different residues need different branches.

## Current-main premise and sources

Remote `main` was verified with `git ls-remote` and fetched at
`d0ba6e344f007f1fb6caaf0e3e29b88322f15455`. The measured checkout is
`92af4237464f98d5943dc677b628592ff3096c9f`. The optimization, configuration,
evaluation, pulse sequences, and the `2st`/`3st` models used here match fetched
main. The only Python-source difference is main's removal of the unused
`2st_rs` model; its documentation/tests are also newer. No branch switch or
production-source edit was needed.

The historical motivating failure **does not reproduce**. The existing CEST
qualification passed: adaptive-Jacobian TRF with finite safety bounds reaches
the L18CD1 basin at chi-square **801.9787003**, rather than 1093.4222258, from
the shipped workflow. The test also evaluates the historical higher-basin
vector on the same objective. Its tolerances are 0.002 for fitted chi-square
and 0.001 for the reference vectors. The complete workflow chi-square is
2879.58; the 801.98 number is one residue component, not the complete fit.
See [the qualification](../tests/test_native_trf_scaling.py) and
[the still-historical issue text](https://github.com/gbouvignies/ChemEx/issues/712).

The prototype reuses these current authorities:

- [direct TRF](../src/chemex/optimize/direct_trf.py): bounded local solver,
  native request budgets, success-only candidates, and fresh materialization;
- [grouped TRF](../src/chemex/optimize/grouped_direct_trf.py): exact components
  and validated reconstruction of the complete root;
- [profiled GRID](../src/chemex/optimize/profiled_grid.py): exact dependency
  factors after holding outer coordinates, factor-local TRF, and fresh
  factor/root objective validation;
- [selected-coordinate DE](https://github.com/gbouvignies/ChemEx/blob/88d8a5bdcb0afef754ebef4f2afd431e867bd4fe/src/chemex/optimize/de_direct_trf.py) and
  [production dispatch](../src/chemex/optimize/native_deterministic.py): declared
  coordinates vary, nuisance coordinates remain at step-start values, and the
  best search vector initializes complete TRF;
- [Method compiler](../src/chemex/optimize/method_compiler.py) and
  [search definitions](../src/chemex/configuration/method_plan.py): the shipped
  GRID/DE scopes and declarations, without new Method syntax;
- [native evaluator](../src/chemex/evaluation/native.py): profile caches, which
  already make complete TRF cheaper than a naive dense-model count suggests.

## Prototype and experiment design

[The runner](basin_discovery.py) uses **singleton GRID** for arbitrary outer
points. No GRID redesign or extraction was necessary to obtain evidence.
A temporary benchmark-only wrapper around `_fit_factor_point` tries independent
factor-local sign starts and retains the best successful local endpoint. GRID
then reconstructs and validates the complete vector. Its internal accepted
GRID record is only a search state here: it is never committed or treated as a
final fit. Every reported final chi-square comes from an accepted complete
grouped/native TRF, with the ordinary FIT scope and scaling policy restored.
No live session values are changed by the research runner.

Alternative nuisance starts are built through the native problem constructor
using detached research snapshots. This matters because `restart_from()` updates
the full TRF start but preserves `independent_items`, while GRID derives its
nuisance starts from those items. Using `restart_from()` alone would silently
benchmark the old nuisance start. Production work should introduce a small,
validated nuisance-start argument to the existing point derivation, rather than
copy the research snapshot technique into a live fit.

Outer designs use scrambled SciPy Sobol with explicit seeds 597/598 and complete
power-of-two designs, plus the supplied starting point kept separately.
This follows [SciPy's Sobol guidance](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html).
Full multistart uses the **same outer starts**, leaves initial nuisance values
at their supplied values, and optimizes the entire normal FIT vector. A second
CPMG comparison adds explicit signed starts so that profiling is not given an
exclusive advantage from knowing about DW branches.

| Scope | Profiles / observations | Full FIT dimension | Outer dimension | Independent nuisance factors | Outer ranges |
| --- | ---: | ---: | ---: | --- | --- |
| CPMG_15N_IP, shipped GRID STEP1 selection | 10 / 230 | 17 | 2 | five factors of 3 coordinates | shipped log KEX 100–600, PB 0.03–0.15 |
| CEST_13C_LABEL_CN, L18CD1/L18CB/K25CA | 4 / 422 | 14 | 2 | three factors of 4 coordinates | exploratory log KEX 100–600, PB 0.03–0.15 |
| DCEST_15N_3States, shipped DE K19N-H selection | 3 / 289 | 10 | 4 | one factor of 6 coordinates | shipped PB/PC 0.001–0.2, KEX_AB/AC 10–5000, all log |

CEST deliberately uses a small research scope that includes the historically
ambiguous residue. It is **not** the shipped eight-residue STEP1 or the complete
15-residue STEP2. Its R2_B=R2_A constraint and supplied values are retained.
The CEST coupling box is explicitly exploratory, borrowed from shipped CPMG
GRID around the CEST start; it is not an established optimum-containing range.
DCEST DE additionally searches both DW coordinates in the shipped linear
range -15 to +15, holding the four fitted relaxation coordinates fixed.
CPMG/CEST DE use the same selected coupling boxes as their Sobol comparisons;
only DCEST has a shipped DE declaration.

**Search ranges are not final fit bounds.** The best CEST endpoint has
KEX≈23.30 and PB≈0.207, outside its outer box. The best DCEST endpoint has
KEX_AB≈5.23 and KEX_AC effectively zero, below its declared DE search ranges.
Releasing the full FIT scope is essential. These results establish attraction
basins, not global optimality or the adequacy of the boxes for exhaustive
discovery. No statistically confident sign assignment follows from a small
chi-square difference alone.

## Measured comparisons

CPMG uses seed 598 below; seed 597 agrees on the best endpoints. CEST/DCEST use
seed 597. CPMG/CEST have nine points (start + eight Sobol); DCEST has five
(start + four Sobol). Profiling rows include **one final full TRF**, launched
from the best profiled candidate. Multistart includes nine/five complete TRFs.
Branch-aware profiling uses 90/54/20 factor-local TRFs, respectively; ordinary
profiling uses 45/27/5. Additional polishes performed to compare selection rules
are excluded from each strategy's cost and retained separately in the traces.

| Problem | Strategy | Best converged final chi-square | Objective requests | Actual profile-kernel calculations | Seconds |
| --- | --- | ---: | ---: | ---: | ---: |
| CPMG | Normal TRF | 434.55519 | 126 | 370 | 0.18 |
| CPMG | Selected DE + TRF | 434.55519 | 486 | 3960 | 1.31 |
| CPMG | Full multistart, 9 starts | 434.55519 | 1285 | 3800 | 1.64 |
| CPMG | Profiled 9 points + TRF | 434.55519 | 2349 | 3960 | 1.80 |
| CPMG | Branch-aware profiled 9 points + TRF | **429.67003** | 4607 | 7470 | 3.14 |
| CEST | Normal TRF | 1829.63856 | 2337 | — | 30.09 |
| CEST | Selected DE + TRF | 1829.63856 | 2259 | — | 30.17 |
| CEST | Full multistart, 9 starts | **1460.03905** | 25477 | — | 328.04 |
| CEST | Profiled 9 points + TRF | **1460.03905** | 20948 | — | 105.64 |
| CEST | Branch-aware profiled 9 points + TRF | **1460.03905** | 50793 | — | 288.28 |
| DCEST | Normal TRF | **1008.65947** | 201 | 609 | 6.13 |
| DCEST | Shipped selected DE + TRF | 1712.39539 | 3358 | 10080 | 104.31 |
| DCEST | Full multistart, 5 starts | **1008.65946** | 1283 | 3879 | 39.56 |
| DCEST | Profiled 5 points + TRF | **1008.65947** | 1876 | 5673 | 57.36 |
| DCEST | Branch-aware profiled 5 points + TRF | **1008.65947** | 8691 | 26157 | 266.97 |

Objective requests include native numerical-Jacobian requests; fresh validation
evaluations are not solver requests. Actual kernel counts instrument
`BoundEvaluator._calculate_unscaled`, after profile caching, and include fresh
validation calculations. The first CEST run preceded this instrumentation, so
its exact kernel counts are unavailable. Traces also retain a conservative
requests-times-profiles work proxy; that proxy is **not** an actual model count.
The CPMG second seed and all bounded DCEST runs have exact kernel counts.

These are single-host exploratory times, not stable performance ratios.
Environment: macOS arm64 27.0.1, 18 logical CPUs, Python 3.13.13, NumPy 2.5.1,
SciPy 1.18.0, all five native thread controls set to 1, and one sequential
solver worker. CEST and DCEST processes overlapped on this multicore host;
timing comparisons have that limitation. Times exclude loading, uncertainty,
plotting, and output, and include measured fitting/validation. Subsecond
differences are not a basis for selecting an algorithm.

An initial larger DCEST run with GRID's root-sized nuisance budget was stopped
after one held point spent several minutes refining. Its incomplete timing is
not in the table. The bounded comparison uses **1000 requests per nuisance
attempt**. One ordinary profiled point exhausts that budget and is rejected;
another sign start makes that point usable in the branch-aware run, but its
full endpoint is still poor. Other cases used the original root-sized GRID
budget (38000 CPMG, 30000 CEST); no final TRF policy was altered.

## What the smaller alternatives establish

| At supplied outer start only, then one full TRF | Final chi-square | Requests | Kernel calculations | Seconds |
| --- | ---: | ---: | ---: | ---: |
| CPMG, ordinary profiling | 434.55519 | 209 | 480 | 0.24 |
| CPMG, branch trials from converged nuisance endpoints | **429.67003** | 273 | 586 | 0.29 |
| CEST, ordinary profiling | **1460.03905** | 3280 | — | 18.85 |
| CEST, branch trials from converged nuisance endpoints | **1460.03905** | 3790 | 4256 | 24.69 |
| DCEST, ordinary profiling | **1008.65947** | 292 | 894 | 9.08 |
| DCEST, branch trials from converged nuisance endpoints | **1008.65948** | 997 | 3027 | 31.24 |

Thus one profiled point suffices on these inputs. CEST's ordinary singleton
profiling changes the nuisance start enough to reach the better full basin;
additional outer sampling is unnecessary here. CPMG's endpoint branch trials
are inexpensive, but two full starts (supplied and all-DW-mirrored) also find
the same lower basin. Signed CPMG full multistart at the same nine outer points
costs about 1.56–1.59 seconds and 1226–1262 requests, and finds 429.67003 in
both seeds. The all-negative first start is intentionally included; the later
Sobol sign patterns are not evidence of a random recovery probability.

### DW initialization and genuine branch behavior

[Spin settings](../src/chemex/parameters/spins.py) initialize independent DW
coordinates to **zero** with finite -100/+100 ppm safety bounds. The
[initial-value builder](../src/chemex/parameters/parameterization.py) fills
missing model-derived quantities; it does not infer a physically justified
nonzero DW magnitude from experimental shifts. Supplied examples override DW:
CPMG_15N_IP uses global 2 ppm, CEST L18CD1 uses 0.01 ppm, and DCEST K19N has
DW_AB=2.51 and DW_AC=-1.97 ppm. Safety bounds are not branch-start magnitudes.

The prototype tests two concrete policies: mirroring supplied nonzero starts,
and mirroring converged factor endpoints while reusing their other nuisance
values. If a supplied DW is zero, the first policy uses the successful normal
factor endpoint magnitude. If that is also zero, it tries **no fabricated
nonzero branch**. In CPMG, all-zero DW starts give full TRF 432.13517, while
factor-local refinement plus endpoint-derived magnitudes recovers 429.67003.
This establishes a usable source in that case, not a guarantee that a zero
start can always escape. Where it cannot, require explicit user ranges/starts
or another validated experimental estimator.

The CPMG sign alternatives are not exact numerical symmetries: at the same
converged positive-branch vector, flipping one DW changes chi-square from
434.55519 to 434.23620, 432.75144, 434.48848, 432.82122, or 434.23210.
Independent branch fits prefer negative DW at the supplied outer point.
These results retain the actual finite pulses/carriers and do not justify
assuming sign symmetry from an idealized CPMG formula.

CEST exhibits a much larger competing-sign difference. In the small scope,
the complete higher basin has L18CD1 DW≈-0.16526 and chi-square 1829.63856;
the lower basin has DW≈+0.21004 and chi-square 1460.03905. At the fixed supplied
coupling point, the L18CD1 nuisance minima are chi-square 1310.54816 and
1729.98990. Mirroring the original 0.01 start requires 2020 requests for the
alternate local fit; mirroring the fitted +0.20913 magnitude and fitted
relaxation/CS values needs only **80**, reaching the same higher local branch.
Negating the original start is therefore neither the only nor the most
efficient branch policy.

DCEST tries all four DW_AB/DW_AC sign combinations **inside its single local
factor**. Some reconverge to the same score; others reach worse minima. Equal
scores alone do not establish an exact sign symmetry or a new distinct basin.
No scientific symmetry deduplication is implemented. A production exemption
from branch trials would need an experiment/model-specific proof, including
finite pulses, carriers, detection, and constraints.

The [focused tests](../tests/test_basin_discovery_prototype.py) demonstrate that
20 proven independent single-DW factors need 20 normal plus 20 alternate fits,
not a 2^20 root product. Multi-DW combinations remain exponential only within
each connected factor; they still need an explicit work budget.

### Recovery, ranking, sample count, and diversity

- Both CPMG coupling-only multistart batches miss the known lower sign basin
  (0/18); signed batches recover it in 1/9 and 2/9 endpoints, respectively.
  The best branch-profiled candidate reaches it in both seeds. Nine normal
  profiled candidates never fix this missing branch information.
- CEST full multistart reaches 1460.03905 in 5/9 runs. The best profiled
  candidate reaches it with and without sign trials. Selected-coordinate DE
  reaches the higher basin. DCEST multistart reaches 1008.65947 in 3/5 runs,
  and both best profiled candidates recover it; DE reaches a worse basin.
- Best-one, best-two, best-four, and best-four with minimum Euclidean distance
  0.25 in normalized outer space produce the **same best final score** for
  every tested prefix (start+2, start+4, start+8 where available). No measured
  benefit justifies compulsory extra polishes or diversity management.
- Ranking is useful but not an endpoint oracle. CEST's second-best ordinary
  profiled point among the first five polishes to the worse basin, while a
  worse-ranked point polishes to the better basin. In CPMG, a distant point
  with a much higher profiled score can polish better than the second-ranked
  point. Scores rank conditional fits, not complete attraction basins.
- All best-one selections choose the original outer start. Consequently these
  traces show that **zero additional samples** suffice on these particular
  inputs; they cannot establish that eight samples generally suffice.
- DCEST's unresolved nuisance basins and failed conditional fit are a concrete
  warning: computed Phi(g) is only a best-found local profile, not a certified
  minimum over nuisance coordinates. Branch-aware profiling improves one
  failed point without improving the final selected fit.

Recovery counts use a descriptive 0.001 chi-square window for these traces,
far below the reported between-basin gaps. This is not a production
basin-equivalence tolerance, statistical confidence estimate, or proof of a
global minimum. Two CPMG seeds and one seed for each expensive case are too
little evidence for general recovery probabilities.

## Smallest production design worth considering

1. Keep normal native full/grouped TRF and its acceptance/commit policy intact.
2. Expose a small internal arbitrary-profiled-point helper from existing GRID,
   sharing exact factor discovery, affine-constraint fallback, nuisance TRF,
   reconstruction, and fresh objective validation. Return a non-authoritative
   search state, not a newly accepted final fit.
3. Allow explicit, factor-local DW branch starts; favor mirroring successful
   local endpoints with their other nuisance values. Validate magnitudes,
   bounds and constraints, and record unsuccessful trials. Restrict this
   initially to the measured constant-DW models; polynomial temperature shifts,
   arbitrary sign symmetries, and large connected multi-DW factors need separate
   scientific qualification.
4. Start with the ordinary coupling point and one complete final TRF. If
   representative future failures require outer discovery, add a bounded
   deterministic 4/8-point Sobol batch in explicitly declared finite ranges.
   Keep simple best-one selection initially; optional best-two is a reasonable
   experiment, not a demonstrated default improvement. Do not add clustering.
5. Commit only a successfully converged complete authoritative TRF endpoint.
   Keep search, nuisance, and final-refinement budgets distinct.

Estimated implementation: **200–350 production LOC** for the reusable point
primitive, validated nuisance starts, branch trials, and final-TRF composition,
plus **200–400 LOC** of focused numerical/integration tests. An optional Sobol
driver and simple top-two policy would add roughly **80–150 LOC**, with Method
integration/provenance documentation additional. These are planning estimates;
no generic optimizer/plugin/scheduling framework is warranted. The research
runner is larger because it includes baselines, diagnostics and trace capture.

### SEARCH.DE

Retain its current compatibility and explicit held-nuisance semantics for now.
It adds no demonstrated advantage here and is actively misleading on the CEST
and DCEST scopes: optimizing with fixed nuisance values chooses a start whose
complete refinement is inferior. This supports discouraging it for strongly
coupled nuisance problems, not silently changing its objective or deleting it
on three examples. It might still help when nuisance values are already good
or basin-defining coordinates are all explicitly searched. Before deprecation
or replacement, qualify that use case and compare a deliberately profiled DE
objective separately. No DE semantics were changed by this task.

Conventional full-dimensional DE was not run. The shipped examples give useful
finite coupling/DW ranges, but no useful search ranges for all fitted
relaxation/CS coordinates. Broad numerical safety bounds are not scientific
global-search priors. A 10–17-dimensional DE comparison using invented ranges
would be expensive and unfair. Its evaluation scaling is also described in
[the SciPy DE documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html).

### Online learning

Do not implement a surrogate yet. The measured two-dimensional cases need
only a supplied-point check, and the four-dimensional case's existing start
already wins. DCEST conditional costs and unresolved branches could justify
future **within-fit** adaptive sampling if many costly outer evaluations are
actually needed. A small GP/Bayesian experiment could consume validated traces,
but must handle failed points and discontinuities from changing local branches.
It should earn its complexity against a modest Sobol batch and ordinary
multistart. These 5/9-point traces do not support a reliable learned model or
training on historical ChemEx datasets. No ML or new dependency was added.

## Reproduction and validation

Set the thread/cache environment in [benchmarks/README.md](README.md), then:

```sh
uv sync --locked
uv run --no-sync python benchmarks/basin_discovery.py cpmg --seed 597 --nuisance-budget 38000 --output /tmp/cpmg-597.json
uv run --no-sync python benchmarks/basin_discovery.py cpmg --seed 598 --nuisance-budget 38000 --output /tmp/cpmg-598.json
uv run --no-sync python benchmarks/basin_discovery.py cest --seed 597 --nuisance-budget 30000 --output /tmp/cest.json
uv run --no-sync python benchmarks/basin_discovery.py dcest --seed 597 --exponent 2 --nuisance-budget 1000 --output /tmp/dcest.json
uv run --no-sync python benchmarks/basin_discovery.py cpmg --diagnostics --seed 597 --output /tmp/cpmg-diagnostics.json
uv run --no-sync python benchmarks/basin_discovery.py cest --diagnostics --output /tmp/cest-diagnostics.json
uv run --no-sync python benchmarks/basin_discovery.py dcest --diagnostics --output /tmp/dcest-diagnostics.json
```

[Retained traces](basin_discovery_results.json) include complete outer points,
factor-local request/terminal evidence, reconstructed/full endpoint vectors,
all tested candidate selections, two CPMG signed-start batches, start
sensitivity, and single-point endpoint-mirroring diagnostics. Timings should
be remeasured rather than treated as stable oracles. Exact kernel counting and
diagnostic packaging were added during the experiment; the final runner
reproduces the numerical policies and additionally records its hash/thread
environment. Original CEST kernel counts cannot be recovered retroactively.

Validation completed:

- `uv sync --locked`: successful; no dependency or lockfile change.
- `uv run pytest -q tests/test_native_trf_scaling.py -k cest`: **1 passed**,
  1 unrelated test deselected; current-main premise qualified.
- `uv run --no-sync pytest -q tests/test_basin_discovery_prototype.py tests/test_profiled_grid.py tests/test_native_de_direct_trf.py`:
  **30 passed**, including a real CPMG factor/reconstruction/final-TRF check
  with absolute chi-square tolerance 0.001.
- Focused Ruff lint and formatting, explicit `ty check` of the runner, and
  `git diff --check`: successful.
- `graphify update .`: completed (AST-only); generated graph outputs are
  ignored repository artifacts.

The complete pytest suite, both Python CI targets, complete scientific
acceptance layer, package build, and website build were **not run**, as this
task changes no production code or public workflow. There were no pre-existing
test failures in the selected checks. An initial prototype-only test assumed
the wrong SciPy exception type for a negative Sobol exponent; that unrelated
assertion was removed before the final passing run.

Files added: `benchmarks/basin_discovery.py`,
`benchmarks/basin_discovery_results.json`, `benchmarks/basin_discovery.md`, and
`tests/test_basin_discovery_prototype.py`. CLI, Method/TOML syntax, parameters,
output/provenance formats, dependencies, public Python APIs, TRF scaling and
scientific calculations are unchanged. No user documentation/examples were
changed because this is a research harness, not a supported fitting feature.
Issue #712 and GitHub were not modified; no PR was opened.

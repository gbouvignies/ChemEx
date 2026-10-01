# Production qualification: optional factor-local start preparation

Production follow-up: [#795](https://github.com/gbouvignies/ChemEx/issues/795).

This follows [the retained #712 research](basin_discovery.md), with original
[traces](basin_discovery_results.json) preserved. Fresh current-main
qualifications and all request-budget measurements are in
[profile_preparation_results.json](profile_preparation_results.json).
The runnable [budget comparison](profile_preparation.py) reuses the research's
real-data construction and measurement helpers; it performs no new outer search.

## Git provenance

The old issue-775 branch was `92af4237464f98d5943dc677b628592ff3096c9f`.
The actual checkout at task entry was already a clean research-only branch at
`1a001ed4cbc219e7413ec1e1e39862e7340d73b2`. Research was preserved on
`codex/issue-712-basin-research` at
`88d8a5bdcb0afef754ebef4f2afd431e867bd4fe`.
After fetching, `origin/main` was
`d0ba6e344f007f1fb6caaf0e3e29b88322f15455`.
`codex/issue-712-profile-preparation` was created directly from that main;
the sole cherry-pick was the four-file research commit, becoming
`0b183045f669474550f5a93632b6cf964dae4dda`.
No issue-775 branch history or unrelated source changes were brought across.

Before production edits, the shipped CEST scaling qualification passed
(1 passed, 1 deselected), confirming the historical 1093-versus-802 motivating
failure still does not reproduce. Endpoint diagnostics reproduced CPMG
429.6700333, CEST subset 1460.0390471, and DCEST 1008.6594758 after complete TRF.
These subset objectives are not the complete shipped CEST workflow's objective.

## Budget comparison

Each policy limits actual ChemEx objective requests per nuisance TRF using
`DirectTrfInvocation.objective_request_budget`; this includes finite-difference
requests. Final complete/grouped TRF retains its existing 2000 × (nvars + 1)
budget. No iteration counter or public tuning option was added.

All comparisons below request the same held coordinates and endpoint mirrors.
The environment is recorded in the JSON: macOS, repository uv Python/runtime,
NumPy/SciPy versions, and one native thread. Times are single observations;
small differences are not performance claims.

| Case | Per-attempt policy | Preparation requests | Final requests | Prep + final seconds | Final χ² | Normal-factor fallback |
|---|---|---:|---:|---:|---:|---|
| CPMG | 100 × (n+1) | 183 | 90 | 0.28 | 429.6700333 | none |
| CPMG | 200 × (n+1) | 183 | 90 | 0.29 | 429.6700333 | none |
| CPMG | 500 × (n+1) | 183 | 90 | 0.27 | 429.6700333 | none |
| CPMG | min(500 × (n+1), 1000) | 183 | 90 | 0.29 | 429.6700333 | none |
| CPMG | min(600 × (n+1), 3000) | 183 | 90 | 0.26 | 429.6700333 | none |
| CEST | 100 × (n+1) | 1012 | 8504 | 118.0 | 1460.0390475 | CD1 |
| CEST | 200 × (n+1) | 1512 | 8504 | 120.5 | 1460.0390475 | CD1 |
| CEST | 500 × (n+1) | 3012 | 8504 | 127.6 | 1460.0390475 | CD1 |
| CEST | min(500 × (n+1), 1000) | 1512 | 8504 | 119.0 | 1460.0390475 | CD1 |
| CEST | min(600 × (n+1), 3000) | 3394 | 396 | 25.0 | 1460.0390471 | none |
| DCEST | all five policies | 799 | 198 | ~31.4 | 1008.6594758 | none |

**Selected policy: min(600 × (nvars + 1), 3000).** The useful CEST normal
nuisance attempt takes 2802 requests. The 100/200/500 policies fail that attempt;
fallback remains scientifically usable and reaches the same final basin, but
requires 8504 final requests. The 600 policy is the smallest tested simple
multiplier admitting that endpoint and saves substantial total work. Mirrors
add no value in this CEST/DCEST qualification; HOLD-only CEST is the appropriate
use. Its three successful nuisance fits total 2884 requests, then ~396 final
requests (the preimplementation qualification measured 1516.1486504 before
1460.0390471 after complete TRF).

The difficult *previously sampled* DCEST state was evaluated solely to stress
the work bound, not to introduce outer sampling into production:

| Request cap | Actual requests | Seconds | Result |
|---|---:|---:|---|
| 1000 | 1000 | 29.9 | unsuccessful normal nuisance attempt; original factor/root start retained |
| 3000 | 3000 | 90.3 | same safe fallback; no alternate trials |

The 3000 cap replaces the former large/root-sized nuisance allowance that could
spend many minutes. It is not a wall-clock timeout: this expensive DCEST attempt
still costs roughly 90 seconds on this machine. A universal 1000 cap is faster
for that attempt but markedly slower overall for the useful CEST workflow.
This is an explicit tradeoff, not a claim that all preparation is cheap.

## Follow-up: bounded MIRROR cardinality and alternate budgets

MIRROR now permits at most **two selected coordinates per connected factor**,
including selected coordinates whose fitted endpoint is zero. The guard runs
before sign enumeration; thus a qualified factor has at most three alternate
trials. Oversized factors retain their successful normal endpoint, record a
failed preparation branch request with the explicit diagnostic, skip all sign
trials, and continue to fresh root validation and the authoritative final TRF.
The runtime reports failure reasons. There is no silent truncation or root sign
product. Twenty independent single-DW factors remain forty local fits.

Normal nuisance preparation retains `min(600*(nvars+1),3000)`. Alternate-only
caps of 100, 200, 500, 1000 and 3000 were compared on the same real-data starts.
The follow-up rows are retained under `alternate_cap_comparison` in the JSON;
normal budgets, local/end-state values, actual objective requests, model counts
and times are retained. Final TRF is run again only when the reconstructed
vector differs; reused final measurements are explicitly marked.

| Case | Alternate cap | Preparation requests | Failed alternate fits | Prep seconds | Final χ² |
|---|---:|---:|---:|---:|---:|
| CPMG | 100 / 200 / 500 / 1000 / 3000 | 183 | 0 | ~0.14 | 429.6700333 |
| CEST | 100 | 3109 | 1 | 15.94 | 1460.0390471 |
| CEST | 200 | 3209 | 1 | 17.17 | 1460.0390471 |
| CEST | 500 | 3394 | 0 | 19.91 | 1460.0390471 |
| CEST | 1000 / 3000 | 3394 | 0 | ~20.0 | 1460.0390471 |
| DCEST | 100 | 389 | 2 | 12.33 | 1008.6594658 |
| DCEST | 200 | 589 | 2 | 18.26 | 1008.6594658 |
| DCEST | 500 | 799 | 0 | 24.87 | 1008.6594758 |
| DCEST | 1000 / 3000 | 799 | 0 | ~25.0 | 1008.6594758 |

**Selected alternate cap: 500 requests.** This is the smallest tested cap
preserving every qualified successful local branch endpoint. CPMG's useful
alternates finish within 100; the longest successful CEST alternate takes 385
requests, and the longest DCEST alternate takes 323. The 100/200 caps save work
in the cases where branches add no measured benefit, and successful-normal
fallback preserves their useful complete basin. The 500 cap conservatively
preserves the qualified branch evidence while reducing the maximum work of
an alternate sixfold compared with the previous 3000 ceiling. It does not save
work on these already-converged branches; this is a worst-case bound reduction.
No small timing difference is treated as meaningful, and the DCEST final
~1e-5 variation is not a distinct basin.

With at most two selected mirrors, each factor costs at most one normal attempt
(up to 3000 requests) plus three alternates (up to 1500 requests together),
excluding fresh materialization evaluations. Final TRF has its existing budget.
The difficult DCEST normal attempt can still take ~90 seconds at 3000 requests;
the smaller alternate cap does not change that normal-fit limitation.

Follow-up validation: **99 tests passed** across preparation (including real
CPMG/CEST/DCEST), GRID, Method/compiler, zero-valued oversized selections,
additive independent factors, rejected connected factors, reporting and branch
budget exhaustion. Ruff/formatting, ty and diff checks pass. The complete suite
was not repeated for this focused follow-up; the initial complete-suite results
and confirmed main-branch failures remain documented below.

Current total production diff from main: **473 added, 1693 removed; net −1220**.
The follow-up itself changes 38 source lines net and updates the existing
benchmark runner/traces, preparation tests, Method docs and changelog.

## Scientific and acceptance semantics

`SEARCH.PROFILE` uses `HOLD` and `MIRROR` selectors, with no ranges or seed.
HOLD coordinates stay at the current committed/root start only during nuisance
preparation. Exact GRID dependency discovery proves factor independence;
affine-coupled problems retain GRID's conservative single-factor treatment.
Explicit child nuisance starts pass through root-owned lineage, bounds, and
feasibility validation, including starts after `restart_from()`.

Each normal successful factor endpoint supplies DW magnitudes and all other
nuisance endpoint values for sign trials. Zero endpoints produce no invented
branch. Combinations occur only within a connected factor. The initial MIRROR
surface supports constant DW_AB/DW_AC; temperature-polynomial coefficients and
other DW parameterizations are rejected. No exact sign symmetry is inferred.

Failed alternate attempts keep successful candidates. Failed normal factors
keep original valid values and skip their branches. Complete reconstruction is
freshly validated; failed reconstruction restores the original root start.
Incomplete preparation is reported. Cancellation/interruption still stops the
user-requested operation. Preparation creates no acceptance or commit authority;
one normal complete/grouped TRF owns final acceptance, commit, covariance and
statistics. No TRF scaling, tolerances, residuals or bounds were changed.

Selected-coordinate DE was removed with its grammar types, compiler/dispatch,
lineage, backend, progress and tests. Legacy input gets a targeted retirement
diagnostic. Historical DE measurements remain in the research traces; rerunning
the old DE comparison requires the dedicated research revision. No replacement
global optimizer or online learning is justified by these cases.

## Validation and change size

Initial production source diff: **435 lines added, 1693 removed; net −1258** (including
blank lines and docstrings, excluding tests, research artifacts and docs).
The initial preparation module was 201 lines. No public solver options or generic
restart/candidate framework were introduced.

Focused checks covered Method/compiler, GRID, native direct/grouped TRF,
production fitting, and the new preparation cases. One larger focused run
reported 349 passes and one new invariant-test failure: the production check
correctly rejected a forged child during construction, earlier than the test
expected. The test was corrected and its targeted recheck passed. An earlier
new integration test used an evaluation-only preparation scope and the wrong
single-step output path; its corrected nuisance-failure scope now passes.
The final fallback-policy recheck passed four tests.

One complete suite was proportionate for removing a product backend and grammar:
`uv run --no-sync pytest -q -n 2` reported **2301 passed, 5 failed** in 196.79s.
Three failures were multiprocessing-manager socket restrictions in the sandbox;
all three passed outside the sandbox (3.94s). The other two reproduce with
identical errors using untouched archived `origin/main` source:

- `tests/test_native_evaluation.py::test_shipped_two_profile_dcest_plan_matches_direct_profiles_completely`:
  normalized calculation differs by 3.35276127e-8 at one of 84 elements.
- `tests/test_native_uncertainty.py::test_multivariate_scaled_svd_correlation_and_joint_propagation_reference`:
  Jacobian differs by 6.87564239e-8 at one of 84 elements (~7.74e-6 relative).

These pre-existing host-envelope failures were not hidden or repaired by changing
unrelated scientific tolerances. The complete suite was not repeated. Ruff,
formatting, `ty check`, and `git diff --check` pass. Graphify AST update was run;
generated maps remain outside the committed changes. No website build or package
build was run because dependencies, packaging and executable website assets
were unchanged.

Exact production-follow-up files (relative to the repository root; `D` deleted):

```text
M	CHANGELOG.md
M	benchmarks/basin_discovery.md
M	benchmarks/basin_discovery.py
A	benchmarks/profile_preparation.md
A	benchmarks/profile_preparation.py
A	benchmarks/profile_preparation_results.json
D	examples/Experiments/DCEST_15N_3States/Methods/method_de.toml
A	examples/Experiments/DCEST_15N_3States/Methods/method_profile.toml
M	src/chemex/configuration/method_expressions.py
M	src/chemex/configuration/method_plan.py
M	src/chemex/configuration/method_v2.py
M	src/chemex/configuration/method_validation.py
M	src/chemex/messages.py
D	src/chemex/optimize/de_direct_trf.py
M	src/chemex/optimize/direct_trf.py
M	src/chemex/optimize/method_compiler.py
M	src/chemex/optimize/native_deterministic.py
A	src/chemex/optimize/profile_preparation.py
M	src/chemex/optimize/profiled_grid.py
M	tests/configuration/test_method_plan.py
M	tests/models/kinetic/test_category_c_binding_models.py
M	tests/models/kinetic/test_nstate_workflows.py
M	tests/test_method_compiler_architecture.py
D	tests/test_native_de_direct_trf.py
M	tests/test_native_production_fitting.py
A	tests/test_profile_preparation.py
M	website/docs/user_guide/fitting/method_files.md
M	website/docs/user_guide/fitting/temperature_dependent_shifts.md
```

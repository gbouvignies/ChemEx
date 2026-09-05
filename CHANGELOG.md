# Changelog

All notable changes to ChemEx will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project uses [Calendar Versioning](https://calver.org/) (YYYY.MM.MICRO).

## [Unreleased]

### Added
- Added explicit `3st_eyring_linear` (A ↔ B ↔ C) and `3st_eyring_fork`
  (B ↔ A ↔ C) temperature-dependent kinetic models. The existing
  `3st_eyring` name remains a compatibility name for the linear topology.

### Changed
- ChemEx now presents known failures and user interruptions once at the CLI
  boundary, reports verified diagnostic paths when available, and exits with
  status 130 for Ctrl-C. Unexpected internal exceptions remain concise and
  private by default; the new global `--debug` option restores their natural
  Python tracebacks and chained causes for diagnosis.
- Standard independent spin, relaxation, chemical-shift, coupling, kinetic,
  binding, oligomerization, and equilibrium parameters now have finite broad
  default safety domains. Scalar parameter-file values inherit these bounds,
  while explicit bounds continue to override them. Current MCMC uses the same
  domains as uniform support; `S2` is now restricted to its physical `[0, 1]`
  domain.
- Three-state Eyring models now represent absent pathways structurally instead
  of using an extreme activation-enthalpy sentinel. All real Eyring `DH_*` and
  `DS_*` state and transition coordinates now have finite broad default bounds;
  explicit parameter-file bounds continue to override them.
- Direct Python callers of `calculate_kij_3st_eyring` must use its linear-model
  signature; the obsolete `DH_AC` and `DS_AC` arguments have been removed.

### Fixed
- Native MCMC automatic burn-in once again uses a valid finite positive
  tentative autocorrelation estimate from a short chain. Such runs publish
  posterior products with explicit tentative-estimate and tentative-burn
  warnings while withholding autocorrelation-derived ESS and MCSE; unavailable,
  invalid, or chain-consuming automatic burn estimates remain fail-closed.
- Native MCMC now uses isolated process workers for genuine CPU-parallel
  likelihood evaluation, while preserving seeded serial/parallel chain
  equivalence and native-library thread coordination.
- Native MCMC capture qualification no longer serially re-evaluates every stored
  walker state through the same authoritative evaluator. Structural, content,
  bounds, transition-mask, and evidence-lineage checks remain fail-closed.
- Native MCMC sampling now reports transition progress interactively and concise
  start/completion status in non-interactive output.

## [2026.09.0] - 2026-09-02

### Added
- Added lightweight schema-2 `run_info`: `parameters_used.toml` now records the
  immutable original invocation start, while one atomically replaced
  `restart.toml` records the latest committed continuation state and remains a
  normal `-p/--parameters` input. `run.toml` records resolved execution/software
  facts and persisted unsigned root seeds before stochastic work begins;
  `outcome.toml` reports both committed and successfully restarted revisions.
- Added canonical version 2 `[STEP.SEARCH.DE]` product execution. Explicitly
  selected linear or logarithmic coordinates are searched with a required seed,
  then the best eligible candidate initializes one normal full-coordinate TRF
  fit. Only that final aggregate TRF result can commit or supply uncertainty and
  statistics; DE failures never fall back to the old committed start.

### Changed
- Relaxation and cross-correlation rates used together in an NMR basis now obey
  the positive-semidefinite domain of their represented relaxation block before
  scientific evaluation. Native fitting uses private feasible coordinates while
  preserving public `ETAXY`, `ETAZ`, `SIGMA`, `MU`, `R1`, and `R2` parameter
  names, restart files, and output. The proton-exchange contribution `KHH` and
  diffusion coefficient `D` are non-negative. Invalid fixed starting states fail
  before the kernel; adaptive Jacobian scaling and fail-closed kernel exceptions
  are unchanged. Exact affine constraints on relaxation diagonals may use finite
  literal offsets, unary signs, constant multiplication, and division by a
  finite non-zero constant; nonlinear controller expressions remain rejected.
- Changed native local TRF refinement to one versioned adaptive inverse-Jacobian-
  column-norm scaling policy across Direct, grouped Direct, GRID nuisance fits,
  DE polishing, and resampling refits. Sensitivity-based trust-region geometry
  removes pathological, still-improving trajectories in shipped mixed-sensitivity
  methyl CPMG fits and is less dependent on parameter magnitudes and physical
  units. Local scaling is not used as an implicit global-search policy: a shipped
  CEST start reaches a different qualified attraction basin. Fixed objective-
  request budgets, convergence tolerances, bounds, residual semantics, and fail-
  closed acceptance remain unchanged. The invocation now binds the qualified
  SciPy 1.18.x dense-TRF compatibility class, and execution evidence fingerprints
  every solver request in received order, including refused and interrupted
  requests.
- Converted all shipped method examples to canonical method format v2 and made
  v2 the primary authoring format in the fitting guide. The guide now documents
  ordered complete-role overrides, explicit `ROLES_FROM`, automatic committed
  value continuity, GRID/selected-coordinate DE/implicit TRF, and the frozen v1
  deprecation window. No fitting or numerical behavior changed.
- Removed the disconnected native fit-manifest, artifact-catalogue, Components,
  PartialEvidence, historical run-reader, canonical-method-envelope, and broad
  workflow-provenance architecture. Current covariance and constrained evidence
  retain their detailed scientific records and a small numerical-runtime
  environment identity.
- Removed the migration-only canonical numerical lane, pinned environment and
  attestation machinery, calibration probes, and historical baseline fixtures.
- Covariance-derived local standard errors are now reported when an otherwise
  valid covariance lies within the existing three-sigma boundary threshold,
  including at a bound. Boundary evidence and the threshold remain unchanged,
  while the terminal retains an aggregate warning. Inline warnings are limited
  to fitted coordinates near their own simple bounds and constrained outputs
  structurally depending on them. Genuinely invalid covariance remains unavailable.
- Unified deterministic fit presentation and output around the authoritative
  aggregate method-step result. Interactive Rich output now shows one transient
  fit-component progress table, followed by one aggregate minimization summary
  and one concise uncertainty status. Multi-component Direct and GRID fits now
  write `Parameters/`, `Data/`, `Plots/`, `Statistics/`, `statistics.toml`, and
  aggregate `Grid/` artifacts directly at the method-step root; legacy `All/`,
  `Groups/`, and `Grid/Groups/` output is no longer generated. Reruns also
  remove stale `All/`, `Groups/`, and `Components/` trees. `run_info/` and
  ordinary scientific fitting and covariance behavior are unchanged.
- Ordinary successful native deterministic fits now derive accepted-point
  covariance evidence from exact retained optimizer Jacobians or an independent
  accepted-point fallback and report covariance-derived standard errors in
  `Parameters/fitted.toml`. Observation uncertainties are treated as absolute
  experimental standard deviations, so covariance is not scaled by reduced
  chi-square. Supported arithmetic constraints receive propagated errors;
  rank, normalization, derivative, or covariance-arithmetic failures withhold
  errors with a concise reason while preserving fitted central values.
  Condition and weak-mode diagnostics do not impose an arbitrary hard gate. Full covariance,
  correlation, constrained, and independent-block diagnostics are written
  under `Statistics/` without populating the legacy generic `stderr` field or
  mixing deterministic errors with MC, BS, BSN, or MCMC summaries.
- Native deterministic minimization once again reports fit progress. The Rich
  display shows objective evaluations, best and reduced chi-square, elapsed
  time, and applicable fit-component or aggregate GRID context; objective
  evaluations are no longer mislabeled as solver iterations.
  Internal MC, BS, and BSN refits remain under the existing statistics-level
  reporting rather than opening a progress stream for every replicate.
- Removed the unreachable legacy optimizer, MCMC, resampling, parameter-container,
  and live legacy-observation compatibility paths. ChemEx deterministic fitting
  now uses bounded SciPy TRF exclusively through the native ChemEx parameter and
  evaluation stack; native MC, BS, BSN, MCMC, and statistics behavior is
  unchanged. Canonical v2 omits `FITMETHOD` because TRF is implicit. Deprecated
  v1 continues to accept `FITMETHOD = "trf"` and its `"least_squares"` alias
  only for the frozen v1 compatibility window; arbitrary optimizer forwarding
  is unavailable. `lmfit`, `asteval`, and `numdifftools` are no longer runtime
  dependencies.
- Renamed misleading statistical summary fields without changing sampling or
  fitted values. MC, BS, and BSN now use `percentile_68_lower`,
  `percentile_68_upper`, and `half_percentile_68_width` instead of
  `lower_1sigma`, `upper_1sigma`, and `stderr`. MCMC now uses
  `credible_interval_68_lower`, `credible_interval_68_upper`, and
  `half_credible_interval_68_width`; genuine `mcse_mean` remains the Monte Carlo
  standard error of the posterior mean. Clarified that
  `run_info/parameters_used.toml` is the original invocation start and
  `run_info/restart.toml` is the latest committed continuation state, while the
  fitted/fixed/constrained parameter files are result reports.
- Fit output now has a success-last `run_info/outcome.toml` lifecycle marker.
  Reruns eagerly invalidate only known ChemEx result locations for every planned
  method step, preventing earlier deterministic or later statistics-family
  results from surviving as current after a failure. MC, BS, and BSN publication
  failures now retain truthful completed samples while removing complete-only
  summaries, correlations, and plots and publishing incomplete diagnostics.
- Routed Direct TRF, grouped Direct, GRID, and grouped GRID method steps through
  the native parameterization, evaluation, aggregate acceptance, and atomic
  commit stack in the production `chemex fit` path. Existing v1 method TOML,
  ordered role inheritance, profile selection, and calculated-data files are
  preserved under the aggregate method-step root. Native deterministic fits do
  not populate the legacy generic `stderr` field. MC, BS, and BSN statistics
  now start from that committed native fit and run wholly native Direct TRF
  refits with stable ordered seeds across serial and worker execution. Existing
  STATISTICS syntax and output directories are preserved; incomplete analyses
  retain explicit samples/failures diagnostics without publishing complete
  summaries or plots. MCMC now also starts from the single committed native
  deterministic fit and samples ChemEx weighted residuals through the native
  evaluation engine, with seeded serial/worker replay and the existing compact
  and nested method forms. Native MCMC requires finite, strictly ordered bounds
  and rejects `UPDATE_PARAMETERS = true`; failed or interrupted chains publish
  incomplete diagnostics without summary, sample, correlation, or plot files.
  The deprecated v1 `FITMETHOD` surface is closed to `trf`, with
  `least_squares` retained as a canonicalized alias for the same frozen
  compatibility window; v2 omits `FITMETHOD`. Other legacy and arbitrary
  optimizer spellings fail during method validation. Ordinary
  filtering, revision-zero model/default resolution, simulation
  back-calculation, and native fit publication no longer construct lmfit
  parameter containers. After successful GRID aggregate acceptance and atomic
  commit, requested statistics analyze the committed deterministic result.

### Fixed
- Canonicalized model-free relaxation derivations by physical N-H and C-H
  orientation, so experiments that reverse their local `i`/`s` basis ordering
  now contribute one unambiguous model-owned expression. Single-spin
  components follow the physical nucleus, symmetric components remain shared,
  and deuteration and proton-exchange contributions retain their established
  semantics.
- Restored explicit method `FIT`/`FIX` overrides for user-overridable spin
  parameters whose ordinary baseline is a state-A equality, such as `J_B` in
  the shipped N15 NH RDC analysis. Scientific estimation capability is now
  declared separately from whether an experiment fits the parameter by
  default; model-owned derivations and unsupported parameters remain protected.
- Restored target-relative constraint resolution across companion spins at the
  same molecular group and atom site. Exact spin context remains preferred,
  companion context outranks global fallback, and unrelated or genuinely
  ambiguous candidates still fail closed without encoded-ID similarity.
- Restored GRID as exact profiled chi-square surface analysis. Declared grid
  coordinates remain held at every point while the remaining independent
  fitted coordinates are optimized with native TRF; exact objective
  factorization avoids Cartesian products of unrelated local axes. ChemEx now
  reconstructs one coherent joint grid and nuisance-coordinate solution,
  freshly validates the complete root state, and commits it once. GRID steps
  withhold deterministic covariance for their discrete selected coordinates,
  while later Direct steps and explicitly requested sampling statistics retain
  their normal behavior. `Grid/` now exposes status-bearing raw factor TSVs,
  reusable exact 1D profiles and factor-supported 2D profile TSVs, matching
  PDFs, and `summary.toml` instead of the native multistart `grid.out` layout.
- Fixed native fit commits rejecting finite expression-derived parameter values
  solely because they fall outside configured numerical-coordinate bounds.
  Independent TRF, GRID, DE, and MCMC coordinates remain bounded at their
  numerical-problem boundaries; derived values continue to follow constraint
  evaluator domains and scientific postconditions without clipping.
- Automatic MCMC burn selection now fails closed when autocorrelation evidence
  cannot establish a defensible retained window. ChemEx preserves the completed
  raw chain and incomplete diagnostics but withholds posterior summaries,
  retained samples, correlations, and plots instead of silently using zero burn.

### Infrastructure
- Versioned the Docusaurus documentation by ChemEx release and changed public
  website deployment to build from the exact latest published stable release
  tag, while keeping prior stable documentation available.
- Removed the unused `B1FieldConfig` model and its `default_distribution_config`
  helper from `b1_config.py`; the flat-table `b1_frq` schema it implemented was
  never wired to any experiment. `B1InhomogeneityMixin` remains the single B1
  configuration path, and the distribution-union parsing it depends on is
  unchanged.
- Removed the obsolete `ParameterStore` compatibility aliases
  (`add_parameters`, `add_parameters_mf`, `sort_parameters`,
  `set_param_values`, `set_param_defaults`, `fix_all_parameters`,
  `reset_parameters`) left over from the `AnalysisSession` migration. Internal
  callers now use the canonical short method names
  (`add_multiple`/`add_multiple_mf`, `sort`, `set_defaults`, `fix_all`,
  `reset`); parameter semantics and fit/simulate behavior are unchanged.
- Kept `AnalysisSession` at the `run_methods` orchestration boundary so native
  deterministic steps can use revisioned analysis values and stale-safe commit;
  experiments now own profiles and their direct required IDs, while runtime and
  output consumers receive sealed metadata plus committed or resolved values
  explicitly. Statistics execution still receives the session's explicit
  `ExecutionSettings`.

## [2026.06.1] - 2026-06-26

### Added
- Added lightweight `run_info/` output for fits, including runtime metadata,
  command and input provenance, optional Git metadata, starting independent
  parameters, and copies of experiment, parameter, and method TOML inputs.
- Added multi-state `observed_state` support for final-magnetization
  experiments. Passing a state list now detects the unweighted sum of the
  selected final magnetization components, while the existing string form keeps
  single-state detection.
- Added `start_state` to control non-equilibrium preparation independently from
  detection. It accepts one state, multiple states, or an empty list to request
  equilibrium preparation when overriding a state-specific default.

### Changed
- Updated CEST, D-CEST, CPMG, and relaxation experiments to derive detection
  expressions and starting terms from normalized state selections.
- Validation now checks `observed_state` and `start_state` against the active
  model states, rejects duplicate state selections, and keeps the first observed
  state as the reference for state-specific defaults and offset filtering.
- Preserved the previous observed-state preparation defaults for `cest_1hn_ap`,
  `cpmg_1hn_ap`, `cpmg_1hn_ap_0013`, and `cpmg_ch3_1h_sq`; use
  `start_state = []` to request equilibrium preparation for these experiments.
- Preserved D-CEST HD model starting terms by default, while allowing
  `start_state` to override them explicitly.

### Removed
- Removed `detect_all_states`; list every active model state in
  `observed_state` to detect all states.
- Removed `cs_evolution_prior`; use `start_state` for non-equilibrium
  preparation or `start_state = []` for equilibrium preparation when overriding
  a state-specific default.

### Fixed
- Fixed `run_info` path resolution so literal tilde-prefixed input paths are not
  expanded before being resolved against the working directory.
- Increased CEST plot output precision so closely spaced PPM offsets stay
  distinct in `.fit` and `.exp` files.

### Documentation
- Updated the fitting guide and experiment reference pages for `observed_state`
  list syntax, `start_state`, replacement syntax for removed options, and the
  limits of multi-state component selection.
- Documented the new `run_info/` output directory.
- Removed duplicate MCMC implementation-plan documentation.

## [2026.06.0] - 2026-06-03

### Added
- Added MCMC posterior sampling as a `STATISTICS` method for fitted parameter
  uncertainty, with text outputs for posterior summaries, samples,
  correlations, diagnostics, and a PDF visual report, including automatic
  autocorrelation-based burn-in reporting, tentative short-chain burn-in
  handling, and 95% equal-tailed credible intervals.
- Added canonical `Statistics/<method>/` output directories for sampling-based
  uncertainty analyses, with TSV sample/correlation tables, summaries, and
  diagnostics and PDF visual reports for Monte Carlo and bootstrap runs.

## [2026.05.1] - 2026-05-13

### Infrastructure
- Resolved all ty 0.0.35 type-checking errors: assert-based narrowing in
  `b1_config.py`, `cast`+`getattr` for dynamic `MinimizerResult.params` in
  `minimizer.py`, and type-correct `default_factory` placeholders in `data.py`.
- Refactored `Data` to use `model_post_init` (idiomatic Pydantic v2) instead of
  a custom `__init__`.
- Bumped ty dev dependency to `>=0.0.35`.

## [2026.05.0] - 2026-05-13

### Fixed
- Fixed Monte Carlo and bootstrap simulations crashing on Python 3.14 with
  `TypeError: cannot pickle 'mappingproxy' object`. `Basis.__deepcopy__` now
  returns `self` (correct for an immutable frozen dataclass), and
  `LiouvillianReadout` stores a plain `dict` copy of its vectors.

### Infrastructure
- Made release uploads idempotent.

## [2026.04.0] - 2026-04-07

### Added
- Added composable kinetic-model suffixes so residue-specific fits can use `.rs` and combine it with `.mf` and `.tc` (for example, `2st.rs.mf`), while keeping `2st_rs` as a legacy alias.

### Changed
- Refactored the NMR Liouvillian internals into reusable `chemex.nmr` engine, pulse, B1, detection, effective-field, magnetization, and tensor components.
- Generalized residue-specific parameter naming across kinetic models and documented the `.rs` suffix in the fitting guide.
- Updated experiment output naming to derive unique stems from input paths so files from different datasets no longer collide when basenames match.

### Fixed
- Hardened constraint and grid parsing with clearer validation errors and correct precedence for more specific grid entries.
- Fixed empty-selection handling so experiment collections with no active profiles are treated as empty.
- Tightened numeric validation around B1 inhomogeneity and shift/R1rho eigenvalue calculations.
- Replaced detection `eval` usage with a parser and handled TOML file read/parse failures consistently.

### Infrastructure
- Expanded regression coverage for NMR engine behavior, output paths, TOML/config parsing, selection handling, and residue-specific models.

## [2026.03.0] - 2026-03-11

### Fixed
- Fixed `simulate` plotting for CEST and CPMG experiments in released packages.
- Fixed profile cache invalidation when data objects are mutated during analysis flows.

### Changed
- Made runtime session state explicit to reduce hidden global state during analysis runs.
- Refined computed-field and typing boundaries across the runtime layer.

### Infrastructure
- Updated assorted build, CI, and website dependencies.

## [2025.10.0] - 2025-10-27

### Added
- **New Experiments**: Added fitting module for 13CH3 13C CPMG [0013] experiment (#374)
- **New Features**:
  - Implemented 4-state Eyring model for chemical exchange kinetics (#384)
  - Added `detect_all_states` option for CEST, DCEST, and CPMG experiments (#385)
- **Documentation**: Added CLAUDE.md for project guidance and development instructions

### Changed
- **Build System**:
  - Migrated to `src/` layout for better package structure
  - Updated build system configuration to use `uv_build`
  - Switched to static versioning (removed hatch-vcs dependency)
- **Type System**: Major refactor of array typing across the codebase
  - Updated type annotations to use `Array` for intensities and magnetization
  - Improved type consistency and clarity across parameters
- **Configuration**: Refactored configuration merging to use `deep_update` for improved readability
- **Documentation**:
  - Updated installation instructions for uv in README and welcome guide
  - Updated Docusaurus config to enable new experimental features

### Fixed
- Fixed initial conditions for two-spin systems in ChemEx (#392)
- Updated FIX parameters in `method.toml` for CPMG_CH3_MQ to fix DW_AB configurations
- Updated suffix handling in CpmgCh313CH2c0013Settings for start_terms and detection properties

### Dependencies
- **Python**: Bumped minimum version to 3.13
- **Core Dependencies**:
  - matplotlib: 3.10.0 → 3.10.1
  - numpy: 2.2.3 → 2.2.5
  - pydantic: 2.10.6 → 2.11.4
  - rapidfuzz: 3.12.1 → 3.14.1
  - lmfit: 1.3.2 → 1.3.3
  - cachetools: 5.5.1 → 5.5.2
  - rich: 13.9.4 → 14.0.0

### Infrastructure
- **CI/CD**:
  - Refactored workflows and fixed dependabot group conflicts
  - Updated GitHub Actions: setup-python (5→6), setup-node (4→6), checkout (4→5)
  - Updated sigstore action (3.0.0→3.1.0)
  - Added Claude Code GitHub Workflow (#383)
  - Updated Node.js version to 20 in deployment workflows
- **Dependabot**: Updated configuration for pip, uv, github-actions, and npm ecosystems

### Website
- Multiple dependency updates for Docusaurus and related packages
- Security updates for various npm packages

## [2025.4.0] - 2025-04-XX

### Added
- Support for 5-state and 6-state kinetic models
- Enhanced CEST plotter with updated line styles

---

For a complete list of changes, see the [GitHub Releases](https://github.com/gbouvignies/chemex/releases) page.

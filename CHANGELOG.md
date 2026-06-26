# Changelog

All notable changes to ChemEx will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project uses [Calendar Versioning](https://calver.org/) (YYYY.MM.MICRO).

## [Unreleased]

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

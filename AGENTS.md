# Working on ChemEx

## Purpose and authority

ChemEx is a production scientific Python package for analysing chemical exchange
detected by biomolecular NMR, including CPMG, CEST, relaxation, and shift
experiments. Treat changes to calculated intensities, residuals, fitted
parameters, uncertainties, and model selection as scientific changes, even when
the code edit looks mechanical.

Use this source-of-truth order:

1. Current source code, tests, configuration models, and shipped examples define
   implemented behaviour.
2. `pyproject.toml`, `uv.lock`, and `.github/workflows/` define supported
   environments and automated checks.
3. `website/docs/`, `README.md`, `CONTRIBUTING.md`, and `CHANGELOG.md` describe
   supported user and contributor workflows.
4. This file is the canonical operational guide for contributors and agents, but
   it is derived from the sources above. Update it when those sources change.

Do not treat generated maps, saved benchmark output, or an old output directory
as authoritative. If documentation and implementation disagree, establish the
current behaviour from code and tests, then correct the documentation or raise
the discrepancy.

## Orient by task

- CLI and top-level orchestration: `src/chemex/cli.py`,
  `src/chemex/chemex.py`, `src/chemex/__main__.py`
- Session state, plugin registration, and worker/native-thread settings:
  `src/chemex/runtime/`
- TOML parsing and Pydantic schemas: `src/chemex/toml.py`,
  `src/chemex/configuration/`
- Experiment construction and registration:
  `src/chemex/experiments/builder.py`,
  `src/chemex/experiments/factories.py`,
  `src/chemex/experiments/loader.py`,
  `src/chemex/experiments/catalog/`
- Experimental data, profiles, and experiment collections:
  `src/chemex/containers/`, `src/chemex/evaluation/`, `src/chemex/filterers.py`
- Kinetic/exchange model definitions and registration:
  `src/chemex/models/model.py`, `src/chemex/models/factory.py`,
  `src/chemex/models/kinetic/`
- Parameter naming, creation, constraints, and lmfit conversion:
  `src/chemex/parameters/`
- NMR basis, Liouvillian calculations, propagators, pulses, distributions, and
  readout: `src/chemex/nmr/`
- Fitting, grouping, grids, resampling, MCMC, and statistics:
  `src/chemex/optimize/`, `src/chemex/uncertainty.py`
- Output tables, fitted parameter files, and plots:
  `src/chemex/printers/`, `src/chemex/plotters/`
- Fit provenance and archived inputs: `src/chemex/run_info.py`
- User workflows and realistic TOML/data: `examples/`, `website/docs/`
- Verification: `tests/`; exploratory performance work: `benchmarks/`

Search for the experiment name, model name, parameter name, or TOML key before
editing. Inspect its schema, registration site, representative example, and
tests together.

## Execution and data flow

The installed `chemex` command calls `chemex.chemex:main`; `python -m chemex`
uses the same function. `fit` and `simulate` follow this path:

1. `cli.py` parses paths, model selection, profile selection, output options,
   and fit execution settings.
2. `AnalysisSession` registers experiment/model plugins and owns the active
   `ModelState`, `ParameterStore`, `ParameterFactory`, and execution settings.
3. `chemex.run()` sets the model, then `build_experiments()` reads every
   experiment TOML.
4. `experiments/builder.py` selects a registered `Creators` bundle, validates
   the experiment-specific configuration, loads data relative to the experiment
   TOML, applies profile selection, and creates a `Profile` for each spin system.
5. Parameter settings from the kinetic model and spin/basis are combined in
   `parameters/factory.py` and stored by stable parameter IDs.
6. A `Profile` updates its `Spectrometer`, calls the pulse sequence, applies the
   established scaling convention, and returns masked, uncertainty-weighted
   residuals.
7. Simulation writes calculated data directly. Fitting applies method sections,
   filtering, grouping, lmfit minimisation or grids, optional resampling/MCMC,
   then writes parameters, data, plots, statistics, and `run_info/`.

Keep CLI concerns in parsing/orchestration, schema concerns in
`configuration/`, experiment composition in catalog modules/builders,
scientific matrix operations in `nmr/`, parameter semantics in `parameters/`
and `models/`, and optimisation policy in `optimize/`. The underscore packages
`nmr/_engine/` and `nmr/_pulses/` are internal numerical machinery; prefer the
`Spectrometer` façade from experiment code and change internals only with
focused numerical tests. Preserve intentional compatibility exports such as
`chemex.nmr.spectrometer.calculate_propagators`.

## Scientific and numerical invariants

Before changing scientific code, identify the implemented equations or reference
behaviour and record representative outputs or a failing regression test.
Preserve unless the change explicitly requires otherwise:

- Units encoded by current schemas and calculations: delays/pulse widths in
  seconds, offsets/B1/couplings in Hz, chemical shifts in ppm, proton Larmor
  frequency in MHz, exchange and relaxation rates in s^-1, populations as
  fractions. `Conditions.temperature` is currently used as degrees Celsius by
  the Eyring models.
- Gyromagnetic-ratio, frequency-sign, coherence-order, carrier/offset, and phase
  conventions in `nmr/constants.py`, `nmr/basis.py`, and `nmr/_engine/`.
- Kinetic state names and ordering. `ModelSpec` derives lowercase states in
  `a`, `b`, `c`, ... order, and basis matrices/vectors and kinetic parameter
  names depend on that order.
- Population normalisation, forward/reverse rate semantics, parameter
  expressions, condition scoping, residue-specific scoping, bounds, and which
  parameters vary.
- Magnetisation component ordering, detection/start-state semantics,
  Liouvillian and propagator multiplication order, array shapes, broadcasting
  axes, real/complex dtype, immutable shared basis/distribution data, and cache
  invalidation.
- Reference-point handling, masks, error shapes, scaling, and residual sign.
  `Data.scale` is the established uncertainty-weighted scale and
  `ProfileEvaluator` uses `(calc - exp) / err` before masking.
- Output naming and formatting, parameter names/IDs, CLI defaults, TOML key
  aliases, multi-file merge behaviour, and paths resolved relative to the
  experiment file.

Use realistic values and `numpy.testing.assert_allclose` or `pytest.approx` with
a tolerance justified by numerical conditioning and the intended change.
Exercise limiting or boundary cases for new scientific behaviour. Do not replace
reference values merely because a rewritten algorithm produces different
numbers; explain and validate any intentional numerical change.

## Extending experiments and models

### Experiment or pulse sequence

Start from the closest stable module in `src/chemex/experiments/catalog/` and
its example files, not from `catalog/wip/`. An experiment module normally
defines:

- a literal experiment name and Pydantic settings/configuration;
- `to_be_fitted` rate/model-free requirements;
- a spectrometer builder that chooses the `Basis`, conditions, carriers, fields,
  distributions, and detection;
- a pulse-sequence object implementing `calculate(spectrometer, data)` and
  `is_reference(metadata)`;
- dataset, filterer, printer, and plotter choices;
- `register()`, which supplies all creators to `experiments.factories`.

Discovery is automatic through `experiments/loader.py`; do not add a second
registry. Add or update schema tests, calculation/regression tests, example TOML
and data when user-facing configuration changes, and website experiment
documentation when the experiment is supported. Validate reference planes,
state selection, pulse ordering, array shape, and at least one realistic
calculation.

### Kinetic or exchange model

Add or change `src/chemex/models/kinetic/settings_<name>.py`. Its `register()`
must register a settings maker in `model_factory`; expression helper functions
must also be registered in `user_function_registry`. Preserve the
`ParamLocalSetting` names, condition-selection tuples, bounds, defaults,
expressions, population constraints, and units because these determine stable
parameter IDs and TOML matching. `ModelSpec.from_name()` derives the state count
from the model name and supports only the extensions declared in
`models/model.py`.

Test settings and expressions directly, then add integration coverage across
relevant conditions. Include conservation/normalisation, non-negative rate,
zero-exchange, equilibrium, or temperature/concentration limits as applicable.

### Configuration, parameters, data, and output

Change shared schemas only when all affected experiment modules and shipped TOML
remain valid or a deliberate migration is documented. Put reusable constrained
types in `configuration/types.py`; keep file parsing in `toml.py` or
`containers/dataset.py`. Parameter naming/search belongs in `parameters/`, not
in pulse sequences. A new output format needs its printer/plotter change plus
path-collision and CLI/integration coverage where relevant.

## Compatibility and repository-owned artifacts

Existing experiment, parameter, and method TOML files are compatibility
fixtures even when not imported by pytest. Preserve documented CLI flags and
defaults, stable output paths, saved input usability, and the fit provenance
schema in `run_info.py`, or provide an explicit migration and changelog entry.
Treat imports exercised by tests, benchmarks, and scripts as quasi-public; do
not casually move or rename them.

- `src/chemex/experiments/catalog/wip/` is experimental but is still
  auto-discovered and registered. It is excluded from Ruff and `ty`; do not use
  that exclusion to weaken stable code, and do not promote or redesign WIP code
  without maintainer intent.
- `dist/`, `website/build/`, caches, example `Output*` directories, and ChemEx
  result directories are generated/local artifacts. Do not edit or commit them.
- `uv.lock` and `website/package-lock.json` are generated lock files; change
  them only with the corresponding dependency change.
- `benchmarks/` and saved benchmark result text are exploratory aids, not
  scientific golden results or CI gates. Inspect a benchmark before relying on
  it and record hardware/software context for performance comparisons.

Update `website/docs/`, examples, `README.md`, `CONTRIBUTING.md`, or
`CHANGELOG.md` only when their audience is affected. Reference authoritative
schemas and examples instead of copying large key lists into prose.

## Validation

Set up the repository with:

```sh
uv sync --all-extras --dev
```

The blocking Python CI target is 3.13; Python 3.14 is currently informational.
Run the narrowest relevant check during development, then the applicable
repository checks:

```sh
uv run pytest -q
uv run ty check
uv run ruff check .
uv run ruff format --check <changed-python-files>
git diff --check
```

Choose additional validation by change type:

- Parser/schema/TOML: focused `tests/configuration/`, `tests/test_toml_io.py`,
  representative shipped TOML, and invalid-input cases.
- Experiment/pulse sequence/NMR kernel: focused `tests/experiments/` and
  `tests/nmr/`, a realistic calculation, reference/boundary cases, and an
  example `simulate` or `fit` when practical.
- Kinetic model/parameters: focused `tests/models/`,
  `tests/test_constraint_expressions.py`, and relevant parameter tests,
  including physical limiting cases.
- CLI/output/provenance: `tests/test_cli.py`, `tests/test_output_paths.py`,
  `tests/test_run_info.py`, related printer/plotter tests, and a concrete CLI
  invocation using explicit TOML paths.
- Optimisation/statistics: focused grid, execution, MCMC, resampling, and
  uncertainty tests; control random seeds where reproducibility matters.
- Performance-sensitive numerical work: correctness tests first, then a
  targeted benchmark with before/after environment and timings.
- Packaging: `uv build`.
- Website: from `website/`, run `npm ci` and `npm run build`.

For bug fixes, first capture the failure in a regression test when practical.
For scientific changes, test both the intended result and important limiting
cases. If the full suite is not run, state exactly what was and was not run.
Never hide a pre-existing failure or attribute it to the current change.

## Scope and workflow

Keep changes focused. Do not perform repository-wide formatting, mass renames,
broad typing migrations, speculative cleanup, or unrelated refactoring. Add
precise types to new or changed interfaces when they clarify shapes or
contracts; do not weaken Ruff or `ty` rules to make a change pass.

Use the proportional workflow:

- Small maintenance: inspect, implement, run focused tests, then relevant
  repository checks.
- Bug fix: reproduce, add a regression test where practical, fix, validate.
- Numerical/scientific change: establish reference behaviour and tolerances
  before implementation, then compare outputs and explain differences.
- Cross-cutting feature: write a short explicit plan before changing multiple
  subsystems.
- Architectural change: document alternatives, dependency impact, compatibility
  and migration strategy before rewriting.

## Completion report

Every completed change should report:

- behaviour changed, including user-visible and scientific effects;
- compatibility implications for CLI, TOML, parameters, outputs, and Python API;
- tests and checks run with outcomes;
- numerical comparisons, reference cases, and tolerances when relevant;
- documentation/examples updated, or why none were needed;
- remaining risks, pre-existing failures, and narrowly justified follow-up.

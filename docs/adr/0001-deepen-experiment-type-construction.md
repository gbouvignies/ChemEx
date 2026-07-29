---
status: accepted
---

# Deepen Experiment Type construction

ChemEx will represent each named acquisition and calculation protocol as an
**Experiment Type** whose adapter supplies the validated configuration type and
the scientific construction of a spectrometer and pulse sequence. A single
deep `chemex.experiments.experiment_types` interface will identify an Experiment
Type from an input file and construct the complete configured **Experiment**,
including data selection, profiles, parameters, presentation objects, and noise
estimation. This concentrates generic assembly behind one seam while keeping
experiment-specific science explicit and strongly typed in catalog modules.

## Considered options

- The existing generic `Creators` bundle was rejected because it exposes every
  construction step and spreads orchestration knowledge across the catalog and
  builder.
- A state-indexed construction DSL was rejected because no current Experiment
  Type needs its added flexibility.
- Four semantic family adapters—CPMG, CEST, relaxation, and shift—were chosen
  because each corresponds to repeated loader, filterer, printer, and plotter
  policy. One-off EXSY support remains explicit composition.

## Consequences

Catalog modules retain their configuration classes, spectrometer builders, and
pulse-sequence classes. Registration remains explicit so automatic discovery
and registry reset behavior stay deterministic; duplicate registration is
idempotent only for the same adapter. Expected input failures and nonfatal
notices cross the interface as typed values, while presentation and process
exit remain outer orchestration concerns. All Experiment Types migrate
atomically, without a compatibility layer for the internal `Creators`,
`Factories`, or single-file `build_experiment` interfaces.

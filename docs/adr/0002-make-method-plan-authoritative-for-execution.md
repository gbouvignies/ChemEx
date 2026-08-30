---
status: accepted
---

# Make Method Plan authoritative for execution

ChemEx will use **Method Plan** as the sole runtime representation of fitting
execution semantics. Version 1 and version 2 TOML, together with the supported
Python `Methods` mapping, remain compatibility inputs. Each supported fitting
façade prepares those inputs through
`chemex.configuration.method_input.prepare_method_plan`, which normalizes and
validates them against the sealed parameter model before side effects occur.

`run_fit` remains the run-level façade and owns provenance, initial filtering,
planned-output invalidation, and outcome publication. The supported
`run_methods` Python façade owns only compatibility preparation. Both call the
single deep procedural seam
`chemex.optimize.method_plan_execution.execute_method_plan`, which returns
`None` or propagates an exception.

The executor owns Method Plan traversal, inherited parameter roles, step-local
profile selection and reset semantics, output paths, deterministic-before-
statistics ordering, and statistics ordering. It delegates numerical fitting,
resampling, and MCMC to the existing specialized modules using canonical
resolved facts: `ProfileSelection`, a compiled active parameterization,
`ResamplingRequest`, and `McmcRequest`. Those modules do not receive a Method
Plan, Method Step, legacy `Method`, or method-format origin.

## Considered options

- An executor object with capability records and protocols was rejected because
  ChemEx has one production implementation of each execution kind. The added
  composition interface would be speculative and shallower than the procedural
  seam.
- A staged migration retaining canonical-to-legacy operational reconstruction
  was rejected because it would preserve two sources of runtime authority and
  weaken the deletion test.
- Repository-wide removal of standalone `Method` and `Selection` interfaces was
  rejected as a separate architectural problem. Existing unrelated supported
  callers retain those interfaces.

## Consequences

Successful Method Steps commit scientific state incrementally. Execution is
fail-fast: a failed step prevents later steps from running while preserving
earlier successful commits. Neither the immutable Method Plan nor an execution
cursor is mutated. Reinvocation interprets the complete Method Plan afresh
against the currently committed scientific state.

The canonical execution path never reconstructs legacy operational method,
selection, statistics, or MCMC configuration objects. Tests establish
compatibility through selected profiles, committed parameter values, artifacts,
specialized-module invocations, statistics non-mutation, failure propagation,
partial commits, and representative native and command-line fits. Existing
numerical regression expectations remain authoritative.

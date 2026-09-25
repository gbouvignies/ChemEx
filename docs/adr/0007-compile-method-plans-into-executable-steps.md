---
status: proposed
---

# Compile Method Plans into executable steps before fitting

This is a proposed successor to ADR 0002. If accepted, it supersedes ADR 0002's
assignment of role inheritance, profile selection, and search resolution to the
executor. ADR 0002 remains the accepted description of current ownership until
this design is implemented. The reason for the change is deletion of the second
method interpreter: execution should consume resolved scientific instructions,
not read Method language again.

## End state

The input `MethodPlan` remains the immutable, renderable description produced by
v1 TOML, v2 TOML, or the supported Python `Methods` mapping. One compilation
interface accepts that plan, a sealed parameter model, and the loaded profile
population. It returns an immutable **Executable Method Plan** or source-located
diagnostics. The global resolver and active compiler form the sole fitting
pipeline that interprets Method selectors, `ROLES_FROM`, ordered role actions,
constraints, profile selection, or GRID/DE declarations. Input adapters only
preserve version-specific syntax and defaults.

An executable fit step owns its name and output identity; its concrete, ordered
profile bindings; effective parameter roles and a restricted constraint/value
program for the active dependency scope; concrete independent FIT IDs; a typed
Direct, GRID, or DE search instruction with resolved coordinate IDs, ordering,
ranges, and seed; and typed statistics requests. The compiler retains source
references while diagnosing invalid input; the archived Method input remains
the provenance record. Source references do not travel into numerical
instructions. The plan is immutable for one sealed model and one loaded profile
population. It does not capture parameter values, observation masks, an
optimizer problem, or a resume cursor. Runtime objects may refer to mutable
profiles, but their membership and order in the plan cannot change.

The executable plan carries no unused plan-level model or population identity
fields. Each compiled parameter program records the sealed model identity and
checks it against the step-time values snapshot. Concrete profile bindings
capture the loaded membership; a persisted population fingerprint remains a
separate, versioned provenance decision if restart verification needs one.

An explicit skipped step owns its name, ordinal, output identity, and
`no_profiles` reason. It has no active parameterization or executable search or
statistics. It still contributes its *global* effective roles and constraints
to compilation of later `ROLES_FROM` steps. The executor reports the step and
its no-data outcome, then continues without a fit, sampling, statistics, seed
record, parameter commit, or step output.

At each non-skipped step, execution binds the already compiled program to the
current committed `AnalysisValues` snapshot, checks model/snapshot identity,
and builds numerical workspaces. It performs fitting, optional statistics,
publication, and incremental commits in existing order. Value-dependent
feasibility and numerical failures remain runtime concerns. Binding cannot
match selectors, choose roles, expand Method constraints, select profiles, or
resolve GRID/DE targets. A later step sees values committed by earlier steps;
compilation does not freeze those values.

## Validation rule, including zero-profile steps

Compilation is whole-plan and fail-fast before the first fit. For **every**
step, including a skipped one, it checks all semantics independent of the
active profile population: syntax and version-specific defaults; unique names
and legal inheritance; every declared selector against the sealed model;
model-owned and estimation authority; every constraint reference in its target
context; each action in declaration order even if a later action overrides it;
the *final* effective constraint graph for cycles; search syntax and numeric
shape; and search targets' global identity and role. Thus a skipped step does
not excuse an unknown selector, invalid reference, protected target, final
cycle, GRID declaration with no global FIT target, or DE target that is not a
unique global FIT coordinate. DE duplicate targets and configured physical
bounds are also profile-independent checks. A transient cycle removed by a
later action is not a final cycle. Inherited roles are resolved even when their
source step is skipped.

The global cycle graph includes final Method constraints and effective model or
baseline derivations, even when no profile is selected. A cycle containing a
Method constraint reports its source location and lists the model derivation
members. This check shares the sealed-expression compiler used by active
parameterization and reads no current numerical values.

For a step with **at least one selected profile**, compilation additionally
computes active parameter requirements and constraint dependency closure. It
checks the active-scope conditions that apply to search declarations, derives
the final independent FIT IDs, and projects each GRID/DE declaration to that
scope. A globally valid role action with no active parameter remains inert, as
it is today; it is not an active-scope error. GRID's later-wins override
applies before bounds checks on the winning concrete axes. The compiler rejects
GRID/DE targets with no active final FIT coordinate and all other active-scope
errors before any step executes. Active profile order is frozen at this point.

For a step with **zero selected profiles**, compilation emits the skipped step
after global checks. It does **not** build an active parameterization, project
GRID axes, or require GRID/DE targets to resolve to active FIT coordinates.
Bounds checks that depend on winning *active* GRID axes do not apply. This is
deliberate: the active set is empty, so those assertions have no subject.

The compiler must reproduce current profile selection and ordering across
steps, including the existing active/filtered reordering, until numerical
compatibility evidence permits an explicit change. Initial filtering can still
change observation masks at runtime; it cannot change which profiles were
loaded or selected. If a future filter can change membership, compilation must
move after that filter or reject such mutation rather than silently using a
stale plan.

## Migration design and deletion gates

1. **Characterization.** Capture current v1/v2 and Python behavior through
   final selected profiles, roles, constraints, search coordinates, diagnostics,
   artifacts, commits, and numerical results. Include skipped steps that
   inherit roles, overridden invalid actions, later-wins GRID axes, and failure
   in a later step. These checks are oracles, not a second implementation.
2. **Integrated cutover.** Make the global Method resolver and active compiler
   one semantic pipeline. Project concrete profiles and typed search
   instructions, compile a value-independent constraint program, then switch
   `run_fit` and `run_methods` to compile before `run_info`, output invalidation,
   or any fit. Keep value binding at each active step. Delete production use of
   `MethodPlan.effective_role_actions`, execution-time
   `Experiments.select_profiles(step.selection)`,
   `AnalysisSession.compile_parameterization_from_actions`, the method-specific
   `_build_action_rules`/`_role_for` path, and execution-time
   `resolve_grid_axes`/`resolve_de_coordinates`. The model-only validator and
   standalone Python previews delegate to the same global semantic resolver;
   they do not retain independent interpretations. Keep simulation and model
   semantics separate. The cutover must pass the deletion test below before it
   is considered complete.

The final deletion test is structural: no production executor or numerical
module imports `StepPlan`, Method selectors/actions, `GridSearch`, or `DeSearch`,
and no fit execution or numerical runtime path accepts a Method expression or
selector string. No second function can independently decide the effective
Method role, constraint reference, selection, or search coordinate. A typed
`ExecutableStep` fails this test if it contains `StepPlan` or forwards to old
interpreters, regardless of its line count.

`MethodPlan`/`StepPlan`, `ProfileSelection`, and `GridSearch`/`DeSearch` remain
input types only. `ConstraintProgram` remains a useful restricted numerical
program, compiled once per active step. `ActiveParameterization` becomes its
per-step, snapshot-bound occurrence rather than a Method interpreter. The
resolved GRID/DE coordinate records can serve the numeric interface if they
contain no source-language behavior. A separate validation result mirroring
the executable plan, or an adapter that reconstructs `Method` for fitting,
has no place in the end state.

`MethodPlan.effective_role_actions()` is removed because it independently
interpreted inheritance and had no production caller. Standalone compatibility
previews such as `resolve_grid_axes()` and
`compile_active_parameterization_from_actions()` remain callable for Python
consumers and tests. They reuse the resolver or projection functions; GRID's
standalone preview assumes the caller has already established global FIT
eligibility. Structural tests forbid their use by fitting modules. The
`Experiments.select_profiles()` compatibility method shares the compiler's
selection projection but is likewise outside the fit path.

## Compatibility and rejected alternatives

v1's implicit selection/role inheritance and duplicate-section behavior, v2's
explicit `ROLES_FROM` and step-local selection, and the Python `Methods` mapping
remain input rules. Their adapters feed the one compiler; diagnostics retain
source locations and, where possible, existing wording. The existing
`MethodPlan.validate(model)` interface can remain as a model-only preflight
through the compiler's global pass; it must not claim to validate active-scope
execution and must not maintain separate semantic rules. The renderable input
plan and archived inputs remain the provenance authority. A restart rebuilds
the executable plan against the sealed model and loaded population, checks any
recorded model/population/plan fingerprint, then binds current saved values when
each step starts. Fingerprint recording would need an additive, versioned
provenance change. An executable plan is not a portable serialized checkpoint.

Keeping whole-plan validation plus a second execution-time interpreter was
rejected: it can accept one interpretation and run another. Prebinding values
for all steps was rejected because earlier fits change the starting state of
later steps. Requiring active FIT coordinates for a skipped step was rejected
because no active scope exists; skipping all its validation was rejected because
invalid global semantics can affect later inheritance and reproducibility.

Compilation before the first fit intentionally changes failure timing: an
invalid later step will prevent earlier valid steps from committing, unlike
today's execution-time active-scope failures. This is a compatibility change to
report explicitly, including effects on partial outputs and restart files.

This migration crosses configuration, parameters, profile selection, and
optimization, so disruption is high. It earns that cost only if the listed
interpretation paths are actually removed and numerical/diagnostic behavior is
verified. The main benefit is correctness and maintainability; earlier error
discovery and simpler restart provenance are secondary. No numerical algorithm
change is intended.

---
status: accepted
---

# Translate exceptions at the CLI boundary

ChemEx will keep internal Python interfaces exception-transparent while one
lightweight terminal boundary owns command-line failure translation.  The
installed `chemex` command and `python -m chemex` share that boundary.  It
covers application imports reached from the boundary as well as failures raised
during command execution.

`SystemExit` remains native so argparse and intentional command exits retain
their established status and stream behavior. The boundary applies exactly
three terminal policies: a known ChemEx failure is rendered concisely and exits
1, a user interruption is rendered without a traceback and exits 130, and an
unexpected internal failure is reported as such and exits 1. User-input,
scientific, numerical, operational, and internal categories are useful for
reasoning about failures, but do not form a parallel exception framework.

A small `ChemExError` marker identifies only failures whose structured
information ChemEx deliberately produced and intends to present to users.
Scientific and domain layers raise typed failures without Rich rendering or
process exits. The console entry point owns final error presentation and process
translation. Intermediate orchestration may add narrow operation context and
verified publication paths without introducing a generic context dictionary or
artifact registry.

Default diagnostics for unexpected internal failures never expose arbitrary
exception contents, notes, causes, or tracebacks. They explicitly direct the
user to rerun with the global `--debug` option. In debug mode an unexpected
exception is not translated: its natural traceback, type, message, notes, and
cause or context chain remain visible. Known ChemEx failures and interruptions
retain their normal semantic classification in debug mode. Rendering failure
falls back to plain stderr without replacing the original exit status.

Workflow layers may catch interruption temporarily to freeze scientifically
valid partial state, serialize already committed values or qualified Evidence,
write truthful incomplete status, and perform required cleanup.  Once
interruption is recognized, finalization does not start further scientific
calculation, sampling, replicates, Method Steps, or optional plotting.  The
interruption remains distinguishable from ordinary failure and propagates to
the terminal boundary after finalization. The final diagnostic includes a stage
or artifact path only when orchestration knows the stage or the artifact was
successfully published; directory existence alone is not evidence of
publication.

## Considered options

- Catching failures in every command or scientific layer was rejected because
  it duplicates presentation policy and risks converting failures into success.
- Immediately re-raising every raw `KeyboardInterrupt` was rejected because it
  can discard qualified Evidence and committed results before truthful
  publication.
- A category-shaped global exception hierarchy was rejected because the CLI
  needs only one safe-to-present marker. Existing domain-specific failures keep
  their useful types, while unclassified exceptions remain unexpected.
- Always exposing tracebacks was rejected because ordinary CLI use must not
  disclose arbitrary exception contents. Permanently suppressing tracebacks was
  also rejected because it prevents useful bug reports; explicit `--debug`
  provides the diagnostic path without changing the safe default.
- Leaving application imports outside the boundary was rejected because common
  dependency and startup failures would still expose tracebacks.

## Consequences

Direct calls to `chemex.chemex.main()` and lower-level APIs continue to raise
their original exceptions. Expected domain failures no longer exit the process
or render terminal errors below the console entry point. Workflow-specific
incomplete exceptions may carry their existing normalized `terminal` value and
verified publication paths through orchestration. Completed
deterministic fits, qualified uncertainty Evidence, complete MCMC transitions,
and completed resampling replicates survive a later interruption without being
promoted to complete scientific conclusions.

The boundary can only translate failures after Python has imported the ChemEx
package and reached it.  Interpreter startup failures, an unavailable package,
bootstrap syntax errors, native process crashes, and forced process termination
remain outside this guarantee.
